"""Alignment engines behind one interface (Dependency-Inversion boundary).

The GUI and the ambiguity tools depend only on ``AlignerEngine``; concrete
engines (the built-in progressive aligner, plus subprocess wrappers for MAFFT /
MUSCLE / Clustal Omega / ProbCons / PRANK / ClustalW) are interchangeable. A
``CodonAware`` decorator turns any engine into a reading-frame-preserving codon
aligner via translate -> align amino acids -> back-translate.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

from . import accel, progressive
from .domain import Alignment, Alphabet, CodingSpec, translate_codon

Record = Tuple[str, str]


def _kind(alphabet: Alphabet) -> str:
    return "nucleotide" if alphabet.is_nucleotide else "protein"


@dataclass
class Param:
    """A tunable parameter an engine exposes to the settings dialog and the CLI.

    A numeric parameter whose ``default`` is ``None`` is left to the tool: CRAIC
    passes no flag at all, and ``blank`` says what the tool does instead. That
    is for defaults CRAIC cannot state as one number — PRANK's and ClustalW's gap
    costs differ between DNA and protein, and the built-in engine estimates its
    gap probabilities from the data. Writing one number in for those would
    override the tool with a value that is right for neither alphabet.

    ``alphabet`` marks a parameter that only means something for one kind of
    data ("protein" or "nucleotide"; empty for both), so the settings dialog can
    show only what applies to the sequences the engine will actually be handed.
    """

    key: str
    label: str
    kind: str                       # "float" | "int" | "choice" | "bool"
    default: object
    choices: Optional[List[str]] = None
    help: str = ""
    #: Inclusive bounds for a numeric parameter, where the tool it is passed to
    #: has them. The dialog refuses values outside these rather than letting the
    #: aligner reject them after the run has started.
    lo: Optional[float] = None
    hi: Optional[float] = None
    #: What leaving a numeric parameter blank means (only when ``default`` is None).
    #: A dict keyed "protein" / "nucleotide" when the tool's default differs.
    blank: Union[str, Dict[str, str]] = ""
    #: "" (any data), "protein" or "nucleotide".
    alphabet: str = ""

    def applies_to(self, alphabet: Optional[Alphabet]) -> bool:
        """Does this parameter mean anything for that data? True when unknown."""
        if not self.alphabet or alphabet is None:
            return True
        return self.alphabet == _kind(alphabet)

    def blank_for(self, alphabet: Optional[Alphabet] = None) -> str:
        """What a blank value means — for this alphabet, when it is known."""
        if isinstance(self.blank, str):
            return self.blank
        if alphabet is not None:
            return self.blank[_kind(alphabet)]
        return "; ".join(f"{kind}: {text}" for kind, text in self.blank.items())

    def coerce(self, raw: str):
        """Read a value written as text (on the command line). ValueError if not allowed."""
        raw = raw.strip()
        if self.kind == "choice":
            choices = self.choices or []
            if raw in choices:
                return raw
            # "L-INS-i" for "L-INS-i (accurate, local)": the first word will do
            short = [c for c in choices if c.split()[0] == raw]
            if len(short) == 1:
                return short[0]
            raise ValueError(f"{self.key} must be one of: {', '.join(choices)}")
        if self.kind == "bool":
            low = raw.lower()
            if low in ("1", "true", "yes", "on"):
                return True
            if low in ("0", "false", "no", "off"):
                return False
            raise ValueError(f"{self.key} must be true or false")
        if raw == "" and self.default is None:
            return None
        value = int(raw) if self.kind == "int" else float(raw)
        if (self.lo is not None and value < self.lo) or (self.hi is not None and value > self.hi):
            raise ValueError(f"{self.key} must be between {self.lo} and {self.hi}")
        return value


class AlignerEngine(ABC):
    key: str = "engine"
    label: str = "engine"

    @abstractmethod
    def available(self) -> bool: ...

    @abstractmethod
    def align(self, records: Sequence[Record], alphabet: Alphabet, **opts) -> Alignment: ...

    def parameters(self) -> List[Param]:
        """Tunable parameters exposed in the GUI settings dialog. Default: none."""
        return []

    def supports(self, alphabet: Alphabet) -> bool:
        """Can this engine align that alphabet? Most can align anything.

        Asked *before* a run is started, so an engine that cannot take the data
        is refused with an explanation rather than by a traceback from inside a
        worker thread.
        """
        return True

    def __repr__(self) -> str:  # pragma: no cover
        return f"<{self.__class__.__name__} {self.key}>"


# --------------------------------------------------------------------------- #
# Built-in
# --------------------------------------------------------------------------- #

class BuiltinProgressive(AlignerEngine):
    key = "builtin"
    label = "CRAIC built-in (progressive)"
    consistency_mem_gb = 1.0   # RAM budget (GB proxy) for the accurate consistency path;
                               # families larger than this stream. Real peak RAM ~2-2.5x this.
    #: Effort used when the caller names none. The GUI and CLI pass one explicitly;
    #: the benchmark harness sets this so that a result set records, and does not
    #: silently choose, the compute level it measured.
    benchmark_effort = "med"

    def available(self) -> bool:
        return True

    def align(self, records, alphabet, effort: Optional[str] = None, consistency_iters=None,
              refine_iters=None, delta=None, epsilon=None, matrix=None,
              estimate: bool = True, progress=None, cancelled=None, **opts):
        effort = self.benchmark_effort if effort is None else effort
        return progressive.align(records, alphabet, effort=effort, delta=delta,
                                 epsilon=epsilon, matrix=matrix, estimate=estimate,
                                 consistency_iters=consistency_iters,
                                 refine_iters=refine_iters,
                                 consistency_mem_gb=self.consistency_mem_gb,
                                 progress=progress, cancelled=cancelled)

    def parameters(self):
        return [
            Param("effort", "Compute level", "choice", "med", ["min", "med", "max"],
                  help="min = fast streaming posterior decoding (scales to large "
                       "alignments); med = consistency transformation + posterior "
                       "guide tree (default, most accurate per unit time); "
                       "max = also runs iterative refinement (slower)."),
            Param("delta", "Gap open probability (δ)", "float", None, lo=0.0001, hi=0.49,
                  blank="estimated from the data",
                  help="Pair-HMM probability of opening a gap. Higher means gaps are "
                       "cheaper, so more and shorter gaps. Blank estimates it from a "
                       "quick pilot alignment."),
            Param("epsilon", "Gap extension probability (ε)", "float", None, lo=0.0, hi=0.99,
                  blank="estimated from the data",
                  help="Pair-HMM probability that a gap continues. Higher means "
                       "longer gaps. Blank estimates it from a quick pilot alignment."),
            Param("matrix", "Protein matrix", "choice", accel.PROTEIN_MATRIX,
                  list(accel.PROTEIN_MATRICES), alphabet="protein",
                  help="Substitution matrix behind the protein emission model. "
                       "Lower numbers suit more divergent sequences."),
            Param("consistency_iters", "Consistency passes", "int", None, lo=0, hi=10,
                  blank="set by the compute level",
                  help="Rounds of the ProbCons consistency transformation. 0 turns "
                       "it off."),
            Param("refine_iters", "Refinement passes", "int", None, lo=0, hi=1000,
                  blank="set by the compute level",
                  help="Rounds of iterative refinement after the progressive "
                       "alignment. Needs consistency to be on."),
        ]


# --------------------------------------------------------------------------- #
# External wrappers
# --------------------------------------------------------------------------- #

_SEARCH_PATH: Optional[str] = None


def _conda_bins() -> List[str]:
    """``bin`` directories of conda installations found in the usual places.

    ``CONDA_PREFIX`` is set only by an *activated* shell, so a Finder-launched
    app sees nothing of a conda installation -- which is exactly where bioconda
    puts MAFFT, MUSCLE, ProbCons, PRANK and Clustal Omega. Look for the
    installations themselves rather than relying on the environment variable.
    """
    import glob

    out: List[str] = []
    for root in ("miniforge3", "mambaforge", "miniconda3", "anaconda3", "opt/anaconda3"):
        base = os.path.expanduser(os.path.join("~", root))
        out.append(os.path.join(base, "bin"))
        out.extend(sorted(glob.glob(os.path.join(base, "envs", "*", "bin"))))
    return out


def _search_path() -> str:
    """PATH to search for external aligners. A GUI / .app launch on macOS inherits
    only a minimal PATH (``/usr/bin:/bin:...``), so tools installed via conda,
    Homebrew, MacPorts or pip go undetected. Augment it with the usual locations."""
    global _SEARCH_PATH
    if _SEARCH_PATH is None:
        extra = [
            os.path.join(sys.prefix, "bin"),                       # active venv / conda env
            os.path.join(os.environ.get("CONDA_PREFIX", ""), "bin"),
            "/opt/homebrew/bin", "/usr/local/bin",                 # Homebrew (arm / intel)
            os.path.expanduser("~/.local/bin"),
            "/opt/local/bin",                                      # MacPorts
        ] + _conda_bins()
        seen, dirs = set(), []
        for d in os.environ.get("PATH", "").split(os.pathsep) + extra:
            if d and d not in seen and os.path.isdir(d):
                seen.add(d)
                dirs.append(d)
        _SEARCH_PATH = os.pathsep.join(dirs)
    return _SEARCH_PATH


def _which(binary: str) -> Optional[str]:
    """Locate an external tool on the augmented search path (absolute path or None)."""
    return shutil.which(binary, path=_search_path())


class ExternalAligner(AlignerEngine):
    binary: str = ""
    timeout = None   # per-call subprocess timeout (s); None = no limit (set by benchmarks)

    def available(self) -> bool:
        return _which(self.binary) is not None

    def _command(self, infile: str, outfile: str, alphabet: Alphabet, **opts) -> Tuple[List[str], str]:
        """Return (argv, output_mode) where mode is 'stdout' or 'file'."""
        raise NotImplementedError

    @classmethod
    def _run(cls, argv: List[str], timeout=None):
        """Run an external aligner on the same augmented PATH that found it.

        ``available()`` resolves the binary against ``_search_path()``, so a
        subprocess launched against the *inherited* PATH can fail to find a tool
        the engine has just reported as available — which is exactly what
        happens in a macOS ``.app`` launch, the case ``_search_path`` exists to
        fix. Every external invocation goes through here so the two cannot drift
        apart again.
        """
        exe = _which(argv[0]) or argv[0]
        return subprocess.run(argv, capture_output=True, text=True, timeout=timeout,
                              executable=exe, env={**os.environ, "PATH": _search_path()})

    def align(self, records, alphabet, **opts):
        from .io import read_records, write_fasta

        opts.pop("progress", None)            # external aligners run as one subprocess;
        opts.pop("cancelled", None)           # they don't report progress or cancel mid-run
        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            outfile = os.path.join(d, "out.fasta")
            safe = [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(records)]
            write_fasta(infile, [x[0] for x in safe], [x[1] for x in safe], wrap=0)
            argv, mode = self._command(infile, outfile, alphabet, **opts)
            proc = self._run(argv, self.timeout)
            if proc.returncode != 0:
                raise RuntimeError(f"{self.label} failed: {proc.stderr[:400]}")
            if mode == "stdout":
                with open(outfile, "w") as fh:
                    fh.write(proc.stdout)
            parsed = read_records(outfile, "fasta")
            # restore original ids by position
            idmap = {f"s{i}": rec[0] for i, rec in enumerate(records)}
            parsed = [(idmap.get(pid, pid), seq.upper()) for pid, seq in parsed]
            order = {rec[0]: i for i, rec in enumerate(records)}
            parsed.sort(key=lambda p: order.get(p[0], 0))
            return Alignment.from_records(parsed, alphabet=alphabet)


class MafftEngine(ExternalAligner):
    key = "mafft"
    label = "MAFFT"
    binary = "mafft"

    _STRATEGY = {
        "auto": ["--auto"],
        "L-INS-i (accurate, local)": ["--localpair", "--maxiterate", "1000"],
        "G-INS-i (accurate, global)": ["--globalpair", "--maxiterate", "1000"],
        "E-INS-i (long gaps)": ["--genafpair", "--maxiterate", "1000"],
        "FFT-NS-2 (fast)": ["--retree", "2"],
    }

    #: Protein scoring matrices, as MAFFT's flags spell them.
    _MATRIX = {
        "BLOSUM62": ["--bl", "62"],
        "BLOSUM30": ["--bl", "30"],
        "BLOSUM45": ["--bl", "45"],
        "BLOSUM80": ["--bl", "80"],
        "JTT (PAM 200)": ["--jtt", "200"],
        "Transmembrane (PAM 200)": ["--tm", "200"],
    }

    def parameters(self):
        return [
            Param("strategy", "Strategy", "choice", "auto", list(self._STRATEGY),
                  help="MAFFT preset / iterative-refinement strategy."),
            Param("op", "Gap open (--op)", "float", 1.53, lo=0.0,
                  help="Gap opening penalty for the progressive and refinement "
                       "stages. MAFFT's default is 1.53."),
            Param("ep", "Gap extend (--ep)", "float", 0.0,
                  help="Offset that works like a gap extension penalty. MAFFT's "
                       "default is 0.0."),
            Param("lop", "Local-pair gap open (--lop)", "float", -2.0,
                  help="Gap opening penalty in the local pairwise stage. Used only "
                       "by L-INS-i and E-INS-i (and by auto when it picks one). "
                       "MAFFT's default is -2.00."),
            Param("lep", "Local-pair gap extend (--lep)", "float", 0.1,
                  help="Offset in the local pairwise stage. Used only by L-INS-i "
                       "and E-INS-i. MAFFT's default is 0.1."),
            Param("maxiterate", "Refinement cycles (--maxiterate)", "int", None,
                  lo=0, hi=1000, blank="set by the strategy",
                  help="Rounds of iterative refinement. The *-INS-i strategies use "
                       "1000, FFT-NS-2 uses 0."),
            Param("matrix", "Protein matrix", "choice", "BLOSUM62", list(self._MATRIX),
                  alphabet="protein",
                  help="Amino-acid scoring matrix. Lower BLOSUM numbers suit more "
                       "divergent sequences."),
            Param("kimura", "DNA model (--kimura)", "choice", "200 PAM",
                  ["200 PAM", "20 PAM", "1 PAM"], alphabet="nucleotide",
                  help="Nucleotide scoring matrix (Kimura, κ = 2). 1 PAM suits very "
                       "close sequences."),
        ]

    def _command(self, infile, outfile, alphabet, **opts):
        flags = list(self._STRATEGY.get(opts.get("strategy", "auto"), ["--auto"]))
        if opts.get("maxiterate") is not None:
            if "--maxiterate" in flags:
                i = flags.index("--maxiterate")
                del flags[i:i + 2]
            flags += ["--maxiterate", str(int(opts["maxiterate"]))]
        argv = ["mafft", *flags,
                "--op", str(float(opts.get("op", 1.53))),
                "--ep", str(float(opts.get("ep", 0.0))),
                "--lop", str(float(opts.get("lop", -2.0))),
                "--lep", str(float(opts.get("lep", 0.1)))]
        if alphabet.is_nucleotide:
            argv += ["--kimura", str(opts.get("kimura", "200 PAM")).split()[0]]
        else:
            argv += self._MATRIX.get(opts.get("matrix", "BLOSUM62"), [])
        return argv + ["--quiet", infile], "stdout"


class MuscleEngine(ExternalAligner):
    key = "muscle"
    label = "MUSCLE"
    binary = "muscle"

    def parameters(self):
        # MUSCLE 5 has no gap penalties on the command line. What it does have
        # is the pair of switches Edgar (2022) uses to build alignment ensembles:
        # a guide-tree permutation and a random perturbation of its HMM.
        return [
            Param("algorithm", "Algorithm", "choice", "align", ["align", "super5"],
                  help="align = full PPP algorithm (default). super5 = the faster "
                       "algorithm for a few hundred sequences or more."),
            Param("perm", "Guide-tree permutation (-perm)", "choice", "none",
                  ["none", "abc", "acb", "bca"],
                  help="Which of the three ways of joining the guide tree's root "
                       "subtrees to use. Changing it gives a different, equally "
                       "valid alignment."),
            Param("perturb", "Perturbation seed (-perturb)", "int", 0, lo=0, hi=1000000,
                  help="0 = no perturbation (default). Any other number perturbs "
                       "the HMM parameters with that random seed, giving one "
                       "replicate from MUSCLE's ensemble."),
        ]

    def align(self, records, alphabet, **opts):
        from .io import read_records, write_fasta

        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            outfile = os.path.join(d, "out.fasta")
            safe = [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(records)]
            write_fasta(infile, [x[0] for x in safe], [x[1] for x in safe], wrap=0)
            algorithm = "-super5" if opts.get("algorithm") == "super5" else "-align"
            v5 = ["muscle", algorithm, infile, "-output", outfile]
            if opts.get("perm", "none") != "none":
                v5 += ["-perm", str(opts["perm"])]
            if int(opts.get("perturb", 0) or 0):
                v5 += ["-perturb", str(int(opts["perturb"]))]
            # try v5 syntax, then v3. A v5 attempt that fails to even start
            # (wrong binary, timeout) must not abort the v3 attempt, so the
            # failure is recorded and the loop continues. v3 has none of the
            # v5 options, so they are dropped there.
            last_error = "no attempt ran"
            for argv in (v5, ["muscle", "-in", infile, "-out", outfile]):
                try:
                    proc = self._run(argv, self.timeout)
                except (OSError, subprocess.SubprocessError) as exc:
                    last_error = str(exc)
                    continue
                if proc.returncode == 0 and os.path.exists(outfile):
                    break
                last_error = proc.stderr[:400]
            else:
                raise RuntimeError(f"MUSCLE failed: {last_error}")
            parsed = read_records(outfile, "fasta")
            idmap = {f"s{i}": rec[0] for i, rec in enumerate(records)}
            parsed = [(idmap.get(pid, pid), seq.upper()) for pid, seq in parsed]
            order = {rec[0]: i for i, rec in enumerate(records)}
            parsed.sort(key=lambda p: order.get(p[0], 0))
            return Alignment.from_records(parsed, alphabet=alphabet)


class ClustalOmegaEngine(ExternalAligner):
    """Clustal Omega aligns profile HMMs to each other (HHalign), so unlike
    ClustalW it has no gap penalties to set — the gap model is part of the HMMs
    it builds. What it exposes is how the guide tree is made and how often the
    tree and the HMMs are re-estimated."""

    key = "clustalo"
    label = "Clustal Omega"
    binary = "clustalo"

    def parameters(self):
        return [
            Param("iterations", "Combined iterations", "int", 0, lo=0, hi=100,
                  help="Guide-tree/HMM refinement iterations; 0 = default."),
            Param("full", "Full distance matrix (--full)", "bool", False,
                  help="Compute every pairwise distance for the guide tree instead "
                       "of the mBed approximation. Slower; matters for small sets."),
            Param("max_guidetree_iterations", "Max guide-tree iterations", "int", None,
                  lo=0, hi=100, blank="same as combined iterations",
                  help="Cap on how many of the combined iterations rebuild the "
                       "guide tree."),
            Param("max_hmm_iterations", "Max HMM iterations", "int", None,
                  lo=0, hi=100, blank="same as combined iterations",
                  help="Cap on how many of the combined iterations rebuild the HMM."),
        ]

    def _command(self, infile, outfile, alphabet, **opts):
        argv = ["clustalo", "-i", infile, "-o", outfile, "--force", "--outfmt=fasta"]
        it = int(opts.get("iterations", 0))
        if it > 0:
            argv += ["--iterations", str(it)]
        if opts.get("full"):
            argv.append("--full")
        for key in ("max_guidetree_iterations", "max_hmm_iterations"):
            if opts.get(key) is not None:
                argv.append(f"--{key.replace('_', '-')}={int(opts[key])}")
        return argv, "file"


class ProbConsEngine(ExternalAligner):
    """ProbCons (Do et al., 2005), the method CRAIC's built-in engine follows.

    Present so that the built-in engine can be compared with the thing it is a
    reimplementation of, rather than only with aligners of a different lineage.
    Without it, a benchmark table reporting the built-in engine below MAFFT
    invites the reader to conclude that posterior-consistency alignment is the
    weaker approach, when the honest question is how much of the gap is the
    method and how much is this implementation of it.

    Protein only: ProbCons has no nucleotide model. Its gap probabilities are
    read from a parameter file rather than set on the command line, so they are
    not offered here.
    """

    key = "probcons"
    label = "ProbCons"
    binary = "probcons"

    def supports(self, alphabet: Alphabet) -> bool:
        return not alphabet.is_nucleotide

    def align(self, records, alphabet, **opts):
        if not self.supports(alphabet):
            raise RuntimeError(
                "ProbCons aligns protein sequences only. For coding nucleotides, "
                "tick 'align as protein' to translate, align and back-translate; "
                "otherwise choose another engine.")
        return super().align(records, alphabet, **opts)

    def parameters(self):
        return [
            Param("consistency", "Consistency passes (-c)", "int", 2, lo=0, hi=5,
                  help="Rounds of the consistency transformation. ProbCons accepts "
                       "0 to 5; the published default is 2."),
            Param("iterations", "Refinement iterations", "int", 100,
                  lo=0, hi=1000,
                  help="Rounds of random-bipartition iterative refinement "
                       "(ProbCons -ir). ProbCons accepts 0 to 1000; the "
                       "published default is 100."),
            Param("pretraining", "Pre-training rounds (-pre)", "int", 0, lo=0, hi=20,
                  help="Rounds of EM re-estimation of the pair-HMM on these "
                       "sequences before aligning. ProbCons accepts 0 to 20; the "
                       "default is 0."),
        ]

    def _command(self, infile, outfile, alphabet, **opts):
        argv = ["probcons"]
        for key, flag, default in (("consistency", "-c", 2), ("iterations", "-ir", 100),
                                   ("pretraining", "-pre", 0)):
            value = int(opts.get(key, default))
            if value != default:
                argv += [flag, str(value)]
        return argv + [infile], "stdout"


class PrankEngine(ExternalAligner):
    key = "prank"
    label = "PRANK (phylogeny-aware)"
    binary = "prank"

    def parameters(self):
        # The gap costs are left blank by default because PRANK's own defaults
        # depend on the data (DNA 0.025 / 0.75, protein 0.005 / 0.5). Passing
        # one pair for both used to override PRANK with values right for neither.
        return [
            Param("gaprate", "Gap opening rate", "float", None, lo=0.0, hi=1.0,
                  blank={"protein": "PRANK default, 0.005",
                         "nucleotide": "PRANK default, 0.025"},
                  help="PRANK -gaprate. A rate, so 0 to 1."),
            Param("gapext", "Gap extension prob.", "float", None, lo=0.0, hi=1.0,
                  blank={"protein": "PRANK default, 0.5",
                         "nucleotide": "PRANK default, 0.75"},
                  help="PRANK -gapext. A probability, so 0 to 1."),
            Param("F", "Keep insertions as insertions (+F)", "bool", False,
                  help="PRANK +F: a residue once inferred to be an insertion is "
                       "never matched again further up the tree. Recommended by "
                       "PRANK's authors for most data."),
            Param("iterate", "Iterations (-iterate)", "int", 5, lo=1, hi=100,
                  help="Rounds of re-alignment, keeping the best. PRANK's default "
                       "is 5."),
            Param("termgap", "Penalise terminal gaps (-termgap)", "bool", False,
                  help="By default PRANK does not penalise gaps at the ends."),
        ]

    def align(self, records, alphabet, **opts):
        from .io import read_records, write_fasta

        with tempfile.TemporaryDirectory() as d:
            infile = os.path.join(d, "in.fasta")
            prefix = os.path.join(d, "out")
            safe = [(f"s{i}", s.replace("-", "")) for i, (_, s) in enumerate(records)]
            write_fasta(infile, [x[0] for x in safe], [x[1] for x in safe], wrap=0)
            argv = self._argv(infile, prefix, alphabet, **opts)
            proc = self._run(argv, self.timeout)
            best = prefix + ".best.fas"
            if not os.path.exists(best):
                best = prefix + ".2.fas"
            if not os.path.exists(best):
                raise RuntimeError(f"PRANK failed: {proc.stderr[:400]}")
            parsed = read_records(best, "fasta")
            idmap = {f"s{i}": rec[0] for i, rec in enumerate(records)}
            parsed = [(idmap.get(pid, pid), seq.upper()) for pid, seq in parsed]
            order = {rec[0]: i for i, rec in enumerate(records)}
            parsed.sort(key=lambda p: order.get(p[0], 0))
            return Alignment.from_records(parsed, alphabet=alphabet)

    @staticmethod
    def _argv(infile, prefix, alphabet, **opts) -> List[str]:
        argv = ["prank", f"-d={infile}", f"-o={prefix}"]
        for key in ("gaprate", "gapext"):
            if opts.get(key) is not None:
                argv.append(f"-{key}={float(opts[key])}")
        if opts.get("F"):
            argv.append("+F")
        if int(opts.get("iterate", 5)) != 5:
            argv.append(f"-iterate={int(opts['iterate'])}")
        if opts.get("termgap"):
            argv.append("-termgap")
        if alphabet.is_nucleotide:
            argv.append("-DNA")
        return argv


class ClustalWEngine(ExternalAligner):
    """ClustalW 2 (Larkin et al., 2007), for teaching gap costs.

    Its successors model gaps inside profile HMMs, where there is no single
    penalty to turn. ClustalW keeps the textbook model — a cost to open a gap and
    a cost to extend it, separately for the pairwise stage that builds the guide
    tree and for the progressive stage — plus the position-specific adjustments
    that made it work in practice. Those are what a student needs to be able to
    change and watch.

    The gap costs are blank by default because ClustalW's defaults depend on the
    data (protein 10 / 0.2, DNA 15 / 6.66); CRAIC then lets ClustalW choose.
    """

    key = "clustalw"
    label = "ClustalW"
    binary = "clustalw2"

    def available(self) -> bool:
        return self._binary() is not None

    def _binary(self) -> Optional[str]:
        """Bioconda and Homebrew install ``clustalw2``; Debian installs ``clustalw``."""
        for name in ("clustalw2", "clustalw"):
            if _which(name):
                return name
        return None

    def parameters(self):
        return [
            Param("gapopen", "Gap open (multiple)", "float", None, lo=0.0, hi=100.0,
                  blank={"protein": "ClustalW default, 10",
                         "nucleotide": "ClustalW default, 15"},
                  help="Cost of opening a gap in the progressive alignment "
                       "(-GAPOPEN)."),
            Param("gapext", "Gap extension (multiple)", "float", None, lo=0.0, hi=10.0,
                  blank={"protein": "ClustalW default, 0.2",
                         "nucleotide": "ClustalW default, 6.66"},
                  help="Cost of each further gap position in the progressive "
                       "alignment (-GAPEXT)."),
            Param("pwgapopen", "Gap open (pairwise)", "float", None, lo=0.0, hi=100.0,
                  blank={"protein": "ClustalW default, 10",
                         "nucleotide": "ClustalW default, 15"},
                  help="Gap opening cost in the pairwise alignments that build the "
                       "guide tree (-PWGAPOPEN)."),
            Param("pwgapext", "Gap extension (pairwise)", "float", None, lo=0.0, hi=10.0,
                  blank={"protein": "ClustalW default, 0.1",
                         "nucleotide": "ClustalW default, 6.66"},
                  help="Gap extension cost in the pairwise alignments that build "
                       "the guide tree (-PWGAPEXT)."),
            Param("matrix", "Protein matrix series", "choice", "GONNET",
                  ["GONNET", "BLOSUM", "PAM", "ID"], alphabet="protein",
                  help="ClustalW picks a matrix from the series to suit each "
                       "alignment step's divergence. ID scores identities only."),
            Param("dnamatrix", "DNA matrix", "choice", "IUB", ["IUB", "CLUSTALW"],
                  alphabet="nucleotide",
                  help="IUB scores ambiguity codes; CLUSTALW is identity."),
            Param("gapdist", "Gap separation distance", "int", 4, lo=0, hi=100,
                  alphabet="protein",
                  help="Gaps closer than this to an existing gap are penalised, "
                       "which spreads gaps apart (-GAPDIST)."),
            Param("nopgap", "Residue-specific gap costs off", "bool", False,
                  alphabet="protein",
                  help="By default a gap costs less next to some residues than "
                       "others (-NOPGAP turns that off)."),
            Param("nohgap", "Hydrophilic gap reduction off", "bool", False,
                  alphabet="protein",
                  help="By default gaps are cheaper in runs of hydrophilic "
                       "residues, which tend to be loops (-NOHGAP turns that off)."),
        ]

    def _command(self, infile, outfile, alphabet, **opts):
        argv = [self._binary() or self.binary, f"-INFILE={infile}", f"-OUTFILE={outfile}",
                "-OUTPUT=FASTA", "-OUTORDER=INPUT", "-QUIET",
                "-TYPE=DNA" if alphabet.is_nucleotide else "-TYPE=PROTEIN"]
        for key in ("gapopen", "gapext", "pwgapopen", "pwgapext"):
            if opts.get(key) is not None:
                argv.append(f"-{key.upper()}={float(opts[key])}")
        if alphabet.is_nucleotide:
            argv.append(f"-DNAMATRIX={opts.get('dnamatrix', 'IUB')}")
        else:
            argv.append(f"-MATRIX={opts.get('matrix', 'GONNET')}")
            argv.append(f"-GAPDIST={int(opts.get('gapdist', 4))}")
            if opts.get("nopgap"):
                argv.append("-NOPGAP")
            if opts.get("nohgap"):
                argv.append("-NOHGAP")
        return argv, "file"


# --------------------------------------------------------------------------- #
# Codon-aware decorator (translate -> align aa -> back-translate)
# --------------------------------------------------------------------------- #

class CodonAware(AlignerEngine):
    """Wrap any engine to align coding nucleotides in amino-acid space.

    Assumes reading frame 0 and intact frames (no internal frameshifts). The
    result is a codon-aware nucleotide alignment carrying a CodingSpec, so the
    nt/codon/aa views are all immediately valid.

    A sequence whose length is not a multiple of three ends in an incomplete
    codon. Its last one or two nucleotides are kept, in columns after the last
    codon: outside the reading frame, so the codon and amino-acid views leave
    them out, but still in the alignment. They used to be dropped, which lost
    data and made the result no longer match a reference of the same sequences.
    """

    def __init__(self, inner: AlignerEngine, table: int = 1):
        self.inner = inner
        self.table = table
        self.key = f"codon:{inner.key}"
        self.label = f"{inner.label} + codon-aware"

    def available(self) -> bool:
        return self.inner.available()

    def parameters(self):
        return self.inner.parameters()

    def align(self, records, alphabet, **opts):
        aa_records: List[Record] = []
        codon_lists, tails = {}, {}
        for name, seq in records:
            s = seq.replace("-", "").upper().replace("U", "T")
            n = len(s) // 3
            codons = [s[3 * i : 3 * i + 3] for i in range(n)]
            aa = "".join(translate_codon(c, self.table) for c in codons)
            aa_records.append((name, aa))
            codon_lists[name] = codons
            tails[name] = s[3 * n:]                  # an incomplete last codon
        tail_width = max((len(t) for t in tails.values()), default=0)

        aa_aln = self.inner.align(aa_records, Alphabet.PROTEIN, **opts)

        nt_rows = []
        for name, aa_row in zip(aa_aln.ids, aa_aln.rows):
            codons = codon_lists[name]
            ci = 0
            buf = []
            for ch in aa_row:
                if ch == "-":
                    buf.append("---")
                else:
                    buf.append(codons[ci] if ci < len(codons) else "---")
                    ci += 1
            nt_rows.append("".join(buf) + tails[name].ljust(tail_width, "-"))
        aln = Alignment(list(aa_aln.ids), nt_rows, Alphabet.DNA,
                        coding=CodingSpec(frame=0, table=self.table))
        order = {rec[0]: i for i, rec in enumerate(records)}
        pairs = sorted(zip(aln.ids, aln.rows), key=lambda p: order.get(p[0], 0))
        return Alignment([p[0] for p in pairs], [p[1] for p in pairs], Alphabet.DNA,
                         coding=CodingSpec(frame=0, table=self.table))


def changed_sequences(records: Sequence[Record], aln: Alignment) -> List[str]:
    """Ids whose residues differ between the input ``records`` and the aligned
    ``aln`` (ignoring gaps, case, and U written as T). An aligner must only
    insert gaps; anything else means lost or altered data, which should be said
    rather than shown as an alignment."""
    def residues(row):
        return row.replace("-", "").upper().replace("U", "T")

    given = {name: residues(seq) for name, seq in records}
    return [i for i, row in zip(aln.ids, aln.rows) if i in given and residues(row) != given[i]]


# --------------------------------------------------------------------------- #
# Registry
# --------------------------------------------------------------------------- #

class _BuiltinVariant(BuiltinProgressive):
    """Built-in aligner pinned to a particular pair-HMM gap-open probability.

    Used so the disagreement map is informative even when no external aligner is
    installed: different gap models give genuinely different MEA alignments that
    disagree exactly in the hard regions.
    """

    def __init__(self, delta: float, label: str):
        self._delta = delta
        self.label = label
        self.key = f"builtin:d{delta}"

    def align(self, records, alphabet, **opts):
        return progressive.align(records, alphabet, delta=self._delta, epsilon=0.5,
                                 estimate=False, consistency_iters=1,
                                 cancelled=opts.get("cancelled"))


def builtin_variants() -> List[AlignerEngine]:
    # Vary the pair-HMM gap-open probability (delta): a rarer-gap vs freer-gap
    # model produces genuinely different posterior-decoded alignments, so the
    # disagreement map stays informative with no external aligner installed.
    return [
        _BuiltinVariant(0.005, "built-in · rare gaps"),
        _BuiltinVariant(0.03, "built-in · default gaps"),
        _BuiltinVariant(0.1, "built-in · frequent gaps"),
    ]


def all_engines() -> List[AlignerEngine]:
    return [BuiltinProgressive(), MafftEngine(), MuscleEngine(),
            ClustalOmegaEngine(), ProbConsEngine(), PrankEngine(), ClustalWEngine()]


def available_engines() -> List[AlignerEngine]:
    return [e for e in all_engines() if e.available()]


def get_engine(key: str) -> Optional[AlignerEngine]:
    for e in all_engines():
        if e.key == key:
            return e
    return None
