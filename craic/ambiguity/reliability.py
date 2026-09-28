"""Per-column / per-residue reliability + live masking.

Two complementary, well-grounded signals:

* **consistency** (TCS / heads-style): for each alignment column, average the
  pair-HMM posterior P(res_i ~ res_j) over the residue pairs that column
  asserts are homologous. High = the unaligned-sequence evidence supports the
  column. Uses the CRAIC core directly.
* **perturbation** (guide-tree / gap-regime sensitivity, in the spirit of
  GUIDANCE but not equivalent to it): re-align the same sequences under an
  ensemble of perturbed guide trees across a few gap regimes, and measure what
  fraction of each residue's asserted homologies survive. See
  :func:`perturbation` for what the ensemble does and does not sample.

Both live in [0, 1]; ``combined`` averages them. A threshold turns either into a
column mask that the viewer previews live before you commit.

Two properties this module is careful about, because getting either wrong
produces a plausible-looking number that is wrong:

* Both scores are computed under parameters **estimated from the alignment
  being scored** (:func:`craic.progressive.estimate_params`), not under a fixed
  default. Scoring a 40%-identity family under a 90%-identity model makes every
  column look unreliable and biases masking toward destroying exactly the
  divergent data that masking decisions are about.
* A failed ensemble is an error, not a score of zero. If no replicate
  completes, :func:`perturbation` raises :class:`EnsembleFailure` rather than
  returning all-nan, which downstream would be indistinguishable from "every
  column is maximally unreliable" and would mask the whole alignment away.
"""

from __future__ import annotations

import warnings
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import accel, domain, progressive
from ..domain import Alignment, Alphabet

ResKey = Tuple[int, int]  # (sequence index, residue index)


class EnsembleFailure(RuntimeError):
    """Every replicate in the perturbation ensemble failed.

    Raised rather than returning an all-nan score, which a threshold would read
    as "mask everything".
    """


def _col_nanmean(mat: np.ndarray) -> np.ndarray:
    """Column-wise nanmean that returns nan (not a warning) for empty columns."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(mat, axis=0)


# --------------------------------------------------------------------------- #
# Consistency (uses the pair-HMM posterior core)
# --------------------------------------------------------------------------- #

#: Sequence pairs the consistency score samples at most (all of them below).
MAX_PAIRS = 300

def consistency(
    aln: Alignment,
    model: Optional[accel.EmissionModel] = None,
    max_pairs: int = MAX_PAIRS,
    seed: int = 0,
    delta: Optional[float] = None,
    epsilon: Optional[float] = None,
    progress=None,
    cancelled=None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (col_score[length], cell_score[n, length] with nan at gaps).

    ``model``, ``delta`` and ``epsilon`` default to the values estimated from
    ``aln`` itself rather than to library defaults; pass them explicitly only
    to score under a model of your own choosing. The pair posteriors are
    computed on several threads and summed in a fixed order, so the scores do
    not depend on the thread count. ``progress(done, total)`` is called once per
    sequence pair; ``cancelled()`` is polled throughout.
    """
    est_model, est_delta, est_epsilon = progressive.estimate_params(aln.rows, aln.alphabet)
    if model is None:
        model = est_model
    if delta is None:
        delta = est_delta
    if epsilon is None:
        epsilon = est_epsilon

    n = aln.n_seqs
    L = aln.length
    maps = [domain.residue_index(r) for r in aln.rows]

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    if len(pairs) > max_pairs:
        rng = np.random.default_rng(seed)
        pairs = [pairs[k] for k in rng.choice(len(pairs), max_pairs, replace=False)]

    # accumulate per-cell sums and counts
    cell_sum = np.zeros((n, L))
    cell_cnt = np.zeros((n, L))
    col_sum = np.zeros(L)
    col_cnt = np.zeros(L)

    def post(p):
        return accel.posterior_matrix(aln.rows[p[0]], aln.rows[p[1]], model, delta, epsilon)

    ahead = accel.in_flight(max(len(domain.ungap(r)) for r in aln.rows))
    for k, ((i, j), P) in enumerate(zip(pairs, accel.imap(post, pairs, cancelled, ahead)), 1):
        if progress is not None:
            progress(k, len(pairs))
        if P.size == 0:
            continue
        mi, mj = maps[i], maps[j]
        for c in range(L):
            ri, rj = mi[c], mj[c]
            if ri < 0 or rj < 0:
                continue
            p = float(P[ri, rj])
            col_sum[c] += p
            col_cnt[c] += 1
            cell_sum[i, c] += p
            cell_cnt[i, c] += 1
            cell_sum[j, c] += p
            cell_cnt[j, c] += 1

    col = np.divide(col_sum, col_cnt, out=np.full(L, np.nan), where=col_cnt > 0)
    cell = np.divide(cell_sum, cell_cnt, out=np.full((n, L), np.nan), where=cell_cnt > 0)
    return col, cell


# --------------------------------------------------------------------------- #
# Perturbation (guide-tree / gap-regime sensitivity)
# --------------------------------------------------------------------------- #

# Pair-HMM gap-open probabilities cycled across the replicates, spanning roughly
# an order of magnitude either side of the usual estimate so that the ensemble
# samples gap-model uncertainty as well as guide-tree uncertainty. These are
# tool defaults, not calibrated quantities: they set the scale of the score, so
# they are exposed as an argument and documented rather than buried.
_PERTURB_DELTAS = (0.01, 0.03, 0.08)


def perturbation(
    aln: Alignment,
    alphabet: Optional[Alphabet] = None,
    n_replicates: int = 16,
    seed: int = 0,
    deltas: Tuple[float, ...] = _PERTURB_DELTAS,
    stats: Optional[Dict[str, int]] = None,
    progress=None,
    cancelled=None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return (col_score[length], cell_score[n, length] with nan at gaps).

    Re-align the same sequences under ``n_replicates`` perturbed guide trees
    (cycling through the gap regimes in ``deltas``) and, for each residue, score
    the fraction of its reference-asserted homologies that survive across the
    ensemble. A residue with no asserted partners (a singleton column) has no
    residue-pair evidence and is left nan rather than scored a free 1.0.

    **What this is.** A sensitivity analysis of one engine: how stable are the
    homologies this alignment asserts, when the guide tree is perturbed and the
    gap model is varied? Low scores mark residues whose placement is an artefact
    of a particular tree or gap cost.

    **What this is not.** It is not GUIDANCE. GUIDANCE bootstraps alignment
    columns to build perturbed guide trees and re-aligns with the *same* aligner
    that produced the reference, over ~100 replicates. Here the replicates come
    from CRAIC's built-in engine under plain posterior decoding (no consistency
    pass), which is a weaker aligner than the one that may have produced the
    reference and a *different* aligner if the reference came from MAFFT or
    PRANK. In that case the score conflates genuine alignment uncertainty with
    systematic between-method difference. Read it as a stability probe, and read
    the disagreement map for between-method evidence.

    Raises :class:`EnsembleFailure` if no replicate completes. If some fail, a
    warning is issued and the surviving replicates are used; ``analyse`` records
    the count on the report.

    ``progress(done, total)`` counts sequence pairs aligned across the whole
    ensemble (``n_replicates`` x n(n-1)/2); ``cancelled()`` is polled
    throughout, and a cancel raises :class:`accel.Cancelled` rather than
    counting as a failed replicate.
    """
    alphabet = alphabet or aln.alphabet
    seqs = [domain.ungap(r) for r in aln.rows]
    records = [(str(i), s) for i, s in enumerate(seqs)]

    model, _est_delta, _est_epsilon = progressive.estimate_params(aln.rows, alphabet)

    wanted = max(1, n_replicates)
    per_rep = len(seqs) * (len(seqs) - 1) // 2
    total = wanted * per_rep
    done = [0]

    def tick():
        done[0] += 1
        if progress is not None:
            progress(done[0], total)

    alts: List[Alignment] = []
    first_error: Optional[BaseException] = None
    for i in range(wanted):
        done[0] = i * per_rep                # a failed replicate still moves the bar on
        try:
            alts.append(progressive.align(
                records, alphabet, model=model,
                delta=deltas[i % len(deltas)], epsilon=0.5,
                estimate=False, consistency_iters=0, guide_seed=seed + i + 1,
                cancelled=cancelled, tick=tick))
        except accel.Cancelled:
            raise
        except Exception as exc:          # noqa: BLE001 - recorded, then re-raised if total
            if first_error is None:
                first_error = exc
    n_failed = wanted - len(alts)
    if not alts:
        raise EnsembleFailure(
            f"all {wanted} perturbation replicates failed; "
            f"first error: {first_error!r}"
        ) from first_error
    if n_failed:
        warnings.warn(
            f"{n_failed} of {wanted} perturbation replicates failed; "
            "scores are based on the remainder",
            RuntimeWarning,
            stacklevel=2,
        )
    if stats is not None:
        stats["n_ok"], stats["n_failed"] = len(alts), n_failed

    n = aln.n_seqs
    L = aln.length

    # reference membership uses *positional* sequence index; alts use str(index)
    def mem_by_index(a: Alignment) -> Dict[ResKey, int]:
        mem: Dict[ResKey, int] = {}
        for row_id, row in zip(a.ids, a.rows):
            si = int(row_id)
            r = 0
            for c, ch in enumerate(row):
                if ch != "-":
                    mem[(si, r)] = c
                    r += 1
        return mem

    ref_mem = {(si, r): c for (si, r), c in domain.membership(aln).items()}
    alt_mems = [mem_by_index(a) for a in alts]

    # partners asserted by the reference, per residue
    partners: Dict[ResKey, List[ResKey]] = defaultdict(list)
    maps = [domain.residue_index(r) for r in aln.rows]
    for c in range(L):
        present = [(si, maps[si][c]) for si in range(n) if maps[si][c] >= 0]
        for a_idx in range(len(present)):
            for b_idx in range(len(present)):
                if a_idx != b_idx:
                    partners[present[a_idx]].append(present[b_idx])

    res_score: Dict[ResKey, float] = {}
    for x, plist in partners.items():
        if not alt_mems or not plist:
            continue  # no ensemble, or a singleton column: no residue-pair evidence -> nan
        tot = ok = 0
        for y in plist:
            for mem in alt_mems:
                tot += 1
                cx, cy = mem.get(x), mem.get(y)
                if cx is not None and cx == cy:
                    ok += 1
        if tot:
            res_score[x] = ok / tot

    cell = np.full((n, L), np.nan)
    for (si, r), score in res_score.items():
        c = ref_mem.get((si, r))
        if c is not None:
            cell[si, c] = score
    col = _col_nanmean(cell)
    return col, cell


# --------------------------------------------------------------------------- #
# Combined report + masking
# --------------------------------------------------------------------------- #

def keep_mask(
    scores: np.ndarray,
    threshold: float,
    unscored: str = "drop",
) -> np.ndarray:
    """Columns to keep: ``scores >= threshold``.

    This is the *one* definition of a reliability mask in CRAIC. It previously
    existed twice with opposite handling of unscoreable (nan) columns — the GUI
    kept them, the benchmark dropped them — so the mask a user exported was not
    the mask the benchmark evaluated.

    A nan column is one with no residue-pair evidence at all (a singleton
    column, or one whose ensemble evidence is missing). ``unscored="drop"`` is
    the default and the conservative choice for downstream phylogenetics: a
    column that could not be assessed is not evidence. ``unscored="keep"``
    retains them.

    Raises ``ValueError`` if *every* column is unscoreable, which means the
    reliability analysis failed rather than that the alignment is worthless.
    """
    if unscored not in ("drop", "keep"):
        raise ValueError(f"unscored must be 'drop' or 'keep', not {unscored!r}")
    scores = np.asarray(scores, dtype=float)
    if scores.size and bool(np.all(np.isnan(scores))):
        raise ValueError(
            "every column is unscored, so no threshold is meaningful; "
            "the reliability analysis did not produce usable scores"
        )
    fill = 0.0 if unscored == "drop" else 1.0
    return np.nan_to_num(scores, nan=fill) >= threshold


@dataclass
class Reliability:
    length: int
    col_consistency: np.ndarray
    col_perturbation: np.ndarray
    col_combined: np.ndarray
    cell_consistency: np.ndarray
    cell_perturbation: np.ndarray
    #: replicates that completed / failed in the perturbation ensemble
    n_replicates_ok: int = 0
    n_replicates_failed: int = 0

    @property
    def cell_combined(self) -> np.ndarray:
        """Per-residue confidence (n x length), nan at gaps — for colouring the
        alignment by reliability."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmean(
                np.stack([self.cell_consistency, self.cell_perturbation]), axis=0
            )

    def column_scores(self, which: str = "combined") -> np.ndarray:
        return {
            "combined": self.col_combined,
            "consistency": self.col_consistency,
            "perturbation": self.col_perturbation,
        }[which]

    def keep_mask(self, threshold: float, which: str = "combined",
                  unscored: str = "drop") -> np.ndarray:
        return keep_mask(self.column_scores(which), threshold, unscored)

    def cell_scores(self, which: str = "combined") -> np.ndarray:
        """Per-residue scores (n x length, nan at gaps) — the residue analogue
        of :meth:`column_scores`."""
        return {
            "combined": self.cell_combined,
            "consistency": self.cell_consistency,
            "perturbation": self.cell_perturbation,
        }[which]


def analyse(
    aln: Alignment,
    do_perturbation: bool = True,
    model: Optional[accel.EmissionModel] = None,
    n_replicates: int = 16,
    progress=None,
    cancelled=None,
) -> Reliability:
    """Both scores and their combination.

    ``progress(done, total)`` counts sequence pairs over both stages together
    (see :func:`work`); ``cancelled()`` is polled throughout and a cancel raises
    :class:`accel.Cancelled`.
    """
    n_c, n_p = work(aln, do_perturbation, n_replicates)

    def stage(offset):
        if progress is None:
            return None
        return lambda done, _total: progress(offset + done, n_c + n_p)

    cons = consistency(aln, model=model, progress=stage(0), cancelled=cancelled)
    pert, stats = None, {}
    if do_perturbation:
        pert = perturbation(aln, n_replicates=n_replicates, stats=stats,
                            progress=stage(n_c), cancelled=cancelled)
    return combine(aln, cons, pert, stats)


def work(aln: Alignment, do_perturbation: bool = True,
         n_replicates: int = 16, max_pairs: Optional[int] = None) -> Tuple[int, int]:
    """Sequence pairs each stage aligns: ``(consistency, perturbation)``.

    Every pair costs about the same (one pair-HMM over the two sequences), so
    the time the consistency stage took predicts the perturbation stage well;
    the perturbation stage is usually by far the larger.
    """
    max_pairs = MAX_PAIRS if max_pairs is None else max_pairs
    n = aln.n_seqs
    pairs = n * (n - 1) // 2
    return min(pairs, max_pairs), (max(1, n_replicates) * pairs if do_perturbation else 0)


def combine(aln: Alignment, cons, pert=None, stats: Optional[Dict[str, int]] = None
            ) -> Reliability:
    """Assemble a report from :func:`consistency`'s and, optionally,
    :func:`perturbation`'s ``(col, cell)`` results. Without the perturbation
    scores, ``combined`` is the consistency score alone and
    ``n_replicates_ok`` is 0."""
    col_c, cell_c = cons
    stats = stats or {}
    if pert is None:
        col_p = np.full(aln.length, np.nan)
        cell_p = np.full((aln.n_seqs, aln.length), np.nan)
    else:
        col_p, cell_p = pert
    combined = _col_nanmean(np.vstack([col_c, col_p]))
    return Reliability(
        aln.length, col_c, col_p, combined, cell_c, cell_p,
        n_replicates_ok=stats.get("n_ok", 0),
        n_replicates_failed=stats.get("n_failed", 0),
    )


def apply_mask(aln: Alignment, keep_mask: np.ndarray) -> Alignment:
    keep = [c for c in range(aln.length) if keep_mask[c]]
    rows = ["".join(r[c] for c in keep) for r in aln.rows]
    meta = dict(aln.meta)
    meta.pop(RESIDUE_MASK_KEY, None)        # residue indices no longer hold once columns go
    return Alignment(list(aln.ids), rows, aln.alphabet, coding=None, meta=meta)


# --------------------------------------------------------------------------- #
# Residue masking
# --------------------------------------------------------------------------- #
#
# Masking a column throws away every residue in it, including the ones that
# are aligned correctly, and on a tree that costs more signal than the error it
# removes (Tan et al. 2015). A residue mask removes only the residues judged
# unreliable — by hand, or below a per-residue score — and keeps the column.
#
# A mask is kept in ``Alignment.meta`` as ``{sequence id: [residue indices]}``,
# indices into the ungapped sequence. Keyed by residue rather than by column,
# it follows its residues through any edit or realignment, since neither
# changes a sequence's residues; and it is JSON-safe, so sessions carry it. It
# changes nothing about the alignment until it is applied on export, where each
# masked residue is written as missing data.

RESIDUE_MASK_KEY = "masked_residues"

ResidueMask = Dict[str, List[int]]


def residue_mask(aln: Alignment) -> ResidueMask:
    """The alignment's residue mask, validated against its sequences.

    Entries for sequences that are not present, or indices beyond a sequence's
    length, are dropped: a mask carried onto an alignment that does not hold
    those residues must not mask something else.
    """
    raw = aln.meta.get(RESIDUE_MASK_KEY) or {}
    lengths = {i: sum(ch != "-" for ch in r) for i, r in zip(aln.ids, aln.rows)}
    out: ResidueMask = {}
    for sid, idx in raw.items():
        if sid in lengths:
            keep = sorted({int(k) for k in idx if 0 <= int(k) < lengths[sid]})
            if keep:
                out[sid] = keep
    return out


def with_residue_mask(aln: Alignment, mask: ResidueMask) -> Alignment:
    """A copy of ``aln`` carrying ``mask`` (an empty mask removes it)."""
    meta = dict(aln.meta)
    clean = {sid: sorted(set(idx)) for sid, idx in mask.items() if idx}
    if clean:
        meta[RESIDUE_MASK_KEY] = clean
    else:
        meta.pop(RESIDUE_MASK_KEY, None)
    return Alignment(list(aln.ids), list(aln.rows), aln.alphabet,
                     coding=aln.coding, meta=meta)


def residues_in(aln: Alignment, rows: Sequence[int], col_start: int,
                col_stop: int) -> ResidueMask:
    """The residues of ``rows`` lying in columns [col_start, col_stop)."""
    out: ResidueMask = {}
    for si in rows:
        idx = domain.residue_index(aln.rows[si])[col_start:col_stop]
        picked = [r for r in idx if r >= 0]
        if picked:
            out[aln.ids[si]] = picked
    return out


def residues_below(aln: Alignment, cell_scores: np.ndarray,
                   threshold: float) -> ResidueMask:
    """Residues whose per-residue score is below ``threshold``.

    An unscored residue (nan: nothing else in its column to be homologous to) is
    not masked. Unlike an unscored column, it asserts no homology, so there is
    nothing in it to be wrong.
    """
    cells = np.asarray(cell_scores, dtype=float)
    if cells.shape != (aln.n_seqs, aln.length):
        raise ValueError(f"cell scores are shaped {cells.shape}, "
                         f"but the alignment is {(aln.n_seqs, aln.length)}")
    out: ResidueMask = {}
    for si, row in enumerate(aln.rows):
        idx = domain.residue_index(row)
        low = [idx[c] for c in range(aln.length)
               if idx[c] >= 0 and cells[si, c] == cells[si, c] and cells[si, c] < threshold]
        if low:
            out[aln.ids[si]] = low
    return out


def merge_masks(a: ResidueMask, b: ResidueMask) -> ResidueMask:
    return {sid: sorted(set(a.get(sid, [])) | set(b.get(sid, []))) for sid in {*a, *b}}


def subtract_masks(a: ResidueMask, b: ResidueMask) -> ResidueMask:
    out = {sid: sorted(set(idx) - set(b.get(sid, []))) for sid, idx in a.items()}
    return {sid: idx for sid, idx in out.items() if idx}


def mask_size(mask: ResidueMask) -> int:
    return sum(len(v) for v in mask.values())


def masked_cells(aln: Alignment, mask: ResidueMask) -> set:
    """``(row, column)`` of every masked residue, for drawing."""
    cells = set()
    for si, (sid, row) in enumerate(zip(aln.ids, aln.rows)):
        wanted = set(mask.get(sid, ()))
        if wanted:
            for c, r in enumerate(domain.residue_index(row)):
                if r in wanted:
                    cells.add((si, c))
    return cells


def missing_symbol(alphabet: Alphabet) -> str:
    """How a masked residue is written: missing data to every tree program."""
    return "N" if alphabet.is_nucleotide else "X"


def apply_residue_mask(aln: Alignment, mask: Optional[ResidueMask] = None) -> Alignment:
    """``aln`` with each masked residue replaced by N (nucleotides) or X (protein).

    Columns are untouched. ``mask`` defaults to the alignment's own.
    """
    mask = residue_mask(aln) if mask is None else mask
    sym = missing_symbol(aln.alphabet)
    rows = []
    for sid, row in zip(aln.ids, aln.rows):
        wanted = set(mask.get(sid, ()))
        if not wanted:
            rows.append(row)
            continue
        out, r = [], 0
        for ch in row:
            if ch == "-":
                out.append(ch)
            else:
                out.append(sym if r in wanted else ch)
                r += 1
        rows.append("".join(out))
    meta = dict(aln.meta)
    meta.pop(RESIDUE_MASK_KEY, None)
    return Alignment(list(aln.ids), rows, aln.alphabet, coding=aln.coding, meta=meta)
