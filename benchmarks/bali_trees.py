#!/usr/bin/env python3
"""How much do column and residue masks change trees on real alignments?

BAliBASE 3 RV11 + RV12 (the families of the paper's Table 1; figure 2). There is
no true tree for these families, so each estimated tree is compared with two
others:

* the tree from the BAliBASE reference (structural) alignment of the same
  sequences, i.e. the tree the correct alignment supports; and
* the tree from the same estimated alignment unmasked, i.e. how far the mask
  moves the tree.

Each family is aligned once (MAFFT or CRAIC's built-in aligner) and the masks
are applied to that alignment; nothing is re-aligned. Masks, all from CRAIC's
consistency score at the default threshold of 0.5, as in the benchmark:
columns below 0.5 removed; the lowest-scoring residues, as many as the column
mask removed, written as missing data; every residue below 0.5; and, as
controls, trimAl gappyout (a rule-based column filter) and random columns /
random residues matched in number (mean of three draws). Trees: IQ-TREE 2 on
the PATH, LG+G4, --fast. Distances are normalised Robinson-Foulds.

    python benchmarks/bali_trees.py BALIBASE_DIR mafft   --out bali_trees_mafft.csv
    python benchmarks/bali_trees.py BALIBASE_DIR builtin --out bali_trees_builtin.csv

BALIBASE_DIR is the unpacked BAliBASE 3 release (holding RV11/ and RV12/).
Alignments are cached in --aligned (default aligned_trees/) and families
already in --out are skipped, so an interrupted run resumes.
"""
import argparse
import csv
import glob
import os
import subprocess
import tempfile

import numpy as np

import metrics as M
import run_benchmark as RB
from craic import engines, io
from craic.ambiguity import reliability as rel_mod
from craic.ambiguity import trimming
from craic.domain import Alignment, Alphabet

MAX_MEM_GB = 3.0        # families above this posterior-memory proxy are skipped, as in Table 1
THRESHOLD = 0.5
N_RANDOM = 3
MIN_COLS = 10           # fewer columns than this left: no tree


def ml_bips(names, rows):
    """Bipartitions of the IQ-TREE LG+G4 tree, or None if no tree can be built."""
    from Bio import Phylo
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, "a.fasta")
        with open(fa, "w") as fh:
            for i, r in enumerate(rows):
                fh.write(f">s{i}\n{r}\n")                 # safe names; mapped back below
        try:
            subprocess.run(["iqtree2", "-s", fa, "-st", "AA", "-m", "LG+G4", "--fast",
                            "-T", "1", "--seed", "1", "-quiet", "-redo"], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=900)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            return None
        tree = Phylo.read(fa + ".treefile", "newick")
    leaves = frozenset(names)
    ref = min(leaves)
    bips = set()
    for clade in tree.get_nonterminals():
        side = frozenset(names[int(t.name[1:])] for t in clade.get_terminals())
        if 2 <= len(side) <= len(leaves) - 2:
            bips.add(side if ref not in side else frozenset(leaves - side))
    return bips


def rows_in(aln, names):
    order = {n: i for i, n in enumerate(aln.ids)}
    return [aln.rows[order[n]] for n in names]


def align(eng, base, seqs, aln_dir):
    path = os.path.join(aln_dir, f"{base}__{eng.key}.fasta")
    if os.path.exists(path):
        names, rows = zip(*io.read_records(path, "fasta"))
        return Alignment(list(names), list(rows), Alphabet.PROTEIN)
    aln = eng.align(seqs, Alphabet.PROTEIN)
    io.write_fasta(path, list(aln.ids), list(aln.rows), wrap=0)
    return aln


def one(path, eng, aln_dir):
    base = os.path.splitext(os.path.basename(path))[0]
    names, true_rows, _core = RB._read_balibase_xml(path)
    seqs = [(n, r.replace("-", "")) for n, r in zip(names, true_rows)]
    maxlen = max(len(s) for _, s in seqs)
    if (len(seqs) * (len(seqs) - 1) // 2) * maxlen ** 2 * 8 / 1e9 > MAX_MEM_GB:
        return None                                         # the benchmark skipped it too
    aln = align(eng, base, seqs, aln_dir)
    ref_by_name = dict(zip(names, true_rows))
    names = list(aln.ids)
    ref_rows = [ref_by_name[n] for n in names]

    t_ref = ml_bips(names, ref_rows)
    t_full = ml_bips(names, rows_in(aln, names))
    rep = rel_mod.analyse(aln, do_perturbation=False)
    cells = rep.cell_consistency
    n_res = sum(ch != "-" for r in aln.rows for ch in r)

    def rf_pair(bips):
        if bips is None or t_ref is None or t_full is None:
            return "", ""
        return round(M.rf_distance(t_ref, bips), 4), round(M.rf_distance(t_full, bips), 4)

    def col_tree(keep):
        if keep.sum() < MIN_COLS:
            return None
        return ml_bips(names, rows_in(rel_mod.apply_mask(aln, keep), names))

    def res_tree(mask):
        return ml_bips(names, rows_in(rel_mod.apply_residue_mask(aln, mask), names))

    def removed(keep):
        return sum(ch != "-" for c in range(aln.length) if not keep[c] for ch in aln.column(c))

    out = dict(family=base, set=os.path.basename(os.path.dirname(path)),
               aligner=eng.key, n_seqs=aln.n_seqs, n_cols=aln.length, n_res=n_res,
               rf_full_ref=round(M.rf_distance(t_ref, t_full), 4) if t_ref and t_full else "")

    keep = rep.keep_mask(THRESHOLD, which="consistency")
    n_col_res = removed(keep)
    out["cols_removed"] = int((~keep).sum())
    out["res_removed_colmask"] = n_col_res
    out["rf_col_ref"], out["rf_col_full"] = rf_pair(col_tree(keep))

    out["rf_resm_ref"], out["rf_resm_full"] = rf_pair(res_tree(M.lowest_residues(aln, cells, n_col_res)))
    below = rel_mod.residues_below(aln, cells, THRESHOLD)
    out["res_below"] = rel_mod.mask_size(below)
    out["rf_resb_ref"], out["rf_resb_full"] = rf_pair(res_tree(below))

    gk = trimming.gappyout_mask(aln)
    out["gappy_cols_removed"] = int((~gk).sum())
    out["rf_gappy_ref"], out["rf_gappy_full"] = rf_pair(col_tree(gk))

    rng = np.random.default_rng(0)
    rc_ref, rc_full, rr_ref, rr_full = [], [], [], []
    for _ in range(N_RANDOM):
        rk = np.ones(aln.length, bool)
        rk[rng.choice(aln.length, int((~keep).sum()), replace=False)] = False
        a, b = rf_pair(col_tree(rk))
        if a != "":
            rc_ref.append(a); rc_full.append(b)
        a, b = rf_pair(res_tree(M.random_residues(aln, n_col_res, rng)))
        if a != "":
            rr_ref.append(a); rr_full.append(b)
    mean = lambda v: round(float(np.mean(v)), 4) if v else ""  # noqa: E731
    out.update(rf_rcol_ref=mean(rc_ref), rf_rcol_full=mean(rc_full),
               rf_rres_ref=mean(rr_ref), rf_rres_full=mean(rr_full))
    return out


FIELDS = ["family", "set", "aligner", "n_seqs", "n_cols", "n_res", "rf_full_ref",
          "cols_removed", "res_removed_colmask", "rf_col_ref", "rf_col_full",
          "rf_resm_ref", "rf_resm_full", "res_below", "rf_resb_ref", "rf_resb_full",
          "gappy_cols_removed", "rf_gappy_ref", "rf_gappy_full",
          "rf_rcol_ref", "rf_rcol_full", "rf_rres_ref", "rf_rres_full"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("balibase", help="BAliBASE 3 release directory (holds RV11/ and RV12/)")
    ap.add_argument("aligner", choices=["mafft", "builtin"])
    ap.add_argument("--out", required=True, help="output CSV (resumed if it exists)")
    ap.add_argument("--aligned", default="aligned_trees", help="alignment cache directory")
    args = ap.parse_args()

    engines.BuiltinProgressive.benchmark_effort = "max"
    engines.BuiltinProgressive.consistency_mem_gb = MAX_MEM_GB
    eng = {"mafft": engines.MafftEngine, "builtin": engines.BuiltinProgressive}[args.aligner]()
    os.makedirs(args.aligned, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.balibase, "RV11", "*.xml"))
                   + glob.glob(os.path.join(args.balibase, "RV12", "*.xml")))
    done = set()
    if os.path.exists(args.out):
        done = {r["family"] for r in csv.DictReader(open(args.out))}
    new = not os.path.exists(args.out)
    with open(args.out, "a", newline="", buffering=1) as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        for f in files:
            base = os.path.splitext(os.path.basename(f))[0]
            if base in done:
                continue
            try:
                r = one(f, eng, args.aligned)
            except Exception as exc:
                print(f"! {base}: {exc}", flush=True)
                continue
            if r is None:
                print(f"~ skip {base}", flush=True)
                continue
            w.writerow(r)
            print(base, r["n_seqs"], r["rf_full_ref"], r["rf_col_ref"], r["rf_resm_ref"],
                  r["rf_gappy_ref"], r["rf_rcol_ref"], flush=True)


if __name__ == "__main__":
    main()
