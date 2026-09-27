# CRAIC benchmarks

A reproducible validation suite for the two things the paper claims:

1. **Alignment accuracy** — how the built-in ProbCons-style aligner scores
   (sum-of-pairs and total-column) against a known-true alignment, alongside any
   external aligner on your `PATH`.
2. **Reliability & masking** — whether CRAIC's per-column reliability score
   actually flags the mis-aligned columns (ROC AUC), whether masking the flagged
   columns raises mean column accuracy, and whether it improves the recovered
   (neighbour-joining) tree.

Ground truth comes from two sources, scored with the same metrics:

- a **simulation** sweep (sequences evolved down random trees with substitutions
  and indels, so every residue's homology and the true topology are known), and
- **real reference alignments** you supply (BAliBASE, HOMSTRAD, PREFAB, ...).

Because reliability is a *per-column* claim, you need per-column ground truth —
which simulation gives exactly and structural reference sets give for their core
blocks. That is why both arms exist.

## Requirements

- CRAIC installed (`pip install -e .` from the repo root), **with the compiled
  Rust core** for speed: `maturin develop --release`. It runs on the NumPy
  fallback too, just far slower — the sweep below is minutes with Rust, much
  longer without.
- Biopython (already a CRAIC dependency) for reading reference formats and NJ.
- Optional: **MAFFT / MUSCLE / Clustal Omega / PRANK** on your `PATH`. The
  harness auto-detects them (via CRAIC's engine registry) and includes each in
  the comparison — no configuration needed.

## Quick smoke test

```bash
python benchmarks/run_benchmark.py --quick
```

Runs a tiny simulation preset and prints per-run and summary metrics.

## Simulation sweep

```bash
python benchmarks/run_benchmark.py --sim \
    --taxa 8 16 --length 400 --reps 20 \
    --divergence 0.15 0.25 0.40 --indel 1.0 2.0 3.0 \
    --out sim_results.csv
```

- `--divergence` are branch-length upper bounds; larger = more substitutions
  (harder). `~0.15` is high identity, `~0.40` approaches saturation.
- `--indel` is the indel load; larger = more gaps to place.
- `--reps` replicates per condition (use >=20 for stable means).
- `--perturbation` uses the combined consistency+perturbation reliability score
  instead of consistency alone (slower; the guide-tree bootstrap re-aligns many
  times per column set).
- `--tree ml` measures tree error with maximum-likelihood trees (IQ-TREE 2 on
  your `PATH`, JC+G4 — the simulating model) instead of neighbour-joining. Much
  slower; the paper runs it at `--divergence 0.25 0.5`.
- `--engines builtin mafft` restricts any run to the named engines — for
  example to re-run one aligner after a fix without re-running the others.

Each `(condition x replicate x aligner)` becomes one CSV row.

## Real reference benchmarks

Point `--reference` at any set of reference alignment files (any format
Biopython reads: FASTA, Clustal, **MSF**, Stockholm, PHYLIP, NEXUS). The
sequences are un-gapped, re-aligned by each available engine, and scored against
the reference.

```bash
# BAliBASE (protein). Recent releases (R9/R10) embed the reference alignment in
# each .xml (there is no .msf); the harness reads it. Older releases with .msf
# also work -- just point the glob at the .msf files instead.
python benchmarks/run_benchmark.py \
    --reference 'RV100/*.xml' --protein \
    --out balibase_results.csv
```

### Official BAliBASE scoring (bali-score)

The `sp`/`tc` in the CSV above are computed against the *full* reference
alignment. BAliBASE's official accuracy scores are computed over the annotated
**core blocks** only, by its own `bali_score` program. To get those numbers,
have the harness emit each aligner's alignment and score them with the tool.

Install the scorer once (Rust reimplementation; you already have `cargo`):

```bash
cargo install --git https://github.com/robinhundt/bali-score
```

Emit the alignments while you run the benchmark:

```bash
python benchmarks/run_benchmark.py --reference 'RV100/*.xml' --protein \
    --write-alignments aligned_out --out balibase.csv
```

That writes one FASTA per (reference, aligner) into `aligned_out/` (with the
reference sequence names, so `bali-score` matches them) plus a tab-delimited
`bali_manifest.tsv`. Score everything against its reference and save the output:

```bash
python -c "import csv,subprocess; [subprocess.run(['bali-score','-t',t,'-r',r]) for t,r,a in list(csv.reader(open('aligned_out/bali_manifest.tsv'),delimiter=chr(9)))[1:]]" | tee bali_scores.txt
```

(The one-liner is used instead of a shell loop because it handles reference
paths that contain spaces.) `bali-score -t <test.fasta> -r <reference.xml>`
prints the core-block SP and TC for each alignment.

Where to get reference sets:

- **BAliBASE** (structure-anchored protein reference alignments) —
  https://www.lbgi.fr/balibase/ . Use the `.msf` reference alignments; add
  `--protein`.
- **HOMSTRAD** — https://mizuguchilab.org/homstrad/ .
- **PREFAB**, **OXBench**, **SABmark** are also usable — anything that ships a
  reference MSA per family.

For nucleotide reference data, drop `--protein` (the default alphabet is DNA).

The `--write-alignments DIR` flag (reference mode) additionally writes each
aligner's alignment to `DIR/` for official `bali_score` scoring (see above).

## Comparing masking methods (`--masking`)

CRAIC's model-based reliability can be compared head-to-head against the
model-free trimmers it now also implements (trimAl-, Gblocks-, MSA_trimmer-style)
on data with a known-true alignment. For each family it aligns the sequences and
asks, for every method: how well does it tell the correct columns from the
mis-aligned ones?

```bash
python benchmarks/run_benchmark.py --reference 'RV100/*.xml' --protein \
    --masking --max-mem-gb 10 --out masking.csv
# or on simulated data:
python benchmarks/run_benchmark.py --sim --masking --reps 20 --out masking_sim.csv
```

Per method it reports **AUC** (how well a per-column score ranks correct columns
above wrong ones; score-based methods only), fraction of columns **kept**,
**acc_kept** (accuracy of the retained columns), and **precision/recall** of the
kept-vs-correct decision. On a BAliBASE family this looks like:

| method | AUC | kept | acc_kept | precision |
|---|---|---|---|---|
| CRAIC reliability | **0.93** | 0.88 | 0.94 | 0.90 |
| gap (MSA_trimmer) | 0.67 | 0.96 | 0.90 | 0.85 |
| similarity (trimAl) | 0.75 | 0.41 | 0.96 | 0.95 |
| trimAl gappyout | – | 0.91 | 0.92 | 0.87 |
| trimAl strict | – | 0.52 | 0.97 | 0.95 |
| Gblocks | – | 0.21 | 1.00 | 1.00 |

The model-based reliability distinguishes correct from incorrect columns best
(highest AUC); the rule-based trimmers (Gblocks, trimAl strict) are very
conservative — near-perfect precision but they discard most columns. These are
faithful reimplementations of the published *algorithms* (`craic.ambiguity.trimming`);
for authoritative numbers run the original trimAl / Gblocks.

## Trees on BAliBASE (`bali_trees.py`)

The paper's figure 2 asks whether masking changes trees on real alignments.
BAliBASE has no true tree, so `bali_trees.py` compares each tree with the tree
from the family's reference alignment. It aligns every RV11 and RV12 family
once, applies each mask to that alignment (nothing is re-aligned) and builds
IQ-TREE 2 trees (LG+G4, `iqtree2` on your `PATH`):

```bash
python benchmarks/bali_trees.py path/to/bb3_release mafft --out bali_trees_mafft.csv
python benchmarks/bali_trees.py path/to/bb3_release builtin --out bali_trees_builtin.csv
```

It resumes if stopped. The columns are described in
[results/README.md](results/README.md).

## Output

A CSV with one row per run and a printed summary of means per aligner. Columns:

| column | meaning |
|---|---|
| `sp`, `tc` | sum-of-pairs / total-column accuracy vs the true (or reference) alignment |
| `score_scope` | `core` if the reference marked core blocks (by case) and SP/TC were restricted to them, else `all-columns` |
| `rel_auc` | ROC AUC of per-column reliability predicting which columns are actually correct (0.5 = chance) |
| `pos_rate` | fraction of scorable columns that are actually correct — the AUC's positive rate |
| `acc_all`, `acc_retained` | mean per-column accuracy over all columns vs. those kept at the reliability threshold |
| `acc_retained_random` | the same, for a mask removing the same *number* of columns at random (mean of 200 draws) |
| `acc_retained_gap` | the same, for a gap-fraction mask removing the same number of columns |
| `mask_z` | (`acc_retained` − `acc_retained_random`) in units of the random-mask standard deviation |
| `n_masked`, `n_cols` | columns dropped by the mask, and total |
| `rf_full`, `rf_masked` | normalised Robinson-Foulds distance of the tree (NJ, or ML with `--tree ml`) to the true topology, before / after masking columns (simulation only) |
| `rf_resid` | the same, with every residue scoring below the threshold masked as missing data and all columns kept |
| `rf_resid_matched` | the same, masking the lowest-scoring residues — exactly as many as the column mask removed |
| `rf_resid_random` | the same, masking that many residues at random (mean of 20 draws for NJ, 5 for ML) |
| `n_res`, `n_res_colmask`, `n_res_resid` | residues in the alignment, removed by the column mask, and masked by the threshold residue mask |

**Reliability works** if `rel_auc > 0.5`.

**Masking helps** if `acc_retained > acc_retained_random` (equivalently `mask_z`
is comfortably positive) and, for the claim to be worth the compute, if
`acc_retained > acc_retained_gap`. It is **not** enough that
`acc_retained > acc_all`: that inequality holds for any score correlated with
column difficulty — including raw gap fraction, including column entropy — simply
because masking removes the hard columns. v0.1 of this script reported that
comparison as the result, which was circular; the random and gap-matched controls
exist to replace it. `rf_masked <= rf_full` remains the downstream test.

**A note on SP and TC.** Published BAliBASE SP/TC figures are computed over the
reference's **core blocks** by the official `bali_score`. The figures this script
produces are restricted to core blocks only when the reference file marks them by
case (`score_scope = core`); otherwise they are full-column scores and are *not*
comparable with any published table. For citable BAliBASE numbers use
`--write-alignments` and the `bali-score` route described above.

Expect the reliability signal to be strongest at intermediate divergence — where
there is real, recoverable error to find. On near-perfect alignments there is
little to flag; near saturation, almost everything is uncertain.

**A note on the simulator.** Its default model — JC69, uniform rates, geometric
indel lengths — is close to what CRAIC's pair-HMM assumes, so CRAIC's inference
model is nearly well-specified on this data while MAFFT and MUSCLE carry
empirical parameters tuned for real sequences. Use `--rate-alpha` (gamma
among-site rate variation) and `--indel-zipf` (heavy-tailed indel lengths) to
break that correspondence, and report both runs.

## What we saw in a small demo run (NumPy fallback, few replicates)

Indicative only — run the full sweep on your machine for paper-grade numbers:

- Aligner: mean SP ~ 0.93, TC ~ 0.77 at intermediate divergence.
- Reliability AUC ~ 0.66-0.85 (predicting mis-aligned columns).
- Masking raised mean column accuracy above a matched random mask in every
  replicate (that is the comparison that counts; see above).
- On a hard, scrambled replicate (SP ~ 0.60) masking the flagged columns took
  the NJ tree from `RF = 0.2` to `RF = 0.0` — i.e. it recovered the true
  topology the raw alignment had lost.

## Notes

- Runtime is dominated by CRAIC's all-pairs posteriors; the Rust core makes the
  sweep practical. Start with `--quick`, then scale `--reps`, `--taxa`,
  `--length`.
- Results are seeded (`seed = replicate index`), so runs are reproducible.
- The simulator is deliberately simple (JC69 + geometric indels). It is a clean
  test of the *reliability* claim, not a model of real protein evolution — which
  is what the reference-benchmark arm is for.
