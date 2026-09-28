# Benchmark results

The raw output behind [docs/validation.md](../../docs/validation.md), so the
tables there can be checked without re-running anything. The commands that
produced them are given at the foot of that page.

Used by `validation.md` and the paper:

| file | what it is |
|---|---|
| `balibase_official.csv` | official `bali_score` SP and TC over the reference core blocks, seven aligners — the figures comparable with published BAliBASE tables. The PRANK rows are the 0.5.10 re-run with PRANK's own default gap costs (RV11 and RV12 only); the ClustalW rows were added in 0.5.10 with `--engines clustalw` |
| `balibase.csv` | the harness's own per-family output: all-column SP and TC, the per-column reliability AUC, and the masking controls (kept / matched-random / gap-fraction, with a z score). Its four PRANK rows are superseded by `prank_percol.csv`; the ClustalW rows (RV11 and RV12) were added in 0.5.10 |
| `prank_percol.csv` | the same per-family output for PRANK, from the 0.5.10 re-run |
| `simulation.csv` | the 360-dataset simulation sweep, five aligners (ClustalW added in 0.5.10): same metrics plus Robinson–Foulds distance of the neighbour-joining tree, full alignment versus masked. Made with the consistency score alone and the built-in engine at `--effort med`; the command is at the foot of `validation.md` |
| `residue_masking_nj.csv` | the 360-dataset simulation sweep for the built-in engine, with tree error (neighbour-joining) when the reliability score masks columns and when it masks residues: below the threshold, the same number as the column mask, and that many at random |
| `residue_masking_ml.csv` | the same comparison with IQ-TREE maximum-likelihood trees (JC+G4), divergence 0.25 and 0.5 — what `--tree ml` produces, with only the tree columns kept |
| `bali_trees_mafft.csv` | the paper's figure 2, from `benchmarks/bali_trees.py`: BAliBASE RV11 and RV12, one MAFFT alignment per family, IQ-TREE 2 trees (LG+G4) before and after each mask. Normalised Robinson–Foulds distance to the tree of the reference alignment (`_ref`) and to the tree of the unmasked alignment (`_full`), for reliability columns (`col`), the matched lowest-scoring residues (`resm`), residues below 0.5 (`resb`), trimAl gappyout (`gappy`), and random columns and residues (`rcol`, `rres`, mean of three draws). Blank where a mask left too little for a tree |
| `bali_trees_builtin.csv` | the same with CRAIC's built-in aligner |

Earlier runs, kept for the record:

| file | what it is |
|---|---|
| `balibase_coreblock.csv`, `balibase_allcolumn.csv` | the 0.5.0 BAliBASE run: four aligners, 181 families |
| `prank_sim.csv` | PRANK on the simulation sweep, made by 0.5.9 with the wrong gap extension (0.5; PRANK's DNA default is 0.75). Not used anywhere |

Environment: MAFFT 7.505, MUSCLE 5.1 (the Debian 5.1.0 package, whose binary
reports itself as 5.2), Clustal Omega 1.2.4, ProbCons 1.12, PRANK v.170427,
ClustalW 2.1, IQ-TREE 2.0.7, BAliBASE 3.0, CRAIC with the compiled Rust core.

The built-in engine's rows in `balibase.csv` (RV11 and RV12), `balibase_official.csv`,
`simulation.csv` and `residue_masking_nj.csv` are the 0.5.11 re-run with the exact
consistency transformation (see the note on `validation.md`); the other aligners'
rows are unchanged. `residue_masking_ml.csv` is from 0.5.10 and has not been re-run.

Two scope notes that matter for reading these numbers:

* The `rel_auc` column is the **consistency** score alone. The runs were made
  without `--perturbation`, which is roughly an order of magnitude more
  expensive; the combined score the interface shows by default adds the
  perturbation ensemble on top.
* Larger BAliBASE families exceed the compute budget — the consistency analysis
  is cubic in the number of sequences — and families that any aligner failed to
  complete are excluded from **every** aligner's mean, so the comparison stays
  like for like.
