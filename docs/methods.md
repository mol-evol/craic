# Methods in full

The complete specification of what CRAIC computes, for anyone who needs to
reproduce it, cite it, or argue with it. [Concepts & methods](concepts.md) is the
short version; this is the one with the parameters in it.

Every quantity below derives from a single object — the matrix of posterior
match probabilities for a pair of sequences — so the page is organised as: build
that matrix, then read it six different ways.

---

## 1. The pair-HMM

Two ungapped sequences **a** (length *m*) and **b** (length *n*) are related by a
three-state hidden Markov model: **M** (the two residues are homologous), **X**
(a residue of **a** is aligned to a gap) and **Y** (a residue of **b** is aligned
to a gap). Gaps are affine, parameterised by two numbers:

| | meaning |
|---|---|
| δ (`delta`) | probability of leaving **M** for **X** or **Y** — gap opening |
| ε (`epsilon`) | probability of staying in **X** or **Y** — gap extension |

so the transitions are

```
M → M   1 − 2δ        X → X   ε         Y → Y   ε
M → X   δ             X → M   1 − ε     Y → M   1 − ε
M → Y   δ
```

This is the ProbCons formulation (Do *et al.*, 2005). Note that δ ≥ 0.5 makes
log(1 − 2δ) undefined; both backends validate the range at the boundary rather
than letting a NaN propagate silently into every downstream score.

### Emission models

State **M** emits a residue pair from a joint distribution *J*; **X** and **Y**
emit a single residue from the background *p*.

**Nucleotides.** *J* is built directly from an estimated identity θ: diagonal
mass θ/4, off-diagonal mass (1 − θ)/12, background uniform at 0.25. θ is read
off the alignment being scored (below), clamped to [0.55, 0.95].

**Protein.** *J* is BLOSUM45 converted to a proper joint. A BLOSUM file holds
*s*(i,j) = *c* · log₂[ *q*(i,j) / *p*(i)*p*(j) ] for a scale *c* that differs
between matrices, so for residues *i*, *j*

*J*(i, j) = *p*(i) *p*(j) · 2^(s(i,j)/c)

normalised to sum to one, with *p* the Robinson & Robinson background
frequencies and *c* = 3 for BLOSUM45. (NCBI ships 45, 50 and 80 in third-bits
and 62 and 90 in half-bits. Using the wrong *c* does not fail — the
normalisation absorbs it — it just yields a matrix of the wrong peakedness, so
the scales are asserted against the matrices in the test suite.)

BLOSUM45 rather than the more usual BLOSUM62 because it measured better: on
BAliBASE core blocks with the official scorer, +0.026 SP on RV11 (n = 38,
Wilcoxon p = 0.005) and +0.007 on RV12 (n = 44, p = 0.047), with BLOSUM80 worse
than both at each level. **The model still does not adapt to divergence** — one
matrix is used whether a family is at 80% identity or 25%. Since both reference
sets tested sit below 40% identity, the choice is supported in the regime CRAIC
is aimed at and untested above it; see [Limitations](limitations.md).

Unknown residues map to a neutral wildcard (`N` / `X`) contributing no homology
signal, and RNA is aliased onto the DNA model so that `U` is not silently
discarded as unknown.

### Parameters are estimated, not assumed

`progressive.estimate_params` is the one place an alignment becomes pair-HMM
parameters, and both the aligner's own pass and every reliability score go
through it. From the committed alignment it reads:

- **θ**, the mean pairwise identity over aligned residue pairs, clamped to
  [0.55, 0.95];
- **δ**, gap openings divided by (residue positions + gap openings), clamped to
  [0.005, 0.2];
- **ε**, 1 − (openings / gapped positions), clamped to [0.1, 0.9].

The point is that a divergent family is not penalised for being divergent. A
tool that scores every alignment under near-identity defaults will report that a
genuinely hard family is unreliable everywhere, which is true and useless. The
library defaults (δ = 0.02, ε = 0.5) apply only when there is no alignment to
estimate from.

### Forward–backward

The posterior is obtained by forward–backward **in log space**, using the
log-sum-exp identity throughout. For sequences of a few hundred residues the
linear-space products underflow to zero well before the matrix is filled, so
this is a correctness requirement rather than a performance choice.

The result is

**P**[i, j] = Pr(residue *i* of **a** is homologous to residue *j* of **b**)

marginalised over *every* alignment of the pair, each weighted by its
probability under the model. This matrix — not any single alignment — is what
CRAIC treats as the state of knowledge about the pair.

**Implementation.** The kernel is Rust via PyO3, with a NumPy mirror that the
test suite cross-validates to ~10⁻⁹. The Rust core releases the GIL for the
duration of the kernel, so the interface stays responsive while an analysis
runs; this was measured rather than assumed (178 time slices delivered to
another Python thread during one 566 ms call). Because the kernel releases the
GIL, independent pairwise posteriors are computed on threads, one per CPU core
less one left for the rest of the computer, and summed in a fixed order, so
results do not depend on the number of cores.
Without the compiled core everything runs on the mirror and returns identical
numbers, tens of times more slowly.

---

## 2. The consistency transformation

A pairwise posterior knows only about its own two sequences. The ProbCons
consistency transformation propagates evidence through the whole family: each
posterior is re-estimated using every third sequence,

**P′**(x, y) = (1/n) Σ_z **P**(x, z) **P**(z, y)

applied for two iterations by default. If *x* and *y* both align residue *k* of
*z* to each other's residue, that agreement raises the posterior for the pair
even when the direct evidence is weak.

This is the most expensive thing CRAIC does: O(n³) products of L×L matrices per
iteration, and it is what sets the largest family the workbench can handle.

**Sparsity.** Posteriors are overwhelmingly near zero, so the transformation can
threshold them at 0.01 and run the products sparse, as ProbCons itself does.
That holds 1.8–6% of the dense memory, but it is an approximation, and not a
small one on divergent sequences, whose posteriors are spread thinly: on
simulated DNA at 0.5 substitutions per site it cost the built-in engine about 0.1
of sum-of-pairs accuracy. Since 0.5.11 CRAIC therefore uses the exact (dense)
transformation whenever it fits in the memory budget, and the sparse one only
when it does not and SciPy is installed. (0.5.10 and earlier switched to it for
any sequence over 350 residues whenever SciPy was present.)

**Compute levels.** The transformation is what distinguishes them:

| level | consistency iterations | guide tree | refinement |
|---|---|---|---|
| min | 0 | k-mer | 0 |
| med (default) | 2 | posterior | 0 |
| max | 2 | posterior | 25 |

---

## 3. Building an alignment

### Guide tree

At *med* and above, guide-tree distances come from the posteriors already
computed:

d(x, y) = 1 − min(1, Σ **P**(x, y) / min(|x|, |y|))

— one minus the expected number of aligned residues over the shorter sequence.
This is a much better signal than k-mer overlap and is essentially free, since
the posteriors exist already. At *min*, where no all-pairs posteriors are held,
the tree falls back to cosine distance on 3-mer counts.

Sequences are merged by nearest-neighbour joining with proportional averaging of
distances to the merged profile.

### Maximum-expected-accuracy decoding

Given the posterior, an alignment still has to be chosen. CRAIC maximises the
**expected number of correctly aligned residue pairs**: over all alignments **A**,

maximise Σ_{(i,j) ∈ A} **P**[i, j]

by a dynamic program over the posterior matrix. Two things follow, and they are
the reason this is worth explaining at length.

**There is no gap penalty.** None. Leaving a residue unpaired simply contributes
nothing to the sum; the cost of a gap is already priced into the posterior by the
HMM's δ and ε. Gap placement is therefore a consequence of the model rather than
of a tuned penalty, which is also why changing δ and ε changes the alignment in a
way that has a probabilistic interpretation.

**Why not the most probable alignment?** Viterbi decoding returns the single
highest-probability path. That optimises the probability of getting the *entire*
alignment exactly right — a quantity that is essentially zero for any real
family, and one an alignment can maximise while still being wrong in most of its
columns. MEA instead optimises the number of individual homology statements
expected to be correct, which is both achievable and the thing the SP score
actually measures.

Profile–profile merges use the summed posterior over the residue pairs the two
profiles would align, so the same criterion applies at every internal node.

### Streaming, for families that do not fit

Holding every pairwise posterior is O(n²L²). Past a memory budget the aligner
switches to computing each posterior on demand and discarding it, so peak memory
is one matrix rather than n². This buys scale and gives up the consistency
transformation, which is where much of the accuracy comes from — the trade is
explicit and reported.

### Codon-aware alignment

Protein-coding nucleotides are translated, aligned in amino-acid space, and the
nucleotides threaded back through that alignment, so gaps land on codon
boundaries and the reading frame survives (the MACSE and PRANK lineage). CRAIC
assumes frame 0 and intact frames; it warns about broken frames and shows
partial-gap codons as `X`, but it does not find the frame, handle frameshifts,
or accept sequences beginning mid-codon.

Translation uses an **NCBI genetic-code table**, carried on the alignment as
`CodingSpec.table` and selectable in the interface, on the command line
(`--code N`) and in a saved session. This is not cosmetic. A stop codon
translates to `*`, which is not one of the twenty residues the emission model
knows, so it encodes to the wildcard and emits at background frequency —
contributing nothing to any posterior. Reading a vertebrate mitochondrial or
Mollicute gene under the standard code therefore does not merely mislabel the
amino-acid view: it deletes the homology signal at every TGA, and TGA is those
genomes' main tryptophan codon. Tryptophan is among the most conserved residues
in a protein, so the positions lost are disproportionately the ones an aligner
would have used as anchors.

Because that failure is silent, `Alignment.internal_stops()` reports every stop
codon that is not a sequence's last, and both front ends warn about them. The
check cannot distinguish a wrong code from a wrong reading frame or a
pseudogene, which is precisely why it is worth making: one test catches all
three.

---

## 4. Reading the posterior as reliability

Two independent signals, both in [0, 1], both computed under parameters
estimated from the alignment being scored.

### Consistency score

For each column, the mean posterior over the residue pairs that column asserts.
A column grouping residues that the pair-HMM considers unlikely to be homologous
scores low. Sampled over at most 300 sequence pairs for large families.

In the spirit of T-Coffee's TCS and of ZORRO, but computed over CRAIC's own
posteriors — the numbers are **not** interchangeable with either method's, and a
threshold tuned on one does not transfer.

Reported per column and per residue. A residue with no partners in its column
asserts nothing and is left undefined rather than scored a free 1.0.

### Perturbation score

Re-align the same sequences 16 times under perturbed guide trees, cycling
through three gap regimes (δ = 0.01, 0.03, 0.08), and for each residue score the
fraction of its originally-asserted homologies that survive the ensemble. Guide
trees are perturbed by bootstrap-resampling the k-mer feature columns, with light
jitter on near-ties.

**This is not GUIDANCE**, and the difference matters enough to state twice.
GUIDANCE bootstraps alignment columns, rebuilds guide trees, and re-runs *the
aligner that produced the reference*, typically ~100 times. CRAIC runs 16
replicates of *its own* engine under plain posterior decoding — a weaker
aligner, and a different one entirely if the reference came from MAFFT or PRANK.
In that case the score conflates genuine alignment uncertainty with systematic
between-method difference, and the disagreement map is the better instrument.
Sixteen replicates is also a small ensemble, so the score is noisy at the level
of an individual column.

If some replicates fail the survivors are used and the count is recorded; if none
completes, that is an error rather than a score of zero.

### Combining, and masking

The combined score is the mean of the two. The mask is a threshold on whichever
track is active, previewed live before anything is exported. Columns with no
score are kept, not dropped — an unscorable column is not a bad column.

Whether masking helps is an empirical question and the answer is mixed; see
[Validation](validation.md), including the finding that masking makes
neighbour-joining tree recovery *worse*.

---

## 5. Reading the posterior other ways

**Posterior explorer** renders **P** directly as a heatmap for a selected region.
Off-diagonal mass is the set of alternative homology hypotheses the aligner
rejected. A residual mode subtracts the committed alignment, leaving only the
alternatives; a profile mode shows one sequence against the rest.

**Disagreement map** re-aligns with every installed engine (or, with none, under
several gap regimes of the built-in one), projects each result onto a common
coordinate by shared residue identity, and colours each column by how often the
methods reproduce it. This is model-free evidence, and where it disagrees with
the model-based score, that disagreement is itself informative.

**Realignment sandbox** re-aligns a selected block only, under alternative
engines and gap regimes, and splices a chosen alternative back without disturbing
the rest.

**Model-free trimmers** — gap fraction, trimAl's gap and similarity scores with
`gappyout` and `strict`, Gblocks-style blocks, and OD-seq-style sequence outlier
detection — are faithful, dependency-free reimplementations offered as
alternative tracks, so that a model-based and a rule-based judgement of the same
alignment can be compared directly.

---

## 6. Measuring accuracy against a known answer

One module (`craic.evaluate`) serves the GUI's truth mode, the command line and
the benchmark harness, so a number shown to a student and a number printed in a
paper cannot diverge.

Given an inferred alignment and a trusted reference of the *same* sequences
(matched on identifier, and refused outright if the residues differ):

| quantity | definition |
|---|---|
| **SP** | fraction of the reference's homologous residue pairs recovered — recall |
| **precision** | fraction of the asserted pairs that are true |
| **TC** | fraction of reference columns reproduced exactly |
| per-column correctness | fraction of a column's asserted pairs that are true; `nan` where it asserts none |
| per-residue correctness | fraction of a residue's column-partners that are truly homologous |
| **reliability AUC** | ROC area of a reliability score predicting which columns are correct |

Per-column and per-residue correctness are exactly consistent: a column's score
is the mean of its residues'. The AUC uses midranks, so heavily tied scores —
which reliability scores always produce — do not make the answer depend on input
order.

The machinery behind all of this is one lookup, `true_column_grid`: entry [s, c]
is the column the reference puts the residue that the inferred alignment placed
at row *s*, column *c*. Two residues are homologous exactly when their entries
are equal.

**Scope caveat.** Published BAliBASE SP and TC figures are restricted to the
reference's core blocks, using the official `bali_score`. An all-column figure is
a different and generally lower quantity. CRAIC reports which scope it used, and
the benchmark produces both.

---

## 7. Complexity

| step | time | memory |
|---|---|---|
| one pairwise posterior | O(L²) | O(L²) |
| all pairs | O(n²L²) | O(n²L²), or O(L²) streaming |
| consistency transformation | O(n³L²) per iteration | O(n²L²), ~2–6% sparse |
| progressive merge | O(n²L²) | O(nL) |
| consistency score | O(n²L²) sampled to 300 pairs | O(nL) |
| perturbation score | 16 × a full alignment | O(nL) |

The cubic term is the binding constraint. Tens of sequences is comfortable, a few
hundred is slow, a few thousand is out of scope — align those with MAFFT and, if
you need per-column reliability on them, use TCS or GUIDANCE.

---

## Lineage

Nothing here is a new statistical idea. The debts, in order of size: ProbCons
(Do *et al.*, 2005) for posterior decoding and the consistency transformation;
T-Coffee's TCS (Chang *et al.*, 2014) and ZORRO (Wu *et al.*, 2012) for
consistency-based column reliability; GUIDANCE (Penn *et al.*, 2010; Sela *et
al.*, 2015) and Heads-or-Tails (Landan & Graur, 2007) for perturbation;
trimAl (Capella-Gutiérrez *et al.*, 2009), Gblocks (Castresana, 2000) and OD-seq
(Jehl *et al.*, 2015) for the model-free filters; PRANK (Löytynoja & Goldman,
2008) and MACSE (Ranwez *et al.*, 2011) for gap placement and codon awareness.

What CRAIC contributes is putting them in one place, computing them where the
alignment lives, keeping them live under editing, and letting them be checked
against a known answer.
