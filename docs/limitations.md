# Limitations

What CRAIC does not do, does not do well, or does not do the way you might
assume. This page exists because a tool that reports uncertainty has no business
being coy about its own.

## The built-in aligner is not competitive

It is a ProbCons-style maximum-expected-accuracy aligner with a consistency
transform, written so that CRAIC works with nothing else installed and so that
the perturbation analysis has an engine to re-align with. On BAliBASE it is
beaten by MAFFT — see [Validation](validation.md) for the measured gap. It has no
sequence weighting, no iterative refinement below the *max* compute level, and no
empirically-tuned parameters.

Use an external engine for anything you intend to publish. CRAIC drives MAFFT,
MUSCLE, Clustal Omega and PRANK if they are on your `PATH`, and the whole
uncertainty layer works identically on their output.

## It does not scale to large families

The consistency transformation is O(N³) in the number of sequences, and holding
every pairwise posterior is O(N²L²) in memory. Past a size threshold the engine
switches to streaming — computing each posterior on demand and discarding it,
which caps peak memory at a single matrix — but that trades away the consistency
transform, which is where much of the accuracy comes from.

Practically: tens of sequences is comfortable, a few hundred is slow, and a few
thousand is out of scope. The reliability analysis has the same shape, so the
same limit applies to scoring an alignment built elsewhere. Large families should
be aligned with MAFFT and, if you want per-column reliability on them, scored
with TCS or GUIDANCE rather than with CRAIC.

The reliability analysis is dominated by the perturbation score, which
re-aligns the sequences sixteen times. For 36 proteins of 1,000 residues the
consistency score takes about a minute of processor time and the perturbation
score about 40 minutes, shared among the cores. The window shows the time left,
asks before a perturbation run of more than two minutes (the consistency score
alone is always available), and can cancel at any point.

## The perturbation score is not GUIDANCE

It is inspired by GUIDANCE and it is not equivalent to it, in two ways that
matter:

- GUIDANCE bootstraps the guide tree and re-runs **the aligner that built the
  alignment**, typically ~100 times. CRAIC runs **16 replicates of its own
  engine**, regardless of which engine produced the alignment being scored. So
  when you score a MAFFT alignment, the perturbation component is asking how
  stable those columns are under *CRAIC's* model, not under MAFFT's.
- 16 replicates is a small ensemble. The score is correspondingly noisy at the
  level of an individual column.

It is a cheap, weak relative of the published method. If you need GUIDANCE
numbers, run GUIDANCE — the server at `guidance.tau.ac.il` is live and
GUIDANCE3 is pip-installable.

## The consistency score is not TCS

Same caveat in a different place. It is in the spirit of T-Coffee's TCS and of
ZORRO — a posterior-based measure of how well the evidence supports each column —
but it is computed over CRAIC's own pair-HMM posteriors, so the numbers are not
interchangeable with either method's, and a threshold tuned on one does not
transfer.

## Masking is reported honestly, and honesty is unflattering

It is easy to make masking look good and hard to show that it helps. Removing
the hardest columns raises mean column accuracy whatever rule picks them, so
"accuracy after masking is higher than before" is not evidence of anything. The
benchmark therefore reports the kept-column accuracy against two controls:
removing the same *number* of columns at random, and removing them by the trivial
gap-fraction rule. Only the margin over those controls is informative, and it is
much smaller than the naive comparison suggests.

Whether masking improves a downstream tree is a separate question again, and the
answer in the literature is genuinely mixed. CRAIC's benchmark measures it
(Robinson–Foulds distance, simulation mode only) rather than assuming it.

## Protein emission modelling is simplistic

The protein emission model is BLOSUM45 converted to a joint distribution via
Robinson background frequencies, and it does not adapt to the divergence
estimated from the data — so a family at 80% identity and one at 25% identity get
the same matrix. Nucleotide data does get data-estimated parameters. Fixing this
for proteins is a known item, not a design decision.

BLOSUM45 was chosen by measurement rather than convention: it beat BLOSUM62 on
both BAliBASE reference sets tested, and BLOSUM80 lost to both. But RV11 and
RV12 are both below 40% identity, so that result says nothing about closely
related sequences, where BLOSUM62 would likely be the better default. A single
matrix chosen for the low-identity case is a reasonable default for a tool aimed
at ambiguous alignments; it is not the right answer for a family at 80%
identity. Keying the matrix to the estimated divergence is the fix, and it would
need the raw identity rather than the clamped ``theta`` the aligner currently
computes.

## Codon-aware alignment assumes a lot

Translate → align → back-translate assumes reading frame 0 and intact frames.
CRAIC warns when frames are broken and shows partial-gap codons as `X`, but it
does not find the frame for you, handle frameshifts, or deal with sequences that
begin mid-codon. MACSE is the tool for data where that matters.

The **genetic code** is chosen, not assumed: the code selector beside *treat as
coding* offers every NCBI translation table, `craic align --codon --code N` does
the same headless, and `craic codes` lists them. It still defaults to the
standard code, so mitochondrial and Mollicute data need the choice made
deliberately. CRAIC warns when a coding alignment contains internal stop codons,
which is the usual symptom of getting it wrong — though that warning cannot
distinguish a wrong code from a wrong frame or a pseudogene.

What it costs to get wrong is smaller than you might expect but not nothing.
Reading Mollicute-like sequences under the standard code turns roughly 1.6% of
codons (the TGA tryptophans) into stops, and a stop is not one of the twenty
residues, so it encodes to the wildcard and contributes no homology signal.
Measured over twelve simulated datasets, correcting the code moved mean SP by
+0.004 and TC by +0.014, and changed the alignment at all in only two of the
twelve. The larger problem is the amino-acid view being simply wrong, and, before
0.5.1, being wrong silently.

## What a reference alignment can and cannot tell you

Truth mode is exactly as good as its reference.

- **Simulated data** gives per-residue truth, which is why the teaching mode uses
  it — but the sequences are generated under a model, and if that model is close
  to CRAIC's own pair-HMM the comparison flatters CRAIC relative to aligners
  carrying empirical parameters. Turning on among-site rate variation makes it a
  harder and fairer test, and is worth doing before drawing conclusions.
- **Structural references** such as BAliBASE are real data, but they are only
  confident about their core blocks. Published BAliBASE SP and TC figures are
  restricted to those blocks; an all-column figure computed over the whole
  alignment is a different quantity and is not comparable with them. CRAIC's own
  scoring reports which scope it used, and the benchmark writes a manifest for
  the official `bali_score` so that the published-comparable numbers can be
  produced.

## Scope

No tree inference, no structure, no annotation transfer, no database search, no
support for very large alignments as a display problem (the canvas is
virtualised, but the analyses are not). CRAIC ends where the alignment ends.

## Platform

Developed and used on macOS and Linux. The GUI is PySide6 and should run on
Windows; it is not routinely tested there. The Rust core is optional everywhere —
without it the NumPy mirror gives identical results more slowly, and the status
bar says which is in use.
