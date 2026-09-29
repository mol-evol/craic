# Changelog

## 0.5.14 — "align as protein" keeps every nucleotide

**Aligning as protein dropped an incomplete last codon.** A sequence whose
length is not a multiple of three lost its last one or two nucleotides when
aligned with "align as protein". Besides losing data, the result no longer
matched a reference alignment of the same sequences, so in teaching mode a
simulated dataset's answer was silently dropped and truth mode then said no
reference was loaded. Those nucleotides are now kept, in columns after the last
codon, outside the reading frame.

**An aligner that changes the sequences is reported.** After aligning, CRAIC
checks that every sequence came back with the same residues, gaps apart, and
warns if not, saying so when the reference no longer applies. The truth track
now explains a missing answer the same way the truth colouring does. 283 tests
pass.

## 0.5.13 — one core left free

Long analyses now use one thread fewer than the computer has cores, leaving a
core for the operating system and other programs, so a run of several minutes
does not make the rest of the computer sluggish. It costs about an eighth of the
speed-up on an eight-core machine; a two-core machine runs single-threaded, as
before 0.5.12. Results are unchanged.

## 0.5.12 — long analyses use every core, show the time left, and can be stopped

A Windows user scored 36 long protein sequences and the window appeared to
hang: the reliability analysis was running, single-threaded, behind nothing but a
status-bar message. It is mostly the perturbation ensemble, which re-aligns the
sequences sixteen times; for 36 sequences of 1,000 residues that is about 40
minutes on one core, and tens of times longer without the compiled core. Four
changes:

- Pairwise posteriors are computed on all CPU cores (`accel.imap`, threads over
  the Rust kernel, which releases the GIL): the consistency score, the aligner's
  all-pairs and streaming passes, and so every perturbation replicate. Results
  are summed in a fixed order, so they are identical whatever the thread count;
  re-running benchmark rows reproduces them exactly. Measured on two cores, the
  36-sequence case is 1.9 times faster.
- The reliability analysis and the aligner-agreement track run behind a progress
  window with the time left and a Cancel button.
- The consistency score comes first. Its time predicts the perturbation
  ensemble's (every sequence pair costs about the same), and if that will be over
  two minutes CRAIC asks before running it; the user can keep the consistency
  score alone, and cancelling the ensemble part-way keeps it too. Choosing the
  Perturbation track later offers the run again. Without the compiled core, the
  progress window and the question say so.
- Exporting the confidence report no longer computes the analysis on the
  interface thread, which froze the window.

`reliability.analyse`, `consistency` and `perturbation` take `progress` and
`cancelled`; `reliability.work` gives the pairs each stage aligns and
`reliability.combine` assembles a report from the two stages. 280 tests pass.
Results are unchanged from 0.5.11.

## 0.5.11 — the exact consistency transformation

**The built-in engine no longer approximates when it does not have to.** With
SciPy installed, 0.5.10 and earlier switched the consistency transformation to
its sparse form, which drops posteriors below 0.01, for any family with a
sequence of 350 residues or more, to save time. On divergent sequences, whose
posteriors are spread thinly, that approximation is not small: on the simulated
DNA it cost the built-in engine about 0.1 of sum-of-pairs accuracy at 0.5
substitutions per site. It also made results depend on whether SciPy happened to
be installed, and SciPy is an optional extra, so the benchmark (with SciPy) and
a default install (without) aligned differently. The exact (dense)
transformation is now used whenever it fits in the memory budget, and the sparse
one only for families too large for it.

Re-run with it, the built-in engine's benchmark numbers change. BAliBASE: SP
0.827 → 0.830, TC 0.674 → 0.676 (161 of the 163 families exact; the two largest
still need the approximation within the 3 GB budget), now slightly ahead of MAFFT
(p = 0.045). Simulations: SP 0.539 → 0.604, level with MAFFT's 0.608, the gain
coming at high divergence (0.598 → 0.697 at 0.5 substitutions per site); the
reliability AUC on its alignments is 0.855 (was 0.875). No other aligner is
affected. `docs/validation.md` and `benchmarks/results/` carry the new rows, and
two new tests check which transformation `align` chooses. 269 tests pass.

**ClustalW in the benchmark.** ClustalW 2.1 joins the BAliBASE and simulation
benchmarks (BAliBASE SP 0.743, TC 0.551); its rows are added to `balibase.csv`,
`balibase_official.csv` and `simulation.csv`.

**The BAliBASE tree comparison behind the paper's figure 2**:
`benchmarks/bali_trees.py` and its results, `bali_trees_mafft.csv` and
`bali_trees_builtin.csv`.

**The simulation command in `docs/validation.md` is corrected.** It gave
`--perturbation` and left the built-in engine at the default effort, whereas the
published figures are the consistency score alone with the built-in engine at
`--effort med`, so running it as written produced different numbers.

## 0.5.10 — every engine's settings, and ClustalW

**PRANK was run with the wrong gap costs.** CRAIC always passed
`-gaprate=0.025 -gapext=0.5`. PRANK's own defaults depend on the data — DNA
0.025 / 0.75, protein 0.005 / 0.5 — so every protein run used five times PRANK's
gap-opening rate and every DNA run a lower extension probability, whether or not
anyone had touched the settings. The BAliBASE PRANK rows in `docs/validation.md`
are re-run with PRANK's own defaults: SP 0.686 → 0.714, TC 0.470 → 0.506. Every
other aligner's numbers, and every conclusion, are unchanged.

The fix is general: a numeric parameter can now default to *blank*, which means
CRAIC passes no flag and the tool chooses. The dialog shows the field empty with
what the tool will do instead ("blank = PRANK default (DNA 0.025, protein
0.005)"), and the History leaves those settings out.

**Settings reached only two of the places an engine runs.** The ⚙ values were
used by **Align** and the **Aligner agreement** track, but not by the
realignment sandbox, **Realign around pinned columns**, the ensemble
co-occurrence figure, or MUSCLE at all. They now are.

**More of each tool's settings are offered.**

- Built-in: gap open and extension probabilities (δ, ε; blank estimates them
  from the data, as before), the protein matrix (BLOSUM45, 50, 62, 80, 90), and
  consistency and refinement passes.
- MAFFT: local-pair gap costs `--lop` / `--lep` (the ones L-INS-i and E-INS-i
  use), `--maxiterate`, the protein matrix (`--bl`, `--jtt`, `--tm`) and the DNA
  model (`--kimura`).
- MUSCLE 5: `-super5`, `-perm` and `-perturb`. MUSCLE 5 has no gap penalties.
- Clustal Omega: `--full` and the guide-tree and HMM iteration caps. It has no
  gap penalties either.
- ProbCons: consistency passes (`-c`) and pre-training (`-pre`).
- PRANK: `+F`, `-iterate` and `-termgap`.

**ClustalW is a new engine**, for teaching gap costs: opening and extension
penalties for both the pairwise and the progressive stage, the protein matrix
series and DNA matrix, gap separation distance, and the residue-specific and
hydrophilic gap adjustments. Its gap costs are blank by default because
ClustalW's defaults also differ between DNA and protein.

**The settings dialog follows the data.** It shows only the settings that apply
to what the engine will be handed — protein or nucleotide, and protein when
*align as protein* is on — and a blank field names the default the tool will
use for that data ("blank = ClustalW default, 10"). Hidden settings keep their
values, and the History lists only the settings that applied. `craic engines`
marks the one-alphabet settings.

**Residue masks.** Masking a column throws away its correctly aligned residues
along with the wrong ones, and on a tree that costs more signal than the error
it removes. A residue mask removes only what you judge wrong: select residues and
**Mask selected residues** (Ctrl/⌘K), or **Mask residues below the threshold**
to mask every residue whose own reliability is below the Mask slider. The
columns are kept; **Export masked…** writes each masked residue as missing data
(`N` or `X`). The mask names residues, not columns, so it follows them through
edits and realignment, is saved in sessions, and every change is on the undo
stack and the History. `craic mask --residues` does the same on the command
line.

**One undo step per edit.** A hand edit pushed the previous alignment onto the
undo stack twice, so the second Undo after one edit appeared to do nothing.

**Does masking residues rescue the tree? No.** On the 360 simulated datasets,
masking the lowest-scoring residues — exactly as many as the column mask removes
— gives worse trees than masking the columns (Robinson–Foulds +0.09,
neighbour-joining and maximum likelihood alike), and only a little better than
masking as many residues at random. The residues the score flags carry signal
as well as error at the default threshold; residue masks are for hand curation,
not a better automatic filter. Details in `docs/validation.md`.

**Benchmark harness.** `--engines` restricts a run to named engines (used to
re-run PRANK alone); `--tree ml` measures tree error with IQ-TREE (JC+G4)
instead of neighbour-joining; simulation runs report tree error when the same
score masks residues instead of columns — at the threshold, matched to the
number of residues the column mask removed, and for that many at random. A
results file written by an older harness is rewritten under the current header
when a run resumes, rather than having rows appended that do not line up.

**Credit and citation.** A **Help** menu opens the website and an About box that
names the author, links to <https://mol-evol.github.io/>, and gives the citation with a
button that copies it; the splash screen shows the author's website, and
`craic --version` and `craic --help` point to the About page. The citation (the bioRxiv
preprint while the paper is under review) is in `CITATION.cff`, the README and the new
**About & citation** page, and a test checks they agree. The website gains a home page,
top-level tabs, and a footer crediting the author on every page.

**A proper app icon.** A grid of alignment cells over a reliability track, with one
weak column in amber, laid out on the macOS icon grid. It is drawn in code, so the
running app and the downloadable apps (which had PyInstaller's default icon) use the
same one. The install guide now gives macOS 15's way of opening an unsigned app.

**`craic make-app`** (macOS) builds a double-clickable CRAIC.app in `~/Applications`
that runs the Python it was made from. Made on the user's own Mac, it is never marked
as downloaded, so Gatekeeper does not block it, which the ready-made download cannot
avoid without a paid Apple signature. It replaces `scripts/make_app.py`, which only
existed in the repository, so pip users could not reach it.

**CI's NumPy-fallback job failed** before running a test: `pytest` (unlike
`python -m pytest`) does not put the project folder on the path, and that job does not
install CRAIC. `pyproject.toml` now sets `pythonpath = ["."]` for pytest.

**Raw results.** `benchmarks/results/` now holds the PRANK re-run (in place of
the old PRANK rows in `balibase_official.csv` and `prank_percol.csv`) and the
residue-masking runs (`residue_masking_nj.csv`, `residue_masking_ml.csv`); its
README says which file is behind which table.

**Command line.** `craic align --param KEY=VALUE` (repeatable) sets any engine
parameter, checked against the engine's own list, and `craic engines` lists them
with ranges and defaults. A new parameter on an engine appears in both without
touching the CLI.

Every engine was run on DNA and protein with each of its settings changed in
turn, against MAFFT 7.505, MUSCLE 5.2, Clustal Omega 1.2.4, ProbCons 1.12,
PRANK v.170427 and ClustalW 2.1. 263 tests pass with the Rust core, 260 + 3 skipped on NumPy.

## 0.5.9 — the benchmark, re-run and re-scoped

`docs/validation.md` is regenerated. Every number on it changed, because the
emission model changed in 0.5.3 and the reliability score reads CRAIC's
posteriors whoever built the alignment.

**Six engines, not four.** ProbCons and PRANK join MAFFT, MUSCLE, Clustal Omega
and the built-in engine. ProbCons because the built-in engine follows its method
and a table without it invites the reader to blame the approach for what is this
implementation; PRANK because CRAIC wraps it and quietly omitting an engine the
workbench offers is a selection the reader cannot see.

**RV11 and RV12 only, complete.** 163 families: RV11 76/76 and RV12 87/88, the
one exclusion being a family whose dense posteriors exceed the memory budget for
every engine alike. RV20–RV50 are dropped rather than partially reported — which
families survive is decided by CRAIC's memory ceiling rather than by the
aligners, so a partial set is a biased sample, and a biased sample is worse than
an absent one.

| aligner | SP | TC |
|---|---|---|
| MUSCLE 5.1 | 0.856 | 0.712 |
| ProbCons 1.12 | 0.850 | 0.705 |
| CRAIC built-in | 0.827 | 0.674 |
| MAFFT 7.505 | 0.818 | 0.667 |
| Clustal Omega 1.2.4 | 0.766 | 0.593 |
| PRANK v.170427 | 0.686 | 0.470 |

The built-in engine and MAFFT are statistically indistinguishable (paired
Wilcoxon p = 0.12); it is behind MUSCLE and ProbCons and ahead of Clustal Omega.
The previous page's “a little behind MAFFT” is no longer supportable in either
direction.

**A new section: no method reconstructs these alignments.** Six aligners, 163
families, and 10 families are reconstructed exactly by any of them — none by all
six. An oracle picking the best method per family reaches TC 0.739 against
MUSCLE's 0.712, and still leaves a quarter of columns misplaced. This is the
premise the workbench rests on, now measured rather than asserted.

**Two tables that were being read wrongly are now labelled.** The reliability AUC
runs almost opposite to accuracy — highest on the two least accurate aligners,
because a worse alignment has more wrong columns to find — so it is not a ranking
of aligners and now says so. And the masking section leads with the margin over
the gap-fraction rule (≈ 0.10) rather than the margin over random (≈ 0.19), since
gap fraction is one line of code and is the bar the model-based score has to
clear.

## 0.5.8 — engine parameters know their own limits

The parameters dialog put the same ``QIntValidator(-1000000, 1000000)`` behind
every integer field, so it happily accepted 10000 refinement passes for ProbCons,
which takes 0 to 1000 — and the run then failed after it had started.

``Param`` now carries optional ``lo`` and ``hi`` bounds, and the dialog uses them:
the validator refuses an out-of-range keystroke, the field shows the range as
placeholder text, and the value is clamped when it is read, since paste can get
past a validator. Bounds are declared where the tool documents them: ProbCons
refinement 0–1000, PRANK gap rate and gap extension 0–1 (they are
probabilities), Clustal Omega iterations 0–100.

## 0.5.7 — align as protein unlocks the protein-only engines

0.5.6 greyed out engines that cannot align the current data, but judged them
against the document's own alphabet. With *align as protein* ticked the
nucleotides are translated before the engine sees them, so ProbCons is perfectly
usable on coding DNA — it is handed amino acids. The check now asks which
alphabet the engine will actually be given, and re-runs when the checkbox is
toggled, so ticking it enables ProbCons rather than leaving it greyed out.

## 0.5.6 — engines say what they can do; failures stop shouting

Adding ProbCons put the first protein-only engine in the list, and nothing in the
interface knew that an engine could have limits. Selecting it on nucleotide data
started the run, raised inside the worker thread, and put a Python traceback on
screen. Three fixes:

* ``AlignerEngine.supports(alphabet)`` — asked before a run starts. Only
  ProbCons answers no, and only for nucleotides.
* Engines that cannot take the current data are greyed out in the **Engine**
  menu, with a tooltip saying why, and stay visible rather than disappearing. If
  the selected engine becomes unusable the selection falls back to the built-in
  one. Choosing one anyway explains the limit and points at *align as protein*
  for coding nucleotides.
* **An alignment failure no longer shows a traceback.** The dialog leads with the
  error's own message, says the alignment on screen is unchanged, and keeps the
  traceback behind *Show Details* — where it is useful to whoever needs it and
  invisible to everyone else.

**Go to next ambiguous region** was the one command that computed reliability
synchronously, freezing the window while the sixteen-replicate perturbation
ensemble ran. It now goes through the same background path as every other
consumer of the analysis, with the progress message and the shared cache.

## 0.5.5 — pinned columns actually exist

**Realign around pinned columns** named a feature that was never built. Nothing
in the package ever wrote ``meta["anchors"]``: the menu item read it, found it
empty every time, showed a dialog advising the user to pin some columns first —
an operation CRAIC did not offer — and then quietly fell back to using the two
ends of the current selection. The documentation described pinning too.

Pinning is now real:

* **Pin selected columns** (`Ctrl/⌘ P`), **Unpin selected columns**
  (`Ctrl/⌘ ⇧ P`) and **Clear all pins**, in *Edit ▸ Edit alignment*.
* Pinned columns are drawn with a gold bar above the grid and a faint wash, above
  whatever colour mode is on — a pin is the user's own assertion and should not
  be mistakable for a computed score.
* Pins are stored on the alignment as a sorted list of column indices, so they
  save and reload with a session with no extra code.
* **Pins follow their columns through an anchored realign.** Re-solving the
  blocks between anchors changes their widths, so the anchors land at new
  indices; ``editing.anchored_realign`` now reports where they went and the
  window adopts those. Any other edit that changes the number of columns drops
  the pins instead of keeping indices that would silently re-pin whatever moved
  into that position.
* The dialog shown when nothing is pinned now says how to pin something.

## 0.5.4 — ProbCons as an engine, and a better default

CRAIC now wraps **ProbCons** (Do et al., 2005) alongside MAFFT, MUSCLE, Clustal
Omega and PRANK. It is protein-only, and like every external engine it appears in
the **Engine** menu only when the ``probcons`` binary is found on the process
PATH.

It is here for a specific reason. CRAIC's built-in engine follows ProbCons's
method — pair-HMM posteriors, consistency transformation, maximum-expected-accuracy
decoding — so a benchmark table that reports the built-in engine against MAFFT
and MUSCLE alone invites the reader to conclude that posterior-consistency
alignment is the weaker approach. Measured on BAliBASE core blocks with the
official scorer, ProbCons scores 0.856 SP against the built-in engine's 0.835, so
most of that difference is this implementation rather than the method: the
built-in engine estimates its pair-HMM parameters per dataset where ProbCons
EM-trains them once, and that is the remaining substantive difference between
them.

The **Engine** menu still starts on the built-in aligner, now as a stated choice
rather than an accident of registry order. It is the only engine guaranteed to be
present, so it is the only one that can be the default on every machine, and it is
the engine whose posteriors the reliability overlay, posterior explorer, sandbox
and perturbation ensemble all read — so starting there makes the workbench
coherent before the user picks anything. On BAliBASE core blocks it is level with
MAFFT (0.835 vs 0.830 SP, n=188, Wilcoxon p=0.51) and ahead of Clustal Omega,
behind MUSCLE and ProbCons; the menu is one click away where that matters.

The benchmark harness gained ``--effort``, defaulting to ``max``. It previously
ran the built-in engine at its ``med`` default, which skips iterative refinement
entirely — measuring the engine with a stage switched off. (For the record,
turning refinement on does not help: 0, 25 and 100 rounds score 0.697, 0.692 and
0.700 SP on a 12-family pilot.)

## 0.5.3 — BLOSUM45

The protein emission model was BLOSUM62, chosen by convention rather than by
measurement. Tested against the alternatives on BAliBASE core blocks with the
official scorer — same guide tree, same gap estimation, same decoder, only the
matrix changing — BLOSUM45 is better at both divergence levels available:

| | RV11 (<20% identity, n=38) | RV12 (20–40%, n=44) |
|---|---|---|
| BLOSUM45 | **0.600** SP / 0.374 TC | **0.922** SP / 0.799 TC |
| BLOSUM62 | 0.574 / 0.347 | 0.915 / 0.780 |
| BLOSUM80 | 0.542 / 0.324 | 0.904 / 0.758 |

BLOSUM45 beats BLOSUM62 by +0.026 SP on RV11 (25 families better, 11 worse, 2
tied; Wilcoxon p=0.005) and +0.007 on RV12 (25/16/3, p=0.047); pooled +0.016,
p=0.0006. The effect is four times larger in the harder set, which is what a
divergence mismatch should look like. A fixed BLOSUM45 also beat a rule that
picks 45 or 62 by estimated identity, so no selection machinery was added.

Both reference sets sit below 40% identity, so this says nothing about closely
related sequences, where BLOSUM62 would likely be the better default. Keying the
matrix to divergence is the real fix and is left for a later version; it would
have to read the raw identity, since ``estimate_params`` clamps ``theta`` to
[0.55, 0.95] and so pins it at 0.55 for every family in RV11 and RV12.

One trap found on the way, now guarded by a test: NCBI ships BLOSUM45, 50 and 80
in third-bits but 62 and 90 in half-bits. Converting log-odds to a joint with the
wrong scale does not raise — the normalisation absorbs it — it just yields a
matrix of the wrong peakedness. A first run of BLOSUM80 under the wrong scale
scored 0.330 on RV11 instead of 0.542, which looks like a result rather than a
bug. ``_BLOSUM_SCALE`` now names each divisor and the test suite checks every one
against the matrix itself.

The trimAl similarity score in ``ambiguity/trimming.py`` still uses BLOSUM62 on
purpose: that matrix is part of the published rule being reimplemented, not a
CRAIC parameter.

Everything in `docs/validation.md` is being regenerated, since changing the
emission model changes every posterior and so every reliability and masking
number, for all four aligners.

## 0.5.2 — the generated answer is scored at once

**Generate dataset with known answer…** attached the simulated truth to the
document but deliberately did not score it, on the reasoning that nothing had
been aligned yet. Truth mode reads its per-residue colours from that score, so
until something else changed the alignment, *Correct placement* refused to turn
on and offered to load a reference alignment — a program that has just generated
an answer telling the user it has none. Attaching the reference now scores it
straight away, against the unaligned sequences, which is a legitimate starting
picture: almost everything is wrong, and the first Align is visibly an
improvement.

The message that appears when there really is no answer named two menus that do
not exist (*Tools ▸ Simulate a dataset…*, *Reference ▸ Open reference
alignment…*); both live under **Teach**. It now names them correctly, and
distinguishes having no reference from holding one that nothing has been scored
against yet.

## 0.5.1 — genetic codes

The codon and amino-acid views, and `align as protein`, always read the standard
genetic code. Mitochondrial and Mollicute genomes do not use it. The most
consequential difference is TGA: tryptophan in vertebrate mitochondria
(table 2), in moulds and protozoa and in *Mycoplasma* and *Spiroplasma*
(table 4), but a stop under the standard code.

That was not merely a display problem. A stop translates to `*`, which is not one
of the twenty residues the emission model knows, so it encodes to the wildcard
and emits at background frequency — contributing nothing to any posterior.
Reading a Mollicute gene with the wrong code therefore deletes the homology
signal at every TGA, and tryptophan is among the most conserved residues in a
protein, so the positions lost are disproportionately the ones an aligner would
have used as anchors.

The data model already carried an NCBI table id on `CodingSpec` and translated
correctly under any table; sessions already round-tripped it. What was missing
was any way to set it. Now:

* a **genetic code** selector sits beside *treat as coding*, listing every table
  Biopython knows (the list is read from Biopython rather than hard-coded, since
  NCBI adds tables and a frozen copy goes stale silently). Changing it re-reads
  the alignment at once.
* `craic align --codon --code N` does the same headless, and `craic codes` lists
  the tables.
* `align as protein` uses the selected code rather than assuming table 1.

### Internal stop codons are now reported

The real defect was silence: nothing said the code was wrong.
`Alignment.internal_stops()` finds every stop that is not a sequence's last
codon, ignoring trailing gaps so a short sequence padded to the alignment width
is not reported as stopping early. Both front ends warn about them, on marking an
alignment as coding, on changing the code, and after aligning.

The check cannot tell a wrong genetic code from a wrong reading frame or a
pseudogene, which is exactly why it earns its place: one test catches all three,
and the message says so rather than guessing.

### How much it was costing

Measured over twelve simulated datasets read as Mollicute data, the wrong code
blanks about 1.6% of codons, and correcting it moves mean SP by +0.004 and TC by
+0.014, changing the alignment at all in two of the twelve. So the alignment cost
is small; the display was simply wrong, and the failure was silent. Both are
fixed. See `docs/limitations.md`.

Nine new tests; the suite is 230 passing with the Rust core.

## 0.5.0 — showing the answer, not just the verdict

Truth mode could tell you that a column was wrong. It could not tell you what the
right answer was, which is most of what a teaching tool is for: a column at 0.8
named neither the residue that broke it nor what should have happened instead.

The two alignments cannot be shown side by side column-for-column — they have
different numbers of columns — but they hold exactly the same residues, so every
residue has a true column. That lookup (`evaluate.true_column_grid`) is what the
three new views are built from.

### Colour by correct placement

A third entry in the **Colour** control paints each residue by the fraction of
its partners in that column that are genuinely homologous
(`evaluate.cell_correct_fraction`). One misplaced residue now shows up red in an
otherwise green column, naming the sequence that broke it.

It is exactly consistent with the correctness track above it — a column's score
is the mean of its residues' scores, which is checked by a test rather than
asserted — and it survives an edit, so dragging a block recolours the residues it
moved. The colour key says which claim the ramp is making, since "unreliable" and
"wrongly placed" are different things painted the same red.

### Column inspector

A third tab in the tools dock states, for the column under the cursor, which
residues the reference groups together, which of them are in place, which
residues do not belong there and how far they would have to move, and which
belong there but were put elsewhere (`evaluate.explain_column`). Every column it
names is a link. It follows the edit cursor and updates after every edit.

It will report a column as incomplete even when that column scores a perfect 1.0,
which is not a contradiction: the score asks whether the homologies a column
*asserts* are true, and a column that quietly dropped a residue asserts nothing
false. Both facts are now visible at once.

### Show the true alignment

**Teach → Show the true alignment** (`Ctrl+T`) puts the reference on screen,
read-only, and `Ctrl+T` returns to the working alignment with the view restored.
The document, the undo stack and the dirty flag are untouched. Editing is refused
at the canvas — the slide keys and drag-to-slide never fire — and again at
`_apply_edit`, because a keystroke aimed at the truth's columns would otherwise
land on the working alignment's.

### Validation

The reliability claim is now measured on real data rather than asserted. 181
BAliBASE 3 families (RV11 and RV12 in full, 17 of RV20) were re-aligned with the
built-in engine, MAFFT 7.505, MUSCLE 5.1 and Clustal Omega 1.2.4 and scored with
the official core-block `bali_score`. The consistency score predicts which
columns are truly correct with a ROC AUC of 0.81–0.86, and does so at least as
well on alignments CRAIC did not build as on its own. Masking beats a matched
random mask by about 0.20 and the trivial gap-fraction rule by about 0.09 in
mean retained-column accuracy. See `docs/validation.md`.

A 360-dataset simulation sweep (six divergence levels, three indel rates, 8 and
16 taxa) adds the gradient the real data cannot give: alignment accuracy falls
off a cliff between divergence 0.25 and 0.5 for every method, while the
reliability AUC stays between 0.83 and 0.90 throughout — the score does not
degrade where the alignment does. It also reproduces the negative result
prominently: masking improves per-column accuracy but makes neighbour-joining
tree recovery *worse* at every divergence above the trivial (Robinson–Foulds
0.16 → 0.47 at divergence 0.5), because it discards the phylogenetic signal
along with the ambiguity. The masking slider is not an automatic filter and the
documentation now says so in several places.

`benchmarks/score_balibase.py` is new: it converts the harness's output to the
GCG/MSF layout the official scorer requires and collects core-block SP and TC,
so the published-comparable numbers can be reproduced.

### Documentation

New pages: background and intent, validation, limitations, FAQ and
troubleshooting, a glossary, and an API reference generated from the docstrings
(`mkdocstrings`). The architecture section gained a walk-through of how the
pieces fit. `scripts/make_screenshots.py` regenerates every image in the docs by
driving the real window offscreen over seeded data, so the screenshots stop
drifting away from the interface.

### Also

* `evaluate.matched_rows` factored out of `compare_to_reference`, so the id
  matching and same-sequences check have one implementation.
* `canvas.scroll_to_nt` — callers working in alignment coordinates no longer need
  to know that a codon view packs three of them into one display column.
* `serve-docs.command` no longer looks for a `.venv` inside a cloud-synced
  folder, where `build_rust.sh` deliberately does not put one.
* Nine new tests. The suite is 221 passing with the Rust core (218 passing and 3
  skipped on the NumPy fallback).

## 0.4.1 — structural clean-up, no behaviour change

No new features and, by design, no different answers: the alignments, reliability
scores, masks, accuracy figures and agreement maps produced by this release are
byte-identical to 0.4.0 on the same inputs. That was checked directly — a
fingerprint of nine outputs across three simulated datasets, compared before and
after — rather than inferred from the tests passing.

The point was to remove the structures that had produced bugs.

### A document, separate from the window

`CraicWindow` was 2,296 lines and about 180 attributes, holding widgets and
working state in one object. Two bugs came out of that directly: background
analyses had no way to tell the alignment had moved on (fixed in 0.2.0 with a
generation counter bolted onto the window), and "everything that replaces the
document" was not a concept anywhere, so the unsaved-work guard was added to one
of the four places that needed it and missed three.

`craic.gui.document.Document` now owns the alignment, the reference and the
accuracy against it, groups, annotations, paths, the dirty flag and the
generation, and emits signals when they move. The window owns widgets. The old
attribute names (`win.aln`, `win.reference`, `win.truth`, `win._dirty`) remain as
thin properties onto the document, so nothing outside had to change at once.

### One description of the tracks

The column overlays were three parallel encodings of the integers 0 to 11: a
label list, a set of which were cheap enough to survive an edit, and an if-chain
dispatching on index. Adding one in the middle broke the other two silently.
`craic.gui.tracks` describes each track once — label, kind, how it is computed —
and the combo box, the cheap-track rule and the dispatch all read from it. There
are now no magic track indices in `app.py`, and a test asserts the three stay
consistent.

### One residue-index implementation

Eight near-identical column-to-residue walks across six modules became
`domain.residue_index`, `domain.column_index` and `domain.membership`. This is
the duplication the 0.2.0 gap-character bug lived in: the copies tested
`ch != "-"` while the definition of a gap said otherwise, and a drifted copy
shifts every residue index and corrupts every score computed from it. It now sits
beside the definition of what a gap is.

### Smaller

- The dialogs (engine parameters, figure export, dataset generator) moved to
  `craic.gui.dialogs`.
- The two copies of the Save / Discard / Cancel prompt, with different wording,
  became one. Quitting is now just another caller of `_confirm_discard`.
- `app.py` is 2,050 lines from 2,296, with 269 lines of dialogs, 172 of document
  and 79 of track definitions extracted from it.

### The one behaviour change

`pair_posteriors` now releases the GIL, as `gotoh_align` and `mea_align` already
did. It is the dominant cost of every reliability run, every perturbation
replicate and the posterior explorer, so it was the one kernel that could freeze
the interface while it ran. Measured: during a 566 ms run another Python thread
received 178 time slices, where before it would have received none.

### Still deferred

The FFI still crosses as Python lists rather than buffers, which discards much of
the Rust speedup at the boundary; that needs a NumPy dependency in the Rust crate
and is optimisation rather than a source of bugs. `figures.py` (609 lines) and
`panels.py` (500) have not been looked at. `sandbox.Alternative.score` is still
raw SP identity.

## 0.4.0 — nothing gets lost

Quitting CRAIC with unsaved hand-curation discarded it silently: there was no
dirty flag, no close prompt, and no record of which file the alignment came from.
For a tool whose whole purpose is human judgement applied to hard regions, that
was the worst possible failure mode — the one thing in the application that
cannot be recomputed was the one thing not protected.

Fixed in three layers, deliberately not by auto-saving to a standard format. A
session is more than its alignment; writing NEXUS or MEGA on the way out would
have preserved the residues and quietly dropped the reference alignment, the
annotations, the groups, the provenance log and the view state — the same class
of silent loss, wearing a respectable format.

### Sessions

`craic.session` is a small, Qt-free document holding the whole working state:
alignment (with coding annotation and provenance), the truth-mode reference and
its label, annotations, sequence groups, the source path, and view state. Plain
readable JSON with a version field, so a file from a newer CRAIC is refused
rather than misread, and a JSON file that is not a session is refused too.

- **File → Save session** (`Ctrl+S`), **Save session as…**, **Open session…**.
  `Ctrl+S` is the session because it is the operation that loses nothing;
  exporting the alignment alone moved to `Ctrl+Shift+S` as **Save alignment as…**.
- A session opens from the ordinary Open dialog as well — a `.craic.json` is
  recognised by content, not by which menu item was used.
- Saving is atomic, via a temporary file and a rename. The autosave rewrites this
  file constantly, and a crash mid-write is precisely when it matters.
- Metadata that will not serialise (NumPy score arrays and the like) is dropped
  rather than treated as fatal: it is all derived and recomputable, and refusing
  to save because of it would turn a cosmetic problem into lost work.
- `craic session FILE` describes one from the command line — sequences, columns,
  reference, groups, annotations and the full provenance log — and `-o` exports
  its alignment to any supported format.

### Unsaved changes

- The document is marked dirty on every edit and clean on save or open, shown as
  a `•` in the title bar alongside the file name.
- Closing with unsaved work offers **Save / Discard / Cancel**, and says plainly
  what each keeps: the session keeps everything, a standard format keeps the
  residues.
- The window now remembers the file it opened, so Save has somewhere to go.
- Saving the alignment to a standard format marks the document clean — the
  residues, groups and annotations are on disk — but if a truth-mode reference is
  loaded it says so, because that is the part a standard format cannot carry.

### Every way of replacing the document, not just quitting

Adding a close prompt fixed one exit and left three others open. Opening another
file, opening a session and generating a dataset each replace the document just
as completely, and each did it without asking — so generating a teaching dataset
straight after aligning discarded the alignment silently. All four now go through
one confirmation.

Realigning is also undoable now. It throws away every hand edit made since the
last alignment, which is as destructive as any single edit, and was the only
destructive action in the application with no way back: `_commit` pushes onto the
undo stack the way `_apply_edit` always did.

### Crash recovery

- The session is autosaved to the application's own state directory, ~4 seconds
  after a change, coalescing bursts of edits into one write. Never beside the
  user's data: opening someone's reference alignment and looking at it must not
  leave files in their directory, and there is a test asserting exactly that.
- On the next launch with no file named on the command line, CRAIC offers to
  recover it, describing what it holds and when it was last touched. Recovered
  work stays marked unsaved, because it is.
- A clean exit and an explicit save both remove the recovery file, so its
  presence always means something actually went wrong. A corrupt one is discarded
  without a dialog.
- Restoring a session does not resume an expensive track: reopening your work
  should not silently start a minutes-long reliability analysis.

### Tests

`tests/test_v04_sessions.py` — 25 tests: session round-trip including coding
annotation and provenance, format-version and foreign-file refusal, atomic
writes, dirty tracking and the title bar, all three close-prompt answers, the
clean-close path, autosave location and content, recovery accepted, declined,
absent and corrupt, a full window round-trip through a saved session, a confirmation before every
document-replacing action, and undoing a realignment.

## 0.3.0 — teaching with a known answer, a headless command line, and a sparse transform

Three additions, chosen because each removes a limit on who can use CRAIC and for
what, rather than adding a feature to an already long list.

### The command line

`craic` and `craic FILE` still open the workbench. Everything else is new and runs
headless, never importing Qt:

    craic simulate --taxa 12 --divergence 0.4 -o truth.fasta --unaligned seqs.fasta
    craic align seqs.fasta -o aln.fasta
    craic score aln.fasta --reference truth.fasta -o scores.tsv
    craic mask  aln.fasta --threshold 0.5 -o masked.fasta
    craic trim  aln.fasta --method gappyout -o trimmed.fasta

`score` writes per-column consistency, perturbation, combined reliability, gap
score and similarity as a table, and with `--reference` adds per-column
correctness and reports SP, TC and the reliability AUC. Output goes to stdout when
`-o` is omitted and progress to stderr, so the subcommands pipe. `craic --core`
reports the active backend. A test blocks the PySide6 import and asserts the
headless path still completes, because "does not need a display" is the whole
point and is easy to break by accident.

### Teaching with a known answer

- **Teach → Generate dataset with known answer…** simulates sequences down a
  random tree (taxa, root length, divergence, indel rate, among-site rate
  variation, seed), loads them *unaligned*, and keeps the true alignment.
- **Teach → Load reference alignment…** does the same for real data with a
  structural reference such as BAliBASE.
- A **Reference correctness** track shows which columns are actually right, beside
  the reliability track that tried to predict it; the status bar reports SP, TC
  and the reliability AUC.
- The reference deliberately survives a realignment, and the accuracy is
  recomputed: watching correctness change as you realign is the point. It is
  dropped, quietly, only when the document becomes a different set of sequences.
- A reference that does not describe the same sequences is refused with an
  explanation rather than scored — otherwise a reference for the wrong data reads
  as an alignment that is entirely wrong.
- **Cheap tracks stay live through an edit.** Every alignment change used to reset
  the track selector to "(none)", so moving a block cleared the overlay you were
  moving it for. Conservation, the gap and similarity scores, the rule-based
  trimmers and Reference correctness are now kept and recomputed after each edit,
  and the status bar reports SP, TC, precision and the number of wrong columns as
  they change. The reliability family and the aligner-agreement map still reset,
  because re-running them on every nudge would make editing unusable.

### Sparse consistency transformation

The ProbCons consistency transformation was dense, which made it the memory and
time ceiling on family size. With SciPy installed (`pip install craic-msa[speed]`)
posteriors are thresholded at 0.01 and the products run sparse, as ProbCons does.

- Memory drops to a few per cent of dense (1.8-6% measured on simulated DNA), and
  the transform re-sparsifies each iteration, without which the fill-in of a sum
  of sparse products makes the second pass dense again.
- Speed depends on size: at ~260 columns dense wins by about 2x, at ~600 columns
  sparse wins by about 4x. The choice is therefore made by size, not merely by
  availability — below `SPARSE_MIN_LEN` (350) the dense product is used.
- This is an approximation, not a refactor. On simulated data the resulting
  alignments are usually identical and never materially less accurate (mean SP
  0.6976 dense vs 0.6957 sparse over four datasets), but they are not guaranteed
  bit-identical. The dense path is kept (`sparse=False`) and the two are compared
  in the test-suite rather than assumed equivalent.
- `_posterior_gb` now counts both generations of posteriors, fixing the roughly
  twofold understatement noted in 0.2.0, and accounts for sparse storage. The
  practical effect is that families previously skipped for memory — 30 sequences
  at 400 columns, for instance — can now use consistency.

### Structure

- Accuracy metrics and the ground-truth simulator moved out of `benchmarks/` and
  into the installed package as `craic.evaluate` and `craic.simulate`. The GUI,
  the command line and the benchmark harness now share one implementation, so a
  number shown to a student and a number printed in a paper cannot drift apart.
  `benchmarks/metrics.py` re-exports them and keeps only what is benchmark-specific.
- Sequence-id matching between two alignments moved from `ambiguity.disagreement`
  to `domain.match_ids`, since truth mode needs the same tolerance for aligners
  that truncate identifiers at whitespace.

### Tests

`tests/test_v03_features.py` — 19 tests covering the CLI end to end (as
subprocesses, including the no-Qt guarantee and that bare `craic` still launches
the workbench), reference scoring and its failure modes, sparse-versus-dense
agreement in both posteriors and final accuracy, the size-based auto choice, the
memory estimate, and the full teaching loop in the GUI.

## 0.2.0 — correctness and benchmark honesty

This release fixes defects that could produce a plausible-looking but wrong
scientific result, and replaces the benchmark's central claim with one that has
a control. Nothing here is a feature; everything is either a correction or a
disclosure.

### Results that could have been wrong

- **Reliability was scored under the wrong model.** `consistency()` used the
  library default emission model (90% identity) and default gap parameters
  rather than the parameters estimated from the alignment being scored, while
  `perturbation()` estimated its own — so the two halves of the combined score
  came from different models. Divergent families were systematically made to
  look unreliable, biasing masking against exactly the data masking decisions
  are about. Both now go through `progressive.estimate_params()`.

- **A failed ensemble masked the whole alignment.** A bare `except: pass` in the
  perturbation replicate loop meant that if every replicate failed, every score
  was `nan`, and `keep_mask` read `nan` as zero — indistinguishable from "every
  column is maximally unreliable". Total failure now raises `EnsembleFailure`;
  partial failure warns and the surviving count is recorded on the report.

- **The exported mask was not the evaluated mask.** The GUI kept unscoreable
  (nan) columns and the benchmark dropped them. There is now one
  `reliability.keep_mask()` used by both, with an explicit `unscored` policy
  (default: drop) and a refusal when no column could be scored at all.

- **`.` and `?` were counted as residues.** `domain._GAP` recognised them but
  every residue-index walk tested `ch != "-"`, so NEXUS and Stockholm input
  shifted every residue index and corrupted every score derived from it. Gap
  characters are now normalised once, in `Alignment.__post_init__`.

- **Truncated identifiers read as total disagreement.** MAFFT truncates FASTA
  ids at whitespace; an unmatched id incremented the pair denominator but not
  the numerator, so every column silently scored zero agreement. Ids are now
  mapped (exact, then whitespace-truncated) and an unmappable alignment is
  excluded with a warning or raises `IdMismatch`.

- **Out-of-range gap parameters produced NaN, silently.** `delta >= 0.5` makes
  `log(1 - 2*delta)` undefined; the Rust core turned that into NaN that
  propagated through every posterior while the NumPy path raised. Both backends
  now validate.

- **`outlier_sequences` flagged everything when the MAD was zero.** Near-identical
  sequences gave `mad = 0`, replaced by `1e-9`, putting the threshold immediately
  below the median. No robust spread now means no outliers.

- **`gappyout_mask` could trim the alignment away.** Cutting at the largest jump
  in the sorted gap-score distribution has no guard when that jump sits at the
  top: `[0.99]*100 + [1.0]` removed 100 of 101 columns. The cut is now bounded
  by `max_removed` (default 0.5).

- **The GUI could paint scores from a previous alignment.** Background
  reliability and agreement analyses were not cancelled or invalidated when the
  alignment changed, and the canvas did not check array shapes. There is now a
  document generation counter (`_stale_guard`) and a shape check at
  `Canvas.set_cell_scores`.

- **MUSCLE and PRANK were launched on the wrong PATH.** `available()` resolves
  against an augmented search path but those two engines called `subprocess.run`
  without it, so they reported available and then failed with `FileNotFoundError`
  in exactly the macOS `.app` case the augmented path exists to fix. All external
  invocations now go through `ExternalAligner._run`.

- **Splicing could re-phase the reading frame.** `sandbox.splice` preserved the
  coding annotation for any block width; a width change that is not a multiple of
  three now drops it rather than silently producing a wrong amino-acid view.

- Smaller guards: empty input and duplicate ids in `progressive.align`, an
  out-of-range residue index in `posterior.probe_residue`, and the shape contract
  of `posterior.ambiguity_profile` on an empty partner.

### Benchmark

- **The masking claim now has a null.** `masking_gain` reported mean accuracy
  over all columns against mean accuracy over retained columns. That inequality
  holds for *any* score correlated with column difficulty — gap fraction
  included — because masking removes hard columns, so it was circular. It now
  reports, at the same number of masked columns, a random-mask null (mean and
  s.d. over 200 draws), a gap-score control, and a z-score. The printed summary
  says so in as many words.

- **SP/TC scope is stated.** The BAliBASE XML reader no longer uppercases the
  sequence data, so core blocks marked by case are honoured and SP/TC can be
  restricted to them; each row carries `score_scope` (`core` or `all-columns`).
  All-column figures are not comparable with published BAliBASE tables, and the
  output says that. The `bali_score` manifest remains the citable route.

- **Every masking method gets an AUC.** `compare_masking` applied one 0.5
  threshold to a posterior, a gap fraction and a BLOSUM-normalised similarity —
  three incommensurable scales — and left the rule-based masks with `auc = nan`,
  so "best column-AUC" was structurally unwinnable by trimAl or Gblocks. Binary
  masks are now scored as the (coarse) classifiers they are, each method is cut
  at its own median, and every method carries a matched random baseline.

- **`--mask-engine`** lets the masking comparison run on a MAFFT or MUSCLE
  alignment instead of CRAIC's own, removing the home-field advantage.

- **Timing measures alignment only.** The clock previously ran across
  `evaluate()` — reliability analysis, NJ tree, RF distance — and reference mode
  recorded no timing at all.

- **Aligner versions are captured** and printed at the start of every run.

- **Skips and failures are reported.** Families skipped by the memory bound (a
  bound derived from CRAIC's dense-posterior cost, applied to every engine) and
  engine failures are now counted and summarised, so the means are not quietly
  taken over unequal subsets.

- **The alphabet is detected from the data** rather than defaulting to DNA, which
  had made it possible to score a protein benchmark under a DNA emission model
  with no warning. `--protein` / `--dna` still force it.

- **`auc` handles ties by midranks.** Reliability scores tie heavily at 1.0, and
  arbitrary tie ranks made the AUC depend on input order.

- **The simulator can break its own model.** `--rate-alpha` adds gamma among-site
  rate variation and `--indel-zipf` a heavy-tailed indel length distribution.
  The defaults are unchanged, and the module now states plainly that the default
  model is close to what CRAIC's pair-HMM assumes, which flatters CRAIC relative
  to MAFFT and MUSCLE.

### Tests

- `tests/test_benchmarks.py` (new): the benchmark layer had no tests at all, and
  it is what produces every number in the paper. SP/TC against hand-computed
  answers, core-mask restriction, AUC including hand-computed ties, the masking
  control (with a regression guard asserting that an *uninformative* score does
  not beat the random null), RF distance, and the simulator's internal
  consistency under all four rate/indel settings.
- `tests/test_v02_fixes.py` (new): one test per defect above, named for the
  failure rather than the function.

### Documentation

- README, `docs/concepts.md`, `docs/index.md` and `benchmarks/README.md` state
  what the perturbation score is and is not, that the consistency score is a
  goodness-of-fit under CRAIC's own model rather than independent evidence, that
  the posterior explorer's window is anchored on the committed alignment, that
  the consistency transform is dense rather than sparsified, that the protein
  emission model ignores the estimated divergence, and that the NumPy fallback's
  cross-validation necessarily skips in the configuration where it matters most.

### Known, deliberately not changed in 0.2.0

- `pair_posteriors` still does not release the GIL and does not honour
  cancellation in the Rust path, so the most expensive kernel is the one that
  can block the UI.
- The consistency transform is still dense; `_posterior_gb` still understates
  peak memory roughly twofold.
- FFI still crosses as Python lists, which discards much of the Rust speedup.
- `sandbox.Alternative.score` is still raw SP identity, which rewards
  over-alignment when ranking alternatives.
- `gui/app.py` is still a single large window class. The generation counter fixes
  the race it caused; the structural fix is a document model, for a later release.

## 0.1.0

Initial release.
