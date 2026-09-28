# FAQ & troubleshooting

## About the tool

**Should I use CRAIC's built-in aligner for real work?**

Generally no. It exists so that CRAIC works with nothing else installed, and so
that the perturbation analysis has an engine to re-align with. It is an honest
ProbCons-style implementation and on BAliBASE it is beaten by MAFFT — see
[Validation](validation.md) for the measured gap. Install MAFFT and select it
from the **Engine** menu; CRAIC finds anything on your `PATH` automatically.

**Then what is CRAIC actually for?**

Looking at an alignment's uncertainty, interrogating it, and curating on the
strength of it. Build the alignment with whatever engine you trust, and use
CRAIC to find out which parts of it you should not rely on. See
[Background & intent](background.md).

**Is the reliability score the same as GUIDANCE / TCS / ZORRO?**

No, and the docs try to be careful about this. The consistency score is in the
spirit of TCS and ZORRO — a posterior-based measure of how well the evidence
supports each column — but it is CRAIC's own implementation over its own
pair-HMM, so the numbers are not interchangeable with theirs. The perturbation
score is inspired by GUIDANCE but is explicitly *not* equivalent: GUIDANCE
bootstraps the guide tree and re-runs the aligner that built the alignment,
typically ~100 times; CRAIC runs 16 replicates of its own engine. It is a
cheaper, weaker relative. If you need GUIDANCE numbers, run GUIDANCE.

**Why does a column of identical residues sometimes score badly?**

Because conservation and reliability are different questions. Conservation asks
how similar the residues in a column are; reliability asks how confident we are
that they belong together. A run of identical residues inside a repeat can be
perfectly conserved and genuinely ambiguous, because it could have been placed
one or two columns over with equal support. Those are exactly the columns worth
knowing about.

**Can CRAIC infer a tree?**

No. It builds a neighbour-joining tree inside the benchmark harness as a
diagnostic — to ask whether masking improves the recovered topology — but tree
inference is not a feature and is not going to be. Export the alignment and use
IQ-TREE or RAxML.

## Using it

**How do I compare two alignments I made elsewhere?**

**File → Compare with alignment…** loads a second alignment of the same
sequences and overlays a per-column agreement track. If one of the two is
trusted — a structural reference, or the known answer from a simulation — use
**Teach → Load reference alignment…** instead, which scores rather than just
compares.

**The Correct placement colour and the Column inspector are greyed out.**

Both need a known answer. Either generate one (**Teach → Generate dataset with
known answer…**) or load a trusted alignment of the same sequences (**Teach →
Load reference alignment…**).

**My reference alignment is refused.**

CRAIC checks that the reference describes the *same residues* as the alignment,
matching on sequence identifier and tolerating the whitespace truncation some
aligners apply. A reference for different sequences is refused rather than
scored, because it would otherwise read as an alignment that is entirely wrong —
a much more misleading outcome than an error message. Check that the names match
and that neither file has been trimmed.

**Can I get the reliability numbers without the GUI?**

Yes. `craic score aln.fasta -o scores.tsv` writes per-column consistency,
perturbation, combined score, gap fraction and similarity; add
`--reference truth.fasta` for per-column correctness alongside. Nothing on the
headless path imports Qt, so it works over SSH and in a job script. See the
[command line reference](reference/command-line.md).

**Where did my work go when I closed the program?**

It should have asked. CRAIC keeps a session document (`*.craic.json`) holding the
alignment, the reference, annotations, groups and the provenance log — `Ctrl+S`
saves it — and autosaves for crash recovery. If you were caught by an older
version, see [Sessions & unsaved work](guide/sessions.md).

## Problems

**"Aligning…" never finishes.**

The built-in engine is pure Python over a Rust kernel and the consistency
transform is O(N³) in the number of sequences. A few hundred sequences will take
a very long time. Cancel, switch the **Engine** to MAFFT, and re-run. If you must
use the built-in engine on a large family, drop the compute level to *min*, which
skips the consistency transform and streams the posteriors.

**"No space" / the program is killed while aligning.**

Same cause. Holding every pairwise posterior for a large family is quadratic in
sequences and in length. The engine switches to streaming automatically past a
size threshold, but a family big enough will still exhaust memory. Use MAFFT.

**The Rust core won't build.**

CRAIC runs without it — `craic/accel.py` is a NumPy mirror of the same kernel and
the test suite cross-validates the two to ~1e-9, so results are identical and
only speed differs, but by tens of times: an analysis that takes minutes with the
core can take hours without it. The status bar says which is in use (`core: rust`
or `core: numpy`), and the progress window warns when the core is missing. If you want the Rust core, you need a Rust toolchain
(`rustup`) and then `./build_rust.sh`.

**`cargo` cannot write `Cargo.lock` (Permission denied) in a cloud-synced folder.**

Dropbox, OneDrive, iCloud Drive and Google Drive folders can leave files in a
state where the Unix permission bits say writable but the filesystem refuses the
write. Copy the file over itself to clear it:

```bash
find . -type f -exec sh -c 'cp -p "$1" "$1.tmp" && mv -f "$1.tmp" "$1"' _ {} \;
```

`chmod` alone will not help, because the mode bits are already correct. Better
still, keep the working copy outside the synced folder — `build_rust.sh` already
puts the virtualenv outside it for the same family of reasons.

**External aligners aren't detected.**

CRAIC looks for `mafft`, `muscle`, `clustalo` and `prank` on the `PATH` of the
process, which on macOS is the `PATH` a GUI application inherits from `launchd`
— not the one your shell has after `.zshrc` has run. Starting CRAIC from a
terminal (`./craic.command`) gets your shell's `PATH`. The status bar lists what
was found.

**My mitochondrial / Mycoplasma sequences are full of stop codons.**

The genetic code is wrong. Pick the right NCBI table from the **genetic code**
selector beside *treat as coding* (vertebrate mitochondria: 2; *Mycoplasma* and
*Spiroplasma*: 4; bacteria and plastids: 11), or pass `--code N` on the command
line — `craic codes` lists them. This is worth fixing rather than ignoring: a
stop codon is not one of the twenty residues, so it encodes to the wildcard and
contributes no homology signal, and TGA is the main tryptophan codon in exactly
those genomes.

If the code is already right and the stops persist, check the reading frame, or
consider whether you are looking at a pseudogene. CRAIC's warning cannot tell
the three apart.

**The amino-acid view is missing for my nucleotide data.**

Tick **treat as coding**. CRAIC assumes reading frame 0 and intact frames; if
the frames are broken it says so, and the codon and amino-acid views will show
`X` for partial-gap codons.

**Everything is very small / very large.**

`Ctrl/⌘` with `+`, `-` or `0`, or `Ctrl/⌘` and scroll over the alignment. The
residue font scales with the cells.
