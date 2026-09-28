# Ambiguity tools

CRAIC's reason for existing: four complementary, evidence-based views of where an alignment is
trustworthy and where it is guessing.

## 1 · Reliability overlay + live masking

Every column and residue gets a confidence score from two independent signals: a **consistency**
score (the pair-HMM posterior that the residues a column groups together really are homologous)
and a **perturbation** score (how many of a residue's asserted homologies survive re-alignment
under perturbed conditions). Choose **Reliability** from the **Track** dropdown (it computes on
first use, with a progress window showing the time left; if the perturbation part will take
more than two minutes CRAIC asks first, and you can keep the consistency score alone), then
read it on the track, as
**Confidence** colouring, or as a live mask preview via the **Mask** slider.
To keep the columns and drop only the unreliable residues in them, use **Mask
residues below the threshold** (Edit ▸ Edit alignment) — see
[residue masks](figures-and-curation.md#masking-residues-instead-of-columns).

## 2 · Multi-aligner disagreement map

Choose **Aligner agreement** from the **Track** dropdown to run several engines (or gap regimes of the built-in
aligner) and colour each column by how often the methods agree. Where they fight is where you
should look — disagreement is the most honest signal of ambiguity.

## 3 · Local realignment sandbox

Select a hard column range, open the **Realignment sandbox**, and press **Generate
alternatives**. CRAIC re-aligns *just that block* under different engines and gap costs, scores
each, and lets you **preview** and **splice** your choice back in without disturbing the rest of
the alignment.

## 4 · Posterior explorer

Instead of one hard column assignment, see the pair-HMM posterior for a region. Because a
confident alignment is mostly a bright, boring diagonal, there are four ways to read it:

**Colour by confidence** — shade residues by reliability, in place.

![Confidence colouring: well-supported residues green, guesswork red](../assets/confidence-colouring.png){ width="100%" }

**Probe a residue** — single-click any residue to paint its *homology cloud*: the other columns
it could plausibly belong to. A tight cloud means solid homology; a spread-out one means it
could slide.

![Probing a residue lights up its candidate columns](../assets/probe.png){ width="100%" }

**Residual mode** — subtract the committed alignment so the diagonal goes dark and only the
genuine alternative-homology mass remains.

![Residual mode — the confident diagonal removed](../assets/residual.png){ width="100%" }

**Profile mode** — one sequence versus the rest, so the picture reflects the whole alignment
rather than an arbitrary pair.

![Profile mode — each row is a residue, smeared across the columns it could occupy](../assets/profile.png){ width="100%" }
