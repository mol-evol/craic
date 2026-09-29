"""CRAIC — Conserved-Region Alignment by Iterative Convergence.

A multiple sequence alignment workbench built for ambiguous regions.

The package is intentionally light to import: the GUI and heavy engines are
imported lazily by their own modules so that headless / scripted use never
pays for PySide6.
"""

__version__ = "0.5.14"
__author__ = "James McInerney"

AUTHOR_URL = "https://mol-evol.github.io/"
WEBSITE = "https://mol-evol.github.io/craic/"
ABOUT_URL = WEBSITE + "about/"

# How to cite CRAIC. While the paper is under review this is the bioRxiv
# preprint; when it is published, replace it with the journal citation. The
# same text is in CITATION.cff, README.md and docs/about.md — change all four
# (docs/editing-docs.md has the checklist).
CITATION = ("McInerney J. 2026. CRAIC: interactive interrogation of multiple "
            "sequence alignment uncertainty. bioRxiv preprint (link to follow).")
