"""Multi-aligner disagreement map.

Run several engines on the same sequences and measure, on a reference's
coordinates, how often the homologies each column asserts are reproduced by the
other alignments. Columns every method agrees on are trustworthy; columns they
fight over are exactly the ambiguous regions worth a human's attention.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import domain
from ..accel import Cancelled
from ..domain import Alignment, Alphabet, IdMismatch, match_ids


def membership_by_id(aln: Alignment) -> Dict[ResKey, int]:
    """``(sequence id, residue) -> column``: comparing two alignments has to
    match on the name, not the row position."""
    return domain.membership(aln, by="id")
from ..engines import AlignerEngine

Record = Tuple[str, str]
ResKey = Tuple[str, int]


@dataclass
class DisagreementResult:
    reference_label: str
    reference: Alignment
    alignments: List[Tuple[str, Alignment]]
    col_agreement: np.ndarray  # per reference column in [0, 1]; nan where undefined

    @property
    def labels(self) -> List[str]:
        return [lbl for lbl, _ in self.alignments]


def compare(reference: Alignment, other: Alignment) -> np.ndarray:
    """Per-column agreement of ``reference`` with another alignment of the same
    sequences: the fraction of each reference column's residue pairs that ``other``
    also places in one column. 1 = identical placement, low = the two disagree.

    Raises :class:`IdMismatch` if the two alignments do not describe the same
    sequences, rather than reporting universal disagreement.
    """
    inv = {oid: rid for rid, oid in match_ids(reference, other).items()}
    memb = {(inv[oid], r): c for (oid, r), c in membership_by_id(other).items() if oid in inv}
    ref_maps = {rid: domain.residue_index(row) for rid, row in zip(reference.ids, reference.rows)}
    L = reference.length
    agree = np.full(L, np.nan)
    for c in range(L):
        present = [(rid, ref_maps[rid][c]) for rid in reference.ids
                   if rid in ref_maps and ref_maps[rid][c] >= 0]
        if len(present) < 2:
            continue
        tot = ok = 0
        for a in range(len(present)):
            for b in range(a + 1, len(present)):
                tot += 1
                cx, cy = memb.get(present[a]), memb.get(present[b])
                if cx is not None and cx == cy:
                    ok += 1
        agree[c] = ok / tot if tot else np.nan
    return agree


def compute(
    records: Sequence[Record],
    alphabet: Alphabet,
    engines: Sequence[AlignerEngine],
    reference: int = 0,
    base_alignment: Optional[Alignment] = None,
    params_by_key: Optional[dict] = None,
    progress=None,
    cancelled=None,
) -> DisagreementResult:
    """Per-column agreement of several aligners.

    If ``base_alignment`` is given, agreement is measured on *its* columns — how
    often each engine reproduces the homologies the current alignment asserts —
    so the result annotates the alignment you already have rather than replacing
    it. Otherwise the first engine's alignment is used as the reference.
    ``progress(done, total, label)`` is called before each engine runs;
    ``cancelled()`` is checked between engines and handed to those that can stop
    part-way (the built-in one).
    """
    others: List[Tuple[str, Alignment]] = []
    engines = [e for e in engines if e.available()]
    for k, eng in enumerate(engines):
        if cancelled is not None and cancelled():
            raise Cancelled()
        if progress is not None:
            progress(k, len(engines), f"{eng.label} ({k + 1} of {len(engines)})")
        try:
            p = params_by_key.get(eng.key, {}) if params_by_key else {}
            others.append((eng.label, eng.align(records, alphabet, cancelled=cancelled, **p)))
        except Cancelled:
            raise
        except Exception:
            continue
    if progress is not None:
        progress(len(engines), len(engines), "")

    if base_alignment is not None:
        ref_label, ref = "current alignment", base_alignment
        against = others
    else:
        if not others:
            raise RuntimeError("no alignment engine succeeded")
        reference = min(reference, len(others) - 1)
        ref_label, ref = others[reference]
        against = [r for i, r in enumerate(others) if i != reference]

    all_runs = ([(ref_label, ref)] if base_alignment is not None else []) + others
    L = ref.length
    agree = np.full(L, np.nan)
    if not against:
        return DisagreementResult(ref_label, ref, all_runs, agree)

    # Re-key every comparison alignment onto the reference's own identifiers, so
    # that an aligner which truncated the ids does not read as total
    # disagreement. An alignment that cannot be mapped is dropped with a warning
    # rather than counted as disagreeing with everything.
    mems = []
    for label, a in against:
        try:
            inv = {oid: rid for rid, oid in match_ids(ref, a).items()}
        except IdMismatch as exc:
            warnings.warn(f"excluding {label} from the disagreement map: {exc}",
                          RuntimeWarning, stacklevel=2)
            continue
        mems.append({(inv[oid], r): c
                     for (oid, r), c in membership_by_id(a).items() if oid in inv})
    if not mems:
        return DisagreementResult(ref_label, ref, all_runs, agree)
    ref_maps = {rid: domain.residue_index(row) for rid, row in zip(ref.ids, ref.rows)}
    for c in range(L):
        present = [(rid, ref_maps[rid][c]) for rid in ref.ids if ref_maps[rid][c] >= 0]
        if len(present) < 2:
            continue
        tot = ok = 0
        for a in range(len(present)):
            for b in range(a + 1, len(present)):
                x, y = present[a], present[b]
                for mem in mems:
                    cx, cy = mem.get(x), mem.get(y)
                    tot += 1
                    if cx is not None and cx == cy:
                        ok += 1
        agree[c] = ok / tot if tot else np.nan

    return DisagreementResult(ref_label, ref, all_runs, agree)
