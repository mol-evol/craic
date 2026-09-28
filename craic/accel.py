"""Acceleration layer: one pair-HMM posterior primitive, Rust or NumPy.

Design notes
------------
* The compiled ``craic_accel`` extension (src/lib.rs) is used when present.
* A NumPy fallback mirrors it line-for-line so the workbench runs anywhere and
  so the two implementations can be cross-validated in the test-suite
  (Liskov: callers cannot tell which one answered).
* Emission models are built here, not in Rust, so the core stays
  biology-agnostic. DNA/RNA use an identity-driven model; protein uses
  BLOSUM45 converted to a proper joint probability via Robinson background
  frequencies.
"""

from __future__ import annotations

import os
import threading
from collections import deque
from concurrent.futures import Future, ThreadPoolExecutor
from dataclasses import dataclass
from math import inf, log, exp
from typing import Callable, Deque, Dict, Iterable, Iterator, List, Optional, Sequence

import numpy as np

try:  # pragma: no cover - exercised by whichever build is installed
    from . import _accel as _rust

    HAVE_RUST = True
except Exception:  # pragma: no cover
    _rust = None
    HAVE_RUST = False


class Cancelled(Exception):
    """Raised when a caller asks an in-progress computation to stop."""


# --------------------------------------------------------------------------- #
# Parallel map over independent pair computations
# --------------------------------------------------------------------------- #

#: Threads used for independent pair-HMM computations: one per core, less one
#: left for the rest of the computer, so a long analysis does not make it
#: sluggish. The Rust core releases the GIL, so the threads run in parallel; the
#: NumPy fallback holds the GIL and gains nothing, but loses nothing either.
WORKERS = max(1, (os.cpu_count() or 1) - 1)

#: Memory (GB) the threads may hold in posterior matrices at once. Streaming
#: exists to keep peak memory at one matrix; running it on many threads must not
#: quietly multiply that for very long sequences.
PARALLEL_MEM_GB = 1.0

_POOL: Optional[ThreadPoolExecutor] = None
_POOL_LOCK = threading.Lock()


def in_flight(max_len: int) -> int:
    """How many ``max_len`` x ``max_len`` posteriors :func:`imap` may hold at once
    within :data:`PARALLEL_MEM_GB`: the float64 matrix, plus the Python list it
    crosses the FFI as (about four times larger) while it is converted."""
    per = max(1, max_len) ** 2 * 8 * 5
    return max(1, min(2 * WORKERS, int(PARALLEL_MEM_GB * 1e9 // per)))


def _pool() -> ThreadPoolExecutor:
    global _POOL
    with _POOL_LOCK:
        if _POOL is None:
            _POOL = ThreadPoolExecutor(WORKERS, thread_name_prefix="craic")
        return _POOL


def imap(fn: Callable, items: Iterable, cancelled=None, ahead: Optional[int] = None
         ) -> Iterator:
    """``fn(item)`` for each item, computed on :data:`WORKERS` threads and
    yielded in input order.

    Results come back in order, so any sum over them is the same, bit for bit,
    as the serial one. At most ``ahead`` results (default ``2 * WORKERS``; see
    :func:`in_flight` for large matrices) are in flight at once, which bounds
    memory. ``fn`` must not call ``imap`` itself: the pool is shared, and nested
    waits could exhaust it. ``cancelled()`` is polled before each submission.
    """
    ahead = 2 * WORKERS if ahead is None else max(1, ahead)
    if WORKERS == 1 or ahead == 1:
        for x in items:
            if cancelled is not None and cancelled():
                raise Cancelled()
            yield fn(x)
        return
    pool = _pool()
    pending: Deque[Future] = deque()
    try:
        for x in items:
            if cancelled is not None and cancelled():
                raise Cancelled()
            pending.append(pool.submit(fn, x))
            if len(pending) >= ahead:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()
    finally:
        for f in pending:          # stopped early (cancel, error, or caller quit)
            f.cancel()


# --------------------------------------------------------------------------- #
# Emission models
# --------------------------------------------------------------------------- #

# Robinson & Robinson (1991) amino-acid background frequencies.
_ROBINSON = {
    "A": 0.07805, "R": 0.05129, "N": 0.04487, "D": 0.05364, "C": 0.01925,
    "Q": 0.04264, "E": 0.06295, "G": 0.07377, "H": 0.02199, "I": 0.05142,
    "L": 0.09019, "K": 0.05744, "M": 0.02243, "F": 0.03856, "P": 0.05203,
    "S": 0.07120, "T": 0.05841, "W": 0.01330, "Y": 0.03216, "V": 0.06441,
}
_AA20 = "ARNDCQEGHILKMFPSTWYV"


@dataclass(frozen=True)
class EmissionModel:
    """Symbols, their index, and the pair-HMM emission distributions.

    The final symbol is always a wildcard ('N' for nucleotides, 'X' for
    protein) whose emissions equal the background product, so unknown residues
    contribute no homology signal rather than a spurious one.
    """

    kind: str
    symbols: str
    joint: np.ndarray  # (k, k) joint match-emission probabilities, sums to 1
    bg: np.ndarray     # (k,) gap-emission (background) probabilities
    index: Dict[str, int]

    @property
    def k(self) -> int:
        return len(self.symbols)

    def encode(self, seq: str) -> List[int]:
        wild = self.k - 1
        idx = self.index
        return [idx.get(c, wild) for c in seq.upper()]


def _finalize(kind: str, syms: str, joint: np.ndarray, bg: np.ndarray, wild: str) -> EmissionModel:
    """Append a neutral wildcard symbol and normalise."""
    k = len(syms)
    bgm = float(bg.mean())
    big = np.zeros((k + 1, k + 1))
    big[:k, :k] = joint
    # wildcard emits like the background: joint(wild, y) = bg[y] * mean_bg
    big[k, :k] = bg * bgm
    big[:k, k] = bg * bgm
    big[k, k] = bgm * bgm
    big /= big.sum()
    bg2 = np.append(bg, bgm)
    bg2 = bg2 / bg2.sum()
    symbols = syms + wild
    return EmissionModel(kind, symbols, big, bg2, {c: i for i, c in enumerate(symbols)})


def _dna_model(identity: float = 0.9) -> EmissionModel:
    syms = "ACGT"
    bg = np.full(4, 0.25)
    off = (1.0 - identity) / 12.0
    joint = np.full((4, 4), off)
    np.fill_diagonal(joint, identity / 4.0)
    joint /= joint.sum()
    model = _finalize("dna", syms, joint, bg, "N")
    # RNA is modelled on the DNA emission model, so uracil is thymine. Without
    # this alias, encode() would send every 'U' to the wildcard and silently
    # strip the homology signal from RNA data.
    model.index["U"] = model.index["T"]
    return model


#: Bit-scale of each BLOSUM matrix as NCBI distributes it. The log-odds in the
#: file are ``s = scale * log2(q_ij / p_i p_j)``, so recovering the joint needs
#: the right divisor: 45, 50 and 80 are in third-bits, 62 and 90 in half-bits —
#: an inconsistency in how NCBI ships them, not a pattern to infer from. Getting
#: it wrong yields a plausible-looking but wrong emission model rather than an
#: error, because ``_protein_model`` normalises the damage away, so every entry
#: is checked against the matrix itself in the test suite.
_BLOSUM_SCALE = {"BLOSUM45": 3.0, "BLOSUM50": 3.0, "BLOSUM62": 2.0,
                 "BLOSUM80": 3.0, "BLOSUM90": 2.0}

#: The shipped protein matrix. BLOSUM45 rather than the more usual BLOSUM62:
#: measured on BAliBASE core blocks with the official scorer, it is better at
#: both divergence levels tested — +0.026 SP on RV11 (<20% identity, n=38,
#: Wilcoxon p=0.005) and +0.007 on RV12 (20–40%, n=44, p=0.047), with BLOSUM80
#: worse than both. One matrix is still used for all protein data; keying the
#: choice to estimated divergence is future work.
PROTEIN_MATRIX = "BLOSUM45"

#: Every protein matrix the emission model can be built from.
PROTEIN_MATRICES = tuple(_BLOSUM_SCALE)


def _protein_model(matrix: str = PROTEIN_MATRIX) -> EmissionModel:
    from Bio.Align import substitution_matrices

    try:
        divisor = _BLOSUM_SCALE[matrix]
    except KeyError:
        raise ValueError(f"no bit-scale known for {matrix!r}; add it to "
                         f"_BLOSUM_SCALE with the scale NCBI distributes it in")
    bl = substitution_matrices.load(matrix)
    bg = np.array([_ROBINSON[a] for a in _AA20])
    joint = np.zeros((20, 20))
    for i, a in enumerate(_AA20):
        for j, b in enumerate(_AA20):
            s = float(bl[a, b])  # log-odds, in units of 1/divisor bits
            joint[i, j] = bg[i] * bg[j] * (2.0 ** (s / divisor))
    joint /= joint.sum()
    return _finalize("protein", _AA20, joint, bg, "X")


_MODEL_CACHE: Dict[str, EmissionModel] = {}


def emission_model(kind: str, identity: float = 0.9,
                   matrix: str = PROTEIN_MATRIX) -> EmissionModel:
    """Return a cached emission model. ``kind`` in {dna, rna, protein}.

    ``matrix`` names the protein substitution matrix and is ignored for DNA.
    """
    kind = kind.lower()
    if kind == "rna":
        kind = "dna"
    key = f"{kind}:{identity}:{matrix}" if kind == "protein" else f"{kind}:{identity}"
    if key not in _MODEL_CACHE:
        _MODEL_CACHE[key] = _protein_model(matrix) if kind == "protein" else _dna_model(identity)
    return _MODEL_CACHE[key]


# --------------------------------------------------------------------------- #
# Posterior decoding
# --------------------------------------------------------------------------- #

def pair_posteriors(
    a_idx: Sequence[int],
    b_idx: Sequence[int],
    k: int,
    joint_flat: Sequence[float],
    bg: Sequence[float],
    delta: float = 0.02,
    epsilon: float = 0.5,
) -> np.ndarray:
    """Dispatch to Rust if available, else the NumPy mirror. Returns (m, n).

    ``delta`` and ``epsilon`` are validated here, in the one place both backends
    go through. Out-of-range values are not a curiosity: ``delta >= 0.5`` makes
    ``log(1 - 2*delta)`` undefined, which the Rust core turns into NaN that
    propagates silently through every posterior and out into every reliability
    score, while the NumPy path raises. A silently wrong number is worse than a
    loud failure, so both now fail loudly.
    """
    _check_gap_params(delta, epsilon)
    if HAVE_RUST:
        post, m, n = _rust.pair_posteriors(
            list(a_idx), list(b_idx), k, list(joint_flat), list(bg), delta, epsilon
        )
        return np.asarray(post, dtype=float).reshape(m, n)
    return _pp_numpy(a_idx, b_idx, k, joint_flat, bg, delta, epsilon)


def _check_gap_params(delta: float, epsilon: float) -> None:
    """Reject pair-HMM gap parameters that do not describe a valid model.

    The three-state chain leaves M with probability ``2*delta`` and leaves a gap
    state with probability ``1 - epsilon``, so ``0 < delta < 0.5`` and
    ``0 <= epsilon < 1`` are the ranges in which the transitions are
    probabilities at all.
    """
    if not (0.0 < delta < 0.5):
        raise ValueError(
            f"delta (gap-open probability) must be in (0, 0.5); got {delta!r}. "
            "At 0.5 or above the M->M transition probability 1 - 2*delta is not positive."
        )
    if not (0.0 <= epsilon < 1.0):
        raise ValueError(
            f"epsilon (gap-extend probability) must be in [0, 1); got {epsilon!r}. "
            "At 1 the chain can never return to the match state."
        )


def posterior_matrix(
    seq_a: str,
    seq_b: str,
    model: EmissionModel,
    delta: float = 0.02,
    epsilon: float = 0.5,
) -> np.ndarray:
    """High-level helper: ungapped strings in, (len_a, len_b) posteriors out."""
    a = model.encode(seq_a.replace("-", ""))
    b = model.encode(seq_b.replace("-", ""))
    return pair_posteriors(a, b, model.k, model.joint.ravel(), model.bg, delta, epsilon)


def gotoh_align(S: np.ndarray, gap_open: float, gap_extend: float, cancelled=None):
    """Affine (Gotoh) profile alignment over the column score matrix ``S`` (la x lb).

    Returns ``(cols_a, cols_b)``: the aligned column index for each output position,
    or -1 for a gap. The Rust core runs the DP with the GIL released; the NumPy
    fallback below mirrors it exactly (cross-validated in the test-suite). The
    NumPy path polls ``cancelled()`` between rows so a slow run can be stopped.
    """
    S = np.ascontiguousarray(S, dtype=float)
    la, lb = S.shape
    # hasattr guard: an older compiled _accel may predate this function, so fall
    # back to NumPy rather than crashing if the build is stale.
    if HAVE_RUST and hasattr(_rust, "gotoh_align"):
        ca, cb = _rust.gotoh_align(S.ravel().tolist(), int(la), int(lb),
                                   float(gap_open), float(gap_extend))
        return list(ca), list(cb)
    return _gotoh_numpy(S, gap_open, gap_extend, cancelled)


def _gotoh_numpy(S: np.ndarray, gap_open: float, gap_extend: float, cancelled=None):
    la, lb = S.shape
    NEG = -1e18
    M = np.full((la + 1, lb + 1), NEG)
    Ix = np.full((la + 1, lb + 1), NEG)   # gap in B (consume A)
    Iy = np.full((la + 1, lb + 1), NEG)   # gap in A (consume B)
    M[0, 0] = 0.0
    for i in range(1, la + 1):
        Ix[i, 0] = gap_open + (i - 1) * gap_extend
    for j in range(1, lb + 1):
        Iy[0, j] = gap_open + (j - 1) * gap_extend
    tbm = np.zeros((la + 1, lb + 1), dtype=np.int8)
    tbx = np.zeros((la + 1, lb + 1), dtype=np.int8)
    tby = np.zeros((la + 1, lb + 1), dtype=np.int8)
    for i in range(1, la + 1):
        if cancelled is not None and cancelled():   # stop mid-run, not just between merges
            raise Cancelled()
        for j in range(1, lb + 1):
            s = S[i - 1, j - 1]
            cand = (M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
            kb = int(np.argmax(cand)); M[i, j] = cand[kb] + s; tbm[i, j] = kb
            cand = (M[i - 1, j] + gap_open, Ix[i - 1, j] + gap_extend, Iy[i - 1, j] + gap_open)
            kb = int(np.argmax(cand)); Ix[i, j] = cand[kb]; tbx[i, j] = kb
            cand = (M[i, j - 1] + gap_open, Ix[i, j - 1] + gap_open, Iy[i, j - 1] + gap_extend)
            kb = int(np.argmax(cand)); Iy[i, j] = cand[kb]; tby[i, j] = kb
    state = int(np.argmax((M[la, lb], Ix[la, lb], Iy[la, lb])))
    i, j = la, lb
    ca, cb = [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and state == 0:
            ca.append(i - 1); cb.append(j - 1); state = int(tbm[i, j]); i -= 1; j -= 1
        elif i > 0 and (state == 1 or j == 0):
            ca.append(i - 1); cb.append(-1); state = int(tbx[i, j]); i -= 1
        else:
            ca.append(-1); cb.append(j - 1); state = int(tby[i, j]); j -= 1
    ca.reverse(); cb.reverse()
    return ca, cb


def mea_align(S: np.ndarray, cancelled=None):
    """Pure maximum-expected-accuracy (gap-free) alignment of the posterior score
    matrix ``S`` (la x lb). Returns ``(cols_a, cols_b)``: the aligned column index
    for each output position, or -1 for a gap. There is no gap penalty -- the
    alignment maximises the expected number of correctly-aligned residue pairs
    (ProbCons decoding). The Rust core runs the DP with the GIL released; the
    NumPy mirror below is identical (cross-validated in the test-suite) and can be
    cancelled between rows.
    """
    S = np.ascontiguousarray(S, dtype=float)
    la, lb = S.shape
    # hasattr guard: an older compiled _accel may predate this function.
    if HAVE_RUST and hasattr(_rust, "mea_align"):
        ca, cb = _rust.mea_align(S.ravel().tolist(), int(la), int(lb))
        return list(ca), list(cb)
    return _mea_numpy(S, cancelled)


def _mea_numpy(S: np.ndarray, cancelled=None):
    la, lb = S.shape
    M = np.zeros((la + 1, lb + 1))
    tb = np.zeros((la + 1, lb + 1), dtype=np.int8)   # 0 diag, 1 up (gap in B), 2 left (gap in A)
    for i in range(1, la + 1):
        if cancelled is not None and cancelled():
            raise Cancelled()
        Mi, Mi1, Si = M[i], M[i - 1], S[i - 1]
        for j in range(1, lb + 1):
            diag = Mi1[j - 1] + Si[j - 1]
            up = Mi1[j]
            left = Mi[j - 1]
            if diag >= up and diag >= left:
                Mi[j] = diag; tb[i, j] = 0
            elif up >= left:
                Mi[j] = up; tb[i, j] = 1
            else:
                Mi[j] = left; tb[i, j] = 2
    i, j = la, lb
    ca, cb = [], []
    while i > 0 or j > 0:
        t = int(tb[i, j]) if (i > 0 and j > 0) else (1 if j == 0 else 2)
        if t == 0:
            ca.append(i - 1); cb.append(j - 1); i -= 1; j -= 1
        elif t == 1:
            ca.append(i - 1); cb.append(-1); i -= 1
        else:
            ca.append(-1); cb.append(j - 1); j -= 1
    ca.reverse(); cb.reverse()
    return ca, cb


def _pp_numpy(a_idx, b_idx, k, joint_flat, bg, delta, epsilon) -> np.ndarray:
    a = list(a_idx)
    b = list(b_idx)
    m, n = len(a), len(b)
    joint = np.asarray(joint_flat, dtype=float).reshape(k, k)
    bgv = np.asarray(bg, dtype=float)
    ljoint = np.log(np.maximum(joint, 1e-300))
    lbg = np.log(np.maximum(bgv, 1e-300))
    ld, le = log(delta), log(epsilon)
    lmm = log(1.0 - 2.0 * delta)
    lgm = log(1.0 - epsilon)

    NEG = -inf

    def lse(*xs):
        mx = max(xs)
        if mx == NEG:
            return NEG
        return mx + log(sum(exp(x - mx) for x in xs))

    fm = [[NEG] * (n + 1) for _ in range(m + 1)]
    fi = [[NEG] * (n + 1) for _ in range(m + 1)]
    fd = [[NEG] * (n + 1) for _ in range(m + 1)]
    fm[0][0] = 0.0
    for i in range(m + 1):
        for j in range(n + 1):
            if i == 0 and j == 0:
                continue
            if i > 0 and j > 0:
                e = ljoint[a[i - 1], b[j - 1]]
                fm[i][j] = e + lse(fm[i - 1][j - 1] + lmm, fi[i - 1][j - 1] + lgm, fd[i - 1][j - 1] + lgm)
            if i > 0:
                fi[i][j] = lbg[a[i - 1]] + lse(fm[i - 1][j] + ld, fi[i - 1][j] + le)
            if j > 0:
                fd[i][j] = lbg[b[j - 1]] + lse(fm[i][j - 1] + ld, fd[i][j - 1] + le)

    z = lse(fm[m][n], fi[m][n], fd[m][n])

    bm = [[NEG] * (n + 1) for _ in range(m + 1)]
    bi = [[NEG] * (n + 1) for _ in range(m + 1)]
    bd = [[NEG] * (n + 1) for _ in range(m + 1)]
    bm[m][n] = bi[m][n] = bd[m][n] = 0.0
    for i in range(m, -1, -1):
        for j in range(n, -1, -1):
            if i == m and j == n:
                continue
            sm = si = sd = NEG
            if i < m and j < n:
                t = bm[i + 1][j + 1] + ljoint[a[i], b[j]]
                sm = lse(sm, t + lmm)
                si = lse(si, t + lgm)
                sd = lse(sd, t + lgm)
            if i < m:
                t = bi[i + 1][j] + lbg[a[i]]
                sm = lse(sm, t + ld)
                si = lse(si, t + le)
            if j < n:
                t = bd[i][j + 1] + lbg[b[j]]
                sm = lse(sm, t + ld)
                sd = lse(sd, t + le)
            bm[i][j] = sm
            bi[i][j] = si
            bd[i][j] = sd

    post = np.zeros((m, n))
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            post[i - 1, j - 1] = min(1.0, exp(fm[i][j] + bm[i][j] - z))
    return post


def backend() -> str:
    return "rust" if HAVE_RUST else "numpy"
