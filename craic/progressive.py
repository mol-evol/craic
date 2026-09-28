"""A small, self-contained progressive aligner (the always-available engine).

ProbCons-style, built entirely on CRAIC's pair-HMM posteriors:

1. estimate the emission divergence and gap parameters from the data (no
   hard-coded assumption that sequences are ~90% identical),
2. compute all pairwise posterior "match" matrices,
3. apply the ProbCons **consistency transformation** (re-estimate each pairwise
   posterior using every third sequence), and
4. align progressively by **pure maximum-expected-accuracy** decoding of the
   transformed posteriors -- gap placement falls out of the model, so there are
   no ad-hoc gap constants.

It is still deliberately modest (no sequence weighting, no iterative refinement)
and is meant to always be available, not to beat MAFFT; but it is genuine
posterior decoding rather than a heuristic. The pair-HMM inner loop runs in the
Rust core when built, with the NumPy fallback otherwise.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import accel
from .accel import Cancelled, EmissionModel, emission_model
from . import domain
from .domain import Alignment, Alphabet

Record = Tuple[str, str]

# Fallback pair-HMM gap probabilities when not estimated from the data.
_DEFAULT_DELTA = 0.02
_DEFAULT_EPSILON = 0.5


# --------------------------------------------------------------------------- #
# Guide tree
# --------------------------------------------------------------------------- #

def _kmer_distance(seqs: Sequence[str], k: int = 3, rng=None) -> np.ndarray:
    """Cosine distance on k-mer counts. With ``rng`` set, the k-mer feature
    columns are bootstrap-resampled (and close distances lightly jittered) to
    sample an alternative guide tree -- the perturbation the reliability score
    re-aligns under (GUIDANCE-style guide-tree uncertainty)."""
    vocab: dict = {}
    counts = []
    for s in seqs:
        s = s.replace("-", "").upper()
        d: dict = {}
        for i in range(len(s) - k + 1):
            kmer = s[i : i + k]
            d[kmer] = d.get(kmer, 0) + 1
            vocab.setdefault(kmer, len(vocab))
        counts.append(d)
    n = len(seqs)
    V = len(vocab) or 1
    mat = np.zeros((n, V))
    for r, d in enumerate(counts):
        for kmer, c in d.items():
            mat[r, vocab[kmer]] = c
    if rng is not None and V > 1:
        mat = mat[:, rng.integers(0, V, size=V)]        # nonparametric bootstrap of features
    norms = np.linalg.norm(mat, axis=1)
    norms[norms == 0] = 1.0
    mat = mat / norms[:, None]
    D = 1.0 - mat @ mat.T
    if rng is not None:
        scale = 0.15 * float(np.mean(D[D > 0])) if np.any(D > 0) else 0.0
        noise = rng.normal(0.0, scale, size=D.shape)
        noise = np.triu(noise, 1)
        D = np.clip(D + noise + noise.T, 0.0, None)
    np.fill_diagonal(D, 0.0)
    return D


# --------------------------------------------------------------------------- #
# Posterior decoding pieces
# --------------------------------------------------------------------------- #

def _all_pairs_posteriors(seqs, model, delta, epsilon,
                          cancelled=None) -> Dict[Tuple[int, int], np.ndarray]:
    n = len(seqs)
    pairs = [(x, y) for x in range(n) for y in range(x + 1, n)]
    posts = accel.imap(lambda p: accel.posterior_matrix(seqs[p[0]], seqs[p[1]], model,
                                                        delta, epsilon),
                       pairs, cancelled, accel.in_flight(max(map(len, seqs))))
    return dict(zip(pairs, posts))


def _post(P, x, y):
    """The posterior for the ordered pair (x, y), always as a dense array.

    The store may hold sparse matrices (see :func:`consistency_transform`);
    profile merging indexes with ``np.ix_`` and so needs dense, and densifying
    one pair at a time keeps peak memory at a single matrix rather than the
    whole set.
    """
    M = P[(x, y)] if x < y else P[(y, x)].T
    return M.toarray() if hasattr(M, "toarray") else M


#: Posterior probabilities below this are dropped when sparsifying. ProbCons
#: uses the same idea and a similar cut-off: the consistency transformation is
#: dominated by a great many near-zero entries that cannot change a decoding.
SPARSE_THRESHOLD = 0.01

def sparse_available() -> bool:
    """Whether SciPy is present, which is what the sparse transform needs."""
    try:
        import scipy.sparse  # noqa: F401
    except Exception:
        return False
    return True


def _sparsify(M, threshold: float):
    from scipy.sparse import csr_matrix

    keep = M >= threshold
    if not keep.any():
        return csr_matrix(M.shape, dtype=float)
    out = np.where(keep, M, 0.0)
    return csr_matrix(out)


def consistency_transform(P, n, iters=2, cancelled=None,
                          sparse: bool = False,
                          threshold: float = SPARSE_THRESHOLD):
    """ProbCons consistency transformation: re-estimate each pairwise posterior
    P(x,y) using every third sequence z, P'(x,y) = (1/n) sum_z P(x,z) P(z,y).

    This is the most expensive thing CRAIC does — ``O(n^3)`` matrix products of
    ``L x L`` matrices per iteration — and it is what sets the largest family the
    workbench can handle with consistency turned on. Posteriors are overwhelmingly
    near-zero, so ``sparse=True`` thresholds the matrices at ``threshold`` and
    runs the products sparse (SciPy), which is how ProbCons itself does it and
    holds a few per cent of the dense memory.

    Dropping entries below ``threshold`` is an approximation, and not a small
    one on divergent sequences, whose posteriors are spread thinly: on simulated
    DNA at 0.5 substitutions per site it cost the built-in engine about 0.1 of
    sum-of-pairs accuracy. So :func:`align` uses it only when the exact (dense)
    transform would not fit in memory.
    """
    if sparse and not sparse_available():
        raise RuntimeError("the sparse consistency transform needs SciPy "
                           "(pip install scipy), or pass sparse=False")
    if sparse:
        P = {k: _sparsify(v, threshold) for k, v in P.items()}

    for _ in range(max(0, iters)):
        new = {}
        for x in range(n):
            for y in range(x + 1, n):
                if cancelled is not None and cancelled():
                    raise Cancelled()
                acc = 2.0 * P[(x, y)]                    # z = x and z = y (identity) terms
                for z in range(n):
                    if z == x or z == y:
                        continue
                    a = P[(x, z)] if x < z else P[(z, x)].T
                    b = P[(z, y)] if z < y else P[(y, z)].T
                    acc = acc + a @ b
                acc = acc / float(n)
                # Re-sparsify each iteration: a sum of sparse products fills in,
                # so without this the second pass is effectively dense again and
                # the memory saving evaporates. ProbCons does the same.
                new[(x, y)] = _sparsify(acc.toarray(), threshold) if sparse else acc
        P = new
    return P


def _merge(rows_a, rows_b, post, cancelled=None, parallel=False, tick=None):
    """Merge two profiles by pure-MEA decoding of their summed posterior scores.

    ``post(x, y)`` returns the pairwise posterior for original sequences x, y --
    either a lookup into stored (consistency-transformed) posteriors, or an
    on-demand computation (streaming). Each profile is a list of
    ``(orig_index, gapped_row)``; returns the merged list. ``parallel`` computes
    the posteriors on several threads (worth it only when ``post`` computes
    rather than looks up); the sum is taken in the same order either way, so
    the result is identical. ``tick()`` is called once per sequence pair."""
    na, nb = len(rows_a), len(rows_b)
    La, Lb = len(rows_a[0][1]), len(rows_b[0][1])

    def cols(rows):
        out = []
        for _, r in rows:
            m = domain.residue_index_array(r)
            c = np.where(m >= 0)[0]
            out.append((c, m[c]))
        return out

    ca_, cb_ = cols(rows_a), cols(rows_b)
    pairs = [(xi, yi) for xi in range(na) for yi in range(nb)]

    def one(p):
        return post(rows_a[p[0]][0], rows_b[p[1]][0])

    longest = max(max((len(c) for c, _ in ca_), default=0), max((len(c) for c, _ in cb_), default=0))
    posts = (accel.imap(one, pairs, cancelled, accel.in_flight(longest)) if parallel
             else map(one, pairs))
    S = np.zeros((La, Lb))
    for (xi, yi), Pm in zip(pairs, posts):
        if cancelled is not None and cancelled():
            raise Cancelled()
        cax, rix = ca_[xi]
        cby, rjy = cb_[yi]
        if cax.size and cby.size:
            S[np.ix_(cax, cby)] += Pm[np.ix_(rix, rjy)]
        if tick is not None:
            tick()
    S /= (na * nb)
    ca, cb = accel.mea_align(S, cancelled)

    def emit(rows, cols):
        return [(g, "".join(r[c] if c >= 0 else "-" for c in cols)) for g, r in rows]

    return emit(rows_a, ca) + emit(rows_b, cb)


def _progressive(seqs, post, dist=None, guide_seed=None, progress=None, cancelled=None,
                 parallel=False, tick=None):
    n = len(seqs)
    if dist is not None:
        D = dist
    else:
        k = min(3, max(1, min(len(s) for s in seqs)))
        rng = None if guide_seed is None else np.random.default_rng(guide_seed)
        D = _kmer_distance(seqs, k=k, rng=rng)
    Dm = {(min(i, j), max(i, j)): D[i, j] for i in range(n) for j in range(i + 1, n)}
    profiles = {i: [(i, seqs[i])] for i in range(n)}
    sizes = {i: 1 for i in range(n)}
    active = list(range(n))
    nxt = n
    total = max(1, n - 1)
    done = 0
    while len(active) > 1:
        if cancelled is not None and cancelled():
            raise Cancelled()
        if progress is not None:
            progress(done, total)
        _, a, b = min((Dm[(min(x, y), max(x, y))], x, y)
                      for ix, x in enumerate(active) for y in active[ix + 1:])
        merged = _merge(profiles[a], profiles[b], post, cancelled, parallel, tick)
        gid = nxt; nxt += 1
        profiles[gid] = merged
        sizes[gid] = sizes[a] + sizes[b]
        for c in active:
            if c in (a, b):
                continue
            da = Dm[(min(a, c), max(a, c))]
            db = Dm[(min(b, c), max(b, c))]
            Dm[(min(gid, c), max(gid, c))] = (da * sizes[a] + db * sizes[b]) / (sizes[a] + sizes[b])
        active = [c for c in active if c not in (a, b)] + [gid]
        done += 1
    if progress is not None:
        progress(total, total)
    return profiles[active[0]]


# --------------------------------------------------------------------------- #
# Guide tree from posteriors + iterative refinement
# --------------------------------------------------------------------------- #

def _posterior_distances(seqs, P) -> np.ndarray:
    """Guide-tree distances from the pairwise posteriors already computed:
    1 - (expected number of aligned residues / length of the shorter sequence).
    A far better signal than k-mer overlap, and essentially free here."""
    n = len(seqs)
    D = np.zeros((n, n))
    for x in range(n):
        for y in range(x + 1, n):
            m = min(len(seqs[x]), len(seqs[y])) or 1
            sim = min(1.0, float(_post(P, x, y).sum()) / m)
            D[x, y] = D[y, x] = 1.0 - sim
    return D


def _drop_allgap_cols(rows):
    """rows = [(gid, gapped_row)]; drop columns that are all-gap within this group."""
    if not rows:
        return rows
    L = len(rows[0][1])
    keep = [c for c in range(L) if any(r[c] != "-" for _, r in rows)]
    return [(g, "".join(r[c] for c in keep)) for g, r in rows]


def _sp_posterior_score(rows, post) -> float:
    """Sum over sequence pairs of the posterior mass on their aligned residues."""
    maps = [(g, domain.residue_index_array(r)) for g, r in rows]
    tot = 0.0
    for i in range(len(maps)):
        gx, mx = maps[i]
        ax = mx >= 0
        for j in range(i + 1, len(maps)):
            gy, my = maps[j]
            both = ax & (my >= 0)
            if both.any():
                Pm = post(gx, gy)
                tot += float(Pm[mx[both], my[both]].sum())
    return tot


def _refine(current, post, n_iters, rng, cancelled=None):
    """Iterative refinement: repeatedly split the sequences into two groups,
    re-align the two sub-profiles, and keep the result if the sum-of-pairs
    posterior score improves. Targets column accuracy the progressive pass misses."""
    best = _sp_posterior_score(current, post)
    gids = [g for g, _ in current]
    for _ in range(n_iters):
        if cancelled is not None and cancelled():
            raise Cancelled()
        rng.shuffle(gids)
        cut = int(rng.integers(1, len(gids)))
        group = set(gids[:cut])
        a = _drop_allgap_cols([(g, r) for g, r in current if g in group])
        b = _drop_allgap_cols([(g, r) for g, r in current if g not in group])
        if not a or not b:
            continue
        merged = _merge(a, b, post, cancelled)
        sc = _sp_posterior_score(merged, post)
        if sc > best + 1e-9:
            current, best = merged, sc
    return current


# --------------------------------------------------------------------------- #
# Parameter estimation
# --------------------------------------------------------------------------- #

def _estimate_identity(rows: Sequence[str]) -> float:
    ids = []
    n = len(rows)
    for a in range(n):
        for b in range(a + 1, n):
            m = [(x, y) for x, y in zip(rows[a], rows[b]) if x != "-" and y != "-"]
            if len(m) > 5:
                ids.append(sum(x == y for x, y in m) / len(m))
    return float(np.mean(ids)) if ids else 0.85


def estimate_params(rows: Sequence[str], alphabet: Alphabet):
    """Emission model and gap parameters read off a committed alignment.

    Returns ``(model, delta, epsilon)``. This is the single place that turns an
    alignment into pair-HMM parameters, so that a column is always *scored*
    under the model the data implies rather than under a fixed default. The
    reliability scores and the aligner's own pilot pass both go through it.

    Note that for protein data ``emission_model`` is currently independent of
    the estimated identity — one matrix, BLOSUM45, is used at every divergence —
    so ``theta`` affects nucleotide models only. The gap parameters are
    estimated for both alphabets. The clamp on ``theta`` below means it carries
    no information under 55% identity in any case, which is where most protein
    reference data sits; keying the protein matrix to divergence would have to
    read the raw identity, not this.
    """
    kind = "protein" if alphabet == Alphabet.PROTEIN else "dna"
    theta = round(min(0.95, max(0.55, _estimate_identity(rows))), 2)
    delta, epsilon = _estimate_gaps(rows)
    return emission_model(kind, theta), delta, epsilon


def _estimate_gaps(rows: Sequence[str]) -> Tuple[float, float]:
    opens = gappos = respos = 0
    for r in rows:
        prev_gap = False
        for ch in r:
            if ch == "-":
                gappos += 1
                if not prev_gap:
                    opens += 1
                prev_gap = True
            else:
                respos += 1
                prev_gap = False
    delta = opens / max(1, respos + opens)
    epsilon = (1.0 - opens / gappos) if gappos else 0.4
    return min(0.2, max(0.005, delta)), min(0.9, max(0.1, epsilon))


# --------------------------------------------------------------------------- #
# Public entry point
# --------------------------------------------------------------------------- #

def _posterior_gb(n: int, maxlen: int, sparse: bool = False) -> float:
    """Approximate peak RAM (GB) for the consistency transformation.

    Counts what is actually live at once, which the previous version did not:
    ``consistency_transform`` builds the new generation of posteriors while the
    old one is still referenced, so the dense path needs room for two. The
    sparse path instead holds two thresholded copies — a few per cent of dense in
    practice, taken here as a conservative 10% — plus one transient dense matrix
    while a pair is re-sparsified.
    """
    one = (maxlen ** 2) * 8 / 1e9
    pairs = n * (n - 1) // 2
    if sparse:
        return 2 * 0.10 * pairs * one + one
    return 2 * pairs * one


_EFFORT = {
    "min": {"consistency_iters": 0, "refine_iters": 0, "guide": "kmer"},
    "med": {"consistency_iters": 2, "refine_iters": 0, "guide": "posterior"},
    "max": {"consistency_iters": 2, "refine_iters": 25, "guide": "posterior"},
}


def align(
    records: Sequence[Record],
    alphabet: Alphabet,
    model: Optional[EmissionModel] = None,
    delta: Optional[float] = None,
    epsilon: Optional[float] = None,
    effort: str = "med",
    consistency_iters: Optional[int] = None,
    refine_iters: Optional[int] = None,
    guide: Optional[str] = None,
    estimate: bool = True,
    consistency_mem_gb: float = 1.0,
    progress=None,
    cancelled=None,
    guide_seed: Optional[int] = None,
    sparse: Optional[bool] = None,
    matrix: Optional[str] = None,
    tick=None,
) -> Alignment:
    """Progressive multiple alignment by posterior decoding, tiered by memory.

    Small/medium families get the full ProbCons **consistency transformation**
    (accurate). Families whose dense pairwise posteriors would exceed
    ``consistency_mem_gb`` get the sparse, thresholded transformation if SciPy is
    installed and that fits (``sparse`` forces one or the other), and otherwise
    fall back automatically to **streaming plain-MEA**: each
    pairwise posterior is computed on demand and discarded, so peak memory is a
    single L x L matrix and the aligner scales to large inputs (less accurate; for
    large jobs an external engine is still preferable). ``progress(done, total)``
    is called before each merge and ``cancelled()`` is polled throughout.
    ``matrix`` names the protein substitution matrix (default
    ``accel.PROTEIN_MATRIX``); it is ignored for nucleotides and when ``model``
    is given. ``tick()`` is called once per sequence pair as the progressive
    pass merges it, n(n-1)/2 times in all, for a caller that needs finer
    progress than one step per merge. Pairwise posteriors are computed on
    :data:`accel.WORKERS` threads; the result does not depend on how many.
    """
    preset = _EFFORT.get(effort, _EFFORT["med"])
    consistency_iters = preset["consistency_iters"] if consistency_iters is None else consistency_iters
    refine_iters = preset["refine_iters"] if refine_iters is None else refine_iters
    guide = preset["guide"] if guide is None else guide

    ids = [r[0] for r in records]
    seqs = [r[1].replace("-", "") for r in records]
    if not seqs:
        raise ValueError("no sequences to align")
    if len(set(ids)) != len(ids):
        # Rows are reordered back onto the input order through a dict keyed on
        # id, and the sandbox splices by id too, so duplicates silently drop
        # sequences rather than failing.
        dup = next(i for i in ids if ids.count(i) > 1)
        raise ValueError(f"duplicate sequence id {dup!r}; ids must be unique")
    if len(seqs) == 1:
        return Alignment(list(ids), list(seqs), alphabet)

    kind = "protein" if alphabet == Alphabet.PROTEIN else "dna"
    matrix = matrix or accel.PROTEIN_MATRIX
    n = len(seqs)
    maxlen = max(len(s) for s in seqs)
    if sparse is None:
        # The sparse transform approximates (see consistency_transform), so it is
        # used only when the exact one would not fit in the memory budget.
        sparse = _posterior_gb(n, maxlen) > consistency_mem_gb and sparse_available()
    use_consistency = (consistency_iters > 0
                       and _posterior_gb(n, maxlen, sparse) <= consistency_mem_gb)

    if model is None and estimate:
        # streaming plain-MEA pilot to read off divergence + gap rates (memory-safe)
        pm = emission_model(kind, matrix=matrix)

        def pilot_post(x, y):
            return accel.posterior_matrix(seqs[x], seqs[y], pm, _DEFAULT_DELTA, _DEFAULT_EPSILON)

        pilot = [r for _, r in _progressive(seqs, pilot_post, cancelled=cancelled,
                                            parallel=True)]
        theta = round(min(0.95, max(0.55, _estimate_identity(pilot))), 2)
        gdelta, gepsilon = _estimate_gaps(pilot)
        model = emission_model(kind, theta, matrix)
        if delta is None:
            delta = gdelta
        if epsilon is None:
            epsilon = gepsilon

    if model is None:
        model = emission_model(kind, matrix=matrix)
    if delta is None:
        delta = _DEFAULT_DELTA
    if epsilon is None:
        epsilon = _DEFAULT_EPSILON

    if use_consistency:
        raw = _all_pairs_posteriors(seqs, model, delta, epsilon, cancelled)
        D = _posterior_distances(seqs, raw) if guide == "posterior" else None
        P = consistency_transform(raw, n, consistency_iters, cancelled, sparse=sparse)

        def post(x, y):
            return _post(P, x, y)
    else:
        D = None                                     # streaming -> cheap k-mer guide tree
        def post(x, y):
            return accel.posterior_matrix(seqs[x], seqs[y], model, delta, epsilon)

    final = _progressive(seqs, post, dist=D, guide_seed=guide_seed,
                         progress=progress, cancelled=cancelled,
                         parallel=not use_consistency, tick=tick)
    if refine_iters > 0 and use_consistency:
        rng = np.random.default_rng(0 if guide_seed is None else guide_seed)
        final = _refine(final, post, refine_iters, rng, cancelled)

    order = {i: k for k, i in enumerate(ids)}
    pairs = sorted(((ids[gi], row) for gi, row in final), key=lambda p: order[p[0]])
    return Alignment([p[0] for p in pairs], [p[1] for p in pairs], alphabet)


def similarity_order(seqs: Sequence[str]) -> List[int]:
    """A reading order that puts similar sequences next to each other:
    a greedy nearest-neighbour chain over k-mer cosine distances."""
    n = len(seqs)
    if n <= 2:
        return list(range(n))
    clean = [s.replace("-", "") for s in seqs]
    k = min(3, max(1, min(len(s) for s in clean)))
    D = _kmer_distance(clean, k=k)
    order = [0]
    remaining = set(range(1, n))
    while remaining:
        last = order[-1]
        nxt = min(remaining, key=lambda j: D[last][j])
        order.append(nxt)
        remaining.discard(nxt)
    return order
