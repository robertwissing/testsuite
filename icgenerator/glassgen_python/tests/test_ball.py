import numpy as np

from glassgen.ball import BallNeighborList, _transpose_csr


def _sets(bn, n):
    ip, ii = bn.indptr, bn.indices
    return [set(ii[ip[i]:ip[i + 1]].tolist()) for i in range(n)]


def test_parallel_transpose_is_exact():
    """The parallel counting-sort CSR transpose must equal the exact
    transpose, ordered with sources ascending within each target."""
    rng = np.random.default_rng(0)
    n = 2000
    counts = rng.integers(0, 40, n)
    indptr = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    m = int(indptr[-1])
    indices = rng.integers(0, n, m).astype(np.int32)
    tip, tii = _transpose_csr(indptr, indices, n)
    assert int(tip[-1]) == m
    fwd = [[] for _ in range(n)]
    for i in range(n):
        for e in range(indptr[i], indptr[i + 1]):
            fwd[indices[e]].append(i)
    for j in range(n):
        assert sorted(fwd[j]) == list(tii[tip[j]:tip[j + 1]])


def test_fused_build_matches_two_pass():
    """The fused single-pass build (store+scatter) must produce the same CSR
    neighbour sets as the robust 2-pass fallback, and a transpose that is the
    exact transpose of the forward graph. Structured h stresses varying ball
    sizes in an anisotropic periodic box."""
    rng = np.random.default_rng(7)
    box = np.array([2.0, 1.0, 0.5])
    n = 5000
    pos = rng.uniform(-0.5, 0.5, (n, 3)) * box
    pos -= box * np.rint(pos / box)
    h = 0.02 + 0.04 * np.abs(pos[:, 0]) / box[0]

    bf = BallNeighborList(box, nsmoothmin=8, skin=0.3, build_transpose=True)
    bf._build(pos.copy(), h.copy())
    sf = _sets(bf, n)

    orig = BallNeighborList._fused_gather
    try:
        # force the 2-pass fallback
        BallNeighborList._fused_gather = lambda self, *a, **k: None
        b2 = BallNeighborList(box, nsmoothmin=8, skin=0.3,
                              build_transpose=True)
        b2._build(pos.copy(), h.copy())
        s2 = _sets(b2, n)
    finally:
        BallNeighborList._fused_gather = orig

    assert int(bf.indptr[-1]) == int(b2.indptr[-1])
    assert all(sf[i] == s2[i] for i in range(n))

    # fused transpose is the exact transpose of the fused forward graph
    ip, ii, tip, tii = bf.indptr, bf.indices, bf.tip, bf.tii
    fwd = set((i, int(j)) for i in range(n) for j in ii[ip[i]:ip[i + 1]])
    rev = set((int(j), i) for i in range(n) for j in tii[tip[i]:tip[i + 1]])
    assert fwd == rev


def _true_ball(pos, box, h2):
    """Brute-force: set of j (!=i) with |x_i - x_j| < 2*h2[i], min-image."""
    d = pos[:, None, :] - pos[None, :, :]
    d -= box * np.rint(d / box)
    dist = np.sqrt((d ** 2).sum(2))
    np.fill_diagonal(dist, np.inf)
    return [set(np.nonzero(dist[i] < 2.0 * h2[i])[0].tolist())
            for i in range(len(pos))]


def test_per_particle_stale_never_drops_a_neighbour():
    """When the per-particle staleness test reports 'no rebuild', the cached
    padded ball MUST still contain every true neighbour within 2*h_now of each
    particle - otherwise the density/force loops silently miss neighbours.
    Drive a small box through random sub-skin steps (with per-particle h drift)
    and assert the superset invariant on every reuse iteration."""
    rng = np.random.default_rng(3)
    box = np.array([1.0, 1.0, 1.0])
    n = 2000
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    h = np.full(n, 0.5 * (64.0 / n) ** (1.0 / 3.0))   # ~nsmooth 64

    bn = BallNeighborList(box, nsmoothmin=8, skin=0.3, build_transpose=False)
    pdisp = np.zeros(n)
    bn.track_displacements(pdisp)
    reuses = 0
    for step in range(40):
        h_in = h * rng.uniform(0.97, 1.04, n)         # gentle per-particle drift
        rebuild = bn.needs_rebuild(h_in)
        _, indices, hout, _, _ = bn.update(pos, h_in, rebuild=rebuild)
        if not rebuild:
            reuses += 1
            ip = bn.indptr
            cached = [set(indices[ip[i]:ip[i + 1]].tolist()) for i in range(n)]
            true = _true_ball(pos, box, hout)
            for i in range(n):
                assert true[i] <= cached[i], (
                    f"step {step} particle {i}: missing {true[i] - cached[i]}")
        # small sub-step: accumulate per-particle displacement as the stepper
        # would, record the global max-move, then advance positions (wrapped)
        s = rng.normal(0, 0.002, (n, 3))
        pdisp += np.abs(s).max(1)
        bn.record_move(float(np.abs(s).max()))
        pos = pos + s
        pos -= box * np.rint(pos / box)
    assert reuses > 0     # the test must actually exercise the reuse path
