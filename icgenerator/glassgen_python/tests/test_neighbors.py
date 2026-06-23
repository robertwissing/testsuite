import numpy as np

from glassgen.neighbors import NeighborList


def brute_force_knn(pos, box, nsmooth):
    """Distance to the (nsmooth-1)-th non-self neighbor -> 2h, per particle."""
    d = pos[:, None, :] - pos[None, :, :]
    if box is not None:
        d -= box * np.rint(d / box)
    dist = np.sqrt((d ** 2).sum(axis=2))
    dist_sorted = np.sort(dist, axis=1)  # column 0 is self (0.0)
    return 0.5 * dist_sorted[:, nsmooth - 1]


def test_h_matches_brute_force_periodic():
    rng = np.random.default_rng(1)
    box = np.array([1.0, 1.0, 1.0])
    pos = rng.uniform(-0.5, 0.5, (500, 3))
    nb = NeighborList(box, nsmooth=32)
    _, h, _, _ = nb.update(pos)
    href = brute_force_knn(pos, box, 32)
    assert np.allclose(h, href)


def test_h_matches_brute_force_open():
    rng = np.random.default_rng(2)
    pos = rng.uniform(-0.5, 0.5, (400, 3))
    nb = NeighborList(None, nsmooth=24)
    _, h, _, _ = nb.update(pos)
    href = brute_force_knn(pos, None, 24)
    assert np.allclose(h, href)


def test_transpose_graph_is_exact_transpose():
    rng = np.random.default_rng(3)
    pos = rng.uniform(-0.5, 0.5, (300, 3))
    nb = NeighborList(np.ones(3), nsmooth=16, margin=4)
    idx, _, indptr, indices = nb.update(pos)
    n, k = idx.shape
    fwd = set((i, int(j)) for i in range(n) for j in idx[i])
    # indices[indptr[i]:indptr[i+1]] lists the j whose ball contains i,
    # i.e. exactly the j with a forward edge (j, i)
    rev = set((int(j), i) for i in range(n)
              for j in indices[indptr[i]:indptr[i + 1]])
    assert fwd == rev


def test_lattice_h():
    # simple cubic lattice spacing a: with nsmooth=33 the ball encloses the
    # 32 nearest non-self sites; the 32nd sorted neighbor is at 2a
    nside = 8
    a = 1.0 / nside
    g = (np.arange(nside) + 0.5) * a - 0.5
    pos = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
    nb = NeighborList(np.ones(3), nsmooth=33)
    _, h, _, _ = nb.update(pos)
    assert np.allclose(2 * h, 2 * a)


def test_grid_build_matches_tree_nonuniform():
    # clustered, anisotropic box: the numba cell-grid build must reproduce
    # the kd-tree exactly (h, neighbor sets, skin) including the fail/retry
    # path for sparse-region particles
    rng = np.random.default_rng(5)
    box = np.array([2.0, 1.0, 0.5])
    n = 4000
    pos = rng.uniform(-0.5, 0.5, (n, 3)) * box
    blob = rng.normal(0, 0.05, (n // 2, 3))  # dense clump at the origin
    pos[:n // 2] = blob - box * np.rint(blob / box)
    nb = NeighborList(box, nsmooth=32, margin=8)
    idx, h, _, _ = nb.update(pos)
    assert nb.last_build == 'grid'
    ref = NeighborList(box, nsmooth=32, margin=8)
    ref._build_tree(pos.astype(np.float64))
    assert np.allclose(h, ref.h, rtol=0, atol=1e-12)
    for i in range(n):
        assert set(idx[i]) == set(ref.idx[i]), i
    assert abs(nb._skin - ref._skin) < 1e-12


def test_skin_reuse_matches_rebuild():
    rng = np.random.default_rng(4)
    box = np.ones(3)
    pos = rng.uniform(-0.5, 0.5, (600, 3))
    nb = NeighborList(box, nsmooth=32, margin=8)
    nb.update(pos)
    # small coherent drift, below the skin trigger
    for _ in range(3):
        step = rng.normal(0, 1e-4, pos.shape)
        pos += step
        pos -= box * np.rint(pos / box)
        nb.record_move(float(np.max(np.abs(step))))
    _, h_reused, _, _ = nb.update(pos)
    nbuilds = nb.nbuilds
    fresh = NeighborList(box, nsmooth=32, margin=8)
    _, h_fresh, _, _ = fresh.update(pos)
    assert nbuilds == 1  # the drift must not have triggered a rebuild
    assert np.allclose(h_reused, h_fresh)
