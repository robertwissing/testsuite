import numpy as np

from glassgen.kernels import dw_wc2, w_wc2
from glassgen.neighbors import NeighborList
from glassgen.sph_core import density_loop, get_force_loop

M_1_PI = 1.0 / np.pi


def make_lattice(nside, jitter=0.0, seed=0):
    a = 1.0 / nside
    g = (np.arange(nside) + 0.5) * a - 0.5
    pos = np.array(np.meshgrid(g, g, g, indexing='ij')).reshape(3, -1).T
    if jitter:
        rng = np.random.default_rng(seed)
        pos += rng.normal(0, jitter * a, pos.shape)
    return np.ascontiguousarray(pos)


def compute_state(pos, box, nsmooth):
    n = len(pos)
    mass = np.full(n, 1.0 / n)
    nb = NeighborList(box, nsmooth)
    idx, h, indptr, indices = nb.update(pos)
    rho = np.empty(n)
    density_loop(pos, mass, h, idx, box, True, rho)
    return mass, idx, h, indptr, indices, rho


def test_density_uniform_lattice():
    pos = make_lattice(10)
    box = np.ones(3)
    _, _, _, _, _, rho = compute_state(pos, box, nsmooth=64)
    # total mass 1 in unit box -> mean density ~1. Wendland C2 has a larger
    # discrete summation-density bias at nSmooth=64 (~1.5% on
    # a lattice - it is more pairing-stable but needs more neighbours for the
    # same density accuracy). The glass only cares about UNIFORMITY (std/mean,
    # still exact on the lattice); the constant scale is absorbed by the mass.
    assert abs(np.mean(rho) - 1.0) < 2e-2
    assert np.std(rho) / np.mean(rho) < 1e-3


def serial_reference(pos, mass, h, rho, pP, box, mode):
    """Literal O(N^2) port of one C ICGenerator sweep (smoothfcn.c:5794):
    visit every ordered pair (p, q in p's 2h_p ball), update both p and q,
    with fNorm1 = 0.5/(pi h_p^5). Also accumulates Q1 and E0."""
    n = len(pos)
    icvel = np.zeros((n, 3))
    q1 = np.zeros(n)
    e0 = np.zeros((n, 3))
    for p in range(n):
        ih2 = 1.0 / (h[p] * h[p])
        fnorm = 0.5 * M_1_PI * ih2 / h[p]
        fnorm1 = fnorm * ih2
        for q in range(n):
            d = pos[p] - pos[q]
            d -= box * np.rint(d / box)
            ar2 = (d @ d) * ih2
            if ar2 >= 4.0:
                continue
            rs = w_wc2(ar2) * fnorm
            drs = dw_wc2(ar2) * fnorm1
            rq = drs * mass[q]
            rp = drs * mass[p]
            if mode == 'gdforce':
                ig2p = 1.0 / (rho[p] * rho[q])
                ig2q = ig2p
            else:
                ig2p = 1.0 / (rho[p] * rho[p])
                ig2q = 1.0 / (rho[q] * rho[q])
            igp = pP[p] * ig2p
            igq = pP[q] * ig2q
            icvel[p] -= (igp + igq) * d * rq
            icvel[q] += (igp + igq) * d * rp
            q1[p] += mass[q] / rho[q] * rs * 2.0
            e0[p] += rho[p] * h[p] * (ig2p + ig2q) * d * rq
            e0[q] -= rho[q] * h[q] * (ig2p + ig2q) * d * rp
    return icvel, q1, e0


def test_force_loop_matches_serial_reference():
    pos = make_lattice(6, jitter=0.15, seed=7)
    box = np.ones(3)
    pos -= box * np.rint(pos / box)
    nsmooth = 32
    mass, idx, h, indptr, indices, rho = compute_state(pos, box, nsmooth)
    rho0 = np.full(len(pos), 1.0)
    pP = (rho / rho0) ** 2.0

    for mode in ('gdforce', 'rho2'):
        loop = get_force_loop(mode)
        n = len(pos)
        icvel = np.empty((n, 3))
        q1d = np.empty(n)
        e0 = np.empty((n, 3))
        loop(pos, mass, h, rho, pP, idx, indptr, indices,
             box, True, icvel, q1d, e0, True)
        ref_icvel, ref_q1, ref_e0 = serial_reference(
            pos, mass, h, rho, pP, box, mode)
        assert np.allclose(icvel, ref_icvel, rtol=1e-10, atol=1e-12), mode
        assert np.allclose(q1d, ref_q1, rtol=1e-10), mode
        assert np.allclose(e0, ref_e0, rtol=1e-10, atol=1e-12), mode


def test_q1_partition_of_unity():
    pos = make_lattice(8, jitter=0.05, seed=8)
    box = np.ones(3)
    pos -= box * np.rint(pos / box)
    mass, idx, h, indptr, indices, rho = compute_state(pos, box, 64)
    pP = np.ones(len(pos))
    loop = get_force_loop('gdforce')
    n = len(pos)
    icvel = np.empty((n, 3))
    q1d = np.empty(n)
    e0 = np.empty((n, 3))
    loop(pos, mass, h, rho, pP, idx, indptr, indices,
         box, True, icvel, q1d, e0, True)
    assert abs(np.mean(q1d) - 1.0) < 0.02
    assert np.max(np.abs(q1d - 1.0)) < 0.1


def test_force_points_from_overdense_to_underdense():
    # uniform lattice, uniform rho0: squeeze x>0 half slightly so it is
    # overdense; net force on squeezed particles must point outward (-x
    # gradient direction reversed: away from the overdensity)
    pos = make_lattice(8)
    box = np.ones(3)
    sel = pos[:, 0] > 0.0
    pos[sel, 0] = 0.25 + (pos[sel, 0] - 0.25) * 0.9
    pos -= box * np.rint(pos / box)
    mass, idx, h, indptr, indices, rho = compute_state(pos, box, 32)
    pP = (rho / 1.0) ** 2.0
    loop = get_force_loop('gdforce')
    n = len(pos)
    icvel = np.empty((n, 3))
    loop(pos, mass, h, rho, pP, idx, indptr, indices,
         box, True, icvel, np.empty(n), np.empty((n, 3)), False)
    # overdense region around x=0.25: particles on its left edge pushed -x,
    # on its right edge +x => correlation between icvel_x and (x - 0.25)
    # (within the squeezed half, away from periodic images)
    inner = sel & (np.abs(pos[:, 0] - 0.25) < 0.15)
    corr = np.corrcoef(pos[inner, 0] - 0.25, icvel[inner, 0])[0, 1]
    assert corr > 0.5
