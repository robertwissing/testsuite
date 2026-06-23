import numpy as np

from glassgen import CallableDensity, TableDensity, relax


def test_uniform_glass_periodic():
    rng = np.random.default_rng(11)
    n = 12 ** 3
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    mass = np.full(n, 1.0 / n)
    # Wendland C2 needs ~64 neighbours for good density accuracy (32 is
    # under-resolved for it), so the e2e quality tests use the nsmooth=64
    # default rather than 32.
    res = relax(pos, mass, lambda p: np.ones(len(p)), box=np.ones(3),
                nsmooth=64, max_iter=1200)
    assert res.converged
    err = np.abs(1.0 - res.rho)
    assert err.mean() < 0.02
    assert res.rho.std() / res.rho.mean() < 0.02
    # positions stayed wrapped
    assert np.all(np.abs(res.pos) <= 0.5)


def sigmoid_profile(x, inner=2.0, outer=1.0, R=0.25, s=0.05):
    # C profile 3 along a coordinate (pkdUpdateDensity)
    return outer + (inner - outer) * 0.5 * (np.tanh((x + R) / s)
                                            - np.tanh((x - R) / s))


def test_sigmoid_profile_along_x():
    rng = np.random.default_rng(12)
    n = 16 ** 3
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    # equal-mass particles: total mass must integrate the target profile
    xq = np.linspace(-0.5, 0.5, 20001)
    mtot = np.trapz(sigmoid_profile(xq), xq)
    mass = np.full(n, mtot / n)
    xt = np.linspace(-0.5, 0.5, 2001)
    dens = TableDensity((xt,), sigmoid_profile(xt), mapping='signed_x')
    res = relax(pos, mass, dens, box=np.ones(3), nsmooth=64, max_iter=1500)
    # binned density vs target in the smooth regions
    bins = np.linspace(-0.5, 0.5, 21)
    ib = np.digitize(res.pos[:, 0], bins) - 1
    prof = np.array([res.rho[ib == b].mean() for b in range(20)])
    centers = 0.5 * (bins[:-1] + bins[1:])
    target = sigmoid_profile(centers)
    # exclude the ramps: kernel smoothing (h ~ ramp width) biases the SPH
    # density estimate there by a few percent, as in the C generator
    smooth = np.abs(np.abs(centers) - 0.25) > 0.1
    relerr = np.abs(prof[smooth] / target[smooth] - 1.0)
    assert relerr.max() < 0.05, relerr


def test_reorder_and_grid_match_tree(monkeypatch):
    # cell-order reordering + grid neighbor build are pure optimizations:
    # a run with both must track a run forced onto the kd-tree path with
    # reordering off (results returned in input order either way)
    from glassgen.neighbors import NeighborList
    rng = np.random.default_rng(14)
    n = 8 ** 3
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    mass = np.full(n, 1.0 / n)
    dens = lambda p: np.ones(len(p))  # noqa: E731
    kw = dict(box=np.ones(3), nsmooth=32, max_iter=8)
    res_fast = relax(pos.copy(), mass, dens, reorder=True, **kw)
    monkeypatch.setattr(NeighborList, '_build_grid', lambda self, p: False)
    res_ref = relax(pos.copy(), mass, dens, reorder=False, **kw)
    assert np.allclose(res_fast.pos, res_ref.pos, rtol=0, atol=1e-9)
    assert np.allclose(res_fast.rho, res_ref.rho, rtol=1e-9)
    assert np.allclose(res_fast.h, res_ref.h, rtol=1e-12)


def test_one_sided_force_glass_quality():
    # the transpose-free walk perturbs forces only for pairs at strong h
    # gradients near the kernel edge; the converged glass must be
    # statistically equivalent (density error, spread), though not
    # particle-identical
    rng = np.random.default_rng(15)
    n = 12 ** 3
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    mass = np.full(n, 1.0 / n)
    dens = lambda p: np.ones(len(p))  # noqa: E731
    kw = dict(box=np.ones(3), nsmooth=64, max_iter=1200)
    res_sym = relax(pos.copy(), mass, dens, one_sided=False, **kw)
    res_one = relax(pos.copy(), mass, dens, one_sided=True, **kw)
    assert res_one.converged
    for res in (res_sym, res_one):
        assert np.abs(1.0 - res.rho).mean() < 0.02
    spread_sym = res_sym.rho.std() / res_sym.rho.mean()
    spread_one = res_one.rho.std() / res_one.rho.mean()
    assert spread_one < 0.02
    # equally good glass: spreads agree to a few 10s of percent relative
    assert abs(spread_one - spread_sym) < 0.3 * max(spread_sym, spread_one)


def test_hmode_density_h_modes_converge():
    # the density-based h options (measured = C DenDVDX h from rho_sph,
    # target = analytic h from rho0) must still produce a valid uniform glass
    rng = np.random.default_rng(21)
    n = 12 ** 3
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    mass = np.full(n, 1.0 / n)
    dens = lambda p: np.ones(len(p))  # noqa: E731
    kw = dict(box=np.ones(3), nsmooth=64, max_iter=1500)
    for hmode in ('measured', 'target'):
        res = relax(pos.copy(), mass, dens, hmode=hmode, **kw)
        assert res.converged, hmode
        assert np.abs(1.0 - res.rho).mean() < 0.02, hmode
        assert res.rho.std() / res.rho.mean() < 0.02, hmode


def test_icr0_stop_banks_deep_polish_on_ball_path():
    # The ball (density-based-h) path enters a long, quality-useless polish
    # grind on gradient targets: denserr saturates early but icr0 keeps random-
    # walking down for thousands more iterations. The default icr0_stop finish
    # (median-smoothed icr0 < threshold) must end the run as soon as icr0 settles
    # onto its critical-step setpoint, producing a glass statistically identical
    # to running to the cap - a direct test that it fires on a real relax and
    # that the banked deep-polish iterations did not change quality.
    rng = np.random.default_rng(7)
    nx = 16
    g = (np.arange(nx) + 0.5) / nx - 0.5
    pos = np.stack(np.meshgrid(g, g, g, indexing='ij'), -1).reshape(-1, 3)
    pos += (rng.random(pos.shape) - 0.5) / nx
    box = np.ones(3)
    pos -= box * np.rint(pos / box)
    rho0 = lambda p: 1.0 + 8.0 / (1.0 + np.exp(-12.0 * p[:, 0]))  # noqa: E731
    mass = np.full(len(pos), rho0(pos).mean() / len(pos))
    kw = dict(box=box, nsmooth=32, hmode='measured', max_iter=2500,
              stop_window=20)

    stopped = relax(pos.copy(), mass, rho0, icr0_stop=3e-3, **kw)
    capped = relax(pos.copy(), mass, rho0, icr0_stop=0.0, **kw)

    # icr0_stop fired, well before the cap (the deep-polish grind is skipped)
    assert stopped.stop_reason == 'icr0_abs'
    assert stopped.niter < capped.niter
    # ...and pstop is only ever fed in polish, so the stop happened in polish
    assert capped.stop_reason == 'max_iter'      # backstop disabled -> cap
    # the skipped iterations were quality-useless: same glass either way
    derr_stop = np.abs(1.0 - stopped.rho / rho0(stopped.pos)).mean()
    derr_cap = np.abs(1.0 - capped.rho / rho0(capped.pos)).mean()
    assert abs(derr_stop - derr_cap) < 0.1 * derr_cap


def test_polish_iter_fixed_budget_stops_after_exactly_n_polish_iters():
    # the fixed polish-iteration budget: with icr0_stop disabled, the ONLY
    # finish is polish_iter, and the run must
    # stop after EXACTLY that many polish-phase iterations.
    rng = np.random.default_rng(3)
    nx = 12
    g = (np.arange(nx) + 0.5) / nx - 0.5
    pos = np.stack(np.meshgrid(g, g, g, indexing='ij'), -1).reshape(-1, 3)
    pos += (rng.random(pos.shape) - 0.5) / nx
    box = np.ones(3)
    pos -= box * np.rint(pos / box)
    mass = np.full(len(pos), 1.0 / len(pos))
    rho0 = lambda p: np.ones(len(p))  # noqa: E731
    npol = 30
    res = relax(pos, mass, rho0, box=box, nsmooth=32, max_iter=4800,
                icr0_stop=0.0, polish_iter=npol)
    assert res.stop_reason == 'polish_iter'
    # history entry = (it, frms, dens_err, icr0, icrloop0, polish); count polish
    n_polish = sum(1 for h in res.history if h[5])
    assert n_polish == npol


def test_polish_iter_composes_icr0_stop_wins_when_earlier():
    # composition: a generous polish_iter plus the default icr0_stop -> icr0_stop
    # converges first on a uniform target (short polish), so it fires, not the
    # fixed budget.
    rng = np.random.default_rng(4)
    nx = 12
    g = (np.arange(nx) + 0.5) / nx - 0.5
    pos = np.stack(np.meshgrid(g, g, g, indexing='ij'), -1).reshape(-1, 3)
    pos += (rng.random(pos.shape) - 0.5) / nx
    box = np.ones(3)
    pos -= box * np.rint(pos / box)
    mass = np.full(len(pos), 1.0 / len(pos))
    rho0 = lambda p: np.ones(len(p))  # noqa: E731
    res = relax(pos, mass, rho0, box=box, nsmooth=32, max_iter=4800,
                icr0_stop=3e-3, polish_iter=100000)
    assert res.stop_reason == 'icr0_abs'


def test_hmode_target_floor_sparse_region():
    # analytic (target) h shrinks where rho0 is high; a particle sitting in a
    # sparse region whose target says "dense" would, unfloored, gather an
    # empty 2h ball -> zero force / NaN. nsmoothmin must keep it finite.
    rng = np.random.default_rng(22)
    n = 14 ** 3
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    mass = np.full(n, 1.0 / n)
    # target 50x denser for x>0 than the (uniform) particles actually are
    dens = lambda p: np.where(p[:, 0] > 0, 50.0, 1.0)  # noqa: E731
    res = relax(pos, mass, dens, box=np.ones(3), nsmooth=32, max_iter=300,
                hmode='target', nsmoothmin=8)
    assert np.all(np.isfinite(res.pos))
    assert np.all(np.isfinite(res.rho))
    assert np.all(res.rho > 0)


def test_open_boundary_smoke():
    # open domains have no confinement (positive pseudo-pressure at a free
    # surface expands without bound), so meaningful glass generation needs
    # a periodic box; here just verify the open code path runs and stays
    # finite over a short window (the surface deficit makes the global
    # density error grow, so no improvement assertion)
    rng = np.random.default_rng(13)
    n = 1000
    pos = rng.uniform(-0.5, 0.5, (n, 3))
    mass = np.full(n, 1.0 / n)
    dens = CallableDensity(lambda r: np.ones(len(r)), mapping='r_sph')
    res = relax(pos, mass, dens, box=None, nsmooth=32, max_iter=10)
    assert np.all(np.isfinite(res.pos))
    assert np.all(np.isfinite(res.rho))
