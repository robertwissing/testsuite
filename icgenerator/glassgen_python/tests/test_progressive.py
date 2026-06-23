import numpy as np

from glassgen import (relax, relax_progressive, default_levels,
                      cascade_levels)


def _uniform(box):
    return lambda p: np.full(len(p), 1.0)


def test_default_levels_coarse_to_fine():
    lv = default_levels(nx_full=16, nx_coarse=8)
    assert [l['name'] for l in lv] == ['coarse', 'fine']
    # factor = nearest power of two to (16/8)^3 = 8
    assert lv[1]['split_factor'] == 8
    assert lv[0]['confine'] is False and lv[0]['reanneal'] is True
    assert lv[1]['confine'] is True and lv[1]['reanneal'] is False


def test_momentum_lands_on_final_fine_only():
    # momentum goes on the FINAL fine relax only - NOT the coarse box-wide
    # stage, and NOT the short cascade intermediate heals (too few iters to
    # settle -> overshoot on a high-force post-split state).
    lv = default_levels(50, 25, momentum=0.3)
    assert 'momentum' not in lv[0] and lv[0]['name'] == 'coarse'
    assert lv[1]['name'] == 'fine' and lv[1]['momentum'] == 0.3
    cv = cascade_levels(50, 12, momentum=0.3)
    assert 'momentum' not in cv[0]            # coarse: no momentum
    for l in cv[1:-1]:                         # intermediate heals: no momentum
        assert l['name'].startswith('cascade') and l['momentum'] == 0.0
    assert cv[-1]['name'] == 'fine' and cv[-1]['momentum'] == 0.3
    # fine_kw without a momentum key must NOT clobber the injected value
    lv2 = default_levels(50, 25, momentum=0.4,
                         fine_kw=dict(confine=True, reanneal=False))
    assert lv2[1]['momentum'] == 0.4


def test_default_levels_coarse_skip():
    # nx_full <= nx_coarse -> single full-res level (no split)
    lv = default_levels(nx_full=8, nx_coarse=16)
    assert len(lv) == 1
    assert lv[0]['split_factor'] == 1


def test_progressive_coarse_to_fine_converges():
    box = np.array([1.0, 1.0, 1.0])
    rng = np.random.default_rng(0)
    nfull = 12
    N = nfull ** 3
    pos = (rng.random((N, 3)) - 0.5) * box
    mass = np.full(N, 1.0 / N)
    levels = default_levels(nfull, 6,
                            coarse_kw=dict(max_iter=300),
                            fine_kw=dict(max_iter=300))
    res = relax_progressive(pos, mass, _uniform(box), box, nx_full=nfull,
                            nx_coarse=6, levels=levels)
    # split factor nearest pow2 to (12/6)^3 = 8 -> final N = 6^3 * 8
    assert len(res.pos) == 6 ** 3 * 8
    err = np.mean(np.abs(1.0 - res.rho / 1.0))
    assert err < 0.05            # uniform glass density well fit
    assert np.isfinite(res.pos).all()
    assert np.all(np.abs(res.pos) <= box / 2 + 1e-9)


def test_cascade_levels_one_split_per_stage():
    # cascade: coarse, then k binary (x2) splits one resolution level at a time,
    # the last being the full relax. k = log2((nx_full/nx_coarse)^3).
    lv = cascade_levels(64, 32, inter_iter=50)
    assert [l['name'] for l in lv] == ['coarse', 'cascade1', 'cascade2', 'fine']
    assert lv[0]['split_factor'] == 1                       # coarse: no split
    assert all(l['split_factor'] == 2 for l in lv[1:])      # every later x2
    # intermediates are short confined smooths; the final is a full confined relax
    assert lv[1]['max_iter'] == 50 and lv[1]['confine'] is True
    assert lv[-1]['max_iter'] == 4800 and lv[-1]['confine'] is True
    assert lv[-1]['reanneal'] is False
    # nx_full <= nx_coarse collapses to a single full-res level
    assert [l['name'] for l in cascade_levels(16, 32)] == ['full']
    # final particle count after k x2 splits == nx_full**3 (matched mass)
    assert 32 ** 3 * 2 ** 3 == 64 ** 3


def test_cascade_progressive_converges():
    # end-to-end: the cascade schedule must still produce a valid glass
    box = np.array([1.0, 1.0, 1.0])
    rng = np.random.default_rng(1)
    nfull = 12
    N = nfull ** 3
    pos = (rng.random((N, 3)) - 0.5) * box
    mass = np.full(N, 1.0 / N)
    levels = cascade_levels(nfull, 6, inter_iter=40,
                            coarse_kw=dict(max_iter=300),
                            fine_kw=dict(max_iter=300))
    res = relax_progressive(pos, mass, _uniform(box), box, nx_full=nfull,
                            nx_coarse=6, levels=levels)
    assert len(res.pos) == 6 ** 3 * 8
    assert [l['name'] for l in res.levels] == \
        ['coarse', 'cascade1.split', 'cascade1', 'cascade2.split', 'cascade2',
         'fine.split', 'fine']
    assert np.mean(np.abs(1.0 - res.rho)) < 0.05
    assert np.isfinite(res.pos).all()


def test_progressive_records_post_split_diagnostic():
    # res.levels must carry a post-split entry between the coarse and fine
    # stages, measuring the glass the fine stage inherits from the split (before
    # any fine relaxation). The split places children off-glass, so its denserr
    # is much worse than the coarse end-state and than the refined fine state.
    box = np.array([1.0, 1.0, 1.0])
    rng = np.random.default_rng(0)
    nfull = 12
    N = nfull ** 3
    pos = (rng.random((N, 3)) - 0.5) * box
    mass = np.full(N, 1.0 / N)
    levels = default_levels(nfull, 6,
                            coarse_kw=dict(max_iter=300),
                            fine_kw=dict(max_iter=300))
    res = relax_progressive(pos, mass, _uniform(box), box, nx_full=nfull,
                            nx_coarse=6, levels=levels)
    names = [l['name'] for l in res.levels]
    assert names == ['coarse', 'fine.split', 'fine']
    coarse, split, fine = res.levels
    assert split['stop_reason'] == 'post-split' and split['niter'] == 0
    # the post-split state is measured at full resolution (already split)
    assert split['n'] == fine['n'] == 6 ** 3 * 8
    for k in ('denserr_mean', 'denserr_p95', 'force_mean'):
        assert np.isfinite(split[k])
    # split degrades the glass; the fine stage then repairs it below the split
    assert split['denserr_mean'] > coarse['denserr_mean']
    assert fine['denserr_mean'] < split['denserr_mean']


def test_progressive_skip_equals_single_relax():
    # nx_full <= nx_coarse: the driver should relax the input directly, giving
    # the same result (same seed) as a plain relax() with matching settings
    box = np.array([1.0, 1.0, 1.0])
    rng = np.random.default_rng(2)
    N = 8 ** 3
    pos = (rng.random((N, 3)) - 0.5) * box
    mass = np.full(N, 1.0 / N)
    lv = default_levels(8, 16, coarse_kw=dict(max_iter=120))
    res_p = relax_progressive(pos.copy(), mass, _uniform(box), box, nx_full=8,
                              nx_coarse=16, levels=lv)
    res_r = relax(pos.copy(), mass, _uniform(box), box=box,
                  max_iter=120, confine=False, reanneal=True)
    assert len(res_p.pos) == N
    assert res_p.niter == res_r.niter
    np.testing.assert_allclose(res_p.pos, res_r.pos)


def test_progressive_merge_seed_restores_input_count():
    # coarse_seed='merge' coarsens the full-res input down by the cascade factor
    # then splits back -> final count == input count, EXACTLY (1:1).
    rng = np.random.default_rng(5)
    box = np.ones(3)
    nx_full, nx_coarse = 16, 8
    N = nx_full ** 3
    pos = rng.uniform(-0.5, 0.5, (N, 3))
    mass = np.full(N, 1.0 / N)
    lv = cascade_levels(nx_full, nx_coarse, inter_iter=5)
    res = relax_progressive(pos, mass, _uniform(box), box, nx_full=nx_full,
                            nx_coarse=nx_coarse, levels=lv, coarse_seed='merge',
                            nsmooth=64, max_iter=200)
    assert len(res.pos) == N


def test_auto_coarse_count_reproduces_mhd_anchor():
    # the density criterion (alpha=2) on the mhdcollapse 360:1 sigmoid must land
    # at the empirical sweet spot ~34k (nx_coarse~25-32): F=8 at the nx50 count,
    # F=1 (no merge) at the nx12 count.
    from glassgen import auto_coarse_count
    from glassgen.progressive import _floor_pow2
    R, Rs, inner, outer = 0.015, 0.0015, 70.7355, 0.196488
    box = np.array([0.15, 0.15, 0.15])

    def rho0(p):
        r = np.sqrt((p ** 2).sum(1))
        return outer + (inner - outer) / (1.0 + np.exp(-(R - r) / Rs))

    ncmin = auto_coarse_count(rho0, box, alpha=2.0, grid=96)
    assert 3.0e4 < ncmin < 3.8e4, ncmin           # ~34k
    assert _floor_pow2(312210 / ncmin) == 8       # nx50 -> F=8
    assert _floor_pow2(4220 / ncmin) == 1         # nx12 -> no merge


def test_progressive_auto_merge_uniform_skips_merge():
    # uniform target -> grad rho=0 -> no large-scale structure to capture, so
    # the auto path picks factor=1 (no merge): a single full-res relax, count
    # unchanged, no coarse/cascade stages.
    rng = np.random.default_rng(7)
    box = np.ones(3)
    N = 16 ** 3
    pos = rng.uniform(-0.5, 0.5, (N, 3))
    mass = np.full(N, 1.0 / N)
    res = relax_progressive(pos, mass, _uniform(box), box, nsmooth=64,
                            max_iter=150)
    assert len(res.pos) == N
    assert [l['name'] for l in res.levels] == ['full']
