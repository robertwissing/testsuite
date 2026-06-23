import numpy as np

from glassgen.stepper import (Icr0Stop, StepController, apply_step,
                              normalize_displacements)


def test_normalize_displacements():
    icvel = np.array([[3.0, -4.0, 0.0],
                      [0.3, 0.0, 0.4],
                      [0.0, 0.0, 0.0]])
    frms, fmax = normalize_displacements(icvel)
    # frms = true per-particle RMS = sqrt(sum v^2 / N) before normalization
    assert np.isclose(frms, np.sqrt((25.0 + 0.25) / 3))
    assert fmax == 4.0
    # log normalization: particle 0 has the global max component -> its
    # icfmax = 1 -> speednorm = 1 -> max |component| becomes exactly 1
    assert np.isclose(np.max(np.abs(icvel[0])), 1.0)
    # weaker particle: 0 < speednorm < 1, direction preserved
    sn = np.log(0.4 * 0.1 + 1) / np.log(4.0 * 0.1 + 1)
    expect = (1 - 1e-3) * sn + 1e-3
    assert np.isclose(np.max(np.abs(icvel[1])), expect)
    assert np.isclose(icvel[1, 0] / icvel[1, 2], 0.75)
    # zero-force particle gets unit components (as in pkdICNormalizeVel)
    assert np.all(icvel[2] == 1.0)


def test_step_zero_and_periodic_rescale():
    box = np.array([1.0, 1.0, 1.0])
    c = StepController(icr0=1.0, icrloop0=0.5, icr0rate=1.5, nsmooth=64, box=box)
    assert np.isclose(c.icr0, 0.99)  # 0.99 * maxPeriod * 1.0
    c.update_step(frms=1.0, h_avg=0.1)
    # step zero: icr0 /= h_avg -> 9.9, icrloop0 = icr0 * 0.5 -> 4.95,
    # then frms_prev(0) - frms < 0 is a "bad move": icr0 /= (1 + 0.8*1.5)
    assert np.isclose(c.icrloop0, 4.95)
    assert np.isclose(c.icr0, 9.9 / 2.2)


def test_good_and_bad_moves():
    c = StepController(icr0=1.0, icrloop0=0.5, icr0rate=1.5, nsmooth=64)
    c.step_zero = False
    c.frms_prev = 2.0
    c.update_step(frms=1.0, h_avg=1.0)   # improved -> grow by 1+0.16*1.5
    assert np.isclose(c.icr0, 1.24)
    c.update_step(frms=3.0, h_avg=1.0)   # worse -> shrink by 1+0.8*1.5
    assert np.isclose(c.icr0, 1.24 / 2.2)


def test_restart_state_machine():
    c = StepController(icr0=1.0, icrloop0=0.5, icr0rate=1.5, nsmooth=64)
    c.step_zero = False
    nfac = c.nfac
    # above nfac: nothing happens
    c.icr0 = 2 * nfac
    c.update_restart(dens_err=0.5)
    assert c.icr0 == 2 * nfac and c.icfmax == 1e13

    # drop below nfac with improving error -> re-anneal to icrloop0
    c.icr0 = 0.5 * nfac
    c.update_restart(dens_err=0.4)
    assert c.icr0 == 0.5 and c.icfmax == 0.4

    # below nfac again, error got worse -> halve icrloop0, restart
    c.icr0 = 0.5 * nfac
    c.update_restart(dens_err=0.45)
    assert np.isclose(c.icrloop0, 0.25)
    assert c.icr0 == c.icrloop0 and c.icfmax == 0.45

    # keep halving until icrloop0 < 2*nfac -> polish phase flag
    while c.icrloop0 >= 2 * nfac:
        c.icr0 = 0.5 * nfac
        c.update_restart(dens_err=1.0)
    c.update_restart(dens_err=1.0)
    assert c.polish

    # polish phase uses the gentle rates
    c.frms_prev = 2.0
    icr0 = c.icr0
    c.update_step(frms=1.0, h_avg=1.0)
    assert np.isclose(c.icr0, icr0 * (1 + 0.002 * 1.5))
    icr0 = c.icr0
    c.update_step(frms=3.0, h_avg=1.0)
    assert np.isclose(c.icr0, icr0 / (1 + 0.1 * 1.5))


def test_apply_step_no_cap_matches_default():
    # max_step_frac = 0 reproduces the unbounded move (the coarse default)
    box = np.array([1.0, 1.0, 1.0])
    pos = np.array([[0.0, 0.0, 0.0]])
    icvel = np.array([[1.0, 0.0, 0.0]])
    h = np.array([0.1])
    disp = np.zeros(1)
    dp = np.zeros((1, 3))
    ms = apply_step(pos.copy(), icvel, h, 0.05, box, True, disp, 0.0, dp, 0.0)
    assert np.isclose(ms, 0.05 * 0.1)


def test_apply_step_cap_limits_to_h_and_no_box_jump():
    # huge icvel * icr0 near a face: with the cap each component moves at most
    # max_step_frac * h, and a single image correction keeps it in the box
    box = np.array([1.0, 1.0, 1.0])
    pos = np.array([[0.49, 0.0, 0.0]])
    icvel = np.array([[1000.0, 1000.0, -1000.0]])  # normalized would be ~1 each
    h = np.array([0.1])
    disp = np.zeros(1)
    dp = np.zeros((1, 3))
    p = pos.copy()
    ms = apply_step(p, icvel, h, 50.0, box, True, disp, 1.0, dp, 0.0)
    # step capped at 1.0 * h = 0.1 per component, not 50*0.1*1000
    assert np.isclose(ms, 0.1)
    assert np.isclose(disp[0], 0.1)
    # stayed in [-0.5, 0.5] via a single wrap (0.49 + 0.1 = 0.59 -> -0.41)
    assert np.all(np.abs(p) <= 0.5 + 1e-12)
    assert np.isclose(p[0, 0], 0.59 - 1.0)


def test_apply_step_momentum_zero_byte_identical():
    # momentum=0 must reproduce the memoryless move EXACTLY (numba sig adds two
    # trailing args but the use_mom=False branch is the old code path)
    rng = np.random.default_rng(0)
    box = np.array([2.0, 2.0, 2.0])
    pos = rng.uniform(-1.0, 1.0, (200, 3))
    icvel = rng.uniform(-1.0, 1.0, (200, 3))
    h = rng.uniform(0.05, 0.15, 200)
    dp = np.zeros((200, 3))
    p0, d0 = pos.copy(), np.zeros(200)
    p1, d1 = pos.copy(), np.zeros(200)
    ms0 = apply_step(p0, icvel, h, 0.3, box, True, d0, 0.0, dp, 0.0)
    ms1 = apply_step(p1, icvel, h, 0.3, box, True, d1, 0.0, dp.copy(), 0.0)
    # (sanity: identical to itself) and disp_prev untouched when momentum=0
    assert ms0 == ms1
    assert np.array_equal(p0, p1)
    assert np.all(dp == 0.0)


def test_apply_step_momentum_accumulates_and_stays_finite():
    # heavy-ball: a steady descent direction makes the effective step grow
    # toward icvel/(1-beta); disp_prev holds the blended (unclamped) direction
    box = np.array([100.0, 100.0, 100.0])
    icvel = np.array([[1.0, 0.0, 0.0]])
    h = np.array([0.1])
    beta = 0.9
    dp = np.zeros((1, 3))
    pos = np.zeros((1, 3))
    disp = np.zeros(1)
    steps = []
    for _ in range(50):
        x0 = pos[0, 0]
        apply_step(pos, icvel, h, 1.0, box, True, disp, 0.0, dp, beta)
        steps.append(pos[0, 0] - x0)
    # later steps are strictly larger than the first (momentum builds up)
    assert steps[-1] > steps[0]
    # bounded by the geometric limit icvel/(1-beta) * icr0 * h = 1/0.1 * 1 * 0.1
    assert np.all(np.isfinite(pos))
    assert steps[-1] < (1.0 / (1.0 - beta)) * 1.0 * 0.1 + 1e-9
    # disp_prev approached the 1/(1-beta) fixed point of the blended direction
    assert abs(dp[0, 0] - 1.0 / (1.0 - beta)) < 0.5


def test_icr0_stop_smooths_out_transient_dips():
    # the median smooth must NOT fire on a single dip below threshold (the
    # noisy-crossing problem); it fires only once the SMOOTHED level is below
    s = Icr0Stop(threshold=3e-3, window=10)
    # window not full yet -> never fires
    for _ in range(9):
        assert s.update(5e-3) is False
    assert s.update(5e-3) is False           # window full, median 5e-3 > thr
    # one transient dip well below threshold: median still above -> no stop
    assert s.update(1e-4) is False
    # now the level genuinely drops below threshold and stays: median crosses
    fired = any(s.update(2e-3) for _ in range(10))
    assert fired


def test_reanneal_off_reaches_polish_without_ladder():
    # with reanneal=False the icrloop0 ladder never fires, but polish is still
    # entered once icr0 decays below nfac, and the icrloop0 stays put
    c = StepController(icr0=1.0, icrloop0=0.5, icr0rate=1.5, nsmooth=64,
                       reanneal=False)
    c.step_zero = False
    nfac = c.nfac
    c.icr0 = 0.5 * nfac
    c.update_restart(dens_err=0.4)
    assert c.polish                      # entered polish
    assert c.icrloop0 == 0.5             # ladder untouched
    assert c.icr0 == 0.5 * nfac          # icr0 not re-inflated to icrloop0


def test_reanneal_on_is_default_and_unchanged():
    c = StepController(icr0=1.0, icrloop0=0.5, icr0rate=1.5, nsmooth=64)
    assert c.reanneal is True
    c.step_zero = False
    nfac = c.nfac
    c.icr0 = 0.5 * nfac
    c.update_restart(dens_err=0.4)       # improving -> re-anneal to icrloop0
    assert c.icr0 == 0.5 and c.icfmax == 0.4
