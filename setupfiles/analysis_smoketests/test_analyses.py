#!/usr/bin/env python3
"""
End-to-end smoke suite for the IC_analysis_* scripts.

For each covered test it stands up a tiny synthetic run (tests/synthdata.py),
then drives the real analysis CLI as a subprocess and asserts:

  1. a plain run exits 0 and writes its plot(s);
  2. `--save-reference` blesses a reference.json (exit 0);
  3. a normal run against that baseline + `--reg-tol` PASSES (exit 0) -- the
     bless->self-compare residual is ~0;
  4. a run whose data was perturbed away from the baseline TRIPS `--reg-tol`
     (exit 2) -- the regression gate the CI was built for.

This is the captured form of the ad-hoc synthetic validation the analyses were
originally written against, and the safety net the framework de-duplication
refactors lean on. Pure stdlib `unittest` (no pytest needed); pytest will also
collect it if installed. Headless (MPLBACKEND=Agg).

Run:  python -m unittest tests.test_analyses        (from the testsuite root)
  or: python tests/test_analyses.py
"""

import glob
import os
import subprocess
import sys
import tempfile
import unittest
import warnings

import numpy as np

# readtipsy's writers open each file with mode 'w' then 'a' without closing the
# first handle, so writing synthetic data leaks file objects -> ResourceWarning.
# That is a quirk of the shared readtipsy module, not the suite; silence it here.
warnings.filterwarnings("ignore", category=ResourceWarning)

# This suite lives in setupfiles/analysis_smoketests/, next to the
# IC_analysis_* scripts it drives (one level up).
HERE = os.path.dirname(os.path.abspath(__file__))
SETUP = os.path.dirname(HERE)
REPO = os.path.dirname(SETUP)
sys.path.insert(0, SETUP)
sys.path.insert(0, HERE)

import synthdata as sd      # noqa: E402


# --------------------------------------------------------------------------- #
#  Shared synthetic-data builders
# --------------------------------------------------------------------------- #

def mhd_snaps(times, bz=1e-3, divb=0.01, m=16):
    """A list of minimal MHD snapshots: particles on an m x m grid in the x-y
    plane (z=0), uniform Bz, a small DivB. Enough for the snapshot fallback of
    the MHD energy/divergence series AND for a non-degenerate x-y render grid;
    the .log is the preferred source for the series."""
    lin = np.linspace(-0.5, 0.5, m)
    X, Y = np.meshgrid(lin, lin, indexing="ij")
    n = X.size
    snaps = []
    for t in times:
        g = sd.gas_table(n)
        g[:, 1] = X.ravel()
        g[:, 2] = Y.ravel()
        B = np.zeros((n, 3))
        B[:, 2] = bz
        snaps.append({"time": t, "gas": g,
                      "aux": {"BField": B, "DivB": np.full(n, divb)}})
    return snaps


def lattice3d(m, half=0.5):
    """An m^3 cubic lattice of positions in [-half, half]^3 and the radius."""
    lin = np.linspace(-half, half, m)
    X, Y, Z = np.meshgrid(lin, lin, lin, indexing="ij")
    pos = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    return pos, np.sqrt(np.sum(pos ** 2, axis=1))


# --------------------------------------------------------------------------- #
#  Subprocess driver
# --------------------------------------------------------------------------- #

def run_analysis(test, run, *args):
    """Invoke setupfiles/IC_analysis_<test>.py with args; return CompletedProcess."""
    env = dict(os.environ)
    env["MPLBACKEND"] = "Agg"
    env["PYTHONPATH"] = SETUP + os.pathsep + env.get("PYTHONPATH", "")
    script = os.path.join(SETUP, "IC_analysis_%s.py" % test)
    return subprocess.run(
        [sys.executable, script, run, *map(str, args)],
        cwd=REPO, capture_output=True, text=True, env=env,
    )


def reference_json(run):
    """Path where a blessed reference.json lands for run dir `run`
    (<parent>/references/<runname>/reference.json)."""
    parent = os.path.dirname(os.path.normpath(run))
    name = os.path.basename(os.path.normpath(run))
    return os.path.join(parent, "references", name, "reference.json")


class AnalysisSmokeBase(unittest.TestCase):
    """Shared bless / gate-pass / gate-fail assertions. Subclasses set `test`
    and implement build_run(parent, diverge=False)."""

    test = None          # IC_analysis_<test>
    extra = ()           # extra CLI flags every invocation gets (e.g. params)

    def setUp(self):
        if self.test is None:
            self.skipTest("base class")
        self._tmp = tempfile.mkdtemp(prefix="smoke_%s_" % self.test)

    def tearDown(self):
        import shutil
        shutil.rmtree(getattr(self, "_tmp", ""), ignore_errors=True)

    # subclasses override -------------------------------------------------- #
    def build_run(self, parent, diverge=False):
        raise NotImplementedError

    def expected_plots(self, save_prefix):
        """A glob the plain run must have produced at least one of."""
        return save_prefix + "*.png"

    # the actual smoke test ------------------------------------------------ #
    def test_smoke(self):
        run = self.build_run(self._tmp)
        save = os.path.join(self._tmp, "out")

        # 1. plain run -> exit 0, plots produced
        r = run_analysis(self.test, run, "--save", save, *self.extra)
        self.assertEqual(r.returncode, 0,
                         msg="plain run failed:\n%s\n%s" % (r.stdout, r.stderr))
        self.assertTrue(glob.glob(self.expected_plots(save)),
                        msg="no plots produced (looked for %s)\n%s"
                            % (self.expected_plots(save), r.stdout))

        # 2. bless -> reference.json exists
        r = run_analysis(self.test, run, "--save-reference", *self.extra)
        self.assertEqual(r.returncode, 0,
                         msg="bless failed:\n%s\n%s" % (r.stdout, r.stderr))
        self.assertTrue(os.path.exists(reference_json(run)),
                        msg="no reference.json at %s" % reference_json(run))

        # 3. self-compare + reg-tol -> PASS (exit 0); residual ~0
        save2 = os.path.join(self._tmp, "cmp")
        r = run_analysis(self.test, run, "--reference", "--reg-tol", "1e-6",
                         "--save", save2, *self.extra)
        self.assertEqual(r.returncode, 0,
                         msg="self-compare gate should PASS:\n%s\n%s"
                             % (r.stdout, r.stderr))

        # 4. perturb the run away from the baseline -> reg-tol TRIPS (exit 2)
        self.build_run(self._tmp, diverge=True)     # overwrite snapshots/log
        save3 = os.path.join(self._tmp, "div")
        r = run_analysis(self.test, run, "--reference", "--reg-tol", "1e-9",
                         "--save", save3, *self.extra)
        self.assertEqual(r.returncode, 2,
                         msg="divergent run should TRIP the gate (exit 2); got "
                             "%d\n%s\n%s" % (r.returncode, r.stdout, r.stderr))

    # shared --convergence smoke check (for tests that support it) ---------- #
    def check_convergence(self, runname_tmpl):
        """Drive `--convergence` on TWO synthetic resolutions (n_x = 32, 64),
        built by build_run with distinct runnames (so parse_resolution reads
        different n_x), and assert exit 0 + a `_convergence.png`. (Identical
        physics at both n_x -> a flat L1; the point is that the convergence code
        path runs and produces output.)"""
        runs = [self.build_run(self._tmp, runname=runname_tmpl % nx)
                for nx in (32, 64)]
        save = os.path.join(self._tmp, "cv")
        r = run_analysis(self.test, runs[0], runs[1], "--convergence",
                         "--save", save, *self.extra)
        self.assertEqual(r.returncode, 0,
                         msg="--convergence failed:\n%s\n%s" % (r.stdout, r.stderr))
        self.assertTrue(glob.glob(save + "_convergence.png"),
                        msg="no _convergence.png produced\n%s" % r.stdout)


# --------------------------------------------------------------------------- #
#  kh -- run_cli metric-vs-time (mode amplitude M, top-5% Ekin)
# --------------------------------------------------------------------------- #

class TestKH(AnalysisSmokeBase):
    test = "kh"

    def build_run(self, parent, diverge=False):
        rng = np.random.default_rng(0)
        n = 1024
        amp = 0.05 * (3.0 if diverge else 1.0)    # divergence: bigger mode
        snaps = []
        for k, t in enumerate((0.0, 0.5, 1.0)):
            g = sd.gas_table(n)
            xy = rng.uniform(-0.5, 0.5, (n, 2))
            g[:, 1] = xy[:, 0]
            g[:, 2] = xy[:, 1]
            g[:, 7] = np.where(np.abs(xy[:, 1]) < 0.25, 2.0, 1.0)   # shear layer
            # seeded KH mode in vy, growing in time
            g[:, 5] = amp * np.exp(0.8 * t) * np.sin(4 * np.pi * xy[:, 0])
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, "kh32_N64latticeGASOLINE", snaps)

    def expected_plots(self, save_prefix):
        return save_prefix + "_M*.png"


# --------------------------------------------------------------------------- #
#  mhdloop -- standalone, E_mag(t)/E_mag(0) from the .log
# --------------------------------------------------------------------------- #

class TestMHDLoop(AnalysisSmokeBase):
    test = "mhdloop"

    def build_run(self, parent, diverge=False):
        n = 512
        # loop B field (z-invariant ring of radius 0.3) for the snapshot fallback
        snaps = []
        for k, t in enumerate((0.0, 10.0, 20.0)):
            g = sd.gas_table(n)
            g[:, 1] = np.linspace(-1, 1, n)
            B = np.zeros((n, 3))
            B[:, 2] = 1e-3 * (1.0 - 0.05 * k)
            snaps.append({"time": t, "gas": g, "aux": {"BField": B}})
        # dense .log Emag series (preferred source); divergence = faster decay
        emag = ([1.0, 0.95, 0.90] if not diverge else [1.0, 0.50, 0.20])
        series = dict(dTime=[0.0, 10.0, 20.0], Emag=emag)
        return sd.make_run(parent, "mhdloop16_N64latticeGASOLINE", snaps,
                           start=0, step=10, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_emag.png"


# --------------------------------------------------------------------------- #
#  sedov -- standalone, binned rho(r) reference + analytic profile + renders
# --------------------------------------------------------------------------- #

class TestSedov(AnalysisSmokeBase):
    test = "sedov"

    def build_run(self, parent, diverge=False):
        # 3-D lattice blast: a dense, fast core inside r0, ambient outside.
        m = 12
        lin = np.linspace(-0.5, 0.5, m)
        X, Y, Z = np.meshgrid(lin, lin, lin, indexing="ij")
        pos = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
        n = len(pos)
        r = np.sqrt(np.sum(pos ** 2, axis=1))
        core = r < 0.2
        g = sd.gas_table(n)
        g[:, 1:4] = pos
        rho = np.where(core, 4.0, 1.0)
        if diverge:
            rho = rho * 2.0                       # shift the binned rho(r) profile
        g[:, 7] = rho
        g[:, 8] = np.where(core, 5.0, 1.0)        # temperature
        with np.errstate(invalid="ignore"):
            vr = np.where(core, 2.0, 0.0)
            g[:, 4:7] = (vr / np.where(r > 0, r, 1.0))[:, None] * pos
        snaps = [{"time": t, "gas": g.copy()} for t in (0.02, 0.04)]
        return sd.make_run(parent, "sedov12_N64latticeselect1.0beta0.0ublast10.0GASOLINE",
                           snaps, start=20, step=20)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"

    def test_render(self):
        """--render exercises the framework render factories (ke_density,
        thermal_pressure) on real synthetic data."""
        run = self.build_run(self._tmp)
        save = os.path.join(self._tmp, "rend")
        r = run_analysis(self.test, run, "--render", "--render-backend", "grid",
                         "--save", save)
        self.assertEqual(r.returncode, 0,
                         msg="render run failed:\n%s\n%s" % (r.stdout, r.stderr))
        self.assertTrue(glob.glob(save + "*render*.png"),
                        msg="no render PNGs produced\n%s" % r.stdout)


# --------------------------------------------------------------------------- #
#  alfven -- standalone profile + energy + convergence (parse_resolution,
#  the _binned_profile wrapper, the energy-block reference, select_at_time)
# --------------------------------------------------------------------------- #

class TestAlfven(AnalysisSmokeBase):
    test = "alfven"

    def build_run(self, parent, diverge=False, runname=None):
        from IC_analysis_alfven import geometry
        runit, e_perp, _xdot0 = geometry(1)        # rotated geometry (default -a 1)
        n = 800
        rng = np.random.default_rng(1)
        pos = rng.uniform(0.0, 1.0, (n, 3))
        x_par = pos @ runit
        amp = 0.1 * (3.0 if diverge else 1.0)      # divergence: bigger wave
        wave = amp * np.sin(2 * np.pi * x_par)
        g = sd.gas_table(n)
        g[:, 1:4] = pos
        g[:, 4:7] = wave[:, None] * e_perp         # v_orth = wave
        B = wave[:, None] * e_perp                 # B_orth = wave
        snaps = [{"time": 5.0, "gas": g, "aux": {"BField": B}}]
        # CP wave: E_mag, E_kin ~ constant (energy-block reference)
        series = dict(dTime=[0.0, 2.5, 5.0], Emag=[1.0, 1.0, 1.0],
                      Ekin=[0.5, 0.5, 0.5])
        return sd.make_run(parent, runname or "alfven32_N64latticerot1GASOLINE", snaps,
                           start=50, step=10, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


# --------------------------------------------------------------------------- #
#  MHD energy / divergence series standalones (.log preferred, like mhdloop)
# --------------------------------------------------------------------------- #

    def test_convergence(self):
        self.check_convergence("alfven%d_N64latticerot1GASOLINE")

class TestOrzag(AnalysisSmokeBase):
    test = "orzag"
    # orzag's --save-reference also blesses render grids; keep the render time on
    # a snapshot we actually dumped (default times include t=2.5, which we omit).
    extra = ("--render-times", "0.5")

    def build_run(self, parent, diverge=False):
        snaps = mhd_snaps((0.0, 0.5, 1.0))
        emag = [1.0, 1.8, 1.2] if not diverge else [1.0, 0.4, 0.1]
        series = dict(dTime=[0.0, 0.5, 1.0], Emag=emag)
        return sd.make_run(parent, "orzag32_N64latticeGASOLINE", snaps,
                           start=0, step=25, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_emag.png"


class TestMHDRotor(AnalysisSmokeBase):
    test = "mhdrotor"

    def build_run(self, parent, diverge=False):
        snaps = mhd_snaps((0.0, 0.15, 0.3))
        emag = [1.0, 1.2, 1.5] if not diverge else [1.0, 3.0, 6.0]
        series = dict(dTime=[0.0, 0.15, 0.3], Emag=emag, Ekin=[1.0, 0.8, 0.6])
        return sd.make_run(parent, "mhdrotor32_N64latticeGASOLINE", snaps,
                           start=0, step=15, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_energy.png"


class TestMHDDivtest(AnalysisSmokeBase):
    test = "mhddivtest"

    def build_run(self, parent, diverge=False):
        snaps = mhd_snaps((0.0, 0.5, 1.0))
        fac = 5.0 if diverge else 1.0
        series = dict(dTime=[0.0, 0.5, 1.0],
                      divBerrAvg=[0.1 * fac, 0.05 * fac, 0.02 * fac],
                      divBerrMax=[1.0 * fac, 0.5 * fac, 0.2 * fac],
                      Emag=[1.0, 0.99, 0.98], Ekin=[0.0, 1e-4, 2e-4],
                      Eth=[1.0, 1.0, 1.0])
        return sd.make_run(parent, "mhddivtest32_N64latticeGASOLINE", snaps,
                           start=0, step=25, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_divberr*.png"


class TestCurrentsheet(AnalysisSmokeBase):
    test = "currentsheet"

    def build_run(self, parent, diverge=False):
        snaps = mhd_snaps((0.0, 2.5, 5.0))
        fac = 6.0 if diverge else 1.0
        series = dict(dTime=[0.0, 2.5, 5.0],
                      Emag=[1.0, 0.99, 0.98], Ekin=[0.0, 0.01, 0.02],
                      Eth=[1.0, 1.0, 1.0],
                      divBerrAvg=[1e-6 * fac, 1e-3 * fac, 2e-3 * fac],
                      divBerrMax=[1e-4 * fac, 1.0 * fac, 1.5 * fac])
        return sd.make_run(parent, "currentsheet32_N64latticeGASOLINE", snaps,
                           start=0, step=25, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_divberr.png"


# --------------------------------------------------------------------------- #
#  Profile-bless standalones (binned profile at --time, like sedov)
# --------------------------------------------------------------------------- #

class TestNoh(AnalysisSmokeBase):
    test = "noh"

    def build_run(self, parent, diverge=False):
        pos, r = lattice3d(12, half=0.5)
        n = len(pos)
        g = sd.gas_table(n)
        g[:, 1:4] = pos
        post = r < 0.2                              # post-shock compressed core
        rho = np.where(post, 16.0, 1.0 + 0.1 / np.where(r > 0, r, 1.0))
        if diverge:
            rho = rho * 2.0
        g[:, 7] = rho
        g[:, 8] = 1e-3                              # cold gas
        with np.errstate(invalid="ignore"):
            g[:, 4:7] = (-1.0 / np.where(r > 0, r, 1.0))[:, None] * pos   # inflow
        snaps = [{"time": t, "gas": g.copy()} for t in (0.6, 0.8)]
        return sd.make_run(parent, "noh16_N64latticed31v1.0GASOLINE", snaps,
                           start=60, step=20)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


class TestZeldovich(AnalysisSmokeBase):
    test = "zeldovich"

    def build_run(self, parent, diverge=False):
        n = 600
        q = np.linspace(-0.5, 0.5, n)
        b = 0.9                                     # pre-caustic (b<1)
        k = 2 * np.pi
        x = q - b * np.sin(k * q) / k              # Zeldovich displacement
        rho = 1.0 / (1.0 - b * np.cos(k * q))
        if diverge:
            rho = rho * 1.5
        g = sd.gas_table(n)
        g[:, 1] = x
        g[:, 7] = rho
        g[:, 4] = -np.sin(k * q) / k              # vx
        snaps = [{"time": a, "gas": g.copy()} for a in (0.3, 0.32)]
        return sd.make_run(parent, "zeldovich16_N64latticeGASOLINE", snaps,
                           start=30, step=2)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


class TestShocktube(AnalysisSmokeBase):
    test = "shocktube"

    def build_run(self, parent, diverge=False):
        n = 400
        x = np.linspace(-0.5, 0.5, n)
        left = x < 0.0                              # Sod-like L/R states
        g = sd.gas_table(n)
        g[:, 1] = x
        rho = np.where(left, 1.0, 0.125)
        if diverge:
            rho = rho * 1.5
        g[:, 7] = rho
        g[:, 4] = np.where(left, 0.0, 0.0) + 0.1 * np.exp(-(x / 0.1) ** 2)  # vx
        P = np.where(left, 1.0, 0.1)
        g[:, 8] = P / (0.4 * rho)                   # temp so P=(gamma-1) rho u
        snaps = [{"time": t, "gas": g.copy()} for t in (0.1, 0.15)]
        return sd.make_run(parent, "shocktube64_N64latticeGASOLINE", snaps,
                           start=10, step=5)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


class TestEvrard(AnalysisSmokeBase):
    test = "evrard"

    def build_run(self, parent, diverge=False):
        pos, r = lattice3d(12, half=1.0)
        keep = (r > 0.01) & (r < 1.2)
        pos, r = pos[keep], r[keep]
        n = len(pos)
        g = sd.gas_table(n)
        g[:, 1:4] = pos
        rho = 1.0 / (2 * np.pi * np.where(r > 0.02, r, 0.02))   # Evrard ~ 1/r
        if diverge:
            rho = rho * 2.0
        g[:, 7] = rho
        g[:, 8] = 0.05
        snaps = [{"time": t, "gas": g.copy()} for t in (0.8, 1.0)]
        series = dict(dTime=[0.0, 0.4, 0.8], Etot=[-0.1, -0.1, -0.1],
                      Ekin=[0.0, 0.05, 0.1], Eth=[0.05, 0.05, 0.05],
                      Epot=[-0.15, -0.2, -0.25])
        return sd.make_run(parent, "evrard16_N64latticeGASOLINE", snaps,
                           start=80, step=20, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


class TestPolytrope(AnalysisSmokeBase):
    test = "polytrope"

    def build_run(self, parent, diverge=False):
        pos, r = lattice3d(12, half=1.0)
        keep = r < 1.0
        pos, r = pos[keep], r[keep]
        n = len(pos)
        g = sd.gas_table(n)
        g[:, 1:4] = pos
        rho = np.maximum(1.0 - r ** 2, 1e-3)        # centrally concentrated star
        if diverge:
            rho = rho * 2.0
        g[:, 7] = rho
        g[:, 8] = 0.1 * rho ** (2.0 / 3.0)          # ~ K rho^(gamma-1)
        snaps = [{"time": t, "gas": g.copy()} for t in (0.0, 35.0)]
        series = dict(dTime=[0.0, 17.5, 35.0], Ekin=[1e-4, 2e-3, 1e-3],
                      Eth=[0.3, 0.3, 0.3], Epot=[-0.6, -0.6, -0.6],
                      Etot=[-0.3, -0.3, -0.3], totentrop=[1.0, 1.01, 1.02])
        return sd.make_run(parent, "polytrope20_N64latticeGASOLINE", snaps,
                           start=0, step=35, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


class TestIsentropicVortex(AnalysisSmokeBase):
    test = "isentropicvortex"

    def build_run(self, parent, diverge=False, runname=None):
        beta_v, vadv, g_ad = 5.0, 1.0, 1.4
        lin = np.linspace(-5.0, 5.0, 40)
        X, Y = np.meshgrid(lin, lin, indexing="ij")
        x, y = X.ravel(), Y.ravel()
        n = x.size
        r2 = x ** 2 + y ** 2
        f = (beta_v / (2 * np.pi)) * np.exp(0.5 * (1.0 - r2))
        T = 1.0 - (g_ad - 1.0) * beta_v ** 2 / (8 * g_ad * np.pi ** 2) * np.exp(1.0 - r2)
        rho = np.maximum(T, 1e-3) ** (1.0 / (g_ad - 1.0))
        snaps = []
        for t in (0.0, 10.0):
            # The reference is E_vort(t)/E_vort(0); to move it, the divergence
            # must change the RATIO -- perturb only the later snapshot's vortex.
            amp = 3.0 if (diverge and t > 0) else 1.0
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 2] = y
            g[:, 4] = vadv - amp * f * y            # vx = vadv - f*y
            g[:, 5] = vadv + amp * f * x            # vy = vadv + f*x
            g[:, 7] = rho
            g[:, 8] = T
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent,
                           runname or "isentropicvortex64_N64latticeb5.0vadv1.0GASOLINE",
                           snaps, start=0, step=100)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"

    def test_convergence(self):
        self.check_convergence("isentropicvortex%d_N64latticeb5.0vadv1.0GASOLINE")


# --------------------------------------------------------------------------- #
#  run_cli metric-vs-time tests (the framework reg-tol gate, like kh)
# --------------------------------------------------------------------------- #

def _grid2d(m, half=0.5):
    lin = np.linspace(-half, half, m)
    X, Y = np.meshgrid(lin, lin, indexing="ij")
    return X.ravel(), Y.ravel()

class TestSquare(AnalysisSmokeBase):
    test = "square"

    def build_run(self, parent, diverge=False):
        x, y = _grid2d(24)
        n = x.size
        v = 0.05 * (4.0 if diverge else 1.0)        # spurious Mach; scales metric
        snaps = []
        for t in (0.0, 1.0):
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 2] = y
            g[:, 7] = np.where((np.abs(x) < 0.25) & (np.abs(y) < 0.25), 4.0, 1.0)
            g[:, 4] = v * np.sin(2 * np.pi * x)     # spurious surface-tension flow
            g[:, 5] = v * np.cos(2 * np.pi * y)
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, "square32_N64latticeGASOLINE", snaps,
                           start=0, step=10)


class TestTaylorGreen(AnalysisSmokeBase):
    test = "taylorgreen"

    def build_run(self, parent, diverge=False):
        x, y = _grid2d(24)
        n = x.size
        v0 = 0.1
        snaps = []
        for t in (0.0, 1.0):
            # Ekin/TGamp are normalised to t=0; perturb only the later snapshot.
            amp = v0 * (3.0 if (diverge and t > 0) else 1.0)
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 2] = y
            g[:, 4] = amp * np.sin(2 * np.pi * x) * np.cos(2 * np.pi * y)
            g[:, 5] = -amp * np.sin(2 * np.pi * y) * np.cos(2 * np.pi * x)
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, "taylorgreen32_N64latticeGASOLINE", snaps,
                           start=0, step=10)


class TestMTI(AnalysisSmokeBase):
    test = "mti"

    def build_run(self, parent, diverge=False):
        lin_x = np.linspace(-0.05, 0.05, 16)
        lin_z = np.linspace(-0.2, 0.2, 24)
        X, Z = np.meshgrid(lin_x, lin_z, indexing="ij")
        x, z = X.ravel(), Z.ravel()
        n = x.size
        snaps = []
        for t in (0.0, 1.7):
            amp = 0.01 * (4.0 if (diverge and t > 0) else 1.0)   # KE norm'd to t=0
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 3] = z
            g[:, 6] = amp * np.sin(4 * np.pi * x / 0.1)          # seeded vz
            g[:, 8] = 1.5                                        # T
            g[:, 10] = np.where(np.abs(z) < 0.15, 0.0, 1.0)      # tracer
            B = np.zeros((n, 3))
            B[:, 0] = 5e-5                                       # weak horizontal Bx
            snaps.append({"time": t, "gas": g, "aux": {"BField": B}})
        return sd.make_run(parent, "mti32_N64latticeGASOLINE", snaps,
                           start=0, step=17)


class TestRT(AnalysisSmokeBase):
    test = "rt"

    def build_run(self, parent, diverge=False):
        lin_x = np.linspace(-0.25, 0.25, 16)
        lin_z = np.linspace(-0.5, 0.5, 32)
        X, Z = np.meshgrid(lin_x, lin_z, indexing="ij")
        x, z = X.ravel(), Z.ravel()
        n = x.size
        amp = 0.1 * (3.0 if diverge else 1.0)       # RT mode amplitude (non-norm.)
        snaps = []
        for t in (0.0, 1.0):
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 3] = z
            g[:, 7] = np.where(z > 0, 2.0, 1.0)     # heavy over light (rhodiff=2)
            g[:, 6] = amp * np.cos(8 * np.pi * x) * np.exp(-8 * np.pi * np.abs(z - 0.5))
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, "rt64_N64latticeGASOLINE", snaps,
                           start=0, step=10)


# --------------------------------------------------------------------------- #
#  More standalones: gresho (Ekin series), coldkeplerian (ring), mhdcollapse
# --------------------------------------------------------------------------- #

class TestGresho(AnalysisSmokeBase):
    test = "gresho"

    def build_run(self, parent, diverge=False, runname=None):
        x, y = _grid2d(40)
        r = np.hypot(x, y)
        vphi = np.where(r <= 0.2, 5 * r,
                        np.where(r <= 0.4, 2 - 5 * r, 0.0))    # Gresho vortex
        with np.errstate(invalid="ignore", divide="ignore"):
            ux = np.where(r > 0, -vphi * y / r, 0.0)
            uy = np.where(r > 0, vphi * x / r, 0.0)
        snaps = []
        for t in (1.0, 3.0):                        # gresho default --times 1 3
            g = sd.gas_table(len(x))
            g[:, 1] = x
            g[:, 2] = y
            g[:, 4] = ux
            g[:, 5] = uy
            snaps.append({"time": t, "gas": g})
        ekin = [1.0, 1.0, 1.0] if not diverge else [1.0, 0.5, 0.3]
        series = dict(dTime=[0.0, 1.0, 3.0], Ekin=ekin)
        return sd.make_run(parent, runname or "gresho32_N64latticeGASOLINE", snaps,
                           start=10, step=20, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


    def test_convergence(self):
        self.check_convergence("gresho%d_N64latticeGASOLINE")

class TestColdKeplerian(AnalysisSmokeBase):
    test = "coldkeplerian"

    def build_run(self, parent, diverge=False, runname=None):
        nr, nth = 12, 60
        rr = np.linspace(1.0, 2.0, nr)
        th = np.linspace(0, 2 * np.pi, nth, endpoint=False)
        R, TH = np.meshgrid(rr, th, indexing="ij")
        snaps = []
        for t in (0.0, 10.0):
            r = R.ravel().copy()
            if diverge and t > 0:                   # AV-driven spreading (later snap)
                r = 1.5 + (r - 1.5) * 1.6
            ang = TH.ravel()
            x, y = r * np.cos(ang), r * np.sin(ang)
            vphi = r ** -0.5                        # exact Keplerian
            g = sd.gas_table(r.size)
            g[:, 1] = x
            g[:, 2] = y
            g[:, 4] = -vphi * np.sin(ang)
            g[:, 5] = vphi * np.cos(ang)
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, runname or "coldkeplerian32_N64latticeGASOLINE", snaps,
                           start=0, step=100)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"


    def test_convergence(self):
        self.check_convergence("coldkeplerian%d_N64latticeGASOLINE")

class TestMHDCollapse(AnalysisSmokeBase):
    test = "mhdcollapse"

    def build_run(self, parent, diverge=False):
        pos, r = lattice3d(14, half=0.5)
        n = len(pos)
        snaps = []
        for k, t in enumerate((0.0, 1.0)):
            g = sd.gas_table(n)
            g[:, 1:4] = pos
            peak = 10.0 ** (k * (2.0 if diverge else 1.0))   # central amplification
            g[:, 7] = 1.0 + peak * np.exp(-(r / 0.1) ** 2)
            g[:, 8] = 0.1
            B = np.zeros((n, 3))
            B[:, 2] = 6e-5
            snaps.append({"time": t, "gas": g, "aux": {"BField": B}})
        series = dict(dTime=[0.0, 0.5, 1.0], Ekin=[0.1, 0.2, 0.3],
                      Epot=[-1.0, -1.5, -2.0], Emag=[0.01, 0.02, 0.03])
        return sd.make_run(parent, "mhdcollapse32_N64latticeGASOLINE", snaps,
                           start=0, step=50, log_series=series)

    def expected_plots(self, save_prefix):
        return save_prefix + "_density.png"


class TestMHDIsowave(AnalysisSmokeBase):
    test = "mhdisowave"

    def build_run(self, parent, diverge=False):
        # Longitudinal sound pulse. The reference bins vx/drho on a FIXED [-2,2]
        # grid keeping empty bins as NaN, so the particles must DENSELY cover
        # [-2,2] (all 161 bins populated) or residuals come out NaN and the gate
        # can't trip. A uniform linspace guarantees that.
        n = 1000
        x = np.linspace(-2.0, 2.0, n)
        h = 0.05
        v0 = 0.01 * (4.0 if diverge else 1.0)      # pulse amplitude -> binned vx
        B0 = np.sqrt(2.0 * 1.0 / 0.1)              # sqrt(2 P/beta), beta=0.1
        snaps = []
        for t in (0.6, 0.7):
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 4] = v0 * np.exp(-(x / (3.0 * h)) ** 2)   # vx pulse
            g[:, 9] = h                                    # SPH smoothing length
            B = np.zeros((n, 3))
            B[:, 0] = B0                                   # Bx must stay = B0
            snaps.append({"time": t, "gas": g, "aux": {"BField": B}})
        return sd.make_run(parent, "mhdisowave128_N64latticebeta0.1GASOLINE",
                           snaps, start=60, step=10)

    def expected_plots(self, save_prefix):
        return save_prefix + "_profile.png"

    def test_convergence(self):
        """--convergence: inputs are the same run at different n_x; L1(v_x) vs
        n_x is plotted (parse_resolution reads the leading '<test><nx>')."""
        runs = []
        for nx, noise in ((32, 0.01), (64, 0.0025)):     # finer -> smaller error
            n = 1000
            x = np.linspace(-2.0, 2.0, n)
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 4] = (0.01 * np.exp(-(x / (3.0 * 0.05)) ** 2)
                       + noise * np.sin(20 * x))          # pulse + n_x-dep noise
            g[:, 9] = 0.05
            B = np.zeros((n, 3))
            B[:, 0] = np.sqrt(20.0)
            runs.append(sd.make_run(
                self._tmp, "mhdisowave%d_N64latticebeta0.1GASOLINE" % nx,
                [{"time": 0.6, "gas": g, "aux": {"BField": B}}],
                start=60, step=10))
        save = os.path.join(self._tmp, "cv")
        # both inputs first (positional nargs+ is greedy), then the flags
        r = run_analysis(self.test, runs[0], runs[1], "--convergence",
                         "--save", save)
        self.assertEqual(r.returncode, 0,
                         msg="convergence run failed:\n%s\n%s" % (r.stdout, r.stderr))
        self.assertTrue(glob.glob(save + "_convergence.png"),
                        msg="no convergence plot\n%s" % r.stdout)


class TestWave(AnalysisSmokeBase):
    test = "wave"

    def build_run(self, parent, diverge=False, runname=None):
        # Sound wave (wavetype 0): only vx projects onto the eigenvector. Keep
        # P = P0 = 1/gamma (dP=0) by setting temp = P0/((gamma-1) rho), tufac=1.
        gamma = 5.0 / 3.0
        temp0 = (1.0 / gamma) / (gamma - 1.0)
        n = 400
        x = np.linspace(0.0, 1.0, n)
        snaps = []
        for t in (0.0, 1.0):
            amp = 1e-3 * (3.0 if (diverge and t > 0) else 1.0)   # rms norm'd to t=0
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 4] = amp * np.sin(2 * np.pi * x)   # vx eigen-perturbation
            g[:, 8] = temp0
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, runname or "wave32_N64latticeGASOLINE", snaps,
                           start=0, step=10)


# --------------------------------------------------------------------------- #
#  blob -- run_cli cloud-survival (mass_cloud / mass_mix, normalised to t=0)
# --------------------------------------------------------------------------- #

    def test_convergence(self):
        self.check_convergence("wave%d_N64latticeGASOLINE")

class TestBlob(AnalysisSmokeBase):
    test = "blob"

    def build_run(self, parent, diverge=False):
        # mass_cloud counts rho > rhodiff/3 (cloud rho=rhodiff=10); mass_mix counts
        # temp in a log band around T_mix = tcloud*sqrt(rhodiff). Both metrics are
        # normalised to mass_cloud(t=0), so the divergence must change the SHAPE:
        # destroy more cloud in the LATER snapshot.
        n = 400
        lin = np.linspace(-0.5, 0.5, 20)
        X, Z = np.meshgrid(lin, lin, indexing="ij")
        snaps = []
        for k, t in enumerate((0.0, 1.0)):
            if k == 0:
                n_cloud = 200                     # full cloud at t=0
            else:
                n_cloud = 30 if diverge else 80   # cloud erodes; divergence erodes more
            g = sd.gas_table(n)
            g[:, 1] = X.ravel()
            g[:, 3] = Z.ravel()
            rho = np.ones(n)                      # wind rho = RHO_OUT = 1
            rho[:n_cloud] = 10.0                  # cloud rho = rhodiff
            temp = np.full(n, 40000.0)            # wind T = tcloud*rhodiff
            temp[:n_cloud] = 4000.0               # cloud T = tcloud
            n_mix = (n - n_cloud) // 4            # some wind cooled into the mix band
            temp[n_cloud:n_cloud + n_mix] = 4000.0 * np.sqrt(10.0)   # T_mix
            g[:, 7] = rho
            g[:, 8] = temp
            snaps.append({"time": t, "gas": g})
        return sd.make_run(parent, "blob32_N64latticeGASOLINE", snaps,
                           start=0, step=10)


# --------------------------------------------------------------------------- #
#  cosmowave -- standalone amplitude-vs-scale-factor (Fourier projection)
# --------------------------------------------------------------------------- #

class TestCosmowave(AnalysisSmokeBase):
    test = "cosmowave"

    def build_run(self, parent, diverge=False):
        # case 3 (default): v_pert=vx, k=2pi, x0=-0.5, dB=dBz. The snapshot TIME is
        # the scale factor a (a_i~0.0143). Amplitudes are raw (not normalised), so
        # a uniform amplitude scale trips the gate.
        k, x0 = 2 * np.pi, -0.5
        n = 300
        x = np.linspace(-0.5, 0.5, n)
        s = np.sin(k * (x - x0))
        amp = 0.01 * (3.0 if diverge else 1.0)
        snaps = []
        for a in (0.02, 0.5, 1.0):                 # scale factors > a_i
            g = sd.gas_table(n)
            g[:, 1] = x
            g[:, 4] = amp * s                      # vx perturbation -> dv
            g[:, 7] = 1.0 + amp * s                # rho perturbation -> drho
            B = np.zeros((n, 3))
            B[:, 2] = amp * s                      # Bz perturbation -> dB (dBz)
            snaps.append({"time": a, "gas": g, "aux": {"BField": B}})
        return sd.make_run(parent, "cosmowave32_N64latticecase3GASOLINE", snaps,
                           start=2, step=50)

    def expected_plots(self, save_prefix):
        return save_prefix + "_compressible.png"


if __name__ == "__main__":
    unittest.main(verbosity=2)
