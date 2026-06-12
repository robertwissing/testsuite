#!/usr/bin/env python3
"""
Analysis for the linear (M)HD wave test (IC_setup_wave).

The canonical Stone et al. (2008, ApJS 178, 137) convergence test: a small-
amplitude eigenmode of the linearised ideal (M)HD equations on a uniform
background. After one wave period the mode returns to the initial condition, so
the deviation measures numerical dissipation and the formal order of
convergence.

`wavetype` (-a) selects the eigenmode (propagation along x):
    0 = sound          (hydro, B=0,  a   = 1)
    1 = fast magnetosonic           (c_f = 2)
    2 = slow magnetosonic           (c_s = 1/2)
    3 = Alfven                      (c_A = 1)
MHD background: rho=1, P=1/gamma, B=(1, sqrt(2), 1/2), gamma=5/3.

The setup perturbs every primitive variable along a unit eigenvector
    R = (drho, dvx, dvy, dvz, dBy, dBz, dP),   |R| = 1,
as   q = q0 + ampl * R_q * sin(k (x - x0)).
The wavetype-agnostic metric is the projection of the simulated perturbation
onto R:
    a(x) = dU . R   with   dU = (rho-rho0, vx, vy, vz, By-By0, Bz-Bz0, P-P0),
which analytically equals  ampl * sin(k(x - x0) - k c_wave t)  (since |R|=1).
This one scalar automatically picks the variables each mode excites (drho,dvx,dP
for sound; dvz,dBz for Alfven; ...) and rejects error orthogonal to the mode.

Modes:
  default       a(x) (and the individual eigen-components) at one wave period vs
                the analytic eigenmode; prints the L1 error mean|a_sim - a_ana|.
  --convergence treat the positional inputs as the SAME wave at different
                resolutions (n_x = the '<test><nx>' number in the run-folder
                name, e.g. 'wave32...' -> 32) and plot L1(a) vs n_x log-log
                against the ideal n_x^-2 line.
The mode RMS-amplitude-vs-time diagnostic (numerical damping) is always produced.

Usage:
    python IC_analysis_wave.py <run-dir> [-a 1] [-b 1e-6] [--time T] [--save out]
    python IC_analysis_wave.py <dir_N32> <dir_N64> <dir_N128> -a 1 --convergence
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    parse_resolution, select_at_time, binned_profile, reg_nbins,
)
from IC_analysis_general import (
    find_files, loaddata, read_bfield, read_tufac,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0               # IC_setup_wave.gamma
RHO0 = 1.0                      # background density
WAVELENGTH = 1.0               # = 2*dx, dx=0.5
WK = 2.0 * np.pi / WAVELENGTH   # wavenumber k
X0 = -0.5                       # IC_setup_wave.xdot0 = -dx

# Variable order of the primitive eigenvector R (matches IC_setup_wave).
VARS = (r"$\delta\rho$", r"$\delta v_x$", r"$\delta v_y$", r"$\delta v_z$",
        r"$\delta B_y$", r"$\delta B_z$", r"$\delta P$")

# runtest.sh setup_wave() parameter order; wavetype and ampl enter the analysis.
PARAM_SPEC = [
    Param("a", "wavetype", 0, int,
          "eigenmode: 0 sound, 1 fast, 2 slow, 3 Alfven"),
    Param("b", "ampl", 1e-4, float, "eigenmode amplitude"),
]


def mhd_primitive_eigvec(rho0, Bx, By0, Bz0, gamma, P0, target):
    """Right eigenvector (real, unit norm) of the 1D ideal-MHD primitive matrix
    for the eigenvalue closest to `target`. Variable order
    (rho, vx, vy, vz, By, Bz, P). Copied verbatim from IC_setup_wave so the
    eigenmode AND its sign convention match the initial condition exactly."""
    A = np.zeros((7, 7))
    A[0, 1] = rho0
    A[1, 4] = By0 / rho0; A[1, 5] = Bz0 / rho0; A[1, 6] = 1. / rho0
    A[2, 4] = -Bx / rho0
    A[3, 5] = -Bx / rho0
    A[4, 1] = By0; A[4, 2] = -Bx
    A[5, 1] = Bz0; A[5, 3] = -Bx
    A[6, 1] = gamma * P0
    vals, vecs = np.linalg.eig(A)
    idx = np.argmin(np.abs(vals.real - target))
    R = vecs[:, idx].real
    nrm = np.linalg.norm(R)
    if nrm > 0:
        R = R / nrm
    for c in R:                                # fix overall sign -> reproducible
        if abs(c) > 1e-12:
            if c < 0:
                R = -R
            break
    return R, vals[idx].real


def background(wavetype):
    """(rho0, P0, B0=(Bx,By0,Bz0), R, c_wave) for `wavetype`, mirroring
    IC_setup_wave.create(): the linearised eigenmode and its phase speed."""
    P0 = 1.0 / GAMMA                           # -> sound speed a = 1
    if int(wavetype) == 0:
        Bx = By0 = Bz0 = 0.0
    else:
        Bx, By0, Bz0 = 1.0, np.sqrt(2.0), 0.5
    a2 = GAMMA * P0 / RHO0
    cA2 = Bx ** 2 / RHO0
    b2 = (Bx ** 2 + By0 ** 2 + Bz0 ** 2) / RHO0
    cf = np.sqrt(0.5 * ((a2 + b2) + np.sqrt((a2 + b2) ** 2 - 4. * a2 * cA2)))
    cs = np.sqrt(0.5 * ((a2 + b2) - np.sqrt(max((a2 + b2) ** 2 - 4. * a2 * cA2, 0.))))
    target = {0: np.sqrt(a2), 1: cf, 2: cs, 3: np.sqrt(cA2)}[int(wavetype)]
    R, cwave = mhd_primitive_eigvec(RHO0, Bx, By0, Bz0, GAMMA, P0, target)
    return RHO0, P0, (Bx, By0, Bz0), R, cwave


def perturbations(fn, params):
    """Per-particle (x, dU, time) for run file fn.

    dU is the (N,7) primitive-perturbation matrix (drho, dvx, dvy, dvz, dBy,
    dBz, dP) relative to the background; x is the propagation coordinate.
    P = (gamma-1) rho u with u = dTuFac * T (col 8). For the MHD modes the B
    perturbation needs the BField aux; if it is absent (e.g. a sound-wave run)
    dBy=dBz=0 (the sound eigenvector has zero B components anyway)."""
    _rho0, P0, (Bx, By0, Bz0), _R, _c = background(params["wavetype"])
    tufac = read_tufac(os.path.dirname(fn) or ".")
    tgdata, _td, _ts, _hdr, time, _N, ngas, _nd, _ns, _h = loaddata(fn)
    x = tgdata[:, 1]
    drho = tgdata[:, 7] - RHO0
    dvx, dvy, dvz = tgdata[:, 4], tgdata[:, 5], tgdata[:, 6]   # background v=0
    P = (GAMMA - 1.0) * tgdata[:, 7] * (tufac * tgdata[:, 8])
    dP = P - P0
    B = read_bfield(fn, ngas) if int(params["wavetype"]) != 0 else None
    if B is None:
        if int(params["wavetype"]) != 0:
            print(f"wave: no BField aux for {fn}; dB set to 0.", file=sys.stderr)
        dBy = np.zeros(len(x))
        dBz = np.zeros(len(x))
    else:
        B = np.asarray(B, dtype=float)
        dBy = B[:, 1] - By0
        dBz = B[:, 2] - Bz0
    dU = np.column_stack([drho, dvx, dvy, dvz, dBy, dBz, dP])
    # A linear eigenmode is zero-mean by construction, but the *theoretical*
    # background we subtract above (RHO0, P0=1/gamma, B0) can differ from the
    # simulation's actual box-mean by a tiny constant -- e.g. the dTuFac unit
    # round-trip stored to ~6 figures in the .log leaves a ~1e-4 offset in P.
    # That offset is fixed in absolute terms, so once normalised by a small ampl
    # (1e-4) it becomes a DC pedestal of order unity that swamps the signal.
    # Remove the measured DC component of each field so the perturbation is
    # genuinely zero-mean and the diagnostic is correct at any amplitude.
    dU = dU - dU.mean(axis=0)
    return x, dU, float(np.asarray(time).ravel()[0])


def eigen_amplitude(dU, R):
    """Project the per-particle perturbation onto the eigenvector: a = dU . R."""
    return dU @ R


def analytic_amp(x, t, ampl, cwave):
    """Analytic eigenmode amplitude ampl*sin(k(x - x0) - k c_wave t).

    Equals dU.R for the exact eigenmode (|R|=1); at integer wave periods
    (t = n/ c_wave) it reduces to the initial condition."""
    return ampl * np.sin(WK * (x - X0) - WK * cwave * t)


def l1_error(x, a_sim, t, ampl, cwave):
    """L1 error of the eigenmode amplitude: mean |a_sim - a_analytic|."""
    return float(np.mean(np.abs(a_sim - analytic_amp(x, t, ampl, cwave))))




# --------------------------------------------------------------------------- #
#  Default mode: eigenmode profile at one period vs analytic
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, params, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    _rho0, _P0, _B0, R, cwave = background(params["wavetype"])
    ampl = float(params["ampl"])
    fig, (axA, axC) = plt.subplots(2, 1, figsize=(8, 9), sharex=True)

    # Reference overlay: the blessed binned eigenmode amplitude a/ampl(x).
    if ref_table is not None and "a_norm" in ref_table:
        axA.plot(ref_table["x_prof"], ref_table["a_norm"], "--", color="0.35",
                 lw=1.6, zorder=5, label="reference")

    l1_txt = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        x, dU, _time = perturbations(fn, params)
        a_sim = eigen_amplitude(dU, R)
        axA.scatter(x, a_sim / ampl, s=3, alpha=0.4, rasterized=True,
                    label=f"{lab} (t={tt:.3g})")
        l1 = l1_error(x, a_sim, tt, ampl, cwave)
        print(f"wave[{lab}]: t={tt:.4g}  c_wave={cwave:.4g}  "
              f"L1(a)={l1:.6g}  L1/ampl={l1/ampl:.4g}")
        l1_txt.append(f"{lab}: {l1:.3g}")

        # Constituent perturbations (only the components this mode excites).
        for j in range(7):
            if abs(R[j]) > 1e-3:
                axC.scatter(x, dU[:, j] / ampl, s=2, alpha=0.3, rasterized=True,
                            label=f"{VARS[j]}")

    # Analytic overlay (eigenmode amplitude and each excited component).
    xg = np.linspace(X0, X0 + WAVELENGTH, 400)
    axA.plot(xg, analytic_amp(xg, t, ampl, cwave) / ampl, "k-", lw=1.8,
             label="analytic")
    for j in range(7):
        if abs(R[j]) > 1e-3:
            axC.plot(xg, R[j] * np.sin(WK * (xg - X0) - WK * cwave * t),
                     lw=1.4, label=f"{VARS[j]} (ana)")

    if l1_txt:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
        axA.text(0.02, 0.04, "L1(a)  " + ",  ".join(l1_txt),
                 transform=axA.transAxes, fontsize=9, va="bottom", bbox=abox)

    axA.set_ylabel(r"$a/\mathrm{ampl}=(\delta U\cdot R)/\mathrm{ampl}$")
    axC.set_ylabel(r"$\delta q/\mathrm{ampl}$")
    axC.set_xlabel(r"$x$ (propagation)")
    wtnames = {0: "sound", 1: "fast", 2: "slow", 3: "Alfven"}
    axA.set_title(f"Linear {wtnames.get(int(params['wavetype']), '?')} wave "
                  f"(t={t:g}, c={cwave:.3g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, (axA, axC),
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Convergence mode: L1(a) vs resolution, against n_x^-2
# --------------------------------------------------------------------------- #

def plot_convergence(inputs, params, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    _rho0, _P0, _B0, R, cwave = background(params["wavetype"])
    ampl = float(params["ampl"])
    nxs, l1s = [], []
    for inp in inputs:
        nx = parse_resolution(inp)
        if nx is None:
            print(f"wave: cannot read '<test><nx>' resolution from '{inp}'; "
                  f"skipping (needed for --convergence).", file=sys.stderr)
            continue
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        x, dU, _time = perturbations(fn, params)
        l1 = l1_error(x, eigen_amplitude(dU, R), tt, ampl, cwave)
        nxs.append(nx)
        l1s.append(l1)
        print(f"wave: n_x={nx:<5d} t={tt:.4g}  L1(a)={l1:.6g}")

    if len(nxs) < 2:
        print("wave: need >=2 resolutions with a valid '<test><nx>' number "
              "for a convergence plot.", file=sys.stderr)
        return

    order = np.argsort(nxs)
    nxs = np.array(nxs)[order]
    l1s = np.array(l1s)[order]

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(nxs, l1s, "o-", label=r"$L_1(a)$")
    ref = l1s[0] * (nxs / nxs[0]) ** (-2.0)     # ideal 2nd order at coarsest
    ax.loglog(nxs, ref, "k--", label=r"$\propto n_x^{-2}$")
    if ref_table is not None and "n_x" in ref_table and "L1_a" in ref_table:
        ax.loglog([float(ref_table["n_x"][0])], [float(ref_table["L1_a"][0])],
                  "D", color="C3", ms=9, mfc="none", label="reference")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$L_1(a)=\langle|a_{sim}-a_{ana}|\rangle$")
    wtnames = {0: "sound", 1: "fast", 2: "slow", 3: "Alfven"}
    ax.set_title(f"Linear {wtnames.get(int(params['wavetype']), '?')} wave "
                 f"convergence (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_convergence.png" if save else None))


# --------------------------------------------------------------------------- #
#  Always-on: mode RMS amplitude vs time (numerical damping)
# --------------------------------------------------------------------------- #

def plot_amplitude_vs_time(inputs, labels, params, save=None, ref_table=None):
    """RMS of the eigenmode amplitude a(x) over each snapshot vs time.

    For the exact traveling eigenmode a(x)=ampl*sin(...), so the spatial RMS is
    ampl/sqrt(2), constant in time. Decay below that line is numerical
    dissipation (wave damping); growth signals an instability. When a reference is
    supplied its blessed (raw) RMS amplitude is overlaid dashed."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    _rho0, _P0, _B0, R, _c = background(params["wavetype"])
    ampl = float(params["ampl"])
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "rms_raw" in ref_table:
        ax.plot(ref_table["t_amp"], ref_table["rms_raw"], "--", color="0.35",
                lw=1.5, zorder=5, label="reference")
    for inp, lab in zip(inputs, labels):
        files = find_files(inp)
        if not files:
            print(f"wave: no snapshots for '{inp}'", file=sys.stderr)
            continue
        ts, rms = [], []
        for fn in files:
            x, dU, t = perturbations(fn, params)
            a = eigen_amplitude(dU, R)
            ts.append(t)
            rms.append(float(np.sqrt(np.mean(a ** 2))))
        o = np.argsort(ts)
        ax.plot(np.array(ts)[o], np.array(rms)[o], "o-", ms=3, label=lab)
    ax.axhline(ampl / np.sqrt(2.0), color="k", ls="--", lw=1.3,
               label=r"analytic $\mathrm{ampl}/\sqrt{2}$")
    ax.set_xlabel("t")
    ax.set_ylabel(r"RMS$_x\,a(x)$")
    ax.set_title("Linear wave mode amplitude vs time (damping)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_amplitude.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the normalised mode-amplitude series)
# --------------------------------------------------------------------------- #

def _resolved_time(params, args):
    """The comparison time used by the profile/convergence plots: an explicit
    --time, else one wave period 1/c_wave (where the analytic equals the IC)."""
    _r, _p, _b, _R, cwave = background(params["wavetype"])
    return args.time if args.time is not None else 1.0 / abs(cwave)


def reference_table(inp, params, args):
    """Multi-block regression table feeding ALL non-render plots:
      gate/amplitude x       : time grid
                     amp     : RMS_x a(x;t)/RMS(0)  (normalised, the gate key)
                     t_amp   : time grid
                     rms_raw : RMS_x a(x;t)         (raw, the amplitude plot's y)
      profile        x_prof  : x bin centres
                     a_norm  : binned a/ampl(x) at the comparison time
                     L1_a    : L1(a) vs analytic (1-elt array)
      convergence    n_x     : the run's '<test><nx>' resolution

    The eigenmode RMS amplitude is analytically constant (ampl/sqrt2), so the
    normalised history is the resolution-independent CI gate. None if neither a
    usable amplitude series nor a profile snapshot."""
    _rho0, _P0, _B0, R, cwave = background(params["wavetype"])
    ampl = float(params["ampl"])
    table = {}
    try:
        files = find_files(inp)
    except FileNotFoundError:
        files = []
    if files:
        ts, rms = [], []
        for fn in files:
            x, dU, t = perturbations(fn, params)
            a = eigen_amplitude(dU, R)
            ts.append(t); rms.append(float(np.sqrt(np.mean(a ** 2))))
        o = np.argsort(ts)
        ts = np.array(ts)[o]; rms = np.array(rms)[o]
        if np.isfinite(rms[0]) and rms[0] != 0:
            table.update({"x": ts, "amp": rms / rms[0],
                          "t_amp": ts, "rms_raw": rms})
    tt = _resolved_time(params, args)
    sel = select_at_time(inp, tt)
    if sel is not None:
        fn, t_at = sel
        x, dU, _time = perturbations(fn, params)
        a = eigen_amplitude(dU, R)
        cx, ab = binned_profile(x, a / ampl, nbins=reg_nbins(inp))
        if cx is not None:
            table["x_prof"] = cx
            table["a_norm"] = ab
            table["L1_a"] = np.array([l1_error(x, a, t_at, ampl, cwave)])
        nx = parse_resolution(inp)
        if nx is not None:
            table["n_x"] = np.array([float(nx)])
    return table or None


def _add_args(parser):
    parser.add_argument("--convergence", action="store_true",
                        help="treat inputs as different resolutions and plot "
                             "L1(a) vs n_x against the n_x^-2 line")


def _plot(inputs, labels, params, args):
    _rho0, _P0, _B0, _R, cwave = background(params["wavetype"])
    t = _resolved_time(params, args)
    print(f"[wave] wavetype={params['wavetype']} ampl={params['ampl']} "
          f"c_wave={cwave:.4g}  requested t={t:.4g}")
    ref = getattr(args, "ref_table", None)
    if args.convergence:
        plot_convergence(inputs, params, t, save=args.save, ref_table=ref)
    else:
        plot_profiles(inputs, labels, params, t, save=args.save, ref_table=ref)
    plot_amplitude_vs_time(inputs, labels, params, save=args.save, ref_table=ref)


def main():
    standalone_cli(
        "wave", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, params,
                                                                  args),
        description="Linear (M)HD wave analysis.",
        reg_keys=("amp",), time_default=None,
        time_help="comparison time (default: one wave period T=1/c_wave, where "
                  "the analytic equals the IC)",
        add_arguments=_add_args)


if __name__ == "__main__":
    main()
