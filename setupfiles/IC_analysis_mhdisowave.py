#!/usr/bin/env python3
"""
Analysis for the isolated MHD wave test (IC_setup_mhdisowave), following
Iwasaki (2015, MNRAS 447, 4021, "The dispersive error of the tensile
instability correction in SPMHD"). A uniform magnetised medium (rho=1, P=1,
isothermal c_s=1, B=(B0,0,0) with B0=sqrt(2 P/beta)) is given a localised
longitudinal velocity pulse at the origin,

    vx(x,0) = v0 exp(-(x/(3h))^2),   v0=0.01,   vy=vz=0,   delta-rho = 0.

Because k || B the longitudinal motion does not bend the field (div B = 0 ->
delta-Bx = 0), so the pulse propagates as a pure (iso)thermal SOUND wave at
c_s = sqrt(P/rho) and does NOT excite the Alfven/fast modes. With zero initial
density perturbation the linear solution is d'Alembert: the pulse splits into
two half-amplitude copies travelling at +/- c_s,

    vx(x,t)      = 1/2 [ f(x - c_s t) + f(x + c_s t) ]
    delta-rho(x,t) = rho0/(2 c_s) [ f(x - c_s t) - f(x + c_s t) ],

with f the analytic initial pulse f(x)=v0 exp(-(x/3h)^2), h the (uniform) SPH
smoothing length read from the snapshot (col 9). The pulse is built analytically
and propagated from t=0 -- NOT taken from a snapshot, because the runs dump no
t=0 output (first dump is already evolved). A perfectly non-dispersive scheme
reproduces this exact solution; the dispersive error of the tensile-instability
correction (Iwasaki's point) shows up as the simulation pulse spreading, breaking
into smaller waves and travelling supersonically (c_emp > c_s) vs this curve.

Produces:
  * profiles at the comparison time vs the d'Alembert solution: vx(x), the
    density perturbation delta-rho(x), and Bx(x) (a robustness check -- it must
    stay = B0), with a volume-weighted L1(vx) error;
  * L1(vx) vs time (numerical dissipation/dispersion as the pulses propagate)
    and the empirical pulse speed vs the sound speed;
  * regression-baseline support (--save-reference / --reg-tol) for CI/CD.

Verified on the run mhdisowave32_N64latticebeta0.1psidec03GASOLINE (t=0.6): the
exact curve is reproduced and the data shows the expected dispersive spreading
(c_emp ~ 1.6 > c_s = 1).

Usage:
    python IC_analysis_mhdisowave.py <run-dir> [-a beta] [--time T] [--save out]
    python IC_analysis_mhdisowave.py <run-dir> --cs 1.0 --save out
    python IC_analysis_mhdisowave.py <run-dir> --save-reference        # bless
    python IC_analysis_mhdisowave.py <run-dir> --reg-tol 1e-6          # CI gate
"""

import os
import sys

import numpy as np

from IC_analysis_framework import (
    Param, standalone_cli, parse_resolution,
    select_at_time, binned_profile as _fw_binned_profile,
)
from IC_analysis_general import (
    find_files, loaddata, read_vector_aux,
    setup_rcparams, finish_figure_with_legend,
)

PRZERO = 1.0                    # IC_setup_mhdisowave.przero
RHOZERO = 1.0                   # background density
V0 = 0.01                       # pulse amplitude (IC)
XBOX = 2.0                      # half box length (dxbound/2 = 2)
REG_NBINS = 161                 # bins on the regression / initial-profile x-grid

# Isothermal run (adi=0) -> isothermal sound speed c_s = sqrt(P/rho).
CS_DEFAULT = float(np.sqrt(PRZERO / RHOZERO))

# runtest.sh setup_mhdisowave() carries only `beta` (sets B0); the analysis is
# independent of B0 (the longitudinal sound wave does not involve B).
PARAM_SPEC = [
    Param("a", "beta", 0.1, float, "plasma beta (sets B0=sqrt(2 P/beta))"),
]


# --------------------------------------------------------------------------- #
#  Per-snapshot fields
# --------------------------------------------------------------------------- #

def profile(fn):
    """Per-particle (x, vx, rho, Bx), the median smoothing length h, and the
    snapshot time. Bx is NaN-filled if the BField aux is absent (then the Bx
    panel/robustness check is skipped). h sets the analytic pulse width (3h)."""
    tgdata, _td, _ts, _hdr, time, _N, ngas, _nd, _ns, _h = loaddata(fn)
    x = tgdata[:, 1]
    vx = tgdata[:, 4]
    rho = tgdata[:, 7]
    h = float(np.median(tgdata[:, 9]))            # col 9 = SPH smoothing length
    B = read_vector_aux(fn, "BField", ngas)
    Bx = np.asarray(B, float)[:, 0] if B is not None else np.full(len(x), np.nan)
    return x, vx, rho, Bx, h, float(np.asarray(time).ravel()[0])


def binned_profile(x, q, nbins=REG_NBINS, xlim=(-XBOX, XBOX)):
    """Per-x-bin mean of `q` on the fixed [xlim] grid via the framework
    binned_profile -> (centres, mean) on the full grid (NaN where empty;
    particles outside xlim excluded)."""
    return _fw_binned_profile(x, q, nbins=nbins, lo=xlim[0], hi=xlim[1],
                              stat="mean", clip_outside=False, drop_empty=False)


def pulse(x, h):
    """Analytic initial velocity pulse f(x) = v0 exp(-(x/3h)^2) (the IC at t=0).

    Built analytically (not from a snapshot) because the runs dump no t=0
    snapshot -- the first output is already partly evolved."""
    return V0 * np.exp(-(x / (3.0 * h))**2)


def dalembert_v(xq, h, cs, t):
    """Exact (non-dispersive) vx(x,t) = 1/2[f(x - cs t) + f(x + cs t)].

    A sound wave has no dispersion, so the pulse keeps its shape and splits into
    two half-amplitude copies at +/- cs t. Iwasaki (2015) Fig. 8 'exact solution'."""
    return 0.5 * (pulse(xq - cs*t, h) + pulse(xq + cs*t, h))


def dalembert_drho(xq, h, cs, t):
    """Exact delta-rho(x,t) = rho0/(2 cs)[f(x - cs t) - f(x + cs t)]."""
    return RHOZERO / (2.0 * cs) * (pulse(xq - cs*t, h) - pulse(xq + cs*t, h))


def l1_vx(x, vx, h, cs, t):
    """Volume-uniform L1 of vx vs the exact solution (mean |.| over particles
    inside the box). Particles are ~equal volume here, so a plain mean suffices."""
    sel = (x >= -XBOX) & (x <= XBOX)
    if not np.any(sel):
        return float("nan")
    return float(np.mean(np.abs(vx[sel] - dalembert_v(x[sel], h, cs, t))))


def empirical_speed(inp):
    """Track the +x pulse peak across snapshots and fit its speed (slope of peak
    position vs time, with the pulse launched at the origin at t=0). A value > cs
    is the dispersive error Iwasaki (2015) reports for xi=1. Returns (c_emp, npts)."""
    files = find_files(inp)
    ts, xs = [], []
    for fn in files:
        x, vx, _r, _b, _h, t = profile(fn)
        if t <= 0:
            continue
        cx, mean = binned_profile(x, vx)
        right = cx > 0
        if not np.any(right) or np.all(np.isnan(mean[right])):
            continue
        peak = cx[right][np.nanargmax(np.abs(mean[right]))]
        ts.append(t); xs.append(peak)
    if len(ts) < 2:
        return float("nan"), len(ts)
    return float(np.polyfit(ts, xs, 1)[0]), len(ts)


# --------------------------------------------------------------------------- #
#  Profiles vs the d'Alembert solution
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, cs, B0, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axV, axR, axB) = plt.subplots(3, 1, figsize=(9, 12), sharex=True)
    # Reference overlay: the blessed binned vx(x) and delta-rho(x) (dashed grey).
    if ref_table is not None and "vx" in ref_table:
        axV.plot(ref_table["x"], ref_table["vx"], "--", color="0.35", lw=1.6,
                 zorder=5, label="reference")
        if "drho" in ref_table:
            axR.plot(ref_table["x"], ref_table["drho"], "--", color="0.35",
                     lw=1.6, zorder=5, label="reference")
    l1_lines = []
    have_B = False
    h_over = None
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        x, vx, rho, Bx, h, _time = profile(fn)
        if h_over is None:
            h_over = h

        # delta-rho: subtract the MEASURED mean density, not the nominal rho0=1,
        # to remove the constant SPH kernel-interpolation bias (the exact d'Alembert
        # delta-rho is zero-mean, so this is the fair comparison).
        axV.scatter(x, vx, s=3, alpha=0.4, rasterized=True, label=lab)
        axR.scatter(x, rho - float(np.mean(rho)), s=3, alpha=0.4, rasterized=True,
                    label=lab)
        if np.any(np.isfinite(Bx)):
            axB.scatter(x, Bx, s=3, alpha=0.4, rasterized=True, label=lab)
            have_B = True

        l1 = l1_vx(x, vx, h, cs, tt)        # exact solution propagated from t=0
        c_emp, npk = empirical_speed(inp)
        print(f"mhdisowave[{lab}]: t={tt:.4g}  h={h:.4g}  c_s={cs:.4g}  "
              f"L1(vx)={l1:.4g}  c_emp={c_emp:.4g} ({npk} pts)")
        l1_lines.append(f"{lab}: L1(vx)={l1:.3g}, c_emp={c_emp:.3g} (c_s={cs:.3g})")

    # Exact (non-dispersive) overlay: the analytic pulse propagated to time tt.
    if inputs and h_over is not None:
        tt = select_at_time(inputs[0], t)[1]
        xq = np.linspace(-XBOX, XBOX, 1000)
        axV.plot(xq, dalembert_v(xq, h_over, cs, tt), "k-", lw=1.6,
                 label="exact (non-dispersive)")
        axR.plot(xq, dalembert_drho(xq, h_over, cs, tt), "k-", lw=1.6,
                 label="exact (non-dispersive)")
    axB.axhline(B0, color="0.5", ls="--", lw=1.2, label=r"$B_0$")  # field must stay B0
    if not have_B:
        axB.text(0.5, 0.5, "no BField aux", ha="center", va="center",
                 transform=axB.transAxes)

    if l1_lines:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
        axV.text(0.02, 0.96, "\n".join(l1_lines), transform=axV.transAxes,
                 fontsize=8, va="top", ha="left", bbox=abox)

    axV.set_ylabel(r"$v_x$")
    axR.set_ylabel(r"$\delta\rho=\rho-\langle\rho\rangle$")
    axB.set_ylabel(r"$B_x$ (should stay $B_0$)")
    axB.set_xlabel(r"$x$")
    axB.set_xlim(-XBOX, XBOX)
    fig.suptitle(f"Isolated MHD wave: longitudinal pulse vs d'Alembert (t={t if t is not None else 'last'})")
    fig.tight_layout()
    finish_figure_with_legend(fig, axV,
                              save=(f"{save}_profile.png" if save else None))


def plot_l1_vs_time(inputs, labels, cs, save=None):
    """L1(vx) vs time -- numerical dissipation/dispersion of the pulse as it
    propagates (should stay small; grows with resolution-limited diffusion)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    for inp, lab in zip(inputs, labels):
        try:
            files = find_files(inp)
        except FileNotFoundError:
            continue
        ts, l1s = [], []
        for fn in files:
            x, vx, _r, _b, h, t = profile(fn)
            ts.append(t)
            l1s.append(l1_vx(x, vx, h, cs, t))   # exact propagated from t=0
        order = np.argsort(ts)
        ts = np.array(ts)[order]; l1s = np.array(l1s)[order]
        ax.plot(ts, l1s, "o-", ms=3, label=lab)
        print(f"mhdisowave[{lab}]: L1(vx) {l1s[0]:.3g}->{l1s[-1]:.3g} "
              f"({len(ts)} snapshots)")
    ax.set_xlabel("t")
    ax.set_ylabel(r"$L_1(v_x)$ vs d'Alembert")
    ax.set_title("Isolated MHD wave: pulse L1 error vs time (dissipation)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax, save=(f"{save}_l1.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (CI/CD)
# --------------------------------------------------------------------------- #

def reference_table(inp, t):
    """Binned vx(x) and delta-rho(x) at the comparison time on the fixed grid,
    as the regression table {x, vx, drho}. None if no snapshot."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _tt = sel
    x, vx, rho, _Bx, _h, _time = profile(fn)
    xg, vxb = binned_profile(x, vx)
    _xg, rb = binned_profile(x, rho - float(np.mean(rho)))   # remove kernel bias
    return {"x": xg, "vx": vxb, "drho": rb}


def plot_convergence(inputs, cs, t, save=None):
    """L1(vx) vs the d'Alembert solution at the comparison time, plotted vs n_x
    (inputs = the SAME run at different resolutions; n_x = the leading
    '<test><nx>' number in the folder name) on log-log axes against the ideal
    n_x^-2 line -- the low-dissipation/dispersion convergence check for the
    longitudinal sound pulse (Iwasaki 2015)."""
    import matplotlib.pyplot as plt
    nxs, l1s = [], []
    for inp in inputs:
        nx = parse_resolution(inp)
        if nx is None:
            print(f"mhdisowave: cannot read '<test><nx>' resolution from '{inp}'; "
                  f"skipping (needed for --convergence).", file=sys.stderr)
            continue
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        x, vx, _rho, _Bx, h, _time = profile(fn)
        l1s.append(l1_vx(x, vx, h, cs, tt))
        nxs.append(nx)
        print(f"mhdisowave: n_x={nx:<5d} t={tt:.4g}  L1(v_x)={l1s[-1]:.6g}")
    if len(nxs) < 2:
        print("mhdisowave: need >=2 resolutions with a valid '<test><nx>' number "
              "for a convergence plot.", file=sys.stderr)
        return
    order = np.argsort(nxs)
    nxs = np.array(nxs)[order]
    l1s = np.array(l1s)[order]
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(nxs, l1s, "o-", label=r"$L_1(v_x)$")
    ref = l1s[0] * (nxs / nxs[0]) ** (-2.0)        # ideal 2nd order at coarsest
    ax.loglog(nxs, ref, "k--", label=r"$\propto n_x^{-2}$")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$L_1(v_x)=\langle|v_{x,sim}-v_{x,exact}|\rangle$")
    ax.set_title(f"Isolated MHD wave convergence (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_convergence.png" if save else None))


def _add_args(parser):
    parser.add_argument("--cs", type=float, default=CS_DEFAULT,
                        help=f"sound speed for the d'Alembert overlay "
                             f"(default isothermal sqrt(P/rho)={CS_DEFAULT:g})")
    parser.add_argument("--convergence", action="store_true",
                        help="treat inputs as different resolutions and plot "
                             "L1(v_x) vs n_x against the n_x^-2 line")


def _plot(inputs, labels, params, args):
    print(f"[mhdisowave] beta={params['beta']}  P={PRZERO:g} rho={RHOZERO:g}  "
          f"c_s={args.cs:g}  comparison "
          f"t={args.time if args.time is not None else 'last'}")
    if args.convergence:
        plot_convergence(inputs, args.cs, args.time, save=args.save)
        return
    B0 = float(np.sqrt(2.0 * PRZERO / params["beta"]))
    plot_profiles(inputs, labels, args.time, args.cs, B0, save=args.save,
                  ref_table=getattr(args, "ref_table", None))
    plot_l1_vs_time(inputs, labels, args.cs, save=args.save)


def main():
    standalone_cli(
        "mhdisowave", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time),
        description="Isolated MHD wave (longitudinal sound pulse) analysis.",
        reg_keys=("vx", "drho"), time_default=0.6,
        time_help="comparison time (default 0.6, as in Iwasaki 2015 Fig. 8)",
        add_arguments=_add_args)


if __name__ == "__main__":
    main()
