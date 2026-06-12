#!/usr/bin/env python3
"""
Analysis for the circularly-polarized Alfven wave test (IC_setup_alfven).

The wave propagates along the unit vector `runit`; for the rotated geometry
(-a/--rotate 1) runit = (1/3, 2/3, 2/3), else (1,0,0). The transverse
perturbation is measured along the wave-frame axis e_perp = (-sinb, cosb, 0):

    B_orth = B . e_perp,   v_orth = v . e_perp,   x_par = pos . runit

with the analytic profile  ampl * sin(k (x_par - x0 - v_wave t))  (ampl=0.1,
k=2pi, v_wave=1, wavelength=1 -> period 1, so t=5 is exactly 5 periods and the
analytic profile equals the initial condition). Mirrors IC_setup_alfven exactly.

Two modes:
  default       compare one (or more) runs at t=5: B_orth and v_orth vs x_par,
                each overlaid with the analytic solution; prints the L1 error of
                B_orth.
  --convergence treat the positional inputs as the SAME test at different
                resolutions (n_x = the '<test><nx>' number in the run-folder
                name, e.g. 'alfven32...' -> 32), and
                plot L1(B_orth) and L1(v_orth) vs n_x on log-log axes against
                the n_x^-2 line (ideal 2nd-order convergence).

Usage:
    python IC_analysis_alfven.py <run-dir> [-a 1] [--time 5] [--save out]
    python IC_analysis_alfven.py <dir_N32> <dir_N64> <dir_N128> --convergence
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, build_parser, params_from_args,
    reference_path_for, save_reference, find_reference, residuals, reg_nbins,
    parse_resolution, select_at_time, binned_profile,
)
from IC_analysis_general import (
    find_files, loaddata, read_bfield, read_log_series,
    setup_rcparams, finish_figure_with_legend,
)

# Wave / box constants, fixed by IC_setup_alfven (not -a..-h parameters).
AMPL = 0.1                      # transverse perturbation amplitude
BZERO = 1.0                     # background field along propagation
WAVELENGTH = 1.0
VWAVE = 1.0                     # Alfven speed B0/sqrt(rho) = 1
WK = 2.0 * np.pi / WAVELENGTH   # wavenumber
DX, DY, DZ = 1.5, 0.75, 0.75    # half-box (distribute.setbound)
X0VEC = np.array([-DX, -DY, -DZ])   # IC_setup_alfven self.x0
VOLUME = (2 * DX) * (2 * DY) * (2 * DZ)              # box volume (rho=1)
# Circular polarization -> |B|^2 = BZERO^2 + ampl^2 and |v|^2 = ampl^2
# everywhere, so the total magnetic and kinetic energies are both analytically
# constant in space AND time (M_total = rho*VOLUME = VOLUME for rho=1).
EMAG_ANALYTIC = 0.5 * (BZERO ** 2 + AMPL ** 2) * VOLUME
EKIN_ANALYTIC = 0.5 * (AMPL ** 2) * VOLUME

DEFAULT_TIME = 5.0              # compare at 5 wave periods by default

# -a matches runtest.sh setup_alfven() (runtest.sh:197) with the same default.
PARAM_SPEC = [
    Param("a", "rotate", 1, int, "rotated (1) vs cartesian (0) wave geometry"),
]


def geometry(rotate):
    """Wave-frame basis from IC_setup_alfven: propagation `runit`, transverse
    `e_perp`, and the phase offset x0.(runit)."""
    if int(rotate) == 1:
        sina, sinb = 2.0 / 3.0, 2.0 / np.sqrt(5.0)
        cosa, cosb = np.sqrt(5.0) / 3.0, 1.0 / np.sqrt(5.0)
    else:
        sina, sinb, cosa, cosb = 0.0, 0.0, 1.0, 1.0
    runit = np.array([cosa * cosb, cosa * sinb, sina])
    e_perp = np.array([-sinb, cosb, 0.0])        # wave-frame 'sin' axis (e2)
    xdot0 = float(np.dot(X0VEC, runit))
    return runit, e_perp, xdot0


def analytic_orth(x_par, t, xdot0):
    """Analytic transverse profile ampl*sin(k(x_par - x0 - v_wave t)).

    Applies to both B_orth and v_orth (rho=1, v_A=1 -> equal perturbations). At
    integer periods (t=5) this is the initial condition."""
    return AMPL * np.sin(WK * ((x_par - xdot0) - VWAVE * t))


def profile(fn, rotate):
    """Per-particle (x_par, B_orth, v_orth) and the snapshot time for run file fn."""
    runit, e_perp, xdot0 = geometry(rotate)
    tgdata, _td, _ts, _hdr, time, _N, ngas, _nd, _ns, _h = loaddata(fn)
    pos = tgdata[:, 1:4]
    vel = tgdata[:, 4:7]
    x_par = pos @ runit
    v_orth = vel @ e_perp
    B = read_bfield(fn, ngas)
    if B is None:
        print(f"alfven: no BField aux for {fn}; B_orth unavailable.",
              file=sys.stderr)
        B_orth = np.full(len(x_par), np.nan)
    else:
        B_orth = np.asarray(B, dtype=float) @ e_perp
    return x_par, B_orth, v_orth, float(np.asarray(time).ravel()[0]), xdot0


def l1_error(x_par, B_orth, t, xdot0):
    """L1 error of B_orth: mean per-particle |B_orth_sim - B_orth_analytic|."""
    return float(np.mean(np.abs(B_orth - analytic_orth(x_par, t, xdot0))))


# --------------------------------------------------------------------------- #
#  Energy diagnostics vs time
# --------------------------------------------------------------------------- #

def energies(fn):
    """(time, E_mag, E_kin) for snapshot fn.

    E_mag = sum 1/2 |B|^2 V_i   (V_i = m_i/rho_i),
    E_kin = sum 1/2 m_i |v_i|^2.
    Both are analytically constant for the circularly-polarized wave; E_mag is
    NaN if the BField aux is missing."""
    tgdata, _td, _ts, hdr, time = tip.readtipsy(fn)
    ngas = int(hdr[2])
    t = float(np.asarray(time).ravel()[0])
    m = tgdata[:, 0]
    v2 = np.sum(tgdata[:, 4:7] ** 2, axis=1)
    e_kin = float(np.sum(0.5 * m * v2))
    B = read_bfield(fn, ngas)
    if B is None:
        return t, np.nan, e_kin
    vol = m / tgdata[:, 7]                       # m / rho
    b2 = np.sum(np.asarray(B, dtype=float) ** 2, axis=1)
    return t, float(np.sum(0.5 * b2 * vol)), e_kin


def energy_series(inp):
    """(t, E_mag, E_kin, source) for run `inp`.

    Prefers the dense gasoline '.log' Emag/Ekin columns (logged every timestep),
    falling back to per-snapshot computation. Returns (None, None, None, None)
    if neither source yields data."""
    log = read_log_series(inp, ["dTime", "Emag", "Ekin"])
    if log is not None and np.any(log["Emag"] > 0):
        return log["dTime"], log["Emag"], log["Ekin"], "log"
    try:
        files = find_files(inp)
    except FileNotFoundError:
        return None, None, None, None
    ts, em, ek = [], [], []
    for fn in files:
        t, e_mag, e_kin = energies(fn)
        ts.append(t); em.append(e_mag); ek.append(e_kin)
    o = np.argsort(ts)
    return np.array(ts)[o], np.array(em)[o], np.array(ek)[o], "snapshots"


def plot_energy_vs_time(inputs, labels, save=None, ref_table=None):
    """E_mag(t) and E_kin(t) for each run, against their analytic constants.

    Uses the dense '.log' energy series when available (see energy_series). When
    a reference is supplied its blessed E_mag(t)/E_kin(t) are overlaid as a line."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axM, axK) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for inp, lab in zip(inputs, labels):
        ts, em, ek, source = energy_series(inp)
        if ts is None:
            print(f"alfven: no snapshots/.log for '{inp}'", file=sys.stderr)
            continue
        style = "-" if source == "log" else "o-"   # dense log -> line
        axM.plot(ts, em, style, ms=3, label=lab)
        axK.plot(ts, ek, style, ms=3, label=lab)
    # Reference overlay: the blessed E_mag(t)/E_kin(t) (dashed line).
    if ref_table is not None and "t_energy" in ref_table:
        te = ref_table["t_energy"]
        if "Emag" in ref_table:
            axM.plot(te, ref_table["Emag"], "--", color="C3", lw=1.6,
                     label="reference")
        if "Ekin" in ref_table:
            axK.plot(te, ref_table["Ekin"], "--", color="C3", lw=1.6,
                     label="reference")
    axM.axhline(EMAG_ANALYTIC, color="k", ls="--", lw=1.3,
                label="analytic (const)")
    axK.axhline(EKIN_ANALYTIC, color="k", ls="--", lw=1.3,
                label="analytic (const)")
    axM.set_ylabel(r"$E_{mag}=\sum\frac{1}{2}|B|^2 V_i$")
    axK.set_ylabel(r"$E_{kin}=\sum\frac{1}{2} m_i|v_i|^2$")
    axK.set_xlabel("t")
    axM.set_title("Alfven wave energy vs time")
    fig.tight_layout()
    finish_figure_with_legend(fig, (axM, axK),
                              save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Default mode: transverse profiles at t (=5) vs analytic
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, rotate, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axB, axv) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    xpar_all = []
    l1B_txt, l1v_txt = [], []
    for i, (inp, lab) in enumerate(zip(inputs, labels)):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        x_par, B_orth, v_orth, _time, xdot0 = profile(fn, rotate)
        xpar_all.append(x_par)
        color = f"C{i % 10}"
        # Faint scatter (behind) + the run's binned-mean line highlighted on top.
        axB.scatter(x_par, B_orth, s=3, alpha=0.12, color=color, rasterized=True,
                    zorder=1, label=f"{lab} (t={tt:.3g})")
        axv.scatter(x_par, v_orth, s=3, alpha=0.12, color=color, rasterized=True,
                    zorder=1, label=f"{lab} (t={tt:.3g})")
        # Binned-mean curve of THIS run (same binning used for the reference);
        # bin count tracks the run's own resolution (reg_nbins ~ n_x).
        cx, Bb, vb = _binned_profile(x_par, B_orth, v_orth,
                                     nbins=reg_nbins(inp))
        if cx is not None:
            axB.plot(cx, Bb, "-", color=color, lw=1.8, zorder=4)
            axv.plot(cx, vb, "-", color=color, lw=1.8, zorder=4)
        l1B = l1_error(x_par, B_orth, tt, xdot0)
        l1v = l1_error(x_par, v_orth, tt, xdot0)
        print(f"alfven[{lab}]: t={tt:.4g}  "
              f"L1(B_orth)={l1B:.6g}  L1(v_orth)={l1v:.6g}")
        l1B_txt.append(f"{lab}: {l1B:.3g}")
        l1v_txt.append(f"{lab}: {l1v:.3g}")

    # Reference overlay: the blessed B_perp/v_perp(x_par) curve (dashed line).
    if ref_table is not None and "B_orth" in ref_table:
        axB.plot(ref_table["x"], ref_table["B_orth"], "--", color="0.35", lw=1.6,
                 zorder=5, label="reference")
        axv.plot(ref_table["x"], ref_table["v_orth"], "--", color="0.35", lw=1.6,
                 zorder=5, label="reference")
        if "L1_B_orth" in ref_table:
            l1B_txt.append(f"ref: {float(ref_table['L1_B_orth'][0]):.3g}")
            l1v_txt.append(f"ref: {float(ref_table['L1_v_orth'][0]):.3g}")

    # Analytic overlay (single line; identical for B_orth and v_orth).
    if xpar_all:
        _, _, xdot0 = geometry(rotate)
        xg = np.linspace(min(x.min() for x in xpar_all),
                         max(x.max() for x in xpar_all), 400)
        ya = analytic_orth(xg, t, xdot0)
        axB.plot(xg, ya, "k-", lw=1.8, zorder=6, label="analytic")
        axv.plot(xg, ya, "k-", lw=1.8, zorder=6, label="analytic")

    # Annotate each panel with its L1 error(s).
    abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8)
    if l1B_txt:
        axB.text(0.02, 0.04, "L1(B_perp)  " + ",  ".join(l1B_txt),
                 transform=axB.transAxes, fontsize=9, va="bottom", bbox=abox)
        axv.text(0.02, 0.04, "L1(v_perp)  " + ",  ".join(l1v_txt),
                 transform=axv.transAxes, fontsize=9, va="bottom", bbox=abox)

    axB.set_ylabel(r"$B_\perp$")
    axv.set_ylabel(r"$v_\perp$")
    axv.set_xlabel(r"$x_\parallel = \mathbf{x}\cdot\hat{\mathbf{n}}$")
    axB.set_title(f"Alfven wave transverse profile (rotate={rotate}, t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, (axB, axv),
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Convergence mode: L1(B_orth) vs resolution, against n_x^-2
# --------------------------------------------------------------------------- #

def plot_convergence(inputs, rotate, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    nxs, l1Bs, l1vs = [], [], []
    for inp in inputs:
        nx = parse_resolution(inp)
        if nx is None:
            print(f"alfven: cannot read '<test><nx>' resolution from '{inp}'; "
                  f"skipping (needed for --convergence).", file=sys.stderr)
            continue
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        x_par, B_orth, v_orth, _time, xdot0 = profile(fn, rotate)
        l1B = l1_error(x_par, B_orth, tt, xdot0)
        l1v = l1_error(x_par, v_orth, tt, xdot0)
        nxs.append(nx)
        l1Bs.append(l1B)
        l1vs.append(l1v)
        print(f"alfven: n_x={nx:<5d} t={tt:.4g}  "
              f"L1(B_orth)={l1B:.6g}  L1(v_orth)={l1v:.6g}")

    if len(nxs) < 2:
        print("alfven: need >=2 resolutions with a valid '<test><nx>' number "
              "for a convergence plot.", file=sys.stderr)
        return

    order = np.argsort(nxs)
    nxs = np.array(nxs)[order]
    l1Bs = np.array(l1Bs)[order]
    l1vs = np.array(l1vs)[order]

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(nxs, l1Bs, "o-", label=r"$L_1(B_\perp)$")
    ax.loglog(nxs, l1vs, "s-", label=r"$L_1(v_\perp)$")
    # Ideal 2nd-order reference, anchored at the coarsest B_orth point.
    ref = l1Bs[0] * (nxs / nxs[0]) ** (-2.0)
    ax.loglog(nxs, ref, "k--", label=r"$\propto n_x^{-2}$")
    # Blessed-baseline L1 point(s) at its own resolution (if available).
    if (ref_table is not None and "n_x" in ref_table
            and "L1_B_orth" in ref_table):
        nxr = float(ref_table["n_x"][0])
        ax.loglog([nxr], [float(ref_table["L1_B_orth"][0])], "D", color="C3",
                  ms=9, mfc="none", label="reference")
        ax.loglog([nxr], [float(ref_table["L1_v_orth"][0])], "D", color="C3",
                  ms=9, mfc="none")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$L_1$")
    ax.set_title(f"Alfven wave convergence (rotate={rotate}, t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_convergence.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the transverse WAVE PROFILE + L1)
# --------------------------------------------------------------------------- #

REG_NBINS = 64                  # bins for the blessed B_perp/v_perp(x_par) curve


def _binned_profile(x_par, B_orth, v_orth, nbins=REG_NBINS):
    """Per-x_par-bin mean of B_orth and v_orth over the FULL [min,max] range, via
    the framework `binned_profile` -> (centres, B_binned, v_binned) or
    (None, None, None). The wave is a function of x_par alone, so the per-bin mean
    recovers a clean sinusoidal reference curve, and the full-extent grid lets the
    reference reach the box edges (no trimmed percentile window)."""
    cx, vals = binned_profile(x_par, {"B": B_orth, "v": v_orth}, nbins=nbins)
    if cx is None:
        return None, None, None
    return cx, vals["B"], vals["v"]


def reference_table(inp, rotate, t):
    """One reference table feeding ALL of alfven's non-render plots:

      profile block   x       : x_par grid       (the profile plot's x-axis)
                      B_orth   : B_perp(x_par)
                      v_orth   : v_perp(x_par)
      energy block    t_energy : time grid        (the energy plot's x-axis)
                      Emag     : E_mag(t)          (raw; volume integral, so
                      Ekin     : E_kin(t)           resolution-independent)
      scalars         L1_B_orth, L1_v_orth : L1 vs analytic (1-element arrays)
                      n_x                  : the run's '<test><nx>' resolution
                                             (the convergence plot's x-axis)

    Each block carries its OWN x-grid (the profile is f(x_par), the energy is
    f(t)), so the overlay/compare code picks the right grid per plot. Blocks are
    included only when computable. None if there is no usable profile snapshot."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, tt = sel
    x_par, B_orth, v_orth, _time, xdot0 = profile(fn, rotate)
    if not np.any(np.isfinite(B_orth)):
        print(f"alfven: no finite B_orth (missing BField aux?) in {inp}; "
              f"cannot bless a profile reference.", file=sys.stderr)
        return None
    cx, Bb, vb = _binned_profile(x_par, B_orth, v_orth, nbins=reg_nbins(inp))
    if cx is None:
        return None
    l1B = l1_error(x_par, B_orth, tt, xdot0)
    l1v = l1_error(x_par, v_orth, tt, xdot0)
    table = {"x": cx, "B_orth": Bb, "v_orth": vb,
             "L1_B_orth": np.array([l1B]), "L1_v_orth": np.array([l1v])}

    # Energy block (its own time grid) for the energy-vs-time overlay.
    ts, em, ek, _src = energy_series(inp)
    if ts is not None:
        table["t_energy"] = ts
        if np.any(np.isfinite(em)):
            table["Emag"] = em
        if np.any(np.isfinite(ek)):
            table["Ekin"] = ek

    nx = parse_resolution(inp)
    if nx is not None:
        table["n_x"] = np.array([float(nx)])
    return table


def _energy_subtable(ts, em, ek):
    """Repackage an energy series as a {x: t, Emag, Ekin} table so the generic
    framework `residuals` (which keys off 'x') can be reused on the time grid."""
    sub = {"x": np.asarray(ts, dtype=float)}
    if em is not None and np.any(np.isfinite(em)):
        sub["Emag"] = np.asarray(em, dtype=float)
    if ek is not None and np.any(np.isfinite(ek)):
        sub["Ekin"] = np.asarray(ek, dtype=float)
    return sub


def compare_to_baseline(inp, rotate, t, explicit=None, ref_table=None):
    """Compare a run against its sibling baseline across ALL non-render diagnostics:
    profile residuals (B_orth, v_orth on the x_par grid), energy residuals
    (Emag, Ekin on the time grid), and the run's L1 next to the baseline's L1.
    Returns the worst RMS (None if no baseline / no profile). `ref_table` may be
    passed pre-loaded to avoid re-discovering it (and re-printing the message)."""
    if ref_table is None:
        ref_table, _p = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = reference_table(inp, rotate, t)
    if table is None:
        return None
    worst = 0.0
    # Profile residuals (keyed off the shared 'x' = x_par grid).
    for key in ("B_orth", "v_orth"):
        res = residuals(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        print(f"alfven[regression] {key} profile vs baseline: "
              f"max|d|={mx:.4g}  rms={rms:.4g}")
    # Energy residuals (own time grid, via the repackaged subtables).
    if "t_energy" in table and "t_energy" in ref_table:
        in_e = _energy_subtable(table["t_energy"], table.get("Emag"),
                                table.get("Ekin"))
        ref_e = _energy_subtable(ref_table["t_energy"], ref_table.get("Emag"),
                                 ref_table.get("Ekin"))
        for key in ("Emag", "Ekin"):
            res = residuals(in_e, ref_e, key)
            if res is None:
                continue
            mx, rms = res
            worst = max(worst, rms)
            print(f"alfven[regression] {key} vs baseline: "
                  f"max|d|={mx:.4g}  rms={rms:.4g}")
    # L1 (vs analytic) of this run next to the blessed baseline's L1.
    for key, lab in (("L1_B_orth", "L1(B_perp)"), ("L1_v_orth", "L1(v_perp)")):
        if key in table and key in ref_table:
            print(f"alfven[regression] {lab}={float(table[key][0]):.6g}  "
                  f"(baseline {float(ref_table[key][0]):.6g})")
    return worst


def main():
    parser = build_parser(PARAM_SPEC,
                          description="Circularly-polarized Alfven wave analysis.")
    parser.add_argument("--times", type=float, nargs="+", default=[DEFAULT_TIME],
                        help=f"comparison time(s) in wave periods (default "
                             f"{DEFAULT_TIME:g}; the analytic = IC at integer t). "
                             f"A single comparison time is just --times <t>.")
    parser.add_argument("--convergence", action="store_true",
                        help="treat inputs as different resolutions and plot "
                             "L1(B_orth) vs n_x against the n_x^-2 line")
    parser.add_argument("--reg-tol", type=float, default=None,
                        help="CI gate: exit non-zero if the worst RMS residual "
                             "vs the regression baseline exceeds this tolerance")
    args = parser.parse_args()
    args.time = args.times[0]      # alfven compares at a single time (the first)
    params = params_from_args(args, PARAM_SPEC)
    rotate = params["rotate"]

    labels = args.labels if args.labels else [
        os.path.basename(os.path.normpath(i)) for i in args.inputs]
    if len(labels) != len(args.inputs):
        parser.error("--labels must match the number of inputs")

    # --save-reference: bless the FIRST input's wave profile + energy + L1.
    if args.save_reference:
        table = reference_table(args.inputs[0], rotate, args.time)
        if table is None:
            print("alfven: no profile to save as reference",
                  file=sys.stderr)
            sys.exit(1)
        sel0 = select_at_time(args.inputs[0], args.time)
        tref = sel0[1] if sel0 else args.time
        save_reference(table, reference_path_for(args.inputs[0]), params,
                       "alfven", time=tref)
        return

    # Auto-discover (or use --reference) the baseline once; it is overlaid on
    # every non-render plot and reused for the residual gate of the first input.
    # Pass the resolved comparison time so a time mismatch vs the baseline warns.
    sel0 = select_at_time(args.inputs[0], args.time)
    tcmp = sel0[1] if sel0 else args.time
    ref_table, _refpath = find_reference(args.inputs[0], explicit=args.reference,
                                         time=tcmp)

    # Primary figure: convergence (inputs = resolutions) or the t=5 profile.
    if args.convergence:
        plot_convergence(args.inputs, rotate, args.time, save=args.save,
                         ref_table=ref_table)
    else:
        plot_profiles(args.inputs, labels, rotate, args.time, save=args.save,
                      ref_table=ref_table)
    # The energy-vs-time diagnostic is always produced.
    plot_energy_vs_time(args.inputs, labels, save=args.save, ref_table=ref_table)

    worst = None
    for i, inp in enumerate(args.inputs):
        # Reuse the already-loaded reference for the first input; others discover
        # their own sibling baseline.
        rt = ref_table if i == 0 else None
        w = compare_to_baseline(inp, rotate, args.time, args.reference,
                                ref_table=rt)
        if w is not None:
            worst = w if worst is None else max(worst, w)
    if args.reg_tol is not None and worst is not None and worst > args.reg_tol:
        print(f"alfven: REGRESSION FAILED -- worst rms {worst:.4g} > tol "
              f"{args.reg_tol:.4g}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()
