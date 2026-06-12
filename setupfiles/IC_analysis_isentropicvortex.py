#!/usr/bin/env python3
"""
Analysis for the isentropic (Yee) vortex test (IC_setup_isentropicvortex).

Yee, Vinokur & Djomehri (2000) JCP 162, 33. A smooth vortex perturbation on a
uniform background (rho = P = T = 1, gamma = 1.4) is an EXACT, steady solution of
the Euler equations: the centrifugal force balances the radial pressure gradient.
With a uniform mean flow (vadv, vadv) the vortex simply advects unchanged, so
after one box crossing it returns to the initial state. Any departure from the
analytic profile is numerical dissipation -> this is the canonical smooth test of
a code's low-dissipation / formal order of convergence (the headline metric is
the L1 density error).

With r the distance from the (moving) vortex centre the analytic profiles are
(IC_setup_isentropicvortex.create / getrhoi):

    v_phi(r) = (beta_v/2pi) r exp((1 - r^2)/2)         (tangential velocity)
    T(r)     = 1 - (gamma-1) beta_v^2/(8 gamma pi^2) exp(1 - r^2)
    rho(r)   = T(r)^(1/(gamma-1))                        (isentropic)
    P(r)     = rho(r)^gamma

The vortex centre starts at the origin and advects to (vadv t, vadv t); on the
periodic box [-5,5]^2 it wraps back, returning to the origin at t = 10/vadv. The
simulated v_phi is recovered after subtracting the mean flow:
    v_phi = (dx (vy - vadv) - dy (vx - vadv)) / r ,  (dx, dy) = pos - centre
(the inverse of the IC's tangential velocity construction).

Three diagnostics (mirroring IC_analysis_gresho):
  default       compare one (or more) runs at --time: radial profiles of rho and
                v_phi vs r, overlaid with the analytic profiles; prints the
                volume-weighted L1 errors of rho (the Yee metric) and v_phi.
  --convergence treat the positional inputs as the SAME test at different
                resolutions (n_x = the '<test><nx>' number in the run-folder
                name, e.g. 'isentropicvortex32...' -> 32) and plot L1(rho) vs n_x
                on log-log axes against the n_x^-2 line (ideal 2nd order).
  energy        the vortex kinetic energy E_vort(t) = sum 1/2 m_i |v_i - vadv|^2
                is always plotted vs time against its analytic (constant) value;
                a decay measures the numerical dissipation of the vortex. (The
                mean-flow KE is subtracted, so unlike gresho this cannot use the
                .log Ekin column -- it is computed per snapshot.)

Usage:
    python IC_analysis_isentropicvortex.py <run-dir> [--time 10] [--save out]
    python IC_analysis_isentropicvortex.py <dir_N32> <dir_N64> ... --convergence
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
    find_files, loaddata,
    setup_rcparams, finish_figure_with_legend,
)

# Constants fixed by IC_setup_isentropicvortex (not -a..-h parameters).
GAMMA = 1.4                     # adiabatic index (setup_isentropicvortex.gamma)
HALF = 5.0                      # domain half-width: box is [-5,5]^2 (dx = dy = 5)
BOX = 2.0 * HALF                # periodic box length (dxbound = 10)

# Default comparison time: one box crossing for vadv = 1 (nsteps*deltastep = 10),
# at which the analytic vortex has returned to the initial state. The analytic is
# steady in the co-moving frame, so any time works (the centre is tracked).
DEFAULT_TIME = 10.0

# -a/-b match runtest.sh setup_isentropicvortex(): beta_v (vortex strength) and
# vadv (mean-flow / advection velocity). Both enter the analytic (beta_v sets the
# amplitude, vadv sets the centre motion). -c/-d are the artificial viscosity
# coefficients, carried for CLI parity but unused by the analytic.
PARAM_SPEC = [
    Param("a", "beta_v", 5.0, float, "vortex strength (Yee standard = 5)"),
    Param("b", "vadv", 1.0, float, "uniform mean-flow (advection) velocity"),
]


def analytic_vphi(r, beta_v):
    """Analytic tangential velocity v_phi(r) = (beta_v/2pi) r exp((1-r^2)/2)."""
    r = np.asarray(r, dtype=float)
    return (beta_v / (2.0 * np.pi)) * r * np.exp(0.5 * (1.0 - r * r))


def analytic_T(r, beta_v):
    """Analytic temperature T(r) = 1 - (g-1) beta_v^2/(8 g pi^2) exp(1-r^2)."""
    r = np.asarray(r, dtype=float)
    dTfac = (GAMMA - 1.0) * beta_v ** 2 / (8.0 * GAMMA * np.pi ** 2)
    return 1.0 - dTfac * np.exp(1.0 - r * r)


def analytic_rho(r, beta_v):
    """Analytic density rho(r) = T(r)^(1/(gamma-1)) (isentropic background)."""
    return analytic_T(r, beta_v) ** (1.0 / (GAMMA - 1.0))


def _wrap(d):
    """Map a coordinate difference into the periodic range [-HALF, HALF)."""
    return (d + HALF) % BOX - HALF


def vortex_centre(t, vadv):
    """The vortex centre at time t: (vadv t, vadv t) wrapped into [-HALF, HALF)."""
    c = _wrap(vadv * t)
    return c, c


def profile(fn, beta_v, vadv):
    """Per-particle (r, rho, v_phi, volume) and the snapshot time for run file fn.

    r is measured from the (moving, periodically-wrapped) vortex centre;
    v_phi = (dx (vy-vadv) - dy (vx-vadv)) / r inverts the IC's tangential velocity
    construction after subtracting the mean flow. The per-particle volume
    V_i = m_i / rho_i weights the spatial L1 norm."""
    tgdata, _td, _ts, _hdr, time, _N, _ngas, _nd, _ns, _h = loaddata(fn)
    t = float(np.asarray(time).ravel()[0])
    cx, cy = vortex_centre(t, vadv)
    dx = _wrap(tgdata[:, 1] - cx)
    dy = _wrap(tgdata[:, 2] - cy)
    r = np.sqrt(dx * dx + dy * dy)
    vpx = tgdata[:, 4] - vadv
    vpy = tgdata[:, 5] - vadv
    with np.errstate(invalid="ignore", divide="ignore"):
        v_phi = np.where(r > 0, (dx * vpy - dy * vpx) / r, 0.0)
    rho = tgdata[:, 7]
    vol = tgdata[:, 0] / rho                         # m / rho (== particle volume)
    return r, rho, v_phi, vol, t


def l1_error(values, analytic, vol):
    """Volume-weighted L1 error: sum V_i |values - analytic| / sum V_i.

    Volume-weighting (V_i = m_i/rho_i) makes this a proper spatial average,
    robust to non-uniform particle masses (the varying-mass IC option)."""
    err = np.abs(values - analytic)
    return float(np.sum(vol * err) / np.sum(vol))




# --------------------------------------------------------------------------- #
#  Vortex kinetic-energy diagnostic vs time
# --------------------------------------------------------------------------- #

def vortex_energy(fn, vadv):
    """(time, E_vort) for snapshot fn, with E_vort = sum 1/2 m_i |v_i - vadv|^2.

    The mean flow (vadv, vadv, 0) is subtracted so only the vortex contributes;
    analytically constant for the steady vortex, so a decay measures numerical
    dissipation. (Because the mean flow is subtracted this cannot use the .log
    Ekin column, which includes it -- it is computed per snapshot.)"""
    tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
    t = float(np.asarray(time).ravel()[0])
    m = tgdata[:, 0]
    vpx = tgdata[:, 4] - vadv
    vpy = tgdata[:, 5] - vadv
    vpz = tgdata[:, 6]
    v2 = vpx * vpx + vpy * vpy + vpz * vpz
    return t, float(np.sum(0.5 * m * v2))


def analytic_vortex_energy(fn, beta_v, vadv):
    """Analytic E_vort on a run's particle set: sum 1/2 m_i v_phi_analytic(r_i)^2.

    Evaluated on the actual particles (positions/masses) so it matches the run's
    thin-slab geometry and resolution exactly -- the constant E_vort is conserved
    against. The centre is tracked from each particle file's own time."""
    tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
    t = float(np.asarray(time).ravel()[0])
    cx, cy = vortex_centre(t, vadv)
    dx = _wrap(tgdata[:, 1] - cx)
    dy = _wrap(tgdata[:, 2] - cy)
    r = np.sqrt(dx * dx + dy * dy)
    return float(np.sum(0.5 * tgdata[:, 0] * analytic_vphi(r, beta_v) ** 2))


def plot_energy_vs_time(inputs, labels, beta_v, vadv, save=None, ref_table=None):
    """E_vort(t) for each run against its analytic (constant) value.

    When a reference is supplied its blessed (raw) E_vort(t) is overlaid dashed."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 5))
    if ref_table is not None and "Evort_raw" in ref_table:
        ax.plot(ref_table["te_raw"], ref_table["Evort_raw"], "--", color="0.35",
                lw=1.5, zorder=5, label="reference")
    for inp, lab in zip(inputs, labels):
        files = find_files(inp)
        if not files:
            print(f"isentropicvortex: no snapshots for '{inp}'", file=sys.stderr)
            continue
        ts, ev = [], []
        for fn in files:
            t, e = vortex_energy(fn, vadv)
            ts.append(t); ev.append(e)
        o = np.argsort(ts)
        ts = np.array(ts)[o]; ev = np.array(ev)[o]
        line, = ax.plot(ts, ev, "o-", ms=3, label=lab)
        e_an = analytic_vortex_energy(files[0], beta_v, vadv)
        ax.axhline(e_an, color=line.get_color(), ls="--", lw=1.0, alpha=0.7)
    ax.set_xlabel("t")
    ax.set_ylabel(r"$E_{vort}=\sum\frac{1}{2} m_i|v_i-v_{adv}|^2$")
    ax.set_title("Yee vortex kinetic energy vs time (dashed = analytic constant)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_energy.png" if save else None))


# --------------------------------------------------------------------------- #
#  Default mode: radial profiles at --time vs analytic
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, beta_v, vadv, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (axr, axv) = plt.subplots(1, 2, figsize=(13, 6))

    # Reference overlay: the blessed binned rho(r) and v_phi(r) (dashed grey).
    if ref_table is not None and "rho" in ref_table:
        axr.plot(ref_table["r_prof"], ref_table["rho"], "--", color="0.35",
                 lw=1.6, zorder=5, label="reference")
        if "v_phi" in ref_table:
            axv.plot(ref_table["r_prof"], ref_table["v_phi"], "--", color="0.35",
                     lw=1.6, zorder=5, label="reference")

    rmax = 5.0
    l1_rho_txt, l1_v_txt = [], []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        r, rho, v_phi, vol, _time = profile(fn, beta_v, vadv)
        rmax = max(rmax, float(np.percentile(r, 99.5)))
        axr.scatter(r, rho, s=3, alpha=0.4, rasterized=True,
                    label=f"{lab} (t={tt:.3g})")
        axv.scatter(r, v_phi, s=3, alpha=0.4, rasterized=True,
                    label=f"{lab} (t={tt:.3g})")
        l1_rho = l1_error(rho, analytic_rho(r, beta_v), vol)
        l1_v = l1_error(v_phi, analytic_vphi(r, beta_v), vol)
        print(f"isentropicvortex[{lab}]: t={tt:.4g}  "
              f"L1(rho)={l1_rho:.6g}  L1(v_phi)={l1_v:.6g}")
        l1_rho_txt.append(f"{lab}: {l1_rho:.3g}")
        l1_v_txt.append(f"{lab}: {l1_v:.3g}")

    rg = np.linspace(0.0, rmax, 400)
    axr.plot(rg, analytic_rho(rg, beta_v), "k-", lw=1.8, label="analytic")
    axv.plot(rg, analytic_vphi(rg, beta_v), "k-", lw=1.8, label="analytic")

    abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8)
    if l1_rho_txt:
        axr.text(0.98, 0.04, "L1(rho)  " + ",  ".join(l1_rho_txt),
                 transform=axr.transAxes, fontsize=9, va="bottom", ha="right",
                 bbox=abox)
    if l1_v_txt:
        axv.text(0.98, 0.96, "L1(v_phi)  " + ",  ".join(l1_v_txt),
                 transform=axv.transAxes, fontsize=9, va="top", ha="right",
                 bbox=abox)

    axr.set_xlabel(r"$r$ (from vortex centre)"); axr.set_ylabel(r"$\rho$")
    axv.set_xlabel(r"$r$ (from vortex centre)"); axv.set_ylabel(r"$v_\phi$")
    fig.suptitle(f"Isentropic (Yee) vortex profiles (t={t:g}, "
                 rf"$\beta_v$={beta_v:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, axv,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Convergence mode: L1(rho) vs resolution, against n_x^-2
# --------------------------------------------------------------------------- #

def plot_convergence(inputs, t, beta_v, vadv, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    nxs, l1s = [], []
    for inp in inputs:
        nx = parse_resolution(inp)
        if nx is None:
            print(f"isentropicvortex: cannot read '<test><nx>' resolution from "
                  f"'{inp}'; skipping (needed for --convergence).", file=sys.stderr)
            continue
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        r, rho, _v_phi, vol, _time = profile(fn, beta_v, vadv)
        l1 = l1_error(rho, analytic_rho(r, beta_v), vol)
        nxs.append(nx)
        l1s.append(l1)
        print(f"isentropicvortex: n_x={nx:<5d} t={tt:.4g}  L1(rho)={l1:.6g}")

    if len(nxs) < 2:
        print("isentropicvortex: need >=2 resolutions with a valid '<test><nx>' "
              "number for a convergence plot.", file=sys.stderr)
        return

    order = np.argsort(nxs)
    nxs = np.array(nxs)[order]
    l1s = np.array(l1s)[order]

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(nxs, l1s, "o-", label=r"$L_1(\rho)$")
    ref = l1s[0] * (nxs / nxs[0]) ** (-2.0)
    ax.loglog(nxs, ref, "k--", label=r"$\propto n_x^{-2}$")
    if (ref_table is not None and "n_x" in ref_table
            and "L1_rho" in ref_table):
        ax.loglog([float(ref_table["n_x"][0])], [float(ref_table["L1_rho"][0])],
                  "D", color="C3", ms=9, mfc="none", label="reference")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$L_1(\rho)$")
    ax.set_title(f"Isentropic (Yee) vortex convergence (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_convergence.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the normalised vortex-energy series)
# --------------------------------------------------------------------------- #

def reference_table(inp, beta_v, vadv, t=None):
    """Multi-block regression table feeding ALL non-render plots:
      gate/energy x        : time grid
                  Evort    : E_vort(t)/E_vort(0)  (normalised, the gate key)
                  te_raw   : time grid
                  Evort_raw: E_vort(t)            (raw, for the energy overlay)
      profile     r_prof   : r bin centres
                  rho/v_phi: binned profiles at the comparison time
                  L1_rho   : volume-weighted L1(rho) vs analytic (1-elt array)
      convergence n_x      : the run's '<test><nx>' resolution

    E_vort (mean flow subtracted) is analytically constant for the advecting Yee
    vortex -> the normalised history is the resolution-independent CI gate. None
    if neither a usable energy series nor a profile snapshot."""
    table = {}
    try:
        files = find_files(inp)
    except FileNotFoundError:
        files = []
    if files:
        ts, ev = [], []
        for fn in files:
            tt, e = vortex_energy(fn, vadv)
            ts.append(tt); ev.append(e)
        o = np.argsort(ts)
        ts = np.array(ts)[o]; ev = np.array(ev)[o]
        if np.isfinite(ev[0]) and ev[0] != 0:
            table.update({"x": ts, "Evort": ev / ev[0],
                          "te_raw": ts, "Evort_raw": ev})
    sel = select_at_time(inp, t)
    if sel is not None:
        fn, _tt = sel
        r, rho, v_phi, vol, _time = profile(fn, beta_v, vadv)
        cx, qb = binned_profile(r, {"rho": rho, "v_phi": v_phi},
                                nbins=reg_nbins(inp), percentile=(0.5, 99.5))
        if cx is not None:
            table["r_prof"] = cx
            table["rho"] = qb["rho"]
            table["v_phi"] = qb["v_phi"]
            table["L1_rho"] = np.array([l1_error(rho, analytic_rho(r, beta_v),
                                                 vol)])
        nx = parse_resolution(inp)
        if nx is not None:
            table["n_x"] = np.array([float(nx)])
    return table or None


def _add_args(parser):
    parser.add_argument("--convergence", action="store_true",
                        help="treat inputs as different resolutions and plot "
                             "L1(rho) vs n_x against the n_x^-2 line")


def _plot(inputs, labels, params, args):
    beta_v, vadv = params["beta_v"], params["vadv"]
    ref = getattr(args, "ref_table", None)
    if args.convergence:
        plot_convergence(inputs, args.time, beta_v, vadv, save=args.save,
                         ref_table=ref)
    else:
        plot_profiles(inputs, labels, args.time, beta_v, vadv, save=args.save,
                      ref_table=ref)
    plot_energy_vs_time(inputs, labels, beta_v, vadv, save=args.save,
                        ref_table=ref)


def main():
    standalone_cli(
        "isentropicvortex", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(
            inp, params["beta_v"], params["vadv"], args.time),
        description="Isentropic (Yee) vortex analysis.",
        reg_keys=("Evort",), time_default=DEFAULT_TIME,
        time_help=f"comparison time (default {DEFAULT_TIME:g} = one box crossing "
                  f"for vadv=1; the vortex is steady and the centre is tracked)",
        add_arguments=_add_args)


if __name__ == "__main__":
    main()
