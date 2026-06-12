#!/usr/bin/env python3
"""
Analysis for the self-gravitating polytrope test (IC_setup_polytrope). A
polytropic sphere P = K rho^gamma in hydrostatic equilibrium, built from the
Lane-Emden solution of index n = 1/(gamma-1) (n = 1.5 for gamma = 5/3), with
M = R = G = 1 and zero initial velocity. It should simply STAY in equilibrium,
so the test measures how well the scheme holds a self-gravitating hydrostatic
structure (the same exact solution SPLASH's exact_polytrope draws).

Exact solution (computed here, matching the IC):
  rho(r) = rho_c theta(xi)^n,  r = alpha xi,  alpha = R/xi_1,
  rho_c  = M / (-4 pi alpha^3 xi_1^2 theta'(xi_1)),
  K      = 4 pi G alpha^2 / ((n+1) rho_c^((1-n)/n)),
and the entropy A = P/rho^gamma = K is UNIFORM (flat) at t=0.

Produces:
  * radial profiles vs the Lane-Emden exact solution, on LINEAR r: rho(r)
    (log y), entropy A(r)=p/rho^gamma (flat = K; rise = spurious heating),
    radial velocity v_r(r) (= 0 in equilibrium; non-zero = spurious
    oscillations), with L1(rho);
  * energies vs time in SEPARATE panels (total / kinetic / thermal /
    potential) -- E_kin on a LOG scale since it is small in equilibrium;
  * total entropy vs time on a LOG scale (the gasoline `.log` totentrop
    column; numerical/AV heating shows as the upward slope);
  * the total system volume V_sys = sum_i m_i/rho_i vs time -- equilibrium
    holds it constant; a drift = the star puffing up / contracting;
  * `--render`: PARTICLE plots in a mid-plane slice (backend="particles",
    project="slice") of density, speed |v| and entropy A;
  * regression-baseline support (--save-reference / --reg-tol) for CI/CD.

Usage:
    python IC_analysis_polytrope.py <run-dir> [--time T] [--save out] [--render]
    python IC_analysis_polytrope.py <lattice> <glass> <changa> --labels ...   # multi
    python IC_analysis_polytrope.py <run-dir> --save-reference        # bless
    python IC_analysis_polytrope.py <run-dir> --reg-tol 1e-6          # CI gate
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli, find_reference,
    RenderSpec, render_panels,
    select_at_time, binned_profile, speed, entropy_function,
)
from IC_analysis_general import (
    find_files, loaddata, read_tufac, read_log_series,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0                       # IC_setup_polytrope.gamma
NPOLY = 1.0 / (GAMMA - 1.0)            # polytropic index n (= 1.5)
M, R, G = 1.0, 1.0, 1.0               # mass, radius, gravitational constant
THETA_CRIT = 0.005                     # surface cut (matches the IC)
REG_NBINS = 40                         # radial bins for the regression baseline

PARAM_SPEC = []                        # setup_polytrope carries no parameters


# --------------------------------------------------------------------------- #
#  Lane-Emden exact solution (matches the IC construction)
# --------------------------------------------------------------------------- #

def lane_emden(n=NPOLY):
    """Solve the Lane-Emden equation; return (r, rho_exact, rho_c, K, Rsurf).

    Reproduces IC_setup_polytrope: alpha=R/xi_1, rho_c from the mass condition,
    K from alpha, rho(r)=rho_c theta^n. The entropy A = K is uniform."""
    import warnings
    from scipy.integrate import odeint

    xi = np.linspace(0.0, 6.0, 6000)      # xi_1(n=1.5) ~ 3.65; 6 is ample

    def ode(y, x):
        th, dth = y
        thn = th**n if th > 0.0 else 0.0   # theta<0 past the surface -> 0 (avoid NaN)
        return [dth, (-2.0*dth/x - thn) if x > 1e-6 else -thn]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol = odeint(ode, [1.0, 0.0], xi, rtol=1e-8, atol=1e-10, mxstep=10000)
    th, dth = sol[:, 0], sol[:, 1]
    thn = np.where(th > 0.0, th, 0.0)**n        # theta<0 past surface -> 0 (no NaN)
    mask = (thn <= THETA_CRIT) | np.isnan(th)
    if np.any(mask):
        i = int(np.argmax(mask))
        xi, th, thn, dth = xi[:i], th[:i], thn[:i], dth[:i]
    xi1, dth1 = xi[-1], dth[-1]
    alpha = R / xi1
    rho_c = M / (-4.0*np.pi * alpha**3 * xi1**2 * dth1)
    K = (alpha**2 * 4.0*np.pi * G) / ((n + 1.0) * rho_c**((1.0 - n)/n))
    r = alpha * xi
    rho = rho_c * thn
    return r, rho, float(rho_c), float(K), float(r[-1])


# --------------------------------------------------------------------------- #
#  Per-snapshot fields
# --------------------------------------------------------------------------- #

def profile(fn, tufac):
    """Per-particle (r, rho, v_r, u, entropy, vol) and snapshot time.

    entropy A = p/rho^gamma = (gamma-1) u rho^(1-gamma); v_r outward radial vel."""
    tgdata, _td, _ts, _hdr, time, _N, _ng, _nd, _ns, _h = loaddata(fn)
    pos, vel = tgdata[:, 1:4], tgdata[:, 4:7]
    rho = tgdata[:, 7]
    u = tufac * tgdata[:, 8]
    r = np.sqrt(np.sum(pos**2, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        v_r = np.where(r > 0, np.sum(pos*vel, axis=1) / r, 0.0)
    entropy = (GAMMA - 1.0) * u * rho ** (1.0 - GAMMA)
    vol = tgdata[:, 0] / rho
    return (r, rho, v_r, u, entropy, vol,
            float(np.asarray(time).ravel()[0]))




def l1_vs_exact(r, field, r_ex, f_ex, vol, rmax, logspace=False):
    """Volume-weighted L1 of `field` vs the exact curve interpolated onto r,
    over particles with r < rmax (log10 space if logspace)."""
    sel = r < rmax
    if not np.any(sel):
        return float("nan")
    if logspace:
        fa = 10.0**np.interp(np.log10(np.clip(r[sel], 1e-30, None)),
                             np.log10(np.clip(r_ex, 1e-30, None)),
                             np.log10(np.clip(f_ex, 1e-30, None)))
        a = np.log10(np.clip(field[sel], 1e-30, None))
        b = np.log10(np.clip(fa, 1e-30, None))
    else:
        a = field[sel]
        b = np.interp(r[sel], r_ex, f_ex)
    return float(np.sum(vol[sel]*np.abs(a - b)) / np.sum(vol[sel]))


# --------------------------------------------------------------------------- #
#  Radial profiles vs the Lane-Emden exact solution
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, exact, save=None, reg_ref=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    r_ex, rho_ex, rho_c, K, Rsurf = exact
    fig, (axRho, axEnt, axVel) = plt.subplots(3, 1, figsize=(8, 13), sharex=True)
    axRho.set_yscale("log")               # LINEAR r (user request), log rho

    # Blessed regression-baseline overlay (the run's OWN binned profile on the
    # fixed grid; grey dashed -- distinct from the black Lane-Emden exact curve).
    if reg_ref is not None and "rho" in reg_ref:
        for ax, key in ((axRho, "rho"), (axEnt, "entropy"), (axVel, "v_r")):
            if key in reg_ref:
                ax.plot(reg_ref["x"], reg_ref[key], "--", color="0.5", lw=1.4,
                        zorder=4,
                        label=("blessed baseline" if ax is axRho else None))

    l1_lines = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        r, rho, v_r, u, entropy, vol, _time = profile(fn, read_tufac(inp))
        axRho.scatter(r, rho, s=3, alpha=0.4, rasterized=True, label=lab)
        axEnt.scatter(r, entropy, s=3, alpha=0.4, rasterized=True, label=lab)
        axVel.scatter(r, v_r, s=3, alpha=0.4, rasterized=True, label=lab)
        l1r = l1_vs_exact(r, rho, r_ex, rho_ex, vol, Rsurf, logspace=True)
        print(f"polytrope[{lab}]: t={tt:.4g}  L1(log rho)={l1r:.4g}  "
              f"<A>/K={np.median(entropy)/K:.4g}  "
              f"max|v_r|={np.nanmax(np.abs(v_r)):.3g}")
        l1_lines.append(f"{lab}: L1(log rho)={l1r:.3g}, <A>/K={np.median(entropy)/K:.3g}")

    # Exact overlays: Lane-Emden rho(r), flat entropy K, zero v_r.
    axRho.plot(r_ex, rho_ex, "k-", lw=1.8, label="Lane-Emden")
    axEnt.axhline(K, color="k", lw=1.8, label=r"$A=K$ (exact)")
    axVel.axhline(0.0, color="k", lw=1.4, label=r"$v_r=0$ (exact)")

    if l1_lines:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
        axRho.text(0.02, 0.04, "\n".join(l1_lines), transform=axRho.transAxes,
                   fontsize=8, va="bottom", ha="left", bbox=abox)

    axRho.set_ylabel(r"$\rho$")
    axEnt.set_ylabel(r"$A=p/\rho^{\gamma}$")
    axVel.set_ylabel(r"$v_r$")
    axVel.set_xlabel(r"$r=\sqrt{x^2+y^2+z^2}$")
    axVel.set_xlim(0.0, Rsurf*1.3)
    ttl = t if t is not None else "last"
    fig.suptitle(f"Polytrope (n={NPOLY:g}) profiles vs Lane-Emden (t={ttl})")
    fig.tight_layout()
    finish_figure_with_legend(fig, axRho,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Time series: energies, entropy, systemic velocity
# --------------------------------------------------------------------------- #

def _snapshot_timeseries(inp, tufac):
    """Per-snapshot (t, Ekin, Eth, Epot, Etot, totA, Vsys) sorted by time, or
    None. totA = sum m A (total entropy); Vsys = sum m/rho (system volume)."""
    try:
        files = find_files(inp)
    except FileNotFoundError:
        return None
    rows = []
    for fn in files:
        tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
        tt = float(np.asarray(time).ravel()[0])
        m, rho = tgdata[:, 0], tgdata[:, 7]
        v = tgdata[:, 4:7]
        u = tufac * tgdata[:, 8]
        ek = float(np.sum(0.5*m*np.sum(v**2, axis=1)))
        eth = float(np.sum(m*u))
        ep = float(0.5*np.sum(m*tgdata[:, 11])) if tgdata.shape[1] > 11 else np.nan
        A = (GAMMA - 1.0) * u * rho**(1.0 - GAMMA)
        totA = float(np.sum(m*A))
        vol_sys = float(np.sum(m / rho))
        rows.append((tt, ek, eth, ep, ek+eth+ep, totA, vol_sys))
    if not rows:
        return None
    return np.array(sorted(rows), float)


def plot_energy_vs_time(inputs, labels, save=None, ref_table=None):
    """Energy budget vs time, one SEPARATE panel per component (the scales
    differ wildly in equilibrium): total, kinetic (LOG y -- it is small and the
    interesting signal is its decay/oscillation), thermal, potential. When a
    reference is supplied its blessed per-component histories are overlaid."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, axarr = plt.subplots(2, 2, figsize=(13, 9), sharex=True)
    (axT, axK), (axH, axP) = axarr
    if ref_table is not None and "Etot" in ref_table:
        ts = ref_table["t_series"]
        for ax, key in ((axT, "Etot"), (axK, "Ekin"), (axH, "Eth"),
                        (axP, "Epot")):
            if key in ref_table:
                ax.plot(ts, ref_table[key], "--", color="0.5", lw=1.3, zorder=5,
                        label=("reference" if ax is axT else None))
    for inp, lab in zip(inputs, labels):
        tufac = read_tufac(inp)
        log = read_log_series(inp, ["dTime", "Ekin", "Eth", "Epot", "Etot"])
        if log is not None and np.any(log["Etot"] != 0):
            t, ek, eth, ep, et, src = (log["dTime"], log["Ekin"], log["Eth"],
                                       log["Epot"], log["Etot"], "log")
        else:
            a = _snapshot_timeseries(inp, tufac)
            if a is None:
                print(f"polytrope: no data for '{inp}'", file=sys.stderr); continue
            t, ek, eth, ep, et, src = (a[:, 0], a[:, 1], a[:, 2], a[:, 3],
                                       a[:, 4], "snapshots")
        style = "-" if src == "log" else "o-"
        axT.plot(t, et, style, ms=3, label=lab)
        axK.plot(t, ek, style, ms=3)
        axH.plot(t, eth, style, ms=3)
        axP.plot(t, ep, style, ms=3)
        drift = (et[-1]-et[0])/abs(et[0]) if et[0] else float("nan")
        print(f"polytrope[{lab}]: E_tot {et[0]:.4g}->{et[-1]:.4g} "
              f"(drift {100*drift:+.2f}%)  E_kin max {np.nanmax(ek):.3g}  ({src})")
    axK.set_yscale("log")
    axT.set_title("total"); axK.set_title("kinetic (log)")
    axH.set_title("thermal"); axP.set_title("potential")
    for ax in (axT, axK, axH, axP):
        ax.set_ylabel("energy")
    for ax in (axH, axP):
        ax.set_xlabel("t")
    fig.suptitle("Polytrope energy budget vs time (equilibrium: E_kin~0)")
    fig.tight_layout()
    finish_figure_with_legend(fig, [axT, axK, axH, axP],
                              save=(f"{save}_energy.png" if save else None))


def plot_entropy_vs_time(inputs, labels, save=None, ref_table=None):
    """Total entropy A vs time on a LOG scale (from the gasoline `.log`
    totentrop column when present). A perfect adiabatic equilibrium holds it
    constant; numerical/AV heating shows as the upward slope. When a reference is
    supplied its blessed total-entropy history is overlaid (dashed)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "totentrop" in ref_table:
        ax.plot(ref_table["t_series"], ref_table["totentrop"], "--", color="0.5",
                lw=1.4, zorder=5, label="reference")
    for inp, lab in zip(inputs, labels):
        log = read_log_series(inp, ["dTime", "totentrop"])
        if log is not None and np.any(log["totentrop"] != 0):
            t, s, src = log["dTime"], log["totentrop"], "log"
        else:
            a = _snapshot_timeseries(inp, read_tufac(inp))
            if a is None:
                continue
            t, s, src = a[:, 0], a[:, 5], "snapshots"
        ax.plot(t, s, "-" if src == "log" else "o-", ms=3, label=lab)
        s0 = s[0] if s[0] else 1.0
        print(f"polytrope[{lab}]: entropy/entropy0 {s[0]/s0:.4g}->{s[-1]/s0:.4g} ({src})")
    ax.set_yscale("log")
    ax.set_xlabel("t"); ax.set_ylabel(r"total entropy $A$ (log)")
    ax.set_title("Polytrope total entropy vs time (numerical heating)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax, save=(f"{save}_entropy.png" if save else None))


def plot_volume_vs_time(inputs, labels, save=None, ref_table=None):
    """Total system volume V_sys = sum_i m_i/rho_i vs time -- the equilibrium
    star keeps it constant; a rise = the star puffing up (spurious heating /
    surface evaporation), a fall = over-contraction. When a reference is supplied
    its blessed V_sys history is overlaid (dashed)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, ax = plt.subplots(figsize=(8, 6))
    if ref_table is not None and "Vsys" in ref_table:
        ax.plot(ref_table["t_series"], ref_table["Vsys"], "--", color="0.5",
                lw=1.4, zorder=5, label="reference")
    for inp, lab in zip(inputs, labels):
        a = _snapshot_timeseries(inp, read_tufac(inp))
        if a is None:
            print(f"polytrope: no snapshots for '{inp}'", file=sys.stderr); continue
        t, vol = a[:, 0], a[:, 6]
        ax.plot(t, vol, "o-", ms=3, label=lab)
        print(f"polytrope[{lab}]: V_sys {vol[0]:.4g}->{vol[-1]:.4g} "
              f"(rel. drift {100*(vol[-1]-vol[0])/vol[0]:+.2f}%)")
    ax.set_xlabel("t"); ax.set_ylabel(r"$V_{\rm sys}=\sum_i m_i/\rho_i$")
    ax.set_title("Polytrope system volume vs time (should stay constant)")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax, save=(f"{save}_volume.png" if save else None))


# --------------------------------------------------------------------------- #
#  Particle renders (--render): mid-plane slice of density, speed, entropy
# --------------------------------------------------------------------------- #

def do_render(inputs, labels, args, save):
    """Particle plots in the mid-plane slice (backend='particles',
    project='slice'): density (log), speed |v|, entropy A. The default extent
    is the star (R=1) with a small margin; --render-* flags forwarded."""
    entries = list(zip(labels, inputs))
    ext = (-1.2, 1.2, -1.2, 1.2)
    specs = [
        (RenderSpec(1, 2, 7, r"$\rho$", True, "slice", None, "inferno",
                    False, ext), "rho"),
        (RenderSpec(1, 2, speed(), r"$|v|$", False, "slice", None, "viridis",
                    False, ext), "speed"),
        (RenderSpec(1, 2, entropy_function(GAMMA), r"$A=p/\rho^\gamma$", False,
                    "slice", None, "magma", False, ext), "entropy"),
    ]
    for sp, slug in specs:
        render_panels(entries, sp, times=args.render_times, n=args.render_n,
                      res=args.render_res,
                      save=(f"{save}_{slug}" if save else None),
                      backend="particles", project="slice",
                      slice_frac=args.render_slice_frac,
                      extent=args.render_extent, aspect=args.render_aspect,
                      nsmooth=args.render_nsmooth)


# --------------------------------------------------------------------------- #
#  Regression baseline (CI/CD)
# --------------------------------------------------------------------------- #

def reference_table(inp, t):
    """Multi-block regression table:
      profile  x        : fixed log-r grid
               rho/entropy/v_r : binned median profiles (the gate keys, NaN-grid)
      series   t_series : time grid
               Etot/Ekin/Eth/Epot : energy budget (the 2x2 energy plot)
               totentrop          : total entropy (the entropy plot)
               Vsys               : system volume (the volume plot)
    None if no profile snapshot."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _tt = sel
    r, rho, v_r, _u, entropy, _vol, _time = profile(fn, read_tufac(inp))
    edges = np.geomspace(1e-2, 1.0, REG_NBINS + 1)
    centres, vals = binned_profile(
        r, {"rho": rho, "entropy": entropy, "v_r": v_r},
        edges=edges, stat="median", clip_outside=False, drop_empty=False,
        geom_centres=True)
    table = {"x": centres, "rho": vals["rho"], "entropy": vals["entropy"],
             "v_r": vals["v_r"]}
    # Time-series block (own time grid) for the energy/entropy/volume overlays.
    a = _snapshot_timeseries(inp, read_tufac(inp))
    if a is not None:
        table["t_series"] = a[:, 0]
        table["Etot"] = a[:, 4]
        table["Ekin"] = a[:, 1]
        table["Eth"] = a[:, 2]
        table["Epot"] = a[:, 3]
        table["totentrop"] = a[:, 5]
        table["Vsys"] = a[:, 6]
    return table


def _reg_resid(table, ref, key):
    """(max_abs, rms) of table vs reference for `key`, over bins finite in BOTH.

    Both use the same fixed radial grid, so compare element-wise -- and a plain
    np.interp residual would be poisoned by the empty (NaN) radial bins."""
    if key not in table or key not in ref:
        return None
    a, b = np.asarray(table[key]), np.asarray(ref[key])
    if a.shape != b.shape:
        return None
    m = np.isfinite(a) & np.isfinite(b)
    if not np.any(m):
        return None
    d = a[m] - b[m]
    return float(np.max(np.abs(d))), float(np.sqrt(np.mean(d**2)))


def compare_to_baseline(inp, t, explicit=None):
    ref_table, _p = find_reference(inp, explicit)
    if ref_table is None:
        return None
    table = reference_table(inp, t)
    if table is None:
        return None
    worst = 0.0
    for key in ("rho", "entropy", "v_r"):
        res = _reg_resid(table, ref_table, key)
        if res is None:
            continue
        mx, rms = res
        worst = max(worst, rms)
        print(f"polytrope[regression] {key:<8s} vs baseline: "
              f"max|d|={mx:.4g}  rms={rms:.4g}")
    return worst


def _plot(inputs, labels, params, args):
    exact = lane_emden()
    _r, _rho, rho_c, K, Rsurf = exact
    print(f"[polytrope] gamma={GAMMA:.4g} n={NPOLY:g}  M={M:g} R={R:g} G={G:g}  "
          f"rho_c={rho_c:.4g} K={K:.4g} R_surf={Rsurf:.4g}")
    ref = getattr(args, "ref_table", None)
    plot_profiles(inputs, labels, args.time, exact, save=args.save, reg_ref=ref)
    plot_energy_vs_time(inputs, labels, save=args.save, ref_table=ref)
    plot_entropy_vs_time(inputs, labels, save=args.save, ref_table=ref)
    plot_volume_vs_time(inputs, labels, save=args.save, ref_table=ref)
    if args.render:
        do_render(inputs, labels, args, save=args.save)


def main():
    standalone_cli(
        "polytrope", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time),
        compare=lambda inp, params, args: compare_to_baseline(inp, args.time,
                                                              args.reference),
        description="Self-gravitating polytrope analysis (vs Lane-Emden).",
        time_default=None,
        time_help="comparison time for the profiles (default: last)")


if __name__ == "__main__":
    main()
