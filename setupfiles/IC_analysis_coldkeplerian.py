#!/usr/bin/env python3
"""
Analysis for the COLD KEPLERIAN DISK / RING test (IC_setup_coldkeplerian).

Cullen & Dehnen (2010) MNRAS 408, 669; Hopkins (2015) MNRAS 450, 53. A cold
(near-pressureless), narrow ring of gas orbits a central point mass (a single
sink, mass = 1, GM = 1) in EXACT Keplerian rotation. The analytic solution is
steady rotation -- the ring should not change. A poorly behaved artificial
viscosity (AV) instead acts as a spurious *shear* viscosity in the differential
rotation, transporting angular momentum so the ring spreads radially and
eventually breaks into clumps. The test therefore isolates the AV switch and the
code's angular-momentum conservation.

The exact velocity field (point mass GM = 1) is purely azimuthal Keplerian,

    v_phi(r) = r^(-1/2),     v_r = 0,

independent of the density power-law -d (q sets the surface-density slope and the
thermal energy, NOT the orbital speed -- IC_setup_coldkeplerian.getveli hardcodes
v_azi = r^(-1/2)). The simulated v_phi is recovered from the Cartesian velocity by
v_phi = (x vy - y vx)/r_cyl and v_r = (x vx + y vy)/r_cyl (the inverse of the
cyl->cart transform in the setup).

Diagnostics:
  default        compare one (or more) runs at --time: v_phi(r_cyl), v_r(r_cyl)
                 and rho(r_cyl) scattered over the ring particles, overlaid with
                 the analytic Keplerian v_phi and the cold v_r = 0; prints the
                 (volume-weighted) L1 error of v_phi and the RMS of v_r (the
                 spurious radial flow driven by the AV).
  always         ring diagnostics vs time (2 panels):
                   * sigma_r(t)/sigma_r(0): mass-weighted radial spread of the
                     ring (ideal = 1; a rise is AV-driven spreading / AM
                     transport -- the headline metric of this test),
                   * L_z(t)/L_z(0): total z angular momentum (ideal = 1; this is
                     conserved exactly by symmetric pairwise forces, so a drift
                     measures the integrator's conservation error).
  --render       face-on (x-y) density map of the ring at --time, the canonical
                 Cullen-Dehnen figure (an intact ring vs a spread / clumped one).

Usage:
    python IC_analysis_coldkeplerian.py <run-dir> [--time 10] [--save out]
    python IC_analysis_coldkeplerian.py <run-dir> --render --save out
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    RenderSpec, render_panels,
    parse_resolution, select_at_time, binned_profile, reg_nbins,
)
from IC_analysis_general import (
    find_files, loaddata,
    setup_rcparams, finish_figure_with_legend,
)

GM = 1.0                # central sink mass (IC_setup_coldkeplerian.create_star)
DEFAULT_TIME = 10.0     # ~one orbit at r=1 every 2*pi; deltastep=2pi/10, default
                        # comparison late enough for AV spreading to show, but the
                        # solution is steady so any time is valid.

# runtest.sh setup_coldkeplerian() parameter order (runtest.sh:819). Only the ring
# radii (rinner/rdisk) bound the analytic window; hr/q are carried for
# labelling -- none of them change the Keplerian v_phi = r^(-1/2).
PARAM_SPEC = [
    Param("a", "hr", 0.05, float,
          "aspect ratio hr (=1/Mach; small=cold; sets thermal/vertical scale)"),
    Param("b", "rinner", 1.0, float, "inner ring radius"),
    Param("c", "rdisk", 2.0, float, "outer ring radius"),
    Param("d", "q", 0.5, float,
          "surface-density power-law (v_phi is Keplerian r^-1/2 regardless)"),
]


def analytic_vphi(r):
    """Analytic Keplerian azimuthal velocity v_phi(r) = sqrt(GM/r), GM = 1."""
    r = np.asarray(r, dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(r > 0, np.sqrt(GM / r), 0.0)


# --------------------------------------------------------------------------- #
#  Per-snapshot kinematics
# --------------------------------------------------------------------------- #

def kinematics(fn):
    """Per-particle (r_cyl, v_phi, v_r, rho, vol, m) and the snapshot time.

    v_phi = (x vy - y vx)/r_cyl and v_r = (x vx + y vy)/r_cyl invert the IC's
    cyl->cart velocity transform; V_i = m_i/rho_i is the per-particle volume that
    weights the spatial L1 norm."""
    tgdata, _td, _ts, _hdr, time, _N, _ngas, _nd, _ns, _h = loaddata(fn)
    x, y = tgdata[:, 1], tgdata[:, 2]
    vx, vy = tgdata[:, 4], tgdata[:, 5]
    m, rho = tgdata[:, 0], tgdata[:, 7]
    rcyl = np.sqrt(x * x + y * y)
    with np.errstate(invalid="ignore", divide="ignore"):
        v_phi = np.where(rcyl > 0, (x * vy - y * vx) / rcyl, 0.0)
        v_r = np.where(rcyl > 0, (x * vx + y * vy) / rcyl, 0.0)
        vol = np.where(rho > 0, m / rho, 0.0)
    return rcyl, v_phi, v_r, rho, vol, m, float(np.asarray(time).ravel()[0])


def ring_mask(rcyl, rinner, rdisk):
    """Particles inside the nominal ring window [0.5*rinner, 2*rdisk].

    Widened beyond [rinner, rdisk] so the L1 error and spread still include
    particles that have migrated out of the original ring (the AV signal),
    while excluding the r->0 region where v_phi diverges."""
    return (rcyl >= 0.5 * rinner) & (rcyl <= 2.0 * rdisk)


def l1_vphi(rcyl, v_phi, vol, mask):
    """Volume-weighted L1 error of v_phi vs the Keplerian analytic, over `mask`."""
    err = np.abs(v_phi[mask] - analytic_vphi(rcyl[mask]))
    w = vol[mask]
    return float(np.sum(w * err) / np.sum(w)) if np.sum(w) > 0 else np.nan


def rms_vr(v_r, vol, mask):
    """Volume-weighted RMS radial velocity over `mask` (analytic v_r = 0)."""
    w = vol[mask]
    return float(np.sqrt(np.sum(w * v_r[mask] ** 2) / np.sum(w))) \
        if np.sum(w) > 0 else np.nan




# --------------------------------------------------------------------------- #
#  Default mode: radial profiles at --time vs the steady analytic
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, params, t, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, axes = plt.subplots(3, 1, figsize=(8, 11), sharex=True)
    ax_vphi, ax_vr, ax_rho = axes

    # Reference overlay: the blessed binned v_phi/v_r/rho(r) (dashed grey).
    if ref_table is not None and "v_phi_prof" in ref_table:
        xr = ref_table["r_prof"]
        for ax, key in ((ax_vphi, "v_phi_prof"), (ax_vr, "v_r_prof"),
                        (ax_rho, "rho_prof")):
            if key in ref_table:
                ax.plot(xr, ref_table[key], "--", color="0.35", lw=1.6, zorder=5,
                        label=("reference" if ax is ax_vphi else None))

    rinner, rdisk = params["rinner"], params["rdisk"]
    rmax = rdisk
    notes = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        rcyl, v_phi, v_r, rho, vol, _m, _time = kinematics(fn)
        mask = ring_mask(rcyl, rinner, rdisk)
        rmax = max(rmax, float(np.percentile(rcyl[mask], 99.5)))
        ax_vphi.scatter(rcyl, v_phi, s=3, alpha=0.4, rasterized=True,
                        label=f"{lab} (t={tt:.3g})")
        ax_vr.scatter(rcyl, v_r, s=3, alpha=0.4, rasterized=True)
        ax_rho.scatter(rcyl, rho, s=3, alpha=0.4, rasterized=True)
        l1 = l1_vphi(rcyl, v_phi, vol, mask)
        vr = rms_vr(v_r, vol, mask)
        print(f"coldkeplerian[{lab}]: t={tt:.4g}  L1(v_phi)={l1:.6g}  "
              f"RMS(v_r)={vr:.6g}")
        notes.append(f"{lab}: L1(v_phi)={l1:.3g}, RMS(v_r)={vr:.3g}")

    # Analytic overlays (time-independent steady solution).
    rg = np.linspace(max(1e-3, 0.5 * rinner), rmax, 400)
    ax_vphi.plot(rg, analytic_vphi(rg), "k-", lw=1.8, label="analytic (Kepler)")
    ax_vr.axhline(0.0, color="k", ls="-", lw=1.5, label="analytic (cold, $v_r=0$)")
    for ax in (ax_vphi, ax_vr):
        ax.axvspan(rinner, rdisk, color="0.85", alpha=0.4, zorder=0)
    ax_rho.axvspan(rinner, rdisk, color="0.85", alpha=0.4, zorder=0,
                   label="initial ring")

    if notes:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.8)
        ax_vphi.text(0.98, 0.96, "\n".join(notes), transform=ax_vphi.transAxes,
                     fontsize=8, va="top", ha="right", bbox=abox)

    ax_vphi.set_ylabel(r"$v_\phi=(x v_y-y v_x)/r_{cyl}$")
    ax_vr.set_ylabel(r"$v_r=(x v_x+y v_y)/r_{cyl}$")
    ax_rho.set_ylabel(r"$\rho$")
    ax_rho.set_xlabel(r"$r_{cyl}=\sqrt{x^2+y^2}$")
    ax_vphi.set_title(f"Cold Keplerian ring profiles (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax_vphi,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  Ring spreading + angular-momentum conservation vs time
# --------------------------------------------------------------------------- #

def ring_series(inp, rinner, rdisk):
    """(t, sigma_r, L_z) arrays for run `inp`, sorted in time.

    sigma_r = mass-weighted std of r_cyl over the ring particles (its growth is
    the AV-driven spreading); L_z = sum m_i (x vy - y vx)_i is the total z angular
    momentum (conserved exactly by symmetric forces). Returns (None,None,None) if
    no snapshots."""
    try:
        files = find_files(inp)
    except FileNotFoundError:
        return None, None, None
    if not files:
        return None, None, None
    ts, sig, lz = [], [], []
    for fn in files:
        tgdata, _td, _ts, _hdr, time = tip.readtipsy(fn)
        x, y = tgdata[:, 1], tgdata[:, 2]
        vx, vy = tgdata[:, 4], tgdata[:, 5]
        m = tgdata[:, 0]
        rcyl = np.sqrt(x * x + y * y)
        mask = ring_mask(rcyl, rinner, rdisk)
        wm = m[mask]
        rbar = np.sum(wm * rcyl[mask]) / np.sum(wm)
        var = np.sum(wm * (rcyl[mask] - rbar) ** 2) / np.sum(wm)
        ts.append(float(np.asarray(time).ravel()[0]))
        sig.append(float(np.sqrt(max(var, 0.0))))
        lz.append(float(np.sum(m * (x * vy - y * vx))))
    o = np.argsort(ts)
    return np.array(ts)[o], np.array(sig)[o], np.array(lz)[o]


def plot_ring_vs_time(inputs, labels, params, save=None, ref_table=None):
    """sigma_r(t)/sigma_r(0) and L_z(t)/L_z(0), both against the ideal line 1.

    When a reference is supplied its blessed sigma_r/L_z histories are overlaid
    (dashed grey)."""
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, (ax_s, ax_l) = plt.subplots(2, 1, figsize=(8, 9), sharex=True)
    rinner, rdisk = params["rinner"], params["rdisk"]
    if ref_table is not None and "x" in ref_table:
        if "sigma_r" in ref_table:
            ax_s.plot(ref_table["x"], ref_table["sigma_r"], "--", color="0.35",
                      lw=1.5, zorder=5, label="reference")
        if "Lz" in ref_table:
            ax_l.plot(ref_table["x"], ref_table["Lz"], "--", color="0.35",
                      lw=1.5, zorder=5, label="reference")
    any_data = False
    for inp, lab in zip(inputs, labels):
        ts, sig, lz = ring_series(inp, rinner, rdisk)
        if ts is None:
            print(f"coldkeplerian: no snapshots for '{inp}'", file=sys.stderr)
            continue
        if sig[0] <= 0 or lz[0] == 0:
            print(f"coldkeplerian[{lab}]: degenerate initial sigma_r/L_z; "
                  f"skipping.", file=sys.stderr)
            continue
        line, = ax_s.plot(ts, sig / sig[0], "o-", ms=3, label=lab)
        ax_l.plot(ts, lz / lz[0], "o-", ms=3, color=line.get_color(), label=lab)
        print(f"coldkeplerian[{lab}]: sigma_r({ts[-1]:.3g})/sigma_r(0)="
              f"{sig[-1] / sig[0]:.4g}  "
              f"L_z drift={abs(lz[-1] / lz[0] - 1.0):.3g}  ({len(ts)} snaps)")
        any_data = True
    for ax in (ax_s, ax_l):
        ax.axhline(1.0, color="k", ls="--", lw=1.3, label="ideal")
    ax_s.set_ylabel(r"$\sigma_r(t)/\sigma_r(0)$  (ring spreading)")
    ax_l.set_ylabel(r"$L_z(t)/L_z(0)$  (AM conservation)")
    ax_l.set_xlabel("t")
    ax_s.set_title("Cold Keplerian ring: spreading & angular-momentum vs time")
    fig.tight_layout()
    if not any_data:
        print("coldkeplerian: no ring series to plot.", file=sys.stderr)
    finish_figure_with_legend(fig, ax_s,
                              save=(f"{save}_ring.png" if save else None))


# --------------------------------------------------------------------------- #
#  Convergence mode: L1(v_phi) vs resolution, against n_x^-2
# --------------------------------------------------------------------------- #

def plot_convergence(inputs, params, t, save=None, ref_table=None):
    """L1(v_phi) vs n_x (log-log) against the ideal 2nd-order n_x^-2 line.

    Inputs are the SAME cold-ring test at different resolutions; n_x is the
    leading '<test><nx>' number in each run-folder name. The Keplerian ring is
    smooth with an exact analytic, so a well-behaved scheme converges as n_x^-2."""
    import matplotlib.pyplot as plt
    rinner, rdisk = params["rinner"], params["rdisk"]
    nxs, l1s = [], []
    for inp in inputs:
        nx = parse_resolution(inp)
        if nx is None:
            print(f"coldkeplerian: cannot read '<test><nx>' resolution from "
                  f"'{inp}'; skipping (needed for --convergence).",
                  file=sys.stderr)
            continue
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        rcyl, v_phi, _vr, _rho, vol, _m, _time = kinematics(fn)
        mask = ring_mask(rcyl, rinner, rdisk)
        l1 = l1_vphi(rcyl, v_phi, vol, mask)
        nxs.append(nx)
        l1s.append(l1)
        print(f"coldkeplerian: n_x={nx:<5d} t={tt:.4g}  L1(v_phi)={l1:.6g}")

    if len(nxs) < 2:
        print("coldkeplerian: need >=2 resolutions with a valid '<test><nx>' "
              "number for a convergence plot.", file=sys.stderr)
        return

    order = np.argsort(nxs)
    nxs = np.array(nxs)[order]
    l1s = np.array(l1s)[order]

    setup_rcparams()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(nxs, l1s, "o-", label=r"$L_1(v_\phi)$")
    ref = l1s[0] * (nxs / nxs[0]) ** (-2.0)   # ideal 2nd order at coarsest n_x
    ax.loglog(nxs, ref, "k--", label=r"$\propto n_x^{-2}$")
    if (ref_table is not None and "n_x" in ref_table
            and "L1_vphi" in ref_table):
        ax.loglog([float(ref_table["n_x"][0])], [float(ref_table["L1_vphi"][0])],
                  "D", color="C3", ms=9, mfc="none", label="reference")
    ax.set_xlabel(r"$n_x$")
    ax.set_ylabel(r"$L_1(v_\phi)$")
    ax.set_title(f"Cold Keplerian ring convergence (t={t:g})")
    fig.tight_layout()
    finish_figure_with_legend(fig, ax,
                              save=(f"{save}_convergence.png" if save else None))


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the normalised ring series)
# --------------------------------------------------------------------------- #

def reference_table(inp, params, t=None):
    """Multi-block regression table feeding ALL non-render plots:
      gate / ring  x        : time grid
                   sigma_r  : sigma_r/sigma_r0  (AV-driven spread, a gate key)
                   Lz       : L_z/L_z0          (AM conservation, a gate key)
      profile      r_prof   : r_cyl bin centres
                   v_phi_prof/v_r_prof/rho_prof : binned profiles at the time t
                   L1_vphi  : volume-weighted L1(v_phi) over the ring (1-elt)
      convergence  n_x      : the run's '<test><nx>' resolution

    The normalised ring spread + angular momentum are the headline diagnostics
    (the resolution-independent gate). None if neither a usable ring series nor a
    profile snapshot."""
    table = {}
    ts, sig, lz = ring_series(inp, params["rinner"], params["rdisk"])
    if ts is not None:
        table["x"] = ts
        if np.isfinite(sig[0]) and sig[0] > 0:
            table["sigma_r"] = sig / sig[0]
        if np.isfinite(lz[0]) and lz[0] != 0:
            table["Lz"] = lz / lz[0]
    sel = select_at_time(inp, t)
    if sel is not None:
        fn, _tt = sel
        rcyl, v_phi, v_r, rho, vol, _m, _time = kinematics(fn)
        mask = ring_mask(rcyl, params["rinner"], params["rdisk"])
        cx, qb = binned_profile(rcyl[mask],
                                {"v_phi": v_phi[mask], "v_r": v_r[mask],
                                 "rho": rho[mask]}, nbins=reg_nbins(inp))
        if cx is not None:
            table["r_prof"] = cx
            table["v_phi_prof"] = qb["v_phi"]
            table["v_r_prof"] = qb["v_r"]
            table["rho_prof"] = qb["rho"]
            table["L1_vphi"] = np.array([l1_vphi(rcyl, v_phi, vol, mask)])
        nx = parse_resolution(inp)
        if nx is not None:
            table["n_x"] = np.array([float(nx)])
    # Need at least one comparable series beyond the bare time grid.
    return table if any(k not in ("x",) for k in table) else None


# --------------------------------------------------------------------------- #
#  Render: face-on density map of the ring
# --------------------------------------------------------------------------- #

def do_render(inputs, labels, t, args, save):
    """Face-on (x-y) density map at time t, forwarding the --render-* flags.

    project defaults to 'average' (the thin disk is z-invariant); the box is
    auto-read from the run .log unless --render-extent is given."""
    spec = RenderSpec(1, 2, 7, "density", True, "average", None,
                      "inferno", False)
    entries = list(zip(labels, inputs))
    render_panels(entries, spec, times=[t], n=1, res=args.render_res,
                  save=(f"{save}_density" if save else None),
                  backend=args.render_backend,
                  project=args.render_project,
                  slice_frac=args.render_slice_frac,
                  extent=args.render_extent, aspect=args.render_aspect,
                  nsmooth=args.render_nsmooth)


def _add_args(parser):
    parser.add_argument("--convergence", action="store_true",
                        help="treat inputs as different resolutions and plot "
                             "L1(v_phi) vs n_x against the n_x^-2 line")


def _plot(inputs, labels, params, args):
    print(f"[coldkeplerian] hr={params['hr']:g}  ring=[{params['rinner']:g},"
          f"{params['rdisk']:g}]  q={params['q']:g}  "
          f"(GM={GM:g}, v_phi=r^-1/2)")
    ref = getattr(args, "ref_table", None)
    if args.convergence:
        plot_convergence(inputs, params, args.time, save=args.save, ref_table=ref)
    else:
        plot_profiles(inputs, labels, params, args.time, save=args.save,
                      ref_table=ref)
    plot_ring_vs_time(inputs, labels, params, save=args.save, ref_table=ref)
    if args.render:
        do_render(inputs, labels, args.time, args, save=args.save)


def main():
    standalone_cli(
        "coldkeplerian", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, params,
                                                                  args.time),
        description="Cold Keplerian disk/ring analysis.",
        reg_keys=("sigma_r", "Lz"), time_default=DEFAULT_TIME,
        time_help=f"comparison/render time (default {DEFAULT_TIME:g}; the "
                  f"analytic ring is steady, so any time works)",
        add_arguments=_add_args)


if __name__ == "__main__":
    main()
