#!/usr/bin/env python3
"""
Analysis for the Noh problem (IC_setup_noh): a cold (P->0) uniform gas with a
uniform radial inflow |v|=v0 directed at the origin. An infinitely strong,
self-similar shock forms at the centre and propagates outward at

    D = 0.5 (gamma-1) v0           (= 1/3 for gamma=5/3, v0=1),  R_shock = D t.

The exact solution (dimensionality d = 3 spherical / 2 cylindrical, rho0=1):
  r < R_shock (post-shock, stationary):
      rho = rho0 ((gamma+1)/(gamma-1))^d   (= 64 in 3D, 16 in 2D),
      v_r = 0,   u = 1/2 v0^2,   P = (gamma-1) rho u
  r > R_shock (pre-shock, still infalling):
      rho = rho0 (1 + v0 t / r)^(d-1),   v_r = -v0,   P ~ pini (cold), u ~ pini.

This script compares the radial profiles of rho, v_r, P and u against that
solution at --time (default 0.6) with volume-weighted L1 errors, and renders the
x-y density map. A stringent test of shock capturing, wall heating (the spurious
density dip at r=0) and symmetry preservation.

Usage:
    python IC_analysis_noh.py <run-dir> [--time 0.6] [--save out] [--render]
"""

import os
import sys

import numpy as np
import readtipsy as tip

from IC_analysis_framework import (
    Param, standalone_cli,
    RenderSpec, render_panels, spec_slug,
    bless_render_specs, render_ref_overlay_path, movie_render_specs,
    select_at_time, binned_profile,
)
from IC_analysis_general import (
    find_files, loaddata, read_tufac,
    setup_rcparams, finish_figure_with_legend,
)

GAMMA = 5.0 / 3.0               # IC_setup_noh.gamma
DEFAULT_TIME = 0.6             # end time = nsteps*deltastep (shock at R=v0 t/3)

# runtest.sh setup_noh() parameter order (runtest.sh:743-761). `is3d` and `v0`
# set the analytic solution; `pini` is the (tiny) ambient pressure; the rest are
# carried for CLI/label parity.
PARAM_SPEC = [
    Param("a", "is3d", 1, int, "3D spherical (1) vs 2D cylindrical (0) geometry"),
    Param("b", "pini", 1e-6, float, "initial (tiny) pressure -> cold gas"),
    Param("c", "v0", 1.0, float, "radial inflow speed"),
]


def noh_profiles(r, t, params):
    """Analytic (rho, v_r, P, u, R_shock) at radii `r` and time t."""
    d = 3 if int(params["is3d"]) == 1 else 2
    v0 = float(params["v0"])
    pini = float(params["pini"])
    rho0 = 1.0
    g = GAMMA
    R = 0.5 * (g - 1.0) * v0 * t                      # shock radius D t
    rho2 = rho0 * ((g + 1.0) / (g - 1.0)) ** d        # post-shock density
    u2 = 0.5 * v0 ** 2                                # all KE -> internal energy
    p2 = (g - 1.0) * rho2 * u2                        # post-shock pressure
    r = np.asarray(r, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rho_pre = rho0 * (1.0 + v0 * t / np.where(r > 0, r, np.nan)) ** (d - 1)
    p_pre = pini * (rho_pre / rho0) ** g              # cold adiabatic (tiny)
    u_pre = p_pre / ((g - 1.0) * rho_pre)
    inside = r < R
    rho = np.where(inside, rho2, rho_pre)
    v = np.where(inside, 0.0, -v0)
    p = np.where(inside, p2, p_pre)
    u = np.where(inside, u2, u_pre)
    return rho, v, p, u, R


def profile(fn, tufac, is3d):
    """Per-particle (r, rho, v_r, P, u, volume) and snapshot time for run file fn.

    r is the spherical (3D) or cylindrical (2D) radius from the origin; v_r is the
    radial velocity; u = tufac*T (col 8 is temperature); P = (gamma-1) rho u."""
    tgdata, _td, _ts, _hdr, time, _N, _ng, _nd, _ns, _h = loaddata(fn)
    pos = tgdata[:, 1:4]
    vel = tgdata[:, 4:7]
    rho = tgdata[:, 7]
    u = tufac * tgdata[:, 8]
    if int(is3d) == 1:
        r = np.sqrt(np.sum(pos ** 2, axis=1))
        rv = np.sum(pos * vel, axis=1)
    else:                                            # 2D cylindrical (x-y only)
        r = np.sqrt(pos[:, 0] ** 2 + pos[:, 1] ** 2)
        rv = pos[:, 0] * vel[:, 0] + pos[:, 1] * vel[:, 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        v_r = np.where(r > 0, rv / r, 0.0)
    P = (GAMMA - 1.0) * rho * u
    vol = tgdata[:, 0] / rho                          # m/rho = particle volume
    return (r, rho, v_r, P, u, vol,
            float(np.asarray(time).ravel()[0]))


def l1_error(r, field, analytic, vol, rmax):
    """Volume-weighted L1 error of `field` vs `analytic`, over particles r<rmax."""
    sel = (r < rmax) & np.isfinite(analytic)
    if not np.any(sel):
        return float("nan")
    err = np.abs(field[sel] - analytic[sel])
    return float(np.sum(vol[sel] * err) / np.sum(vol[sel]))


# --------------------------------------------------------------------------- #
#  Profiles: rho, v_r, P, u vs r with analytic overlay + L1
# --------------------------------------------------------------------------- #

def plot_profiles(inputs, labels, t, params, save=None, ref_table=None):
    import matplotlib.pyplot as plt
    setup_rcparams()
    fig, axarr = plt.subplots(2, 2, figsize=(14, 9), sharex=True)
    axRho, axVel = axarr[0]
    axP, axU = axarr[1]
    axes = (axRho, axVel, axP, axU)
    rmax_box = 1.0                           # inscribed sphere of the [-1,1]^3 box

    # Reference overlay: the blessed binned profile of each panel (dashed grey).
    if ref_table is not None and "rho" in ref_table:
        xr = ref_table["x"]
        for ax, key in ((axRho, "rho"), (axVel, "v_r"), (axP, "P"), (axU, "u")):
            if key in ref_table:
                ax.plot(xr, ref_table[key], "--", color="0.35", lw=1.6, zorder=5,
                        label=("reference" if ax is axRho else None))

    l1_lines = []
    for inp, lab in zip(inputs, labels):
        sel = select_at_time(inp, t)
        if sel is None:
            continue
        fn, tt = sel
        tufac = read_tufac(inp)
        r, rho, v_r, P, u, vol, _time = profile(fn, tufac, params["is3d"])

        axRho.scatter(r, rho, s=3, alpha=0.4, rasterized=True, label=lab)
        axVel.scatter(r, v_r, s=3, alpha=0.4, rasterized=True, label=lab)
        axP.scatter(r, P, s=3, alpha=0.4, rasterized=True, label=lab)
        axU.scatter(r, u, s=3, alpha=0.4, rasterized=True, label=lab)

        rho_a, v_a, p_a, u_a, R = noh_profiles(r, tt, params)
        l1r = l1_error(r, rho, rho_a, vol, rmax_box)
        l1v = l1_error(r, v_r, v_a, vol, rmax_box)
        l1p = l1_error(r, P, p_a, vol, rmax_box)
        l1u = l1_error(r, u, u_a, vol, rmax_box)
        print(f"noh[{lab}]: t={tt:.4g}  R_shock={R:.4g}  L1(rho)={l1r:.4g}  "
              f"L1(v_r)={l1v:.4g}  L1(P)={l1p:.4g}  L1(u)={l1u:.4g}")
        l1_lines.append(f"{lab}: rho {l1r:.3g}, v {l1v:.3g}, "
                        f"P {l1p:.3g}, u {l1u:.3g}")

    # Analytic overlay (step profile) using the first input's time.
    if inputs:
        tt = select_at_time(inputs[0], t)[1]
        rg = np.linspace(1e-3, rmax_box, 1000)
        rho_a, v_a, p_a, u_a, R = noh_profiles(rg, tt, params)
        axRho.plot(rg, rho_a, "k-", lw=1.8, label="analytic")
        axVel.plot(rg, v_a, "k-", lw=1.8, label="analytic")
        axP.plot(rg, p_a, "k-", lw=1.8, label="analytic")
        axU.plot(rg, u_a, "k-", lw=1.8, label="analytic")
        for ax in axes:
            ax.axvline(R, color="0.5", ls=":", lw=1.2)

    if l1_lines:
        abox = dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85)
        axRho.text(0.98, 0.96, "\n".join(l1_lines), transform=axRho.transAxes,
                   fontsize=8, va="top", ha="right", bbox=abox)

    axRho.set_ylabel(r"$\rho$")
    axVel.set_ylabel(r"$v_r$")
    axP.set_ylabel(r"$P$")
    axU.set_ylabel(r"$u=p/((\gamma-1)\rho)$")
    rlabel = (r"$r=\sqrt{x^2+y^2+z^2}$" if int(params["is3d"]) == 1
              else r"$r=\sqrt{x^2+y^2}$")
    axP.set_xlabel(rlabel)
    axU.set_xlabel(rlabel)
    axP.set_xlim(0.0, rmax_box)              # shared x
    geo = "3D" if int(params["is3d"]) == 1 else "2D"
    fig.suptitle(f"Noh problem profiles (t={t:g}, {geo})")
    fig.tight_layout()
    finish_figure_with_legend(fig, axRho,
                              save=(f"{save}_profile.png" if save else None))


# --------------------------------------------------------------------------- #
#  x-y density render
# --------------------------------------------------------------------------- #

IC_BOX = [-1.0, 1.0, -1.0, 1.0]    # default render domain (the [-1,1]^2 IC box)


def render_spec():
    """The x-y density map (slug 'rho')."""
    return RenderSpec(1, 2, 7, r"$\rho$", True, "slice", None, "inferno", False,
                      slug="rho")


def do_render(inputs, labels, t, args, save, reference=None):
    """Render the x-y density map, forwarding the --render-* CLI flags.

    The CLI overrides win: --render-backend, --render-project (None -> the
    RenderSpec's "slice" default), --render-slice-frac, --render-res, and
    --render-extent (defaults to the IC box [-1,1]^2 when not given). When
    `reference` is given and there is a single input, the blessed density render
    grid is overlaid (sim | reference | difference)."""
    spec = render_spec()
    entries = list(zip(labels, inputs))
    extent = args.render_extent if args.render_extent else IC_BOX
    overlay = reference is not None and len(inputs) == 1
    ref_path = (render_ref_overlay_path(inputs[0], reference, spec_slug(spec))
                if overlay else None)
    render_panels(entries, spec, times=[t], n=1, res=args.render_res, save=save,
                  backend=args.render_backend, project=args.render_project,
                  slice_frac=args.render_slice_frac, extent=extent,
                  aspect=args.render_aspect, render_reference=ref_path,
                  residual=getattr(args, "residual", False))


def _save_extra(args, params):
    """Bless the density render-grid baseline (at the comparison time) for
    --save-reference, on the IC box so the overlay aligns."""
    extent = args.render_extent if args.render_extent else IC_BOX
    bless_render_specs(args.inputs[0], [render_spec()], args, "noh",
                       times=[args.time], extent=extent)


# --------------------------------------------------------------------------- #
#  Regression baseline (bless / compare the binned density profile at --time)
# --------------------------------------------------------------------------- #

REG_NBINS = 64


def reference_table(inp, t, params):
    """Regression table {x: r-grid, rho/v_r/P/u: binned profiles} at the comparison
    time t -- one blessed curve per profile panel (overlaid on each), all sharing
    the radial grid. The radial density profile (post-shock plateau + outer infall)
    is the headline Noh diagnostic and the gate key. None if no snapshot / too few
    bins."""
    sel = select_at_time(inp, t)
    if sel is None:
        return None
    fn, _tt = sel
    r, rho, v_r, P, u, _vol, _time = profile(fn, read_tufac(inp), params["is3d"])
    cx, qb = binned_profile(r, {"rho": rho, "v_r": v_r, "P": P, "u": u},
                            nbins=REG_NBINS, percentile=(0.5, 99.5))
    if cx is None:
        return None
    return {"x": cx, **qb}


def _plot(inputs, labels, params, args):
    plot_profiles(inputs, labels, args.time, params, save=args.save,
                  ref_table=getattr(args, "ref_table", None))
    if args.render:
        do_render(inputs, labels, args.time, args, save=args.save,
                  reference=args.reference)
    if args.movie:
        extent = args.render_extent if args.render_extent else IC_BOX
        movie_render_specs(inputs[0], labels[0], [render_spec()], args, args.save,
                           "noh", extent=extent)


def main():
    standalone_cli(
        "noh", PARAM_SPEC, _plot,
        reference_table=lambda inp, params, args: reference_table(inp, args.time,
                                                                  params),
        description="Noh problem analysis.",
        reg_keys=("rho",), time_default=DEFAULT_TIME,
        time_help=f"comparison / render time (default {DEFAULT_TIME:g})",
        save_extra=_save_extra)


if __name__ == "__main__":
    main()
