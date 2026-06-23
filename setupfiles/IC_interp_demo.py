#!/usr/bin/env python
"""Showcase the `sph_interp` interpolation across every target geometry.

Takes ONE tipsy snapshot, deposits the SAME field with the SAME method onto each
target geometry — UniformGrid, AMRGrid, VoronoiGrid, Projection2D — and renders
them side by side in one figure, next to the raw SPH particles, so the
per-geometry deposit can be eyeballed:

  particles | UniformGrid | AMRGrid (octree) | VoronoiGrid | Projection2D

A high density-contrast snapshot (default: mhdcollapse) makes the geometry
differences obvious — the AMR octree refines on the collapsing core and the
Voronoi cells shrink where the particles crowd, both tracking the structure the
uniform grid resolves bluntly.

  interp_demo.sh [SNAPSHOT] [--method sph|petkova] [--field rho] [--npx N] ...

Each geometry's panel shares a colour scale (the 3D volume grids show a
mid-plane slab slice; Projection2D shows its column map with its own scale).
Use --per-geometry to ALSO write one full SPH-vs-grid diagnostic figure per
geometry (the `--visualise` panel set) for a closer look.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# A sensible high-contrast default if the user names no snapshot.
DEFAULT_SNAP = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "test_cases/mhdcollapse/"
    "mhdcollapse32_N64latticemu10rd360Erat0.045facQ1GASOLINE/"
    "mhdcollapse32_N64latticemu10rd360Erat0.045facQ1GASOLINE.00100")


def show_figure(path, enabled):
    """Pop `path` open in an image viewer (non-blocking, survives this process)."""
    import shutil
    import subprocess
    if not enabled:
        return
    for viewer in ("display", "feh", "eog", "xdg-open"):
        if shutil.which(viewer):
            try:
                subprocess.Popen([viewer, path], start_new_session=True,
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL)
                print(f"  displaying {path} ({viewer})")
            except Exception as e:
                print(f"  could not display ({viewer}): {e}")
            return
    print("  no image viewer found (display/feh/eog/xdg-open) — PNG written only")


def geom_cells(res, fld):
    """(centres (N,3), density (N,), volume (N,)) for any 3D target, else None."""
    from sph_interp import UniformGrid, AMRGrid, VoronoiGrid
    grid = res.target
    if isinstance(grid, UniformGrid):
        cx, cy, cz = grid.centres()
        C = np.stack(np.meshgrid(cx, cy, cz, indexing="ij"), axis=-1).reshape(-1, 3)
        return (C, np.asarray(res.data[fld]).ravel(),
                np.full(C.shape[0], grid.cell_volume))
    if isinstance(grid, (AMRGrid, VoronoiGrid)):
        return np.asarray(grid.centres()), np.asarray(res.data[fld]), grid.volumes()
    return None                                      # Projection2D: no 3D density


def report_budget(results, coarse, p, fld):
    """Print total mass (Sum rho_cell*V_cell) and total volume (Sum V_cell) for
    every grid, against the SPH particles (mass Sum m_i, volume Sum m_i/rho_i)."""
    sph_m = float(np.sum(p.mass))
    sph_v = float(np.sum(p.mass / np.maximum(p.rho, 1e-30)))
    box_v = float(np.prod(p.box)) if p.box is not None else float("nan")
    print("=" * 70)
    print(f"  mass / volume budget (field '{fld}')")
    print(f"  SPH particles : mass = {sph_m:.6g}   "
          f"volume = {sph_v:.6g}  (Sum m_i/rho_i)")
    print(f"  domain (box)  : volume = {box_v:.6g}")
    for label, meshes in (("full-resolution", results), ("coarse", coarse)):
        if not meshes:
            continue
        print(f"  {label} meshes:")
        for name in ("UniformGrid", "AMRGrid", "VoronoiGrid"):
            gc = geom_cells(meshes[name], fld)
            if gc is None:
                continue
            _, dens, vol = gc
            m = float(np.sum(dens * vol)); v = float(np.sum(vol))
            print(f"    {name:<12s} mass = {m:.6g} ({100*(m/sph_m-1):+6.2f}%)   "
                  f"volume = {v:.6g} ({100*(v/sph_v-1):+6.2f}%)")
    print("=" * 70)


def shell_density(r, mass, vol, edges):
    """Volume-weighted mean density per radial shell = sum(mass)/sum(volume)."""
    m, _ = np.histogram(r, bins=edges, weights=mass)
    v, _ = np.histogram(r, bins=edges, weights=vol)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(v > 0, m / v, np.nan)
    return d


def sample_field_at_points(res, fld, pts):
    """Look up the mesh field at arbitrary 3D points (which cell contains each
    point) for any 3D target; returns (Npts,) with NaN where a point is outside
    the mesh. Lets the radial profile evaluate the ACTUAL field, so a coarse /
    single-cell mesh reads as its constant density at every radius."""
    from sph_interp import UniformGrid, AMRGrid, VoronoiGrid
    grid = res.target
    vals = np.asarray(res.data[fld])
    n = pts.shape[0]
    out = np.full(n, np.nan)
    if isinstance(grid, UniformGrid):
        npx = np.asarray(grid.npx)
        idx = np.floor((pts - grid.lo) / grid.pixwidth).astype(np.int64)
        inside = np.all((idx >= 0) & (idx < npx), axis=1)
        ic = np.clip(idx, 0, npx - 1)
        out[inside] = vals[ic[inside, 0], ic[inside, 1], ic[inside, 2]]
        return out
    if isinstance(grid, VoronoiGrid):
        from sklearn.neighbors import BallTree           # membership = nearest gen
        j = BallTree(grid.centres()).query(pts, k=1, return_distance=False).ravel()
        return vals[j]
    if isinstance(grid, AMRGrid):
        nb = np.asarray(grid.base_npx)
        bidx = np.floor((pts - grid.lo) / grid.base_pw).astype(np.int64)
        inside = np.all((bidx >= 0) & (bidx < nb), axis=1)
        bc = np.clip(bidx, 0, nb - 1)
        bflat = (bc[:, 0] * nb[1] + bc[:, 1]) * nb[2] + bc[:, 2]
        clo = grid.cell_lo; chi = grid.cell_lo + grid.cell_size
        for k in range(n):                              # few leaves per base cell
            if not inside[k]:
                continue
            b = bflat[k]; P = pts[k]
            for li in range(grid.base_start[b], grid.base_start[b] + grid.base_count[b]):
                if (clo[li, 0] <= P[0] < chi[li, 0] and clo[li, 1] <= P[1] < chi[li, 1]
                        and clo[li, 2] <= P[2] < chi[li, 2]):
                    out[k] = vals[li]
                    break
        return out
    return out                                          # Projection2D: no 3D field


def sampled_radial_profile(res, fld, cc, rmid, box, periodic, rng, nsphere=600):
    """Radial profile by SAMPLING the mesh field: at each shell radius, average
    the field over `nsphere` random directions on the sphere of that radius about
    the core `cc`. Unlike binning cell centres this reads the actual field, so a
    coarse / single-cell mesh shows a flat line at its density (periodic axes are
    wrapped into the centered box for the lookup)."""
    prof = np.full(len(rmid), np.nan)
    half = 0.5 * np.asarray(box) if box is not None else None
    per = np.atleast_1d(periodic)
    for i, rr in enumerate(rmid):
        d = rng.normal(size=(nsphere, 3))
        d /= np.linalg.norm(d, axis=1, keepdims=True)
        pts = cc + rr * d
        if half is not None:
            for a in range(3):
                if (per[a] if per.size == 3 else per[0]):
                    pts[:, a] = (pts[:, a] + half[a]) % box[a] - half[a]
        v = sample_field_at_points(res, fld, pts)
        if np.isfinite(v).any():
            prof[i] = np.nanmean(v)
    return prof


def build_parser():
    pa = argparse.ArgumentParser(
        prog="interp_demo",
        description="Showcase sph_interp across all target geometries.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("snapshot", nargs="?", default=DEFAULT_SNAP,
                    help="tipsy snapshot (default: an mhdcollapse snapshot)")
    pa.add_argument("--method", choices=["sph", "petkova"], default="sph")
    pa.add_argument("--field", default="rho", help="field to interpolate/show")
    pa.add_argument("--npx", type=int, default=128,
                    help="resolution: uniform/projection cells per axis; AMR "
                         "finest level; image raster for AMR/Voronoi slices")
    pa.add_argument("--axes", choices=["xy", "xz", "yz"], default="xy",
                    help="slice / projection plane")
    pa.add_argument("--cmap", default="turbo",
                    help="colormap for the density / field panels. Default turbo "
                         "(no black — voids read blue, not black). Other no-black "
                         "options: viridis, plasma, cividis, Spectral_r, jet, "
                         "rainbow. (inferno/magma/gist_heat are black at the low "
                         "end.)")
    pa.add_argument("--base-npx", type=int, default=16, help="AMR base grid/axis")
    pa.add_argument("--voronoi-ncells", type=int, default=40000,
                    help="Voronoi generators (subsampled particles); 0 = use all")
    pa.add_argument("--project", default="column",
                    choices=["column", "average", "rhocolumn", "slice"],
                    help="Projection2D mode")
    pa.add_argument("--mesh-zoom", type=float, default=0.25,
                    help="fraction of the box (centred) for the actual-cells "
                         "row — cells are only legible zoomed in")
    pa.add_argument("--out", default="interp_demo.png",
                    help="output figure (png/pdf)")
    pa.add_argument("--no-mesh3d", dest="mesh3d", action="store_false",
                    help="skip the second figure (true 3D cells: Uniform/AMR/"
                         "Voronoi polyhedra in a sub-box)")
    pa.set_defaults(mesh3d=True)
    pa.add_argument("--no-profile", dest="profile", action="store_false",
                    help="skip the density-vs-radius profile figure")
    pa.set_defaults(profile=True)
    pa.add_argument("--mesh3d-cells", type=int, default=1500,
                    help="max cells to draw in the 3D figure (per geometry; "
                         "keeps the highest-value cells, never random)")
    pa.add_argument("--mesh3d-coarse", type=int, default=12,
                    help="3D figure uses a COARSE dedicated mesh — this is its "
                         "uniform cells/axis (also sets the coarse AMR finest "
                         "level, derived from it and --mesh3d-base-npx)")
    pa.add_argument("--mesh3d-base-npx", type=int, default=6,
                    help="base grid/axis for the COARSE AMR mesh (3D figure + "
                         "profile + coarse budget); finest level is derived so "
                         "base*2^level >= --mesh3d-coarse")
    pa.add_argument("--mesh3d-voronoi", type=int, default=1500,
                    help="generators for the coarse 3D Voronoi mesh "
                         "(~matches the AMR leaf count shown)")
    pa.add_argument("--mesh3d-cull", type=float, default=None,
                    help="drop 3D cells below this value percentile (declutter "
                         "voids); default None — voids fade via value-alpha")
    pa.add_argument("--mesh3d-clip",
                    choices=["none", "half", "quarter", "octant"],
                    default="quarter",
                    help="cutaway for the 3D figure: remove the near "
                         "half/quarter/octant (centred on the densest cell) to "
                         "expose the core and central AMR refinement")
    g = pa.add_mutually_exclusive_group()
    g.add_argument("--show", dest="show", action="store_true",
                   help="pop the figures open in an X window (default: on when "
                        "$DISPLAY is set)")
    g.add_argument("--no-show", dest="show", action="store_false",
                   help="only write the PNGs, do not open a viewer")
    pa.set_defaults(show=None)                       # None -> auto from $DISPLAY
    pa.add_argument("--per-geometry", action="store_true",
                    help="also write one SPH-vs-grid diagnostic figure per "
                         "geometry (<out-stem>_<geometry>.pdf)")
    pa.add_argument("--nsmooth", type=int, default=64)
    return pa


def main(argv=None):
    args = build_parser().parse_args(argv)

    # a 3D Voronoi tessellation needs at least a few generators; below 3 it is
    # degenerate (and pyvoro can fail). Warn and floor at 3. (0 = "use all".)
    for attr in ("voronoi_ncells", "mesh3d_voronoi"):
        v = getattr(args, attr)
        if 0 < v < 3:
            print(f"  WARNING: --{attr.replace('_', '-')}={v} is too few Voronoi "
                  f"generators; using 3.")
            setattr(args, attr, 3)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    from sph_interp import (from_tipsy, interpolate, UniformGrid, AMRGrid,
                            VoronoiGrid, Projection2D, refine_smoothing_length,
                            refine_count, field_image, visualise, draw_cells,
                            draw_cells_3d)

    a1, a2 = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}[args.axes]
    fld = args.field
    show = args.show if args.show is not None else bool(os.environ.get("DISPLAY"))

    print(f"reading {args.snapshot}  field={fld}")
    p = from_tipsy(args.snapshot, fields=(fld,), nsmooth=args.nsmooth)
    print(f"  {p.n} gas particles, box {np.round(p.box, 4)}")
    box = p.box

    # --- build + interpolate each target geometry --------------------------
    results = {}

    print("  UniformGrid ...")
    g = UniformGrid.centered(box, args.npx)
    results["UniformGrid"] = interpolate(p, g, method=args.method, values=[fld])

    print("  AMRGrid (octree) ...")
    import math
    max_level = max(0, math.ceil(math.log2(max(args.npx, 1) / args.base_npx)))
    g = AMRGrid.build(p, base_npx=args.base_npx, max_level=max_level,
                      criterion=refine_smoothing_length())
    print(f"    {g.ncell} leaf cells (base {args.base_npx}^3, max_level {max_level})")
    results["AMRGrid"] = interpolate(p, g, method=args.method, values=[fld])

    print("  VoronoiGrid ...")
    half = 0.5 * np.asarray(box)
    # Periodic box -> bounds ARE the period (ghost images handle the boundary, so
    # Sum cell volume = box volume exactly). Non-periodic -> a small pad.
    vper = p.periodic if any(np.atleast_1d(p.periodic)) else False
    if vper:
        bounds = [[-half[a], half[a]] for a in range(3)]
    else:
        bounds = [[-half[a] - 0.02 * box[a], half[a] + 0.02 * box[a]]
                  for a in range(3)]
    if args.voronoi_ncells and args.voronoi_ncells < p.n:
        idx = np.random.default_rng(0).choice(p.n, args.voronoi_ncells,
                                              replace=False)
        pts = p.pos[idx]
        print(f"    {args.voronoi_ncells} generators (subsampled particles)"
              f"{' periodic' if any(np.atleast_1d(vper)) else ''}")
    else:
        pts = p.pos
        print(f"    {p.n} generators (all particles)")
    g = VoronoiGrid.from_points(pts, bounds=bounds, periodic=vper)
    results["VoronoiGrid"] = interpolate(p, g, method=args.method, values=[fld])

    print(f"  Projection2D [{args.project}] ...")
    g = Projection2D.centered(box, a1, a2, npx=args.npx, mode=args.project)
    results["Projection2D"] = interpolate(p, g, method=args.method, values=[fld])

    # --- assemble images; share a colour scale across the 3D volume grids ---
    order = ["UniformGrid", "AMRGrid", "VoronoiGrid", "Projection2D"]
    imgs = {}
    for name in order:
        img, extent, _ = field_image(results[name], fld, plane=(a1, a2),
                                     slab_frac=0.05, img_npx=args.npx)
        imgs[name] = (img, extent)

    vol = np.concatenate([imgs[n][0][np.isfinite(imgs[n][0])].ravel()
                          for n in ("UniformGrid", "AMRGrid", "VoronoiGrid")])
    vpos = vol[vol > 0]
    use_log = vpos.size and (vpos.max() / vpos.min() > 50)
    vnorm = (LogNorm(vmin=vpos.min(), vmax=vpos.max()) if use_log
             else Normalize(vmin=vol.min(), vmax=vol.max()))

    # --- mesh-resolution images: paint each geometry by its CELL SIZE, so the
    #     adaptivity is visible (AMR refines / Voronoi cells shrink on the core;
    #     uniform & projection are flat by construction). Reuse field_image by
    #     injecting a synthetic per-cell "cell width" field into each result.
    def cell_width_field(res):
        grid = res.target
        if isinstance(grid, UniformGrid):
            return np.full(tuple(grid.npx), float(grid.pixwidth.min()))
        if isinstance(grid, Projection2D):
            return np.full(tuple(grid.npx), float(grid.pixwidth.min()))
        return grid.volumes() ** (1.0 / 3.0)          # AMR / Voronoi: leaf width

    struct = {}
    for name in order:
        res = results[name]
        res.data["_cellw"] = cell_width_field(res)
        img, extent, _ = field_image(res, "_cellw", plane=(a1, a2),
                                     slab_frac=0.05, img_npx=args.npx)
        struct[name] = (img, extent)
    sall = np.concatenate([struct[n][0][np.isfinite(struct[n][0])].ravel()
                           for n in ("UniformGrid", "AMRGrid", "VoronoiGrid")])
    snorm = LogNorm(vmin=sall[sall > 0].min(), vmax=sall.max())

    # --- zoom window (centred) for the actual-cells row --------------------
    cen = 0.5 * (p.pos.min(0) + p.pos.max(0))
    qh = 0.5 * args.mesh_zoom * box
    zoom = (cen[a1] - qh[a1], cen[a1] + qh[a1],
            cen[a2] - qh[a2], cen[a2] + qh[a2])

    # --- the showcase figure: 3 rows x 5 cols ------------------------------
    #     row 0 = interpolated field; row 1 = mesh resolution (cell width);
    #     row 2 = the ACTUAL cells (boundaries), zoomed in so they are legible.
    axlabel = "xyz"
    fig, axes = plt.subplots(3, 5, figsize=(26, 15.8))
    fig.suptitle(f"sph_interp showcase — {os.path.basename(args.snapshot)} — "
                 f"field '{fld}', method '{args.method}', plane {args.axes}",
                 fontsize=16, fontweight="bold")
    cmap = args.cmap
    scmap = "viridis_r"                              # small cells -> bright

    los = ({0, 1, 2} - {a1, a2}).pop()
    pos = p.pos
    zm = 0.5 * (pos[:, los].min() + pos[:, los].max())
    half_s = 0.5 * 0.05 * (pos[:, los].max() - pos[:, los].min())
    sel = np.abs(pos[:, los] - zm) <= half_s
    cval = p.values[fld] if fld in p.values else p.rho

    # row 0, col 0: particle field
    ax = axes[0, 0]
    ax.scatter(pos[sel, a1], pos[sel, a2], c=cval[sel], s=2, cmap=cmap,
               norm=vnorm, linewidths=0)
    ax.set_title(f"SPH particles\n(slab, {sel.sum()} pts)")
    ax.set_ylabel("interpolated field", fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    # row 1, col 0: particle smoothing length (the SPH resolution)
    ax = axes[1, 0]
    ax.scatter(pos[sel, a1], pos[sel, a2], c=p.h[sel], s=2, cmap=scmap,
               norm=LogNorm(), linewidths=0)
    ax.set_title("smoothing length h")
    ax.set_ylabel("mesh resolution (cell width)", fontsize=12, fontweight="bold")
    ax.set_aspect("equal")
    # row 2, col 0: particles inside the zoom window (what the cells represent)
    ax = axes[2, 0]
    selz = (sel & (pos[:, a1] >= zoom[0]) & (pos[:, a1] <= zoom[1])
            & (pos[:, a2] >= zoom[2]) & (pos[:, a2] <= zoom[3]))
    ax.scatter(pos[selz, a1], pos[selz, a2], c=cval[selz], s=6, cmap=cmap,
               norm=vnorm, linewidths=0)
    ax.set_xlim(zoom[0], zoom[1]); ax.set_ylim(zoom[2], zoom[3])
    ax.set_title(f"particles (zoom, {selz.sum()} pts)")
    ax.set_ylabel(f"actual cells (zoom x{1/args.mesh_zoom:.0f})",
                  fontsize=12, fontweight="bold")
    ax.set_xlabel(axlabel[a1]); ax.set_aspect("equal")

    for k, name in enumerate(order, start=1):
        # --- row 0: field ---
        ax = axes[0, k]
        img, extent = imgs[name]
        if name == "Projection2D":                # own scale (column units)
            d = img[np.isfinite(img)]; dp = d[d > 0]
            pn = (LogNorm(vmin=dp.min(), vmax=dp.max())
                  if dp.size and dp.max() / dp.min() > 50
                  else Normalize(vmin=d.min(), vmax=d.max()))
            im = ax.imshow(np.ma.masked_invalid(img.T), origin="lower",
                           extent=extent, cmap=cmap, norm=pn, aspect="equal",
                           interpolation="nearest")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            ax.set_title(f"Projection2D\n[{args.project}] (own scale)")
        else:
            ax.imshow(np.ma.masked_invalid(img.T), origin="lower", extent=extent,
                      cmap=cmap, norm=vnorm, aspect="equal",
                      interpolation="nearest")
            extra = (f"\n{results[name].target.ncell} leaves" if name == "AMRGrid"
                     else f"\n{results[name].target.ncell} cells"
                     if name == "VoronoiGrid" else "")
            ax.set_title(f"{name}{extra}")
        # --- row 1: mesh resolution ---
        ax = axes[1, k]
        simg, sext = struct[name]
        flat = name in ("UniformGrid", "Projection2D")
        im = ax.imshow(np.ma.masked_invalid(simg.T), origin="lower", extent=sext,
                       cmap=scmap, norm=snorm, aspect="equal",
                       interpolation="nearest")
        ax.set_title("regular (constant cell)" if flat
                     else "adaptive (cell tracks density)")

        # --- row 2: the actual cells (boundaries), zoomed in ---
        ax = axes[2, k]
        cnorm = vnorm if name != "Projection2D" else None
        _, ncelltxt = draw_cells(ax, results[name], fld, plane=(a1, a2),
                                 zoom=zoom, norm=cnorm, cmap=cmap)
        ax.set_title(ncelltxt)
        ax.set_xlabel(axlabel[a1])

    # shared colorbars
    sm = plt.cm.ScalarMappable(norm=vnorm, cmap=cmap)
    fig.colorbar(sm, ax=axes[0, 1:4], fraction=0.025, pad=0.01,
                 label=f"{fld} (shared)")
    sm2 = plt.cm.ScalarMappable(norm=snorm, cmap=scmap)
    fig.colorbar(sm2, ax=axes[1, 1:4], fraction=0.025, pad=0.01,
                 label="cell width (shared)")
    for name in order:                                # tidy the injected field
        results[name].data.pop("_cellw", None)

    fig.savefig(args.out, dpi=130, bbox_inches="tight")
    print(f"wrote showcase -> {args.out}")
    show_figure(args.out, show)

    # --- coarse dedicated meshes (used by BOTH the 3D figure and the profile) -
    #     The interpolation meshes are far too fine for a 3D cell picture (a
    #     sub-box slices into thousands of tiny fragments). Build a coarse,
    #     COMPLETE mesh per geometry — few, large cells, like the original
    #     interpolate_petkova.py polyhedra figures.
    cc = p.pos[np.argmax(cval)]                      # densest cell: core centre
    coarse = None
    if args.mesh3d or args.profile:
        print("  building coarse dedicated meshes ...")
        ncu = args.mesh3d_coarse
        nab = max(1, args.mesh3d_base_npx)
        # coarse AMR: base nab + a few levels (finest ~ mesh3d_coarse), refined on
        # PARTICLE COUNT so it stays coarse in voids and nests on the crowded core.
        # By DEFAULT floor the depth at 3 levels so the octree nesting is visible
        # in the 3D figure; but if the user explicitly sets either coarse-AMR knob,
        # honour it EXACTLY with no floor (so --mesh3d-coarse 1 --mesh3d-base-npx 1
        # really is a single AMR cell, not 3 forced levels -> ~27 leaves).
        _d = build_parser()
        user_amr = (ncu != _d.get_default("mesh3d_coarse")
                    or args.mesh3d_base_npx != _d.get_default("mesh3d_base_npx"))
        mlc = max(0 if user_amr else 3,
                  math.ceil(math.log2(max(ncu, 1) / nab)))
        ncrit = max(32, int(p.n / (nab ** 3 * 8)))   # ~ a few-per-base-cell target
        coarse = {}
        gu = UniformGrid.centered(box, ncu)
        coarse["UniformGrid"] = interpolate(p, gu, method=args.method, values=[fld])
        ga = AMRGrid.build(p, base_npx=nab, max_level=mlc,
                           criterion=refine_count(ncrit))
        coarse["AMRGrid"] = interpolate(p, ga, method=args.method, values=[fld])
        nvg = min(args.mesh3d_voronoi, p.n)
        iv = np.random.default_rng(1).choice(p.n, nvg, replace=False)
        gv = VoronoiGrid.from_points(p.pos[iv], bounds=bounds, periodic=vper)
        coarse["VoronoiGrid"] = interpolate(p, gv, method=args.method, values=[fld])
        print(f"    coarse meshes: uniform {ncu}^3, AMR base{nab}/level{mlc} "
              f"({ga.ncell} leaves), Voronoi {nvg} cells")

    # --- mass / volume budget (printed text) --------------------------------
    report_budget(results, coarse, p, fld)

    # --- second figure: the TRUE 3D cells on the coarse meshes --------------
    if args.mesh3d:
        print("  rendering true 3D cells ...")
        cvall = np.concatenate([np.asarray(coarse[n].data[fld]).ravel()
                                for n in coarse])
        cvp = cvall[cvall > 0]
        cnorm3 = (LogNorm(vmin=cvp.min(), vmax=cvp.max())
                  if cvp.size and cvp.max() / cvp.min() > 50
                  else Normalize(vmin=cvall.min(), vmax=cvall.max()))

        # cutaway centred on the densest cell -> the cut faces pass through the
        # core, exposing it (and the central AMR refinement). View is +x+y+z, so
        # remove the HIGH (near) side.
        clipdef = {"none": None,
                   "half": [(1, "high", cc[1])],
                   "quarter": [(0, "high", cc[0]), (1, "high", cc[1])],
                   "octant": [(0, "high", cc[0]), (1, "high", cc[1]),
                              (2, "high", cc[2])]}[args.mesh3d_clip]

        fig3 = plt.figure(figsize=(21, 7))
        fig3.suptitle(f"sph_interp — actual 3D cells (coarse meshes, "
                      f"{args.mesh3d_clip} cutaway) — "
                      f"{os.path.basename(args.snapshot)}, field '{fld}', "
                      f"method '{args.method}'", fontsize=15, fontweight="bold")
        for j, name in enumerate(["UniformGrid", "AMRGrid", "VoronoiGrid"], 1):
            ax3 = fig3.add_subplot(1, 3, j, projection="3d")
            _, lab = draw_cells_3d(ax3, coarse[name], fld, norm=cnorm3, cmap=cmap,
                                   max_cells=args.mesh3d_cells,
                                   cull_percentile=args.mesh3d_cull,
                                   clip=clipdef)
            ax3.set_title(f"{name}\n{lab}")
        sm3 = plt.cm.ScalarMappable(norm=cnorm3, cmap=cmap)
        fig3.colorbar(sm3, ax=fig3.axes, fraction=0.015, pad=0.02, label=fld)
        out3 = os.path.splitext(args.out)[0] + "_mesh3d.png"
        fig3.savefig(out3, dpi=130, bbox_inches="tight")
        print(f"wrote 3D cells -> {out3}")
        show_figure(out3, show)

    # --- third figure: density vs radius, full-res and coarse meshes --------
    if args.profile:
        print("  building density-vs-radius profile ...")
        # radial bins centred on the core (densest cell)
        rpart = np.linalg.norm(p.pos - cc, axis=1)
        rmax = np.percentile(rpart, 99.5)
        redges = np.linspace(0.0, rmax, 40)         # linear radial bins
        rmid = 0.5 * (redges[:-1] + redges[1:])
        # SPH particle reference: volume-weighted mean density per shell
        pvol = p.mass / np.maximum(p.rho, 1e-30)
        sph_prof = shell_density(rpart, p.mass, pvol, redges)

        gcol = {"UniformGrid": "tab:blue", "AMRGrid": "tab:green",
                "VoronoiGrid": "tab:red"}
        panels = [("full-resolution meshes", results),
                  ("coarse meshes", coarse)]
        figp, axp = plt.subplots(1, 2, figsize=(15, 6.2), sharex=True, sharey=True)
        for ax, (label, meshes) in zip(axp, panels):
            # light particle scatter for context (subsampled)
            sub = (np.random.default_rng(0).choice(p.n, 20000, replace=False)
                   if p.n > 20000 else np.arange(p.n))
            ax.scatter(rpart[sub], p.rho[sub], s=1, alpha=0.08, color="0.6",
                       linewidths=0, zorder=1)
            ax.plot(rmid, sph_prof, "k-", lw=2.6, label="SPH particles", zorder=5)
            rng_prof = np.random.default_rng(7)
            for name in ("UniformGrid", "AMRGrid", "VoronoiGrid"):
                if geom_cells(meshes[name], fld) is None:
                    continue
                # sample the mesh field on each radial shell (not bin cell centres),
                # so a coarse / single-cell mesh reads as a flat line at its density
                prof = sampled_radial_profile(meshes[name], fld, cc, rmid,
                                              p.box, p.periodic, rng_prof)
                ax.plot(rmid, prof, lw=2.0, color=gcol[name], label=name, zorder=4)
            ax.set_yscale("log")                    # density spans decades
            ax.set_xlim(0, rmax)
            ax.set_xlabel("radius from core")
            ax.set_title(label)
            ax.grid(True, which="both", ls=":", alpha=0.3)
            ax.legend(fontsize=9)
        axp[0].set_ylabel(f"{fld} (radial mean: mesh field-sampled, SPH shell-binned)")
        figp.suptitle(f"sph_interp — {fld} vs radius — "
                      f"{os.path.basename(args.snapshot)}, method '{args.method}'",
                      fontsize=14, fontweight="bold")
        figp.tight_layout(rect=[0, 0, 1, 0.96])
        outp = os.path.splitext(args.out)[0] + "_profile.png"
        figp.savefig(outp, dpi=130, bbox_inches="tight")
        print(f"wrote profile -> {outp}")
        show_figure(outp, show)

    if args.per_geometry:
        stem = os.path.splitext(args.out)[0]
        for name in order:
            out = f"{stem}_{name}.pdf"
            visualise(results[name], p, field=fld, out=out, plane=(a1, a2),
                      title=f"{name} / {args.method} — {fld}")

    print("done.")


if __name__ == "__main__":
    main()
