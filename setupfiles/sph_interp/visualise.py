"""Visual diagnostics for `sph_interp` — see the interpolated grid and judge how
faithfully it reproduces the SPH particle field.

Two entry points:

  field_image(res, field, ...)        -> a 2D slab image of ANY target geometry
                                          (UniformGrid / AMRGrid / VoronoiGrid /
                                          Projection2D), so every geometry can be
                                          eyeballed the same way.
  visualise(res, particles, field,    -> a multi-panel figure comparing the SPH
            out=...)                     particles to the interpolated grid:
                                          particle truth, interpolated grid,
                                          density PDFs, and a per-particle
                                          grid-vs-SPH scatter + mass-conservation
                                          line. Saved to PDF/PNG.

Geometry-aware rendering of the SAME plane for each target:
  UniformGrid   slab-averaged slice (LOS cells within the slab averaged);
  AMRGrid       leaf mosaic (each leaf intersecting the slab painted as its box,
                fine cells over coarse) -> the octree refinement is visible;
  VoronoiGrid   nearest-generator raster at the slice plane -> the true Voronoi
                cells; markers fall back to a generator scatter if asked;
  Projection2D  the 2D map as-is.
"""

import numpy as np

# Headless: pick a non-interactive backend before pyplot is imported.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize


# --------------------------------------------------------------------------- #
#  Geometry-aware 2D image of an interpolated field
# --------------------------------------------------------------------------- #

def _norm_for(values, log=None):
    """A LogNorm for strictly-positive wide-range data, else a linear Normalize."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return Normalize()
    pos = v[v > 0]
    if log is None:
        log = pos.size > 0 and (pos.max() / max(pos.min(), pos.max() * 1e-12) > 50)
    if log and pos.size:
        return LogNorm(vmin=pos.min(), vmax=pos.max())
    return Normalize(vmin=float(v.min()), vmax=float(v.max()))


def field_image(res, field, plane=(0, 1), slab_frac=0.05, img_npx=None,
                zmid=None):
    """Rasterize `res.data[field]` onto a 2D image in the `plane` (a1, a2).

    Returns (image2D, extent, (a1, a2)) where image2D is indexed [col=a1, row=a2]
    and `extent = [lo_a1, hi_a1, lo_a2, hi_a2]` (display with image.T,
    origin='lower'). For 3D targets a slab of half-thickness 0.5*slab_frac*box_los
    centred on `zmid` (default mid-plane) is collapsed onto the plane.
    """
    from .targets import UniformGrid, AMRGrid, VoronoiGrid, Projection2D
    grid = res.target
    a1, a2 = plane
    if field not in res.data:
        raise KeyError(f"field {field!r} not in result ({sorted(res.data)})")
    data = np.asarray(res.data[field])

    if isinstance(grid, Projection2D):
        extent = [grid.lo[0], grid.hi[0], grid.lo[1], grid.hi[1]]
        return data, extent, (grid.axis1, grid.axis2)

    if isinstance(grid, UniformGrid):
        los = ({0, 1, 2} - {a1, a2}).pop()
        centres = grid.centres()
        zc = centres[los]
        zm = 0.5 * (grid.lo[los] + grid.hi[los]) if zmid is None else zmid
        half = 0.5 * slab_frac * (grid.hi[los] - grid.lo[los])
        sel = np.abs(zc - zm) <= max(half, grid.pixwidth[los])
        if not sel.any():
            sel[np.argmin(np.abs(zc - zm))] = True
        img = data.take(np.where(sel)[0], axis=los).mean(axis=los)
        # img is now indexed by the two remaining axes in ascending order
        rem = sorted({0, 1, 2} - {los})
        if (a1, a2) != tuple(rem):           # caller wants a2 first
            img = img.T
        extent = [grid.lo[a1], grid.hi[a1], grid.lo[a2], grid.hi[a2]]
        return img, extent, (a1, a2)

    # AMR / Voronoi: rasterize onto an image grid in the plane.
    lo = grid.lo
    hi = grid.hi
    los = ({0, 1, 2} - {a1, a2}).pop()
    zm = 0.5 * (lo[los] + hi[los]) if zmid is None else zmid
    npx = img_npx or 512
    e1 = np.linspace(lo[a1], hi[a1], npx + 1)
    e2 = np.linspace(lo[a2], hi[a2], npx + 1)
    img = np.full((npx, npx), np.nan)

    if isinstance(grid, AMRGrid):
        clo = grid.cell_lo
        csz = grid.cell_size
        chi = clo + csz
        inslab = (clo[:, los] <= zm) & (chi[:, los] >= zm)
        order = np.argsort(grid.cell_level[inslab])      # coarse first
        idx = np.where(inslab)[0][order]
        pw1 = (hi[a1] - lo[a1]) / npx
        pw2 = (hi[a2] - lo[a2]) / npx
        for c in idx:
            i0 = int(np.floor((clo[c, a1] - lo[a1]) / pw1))
            i1 = int(np.ceil((chi[c, a1] - lo[a1]) / pw1))
            j0 = int(np.floor((clo[c, a2] - lo[a2]) / pw2))
            j1 = int(np.ceil((chi[c, a2] - lo[a2]) / pw2))
            img[max(i0, 0):min(i1, npx), max(j0, 0):min(j1, npx)] = data[c]
        extent = [lo[a1], hi[a1], lo[a2], hi[a2]]
        return img, extent, (a1, a2)

    if isinstance(grid, VoronoiGrid):
        from scipy.spatial import cKDTree
        c1 = 0.5 * (e1[:-1] + e1[1:])
        c2 = 0.5 * (e2[:-1] + e2[1:])
        G1, G2 = np.meshgrid(c1, c2, indexing="ij")
        pts = np.empty((G1.size, 3))
        pts[:, a1] = G1.ravel()
        pts[:, a2] = G2.ravel()
        pts[:, los] = zm
        _, nn = cKDTree(grid.centres()).query(pts)   # central cells only
        img = data[nn].reshape(npx, npx)
        extent = [lo[a1], hi[a1], lo[a2], hi[a2]]
        return img, extent, (a1, a2)

    raise NotImplementedError(f"no image for target {type(grid).__name__}")


# --------------------------------------------------------------------------- #
#  Draw the ACTUAL cells (mesh boundaries), filled by the field value
# --------------------------------------------------------------------------- #

def draw_cells(ax, res, field, plane=(0, 1), zoom=None, norm=None,
               cmap="inferno", edgecolor="0.25", lw=0.3, slab=None):
    """Draw the real cell boundaries of `res.target` in the `plane`, each cell
    filled by `res.data[field]`, onto Matplotlib axis `ax`.

    This shows the MESH itself (not a raster): UniformGrid -> regular squares;
    AMRGrid -> the octree leaf boxes that cross the slice plane; VoronoiGrid ->
    the Voronoi polygons (2D tessellation of the generators within a thin slab,
    the faithful slice-plane depiction); Projection2D -> the pixel grid.

    zoom : (x0, x1, y0, y1) in world coords to crop to — cells are only legible
           zoomed in; default is the central quarter of the domain.
    slab : LOS half-thickness for the Voronoi generator selection (default a few
           mean spacings). Ignored by the other geometries.
    """
    from matplotlib.collections import PatchCollection, PolyCollection
    from matplotlib.patches import Rectangle
    from .targets import UniformGrid, AMRGrid, VoronoiGrid, Projection2D
    grid = res.target
    a1, a2 = (grid.axis1, grid.axis2) if isinstance(grid, Projection2D) else plane
    data = np.asarray(res.data[field])
    if norm is None:
        norm = _norm_for(data)

    # default zoom: central quarter of the plane
    if isinstance(grid, Projection2D):
        lo = np.array([grid.lo[0], grid.lo[1]]); hi = np.array([grid.hi[0], grid.hi[1]])
    else:
        lo = np.array([grid.lo[a1], grid.lo[a2]]); hi = np.array([grid.hi[a1], grid.hi[a2]])
    if zoom is None:
        c = 0.5 * (lo + hi); q = 0.25 * (hi - lo)
        zoom = (c[0] - q[0], c[0] + q[0], c[1] - q[1], c[1] + q[1])
    zx0, zx1, zy0, zy1 = zoom

    if isinstance(grid, (UniformGrid, Projection2D)):
        if isinstance(grid, Projection2D):
            e1, e2 = grid.edges(); img = data
        else:
            ed = grid.edges(); e1, e2 = ed[a1], ed[a2]
            img, _, _ = field_image(res, field, plane=(a1, a2), slab_frac=0.05)
        ax.pcolormesh(e1, e2, np.ma.masked_invalid(img.T), cmap=cmap, norm=norm,
                      edgecolors=edgecolor, linewidth=lw, shading="flat")
        ncell = "regular pixels"

    elif isinstance(grid, AMRGrid):
        los = ({0, 1, 2} - {a1, a2}).pop()
        zm = 0.5 * (grid.lo[los] + grid.hi[los])
        clo, csz = grid.cell_lo, grid.cell_size
        chi = clo + csz
        inslab = ((clo[:, los] <= zm) & (chi[:, los] >= zm)
                  & (chi[:, a1] >= zx0) & (clo[:, a1] <= zx1)
                  & (chi[:, a2] >= zy0) & (clo[:, a2] <= zy1))
        idx = np.where(inslab)[0]
        rects = [Rectangle((clo[c, a1], clo[c, a2]), csz[c, a1], csz[c, a2])
                 for c in idx]
        pc = PatchCollection(rects, cmap=cmap, norm=norm,
                             edgecolor=edgecolor, linewidth=lw)
        pc.set_array(data[idx])
        ax.add_collection(pc)
        ncell = f"{idx.size} leaf cells in slice"

    elif isinstance(grid, VoronoiGrid):
        from scipy.spatial import Voronoi
        los = ({0, 1, 2} - {a1, a2}).pop()
        gen = grid.centres()                         # central cells only
        zm = 0.5 * (grid.lo[los] + grid.hi[los])
        if slab is None:
            slab = 2.0 * (np.prod(grid.hi - grid.lo) / grid.ncell) ** (1.0 / 3.0)
        sel = (np.abs(gen[:, los] - zm) <= slab)
        # keep a margin around the zoom so boundary polygons close
        mx = 0.15 * (zx1 - zx0); my = 0.15 * (zy1 - zy0)
        sel &= ((gen[:, a1] >= zx0 - mx) & (gen[:, a1] <= zx1 + mx)
                & (gen[:, a2] >= zy0 - my) & (gen[:, a2] <= zy1 + my))
        idx = np.where(sel)[0]
        pts2 = np.column_stack([gen[idx, a1], gen[idx, a2]])
        vor = Voronoi(pts2)
        polys, vals = [], []
        for k, reg_i in enumerate(vor.point_region):
            reg = vor.regions[reg_i]
            if not reg or -1 in reg:                 # skip open (boundary) cells
                continue
            polys.append(vor.vertices[reg])
            vals.append(data[idx[k]])
        pcoll = PolyCollection(polys, cmap=cmap, norm=norm,
                               edgecolor=edgecolor, linewidth=lw)
        pcoll.set_array(np.asarray(vals))
        ax.add_collection(pcoll)
        ncell = f"{len(polys)} Voronoi cells in slab"
    else:
        raise NotImplementedError(type(grid).__name__)

    ax.set_xlim(zx0, zx1); ax.set_ylim(zy0, zy1)
    ax.set_aspect("equal")
    return norm, ncell


def _cube_faces(lo, size):
    """The 6 quad faces of an axis-aligned box, as a list of (4,3) vertex arrays."""
    x0, y0, z0 = lo
    x1, y1, z1 = lo + size
    v = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                  [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]])
    f = [(0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6),
         (1, 2, 6, 5), (0, 3, 7, 4)]
    return [v[list(q)] for q in f]


def draw_cells_3d(ax, res, field, subbox=None, norm=None, cmap="inferno",
                  alpha=0.55, alpha_min=0.04, edgecolor="k", lw=0.3,
                  max_cells=1500, cull_percentile=None, value_alpha=True,
                  clip=None):
    """Render the ACTUAL 3D cells of `res.target` as polyhedra/boxes on a 3D axis
    `ax`, each face coloured by `res.data[field]`.

    The faithful counterpart to `draw_cells`: VoronoiGrid draws the true polyhedra
    from the stored face/vertex CSR (not a 2D slab approximation); AMRGrid and
    UniformGrid draw their leaf/cell cubes.

    This is only legible for a SMALL number of LARGE cells — feed it a COARSE
    mesh (a few hundred cells), not the interpolation mesh. Decluttering, like
    the original Voronoi-density figure:
      value_alpha       per-cell opacity scales `alpha_min`..`alpha` with the
                        normalised field value, so voids fade out;
      cull_percentile   drop cells whose value is below this percentile;
      max_cells         if still too many, keep the HIGHEST-value cells (a
                        principled cull — never a random subsample, which would
                        punch holes in the tessellation).
    subbox (xlo,xhi,...,zhi) crops to a region; default = the whole domain.
    clip : a CUTAWAY — list of (axis, 'low'|'high', pos). A cell is removed only
           if its centre is on the named side of EVERY plane (AND), so one entry
           cuts a half, two a quarter wedge, three the near octant — exposing the
           core / central refinement behind the cut faces.
    """
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    from .targets import UniformGrid, AMRGrid, VoronoiGrid, Projection2D
    grid = res.target
    if isinstance(grid, Projection2D):
        raise NotImplementedError("Projection2D is a 2D target — use draw_cells")
    data = np.asarray(res.data[field])
    if norm is None:
        norm = _norm_for(data)
    cm = plt.get_cmap(cmap)

    # per-cell value + geometry (flatten UniformGrid to (Ncell,) C-order)
    if isinstance(grid, UniformGrid):
        cx, cy, cz = grid.centres()
        C = np.stack(np.meshgrid(cx, cy, cz, indexing="ij"), axis=-1)
        centres = C.reshape(-1, 3)
        cellval = data.ravel()
    else:
        centres = np.asarray(grid.centres())
        cellval = data
    lo = np.asarray(grid.lo); hi = np.asarray(grid.hi)
    if subbox is None:
        subbox = (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])
    bx = np.array(subbox).reshape(3, 2)
    idx = np.where(np.all((centres >= bx[:, 0]) & (centres <= bx[:, 1]),
                          axis=1))[0]

    if clip:                                              # cutaway: drop the wedge
        cen = centres[idx]
        rm = np.ones(idx.size, dtype=bool)
        for ax_, side, pos in clip:
            rm &= (cen[:, ax_] > pos) if side == "high" else (cen[:, ax_] < pos)
        idx = idx[~rm]

    v_in = cellval[idx]
    if cull_percentile is not None and idx.size:          # drop the void cells
        keep = v_in >= np.percentile(v_in, cull_percentile)
        idx, v_in = idx[keep], v_in[keep]
    if idx.size > max_cells:                               # keep densest cells
        top = np.argsort(v_in)[-max_cells:]
        idx, v_in = idx[top], v_in[top]

    # per-cell opacity (voids fade); repeated per face below
    a01 = np.clip(np.asarray(norm(v_in)), 0.0, 1.0)
    cell_alpha = (alpha_min + (alpha - alpha_min) * a01 if value_alpha
                  else np.full(idx.size, alpha))

    faces, facevals, facealpha = [], [], []
    if isinstance(grid, VoronoiGrid):
        for c, ci in enumerate(idx):
            for f in range(grid.cf_start[ci], grid.cf_start[ci + 1]):
                s, e = grid.fv_start[f], grid.fv_start[f + 1]
                faces.append(np.column_stack([grid.fvx[s:e], grid.fvy[s:e],
                                              grid.fvz[s:e]]))
                facevals.append(cellval[ci]); facealpha.append(cell_alpha[c])
        label = f"{idx.size} Voronoi polyhedra"
    else:
        if isinstance(grid, UniformGrid):
            i, j, k = np.unravel_index(idx, grid.npx)
            pw = grid.pixwidth
            cell_lo = np.column_stack([grid.lo[a] + ijk * pw[a]
                                       for a, ijk in enumerate((i, j, k))])
            cell_sz = np.tile(pw, (idx.size, 1))
        else:                                             # AMRGrid
            cell_lo, cell_sz = grid.cell_lo[idx], grid.cell_size[idx]
        for c in range(idx.size):
            for face in _cube_faces(cell_lo[c], cell_sz[c]):
                faces.append(face)
                facevals.append(v_in[c]); facealpha.append(cell_alpha[c])
        label = f"{idx.size} cells"

    fa = np.asarray(facealpha)
    colors = cm(norm(np.asarray(facevals)))
    colors[:, 3] = fa
    # fade the edges with the cell too, else void cells draw a solid black
    # wireframe that swamps the structure.
    ergba = np.array(matplotlib.colors.to_rgba(edgecolor))
    ecolors = np.tile(ergba, (fa.size, 1))
    ecolors[:, 3] = np.clip(fa * 1.1, 0.0, 1.0)
    pc = Poly3DCollection(faces, facecolors=colors, edgecolors=ecolors,
                          linewidths=lw)
    ax.add_collection3d(pc)
    ax.set_xlim(subbox[0], subbox[1]); ax.set_ylim(subbox[2], subbox[3])
    ax.set_zlim(subbox[4], subbox[5])
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.view_init(elev=25, azim=45)
    return norm, label


# --------------------------------------------------------------------------- #
#  Sample the interpolated grid back at the particle positions
# --------------------------------------------------------------------------- #

def sample_at_particles(res, particles, field):
    """Interpolated value at each particle position, or None for Projection2D.

    UniformGrid -> exact cell lookup; VoronoiGrid -> nearest generator (which IS
    the containing cell, by definition); AMRGrid -> nearest leaf centre (an
    approximation across refinement boundaries, adequate for a scatter).
    """
    from .targets import UniformGrid, AMRGrid, VoronoiGrid, Projection2D
    grid = res.target
    data = np.asarray(res.data[field])
    pos = np.asarray(particles.pos, dtype=np.float64)

    if isinstance(grid, Projection2D):
        return None                          # a LOS integral, no per-particle value

    if isinstance(grid, UniformGrid):
        idx = []
        for a in range(3):
            i = np.floor((pos[:, a] - grid.lo[a]) / grid.pixwidth[a]).astype(int)
            idx.append(np.clip(i, 0, grid.npx[a] - 1))
        return data[idx[0], idx[1], idx[2]]

    from scipy.spatial import cKDTree
    if isinstance(grid, VoronoiGrid):
        _, nn = cKDTree(grid.centres()).query(pos)   # central cells only
        return data[nn]
    if isinstance(grid, AMRGrid):
        _, nn = cKDTree(grid.centres()).query(pos)
        return data[nn]
    return None


def _cell_volumes(grid, field_shape):
    """Per-cell volume array matching the flattened field, or None (2D)."""
    from .targets import UniformGrid, AMRGrid, VoronoiGrid
    if isinstance(grid, UniformGrid):
        return np.full(int(np.prod(grid.npx)), grid.cell_volume)
    if isinstance(grid, (AMRGrid, VoronoiGrid)):
        return grid.volumes()
    return None


def _mass_check(res, particles):
    """(grid_mass, particle_mass) for a 'rho' field, or None."""
    from .targets import UniformGrid, AMRGrid, VoronoiGrid
    grid = res.target
    pm = float(np.asarray(particles.mass).sum())
    if "_mass" in res.data:
        return float(np.asarray(res.data["_mass"]).sum()), pm
    if "rho" not in res.data:
        return None
    if isinstance(grid, UniformGrid):
        return float(res.data["rho"].sum() * grid.cell_volume), pm
    if isinstance(grid, (AMRGrid, VoronoiGrid)):
        return float((res.data["rho"] * grid.volumes()).sum()), pm
    return None


# --------------------------------------------------------------------------- #
#  The full diagnostic figure
# --------------------------------------------------------------------------- #

def visualise(res, particles, field="rho", out=None, plane=(0, 1),
              slab_frac=0.05, img_npx=None, title=None, log=None,
              cmap="turbo"):
    """Compare the SPH particles to the interpolated grid; save to `out`.

    Panels (3D targets): particle truth (slab scatter) | interpolated grid
    (geometry-aware slice) | field PDFs (particles vs grid cells) | per-particle
    grid-vs-SPH scatter with the 1:1 line + a mass-conservation line.
    For Projection2D the per-particle scatter is replaced by an info panel.

    Returns the Matplotlib Figure.
    """
    from .targets import Projection2D
    grid = res.target
    gname = type(grid).__name__
    is_proj = isinstance(grid, Projection2D)
    a1, a2 = (grid.axis1, grid.axis2) if is_proj else plane
    axlabel = "xyz"

    img, extent, (pa1, pa2) = field_image(res, field, plane=(a1, a2),
                                           slab_frac=slab_frac, img_npx=img_npx)
    pvals = np.asarray(particles.values[field]) if field in particles.values \
        else (np.asarray(particles.rho) if field == "rho" else None)

    norm = _norm_for(img[np.isfinite(img)], log=log)

    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    fig.suptitle(title or f"{gname} / {res.method} — field '{field}'",
                 fontsize=14, fontweight="bold")

    # --- panel A: particle truth (slab scatter) -----------------------------
    ax = axes[0, 0]
    pos = np.asarray(particles.pos)
    if is_proj:
        sel = np.ones(particles.n, dtype=bool)
        slabtxt = "all particles (projected)"
    else:
        los = ({0, 1, 2} - {a1, a2}).pop()
        zm = 0.5 * (pos[:, los].min() + pos[:, los].max())
        half = 0.5 * slab_frac * (pos[:, los].max() - pos[:, los].min())
        sel = np.abs(pos[:, los] - zm) <= half
        slabtxt = f"|{axlabel[los]}-mid| < {half:.3g}  ({sel.sum()} pts)"
    csc = (pvals[sel] if pvals is not None else None)
    sc = ax.scatter(pos[sel, a1], pos[sel, a2], c=csc, s=3, cmap=cmap,
                    norm=norm if csc is not None else None, linewidths=0)
    ax.set_title(f"SPH particles — {slabtxt}")
    ax.set_xlabel(axlabel[a1]); ax.set_ylabel(axlabel[a2])
    ax.set_aspect("equal")
    if csc is not None:
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label=field)

    # --- panel B: interpolated grid (geometry-aware) ------------------------
    ax = axes[0, 1]
    im = ax.imshow(np.ma.masked_invalid(img.T), origin="lower", extent=extent,
                   cmap=cmap, norm=norm, aspect="equal", interpolation="nearest")
    modetxt = f" [{grid.mode}]" if is_proj else " slab-slice"
    ax.set_title(f"{gname} interpolation{modetxt}")
    ax.set_xlabel(axlabel[a1]); ax.set_ylabel(axlabel[a2])
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label=field)

    # --- panel C: mass-weighted field PDFs (particles vs grid cells) --------
    # Mass-weighting makes the two directly comparable: empty grid cells (voids)
    # carry ~0 mass and drop out, so this is "how is the material distributed in
    # `field`" for particles vs grid, not biased by the void volume.
    ax = axes[1, 0]
    cells = np.asarray(res.data[field]).ravel()
    cfin = np.isfinite(cells)
    cvol = _cell_volumes(grid, cells.shape)
    if "rho" in res.data and cvol is not None:           # cell mass = rho * vol
        cwt = (np.asarray(res.data["rho"]).ravel() * cvol)[cfin]
    elif cvol is not None:
        cwt = cvol[cfin]
    else:
        cwt = None
    cells = cells[cfin]
    parts = pvals.ravel() if pvals is not None else None
    pwt = np.asarray(particles.mass) if parts is not None else None
    use_log = isinstance(norm, LogNorm)
    if use_log:
        cpos = cells[cells > 0]
        bins = np.logspace(np.log10(cpos.min()), np.log10(cpos.max()), 60) \
            if cpos.size else 60
        ax.set_xscale("log")
    else:
        bins = np.linspace(np.nanmin(cells), np.nanmax(cells), 60)
    wlab = "mass-weighted " if cwt is not None else ""
    if parts is not None:
        ax.hist(parts, bins=bins, weights=pwt, density=True, histtype="step",
                lw=2, color="tab:blue", label="SPH particles")
    ax.hist(cells, bins=bins, weights=cwt, density=True, histtype="stepfilled",
            alpha=0.4, color="tab:orange", label=f"{gname} cells")
    ax.set_title(f"{field} {wlab}distribution (PDF)")
    ax.set_xlabel(field); ax.set_ylabel("probability density")
    ax.legend(fontsize=9)

    # --- panel D: per-particle grid-vs-SPH scatter or info ------------------
    ax = axes[1, 1]
    gsamp = None if is_proj else sample_at_particles(res, particles, field)
    mass = _mass_check(res, particles)
    if gsamp is not None and parts is not None:
        m = np.isfinite(gsamp) & np.isfinite(parts)
        x = parts[m]; y = gsamp[m]
        if m.sum() > 50000:                  # subsample for a readable scatter
            sub = np.random.default_rng(0).choice(m.sum(), 50000, replace=False)
            x, y = x[sub], y[sub]
        ax.scatter(x, y, s=2, alpha=0.2, color="k", linewidths=0)
        good = (x > 0) & (y > 0) if use_log else np.ones(x.shape, bool)
        if use_log:
            ax.set_xscale("log"); ax.set_yscale("log")
        lim_lo = min(x[good].min(), y[good].min())
        lim_hi = max(x[good].max(), y[good].max())
        ax.plot([lim_lo, lim_hi], [lim_lo, lim_hi], "r-", lw=1.5, label="1:1")
        xx = x[good]; yy = y[good]
        if use_log:
            xx, yy = np.log10(xx), np.log10(yy)
        r = np.corrcoef(xx, yy)[0, 1] if xx.size > 1 else float("nan")
        med = np.median(y[good] / x[good]) if good.any() else float("nan")
        ax.set_title(f"grid sampled at particles vs SPH\n"
                     f"r={r:.4f}, median(grid/SPH)={med:.3f}")
        ax.set_xlabel(f"SPH particle {field}")
        ax.set_ylabel(f"interpolated {field}")
        ax.legend(fontsize=9)
    else:
        ax.axis("off")
        lines = [f"target  : {gname}",
                 f"method  : {res.method}",
                 f"field   : {field}",
                 f"plane   : {axlabel[a1]}{axlabel[a2]}"]
        if is_proj:
            lines.append(f"mode    : {grid.mode}")
            lines.append("(projection is a LOS integral —")
            lines.append(" no per-particle 1:1 comparison)")
        ax.text(0.05, 0.95, "\n".join(lines), va="top", ha="left",
                family="monospace", fontsize=12, transform=ax.transAxes)
    if mass is not None:
        gm, pm = mass
        tag = ("exact" if "_mass" in res.data else "approx (SPH)")
        fig.text(0.5, 0.005,
                 f"mass: grid {gm:.6g} vs particles {pm:.6g}  "
                 f"({100*(gm/pm - 1):+.3f}%, {tag})",
                 ha="center", fontsize=11)

    fig.tight_layout(rect=[0, 0.02, 1, 0.97])
    if out:
        fig.savefig(out, dpi=130)
        print(f"  wrote visualisation -> {out}")
    return fig
