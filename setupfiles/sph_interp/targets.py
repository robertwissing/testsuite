"""Target geometries for interpolation.

A *target* defines the cells that particle quantities are deposited onto. The
interpolation *method* (regular SPH or Petkova exact) is orthogonal: any method
can in principle deposit onto any target, because a cell is a polyhedron either
way.

Implemented: `UniformGrid` (3D regular cubic cells), `Projection2D` (2D pixel
grid), `AMRGrid` (octree-refined cubic cells over a base grid), and `VoronoiGrid`
(unstructured polyhedral cells, the native Petkova case).
"""

from dataclasses import dataclass
import numpy as np


@dataclass
class UniformGrid:
    """A regular 3D grid of cubic cells spanning [lo, hi] with `npx` cells/axis.

    lo, hi : (3,) lower/upper corners of the domain.
    npx    : int or (3,) number of cells per axis.

    Cell centres are lo + (i + 0.5) * pixwidth. Use `centered(box, npx)` to get a
    grid spanning [-box/2, +box/2] (the convention of the periodic test boxes).
    """
    lo: np.ndarray
    hi: np.ndarray
    npx: object = 64

    def __post_init__(self):
        self.lo = np.asarray(self.lo, dtype=np.float64).reshape(3)
        self.hi = np.asarray(self.hi, dtype=np.float64).reshape(3)
        npx = self.npx
        if np.isscalar(npx):
            npx = (int(npx), int(npx), int(npx))
        self.npx = tuple(int(x) for x in npx)
        if any(n <= 0 for n in self.npx):
            raise ValueError("npx must be positive")
        if np.any(self.hi <= self.lo):
            raise ValueError("hi must exceed lo on every axis")

    @classmethod
    def centered(cls, box, npx=64):
        """Grid spanning [-box/2, +box/2] for a box size (scalar or (3,))."""
        box = np.asarray(box, dtype=np.float64).reshape(-1)
        if box.size == 1:
            box = np.repeat(box, 3)
        half = 0.5 * box
        return cls(lo=-half, hi=half, npx=npx)

    @property
    def box(self):
        return self.hi - self.lo

    @property
    def pixwidth(self):
        return self.box / np.asarray(self.npx, dtype=np.float64)

    @property
    def cell_volume(self):
        pw = self.pixwidth
        return float(pw[0] * pw[1] * pw[2])

    def edges(self):
        """Per-axis cell edge coordinates: list of three (npx[a]+1,) arrays."""
        return [np.linspace(self.lo[a], self.hi[a], self.npx[a] + 1)
                for a in range(3)]

    def centres(self):
        """Per-axis cell centre coordinates: list of three (npx[a],) arrays."""
        pw = self.pixwidth
        return [self.lo[a] + (np.arange(self.npx[a]) + 0.5) * pw[a]
                for a in range(3)]


@dataclass
class AMRGrid:
    """Adaptive (octree-refined) cubic cells over a base grid.

    A forest of octrees: the domain [lo,hi] is tiled by a `base_npx` base grid
    (level 0); a refinement criterion recursively splits cells into 8 children
    (down to `max_level`) where the data warrants it. The result is a flat list
    of LEAF cells (true memory-efficient AMR — storage scales with the number of
    leaves, not base_npx**3 * 8**max_level). Build one with `AMRGrid.build`.

    Geometry arrays (one row per leaf, sorted by the leaf's base-grid cell):
      cell_lo    (Ncell,3)  lower corner of each leaf
      cell_size  (Ncell,3)  edge lengths of each leaf
      cell_level (Ncell,)   refinement level (0 = base)
      cell_base  (Ncell,)   flat base-grid index the leaf descends from
    Lookup for the deposit (which leaves live in base cell b):
      base_start (nbase,)   first leaf index of base cell b (leaves are sorted)
      base_count (nbase,)   number of leaves in base cell b
    Every method deposits by looping the base cells a particle's kernel overlaps
    (O(1) index arithmetic on the base grid), then that base cell's few leaves.
    """
    lo: np.ndarray
    hi: np.ndarray
    base_npx: object
    cell_lo: np.ndarray
    cell_size: np.ndarray
    cell_level: np.ndarray
    cell_base: np.ndarray
    base_start: np.ndarray
    base_count: np.ndarray
    periodic: object = (False, False, False)

    def __post_init__(self):
        self.lo = np.asarray(self.lo, dtype=np.float64).reshape(3)
        self.hi = np.asarray(self.hi, dtype=np.float64).reshape(3)
        nb = self.base_npx
        if np.isscalar(nb):
            nb = (int(nb), int(nb), int(nb))
        self.base_npx = tuple(int(x) for x in nb)
        p = self.periodic
        if np.isscalar(p):
            p = (bool(p), bool(p), bool(p))
        self.periodic = tuple(bool(x) for x in p)

    @property
    def ncell(self):
        return self.cell_lo.shape[0]

    @property
    def box(self):
        return self.hi - self.lo

    @property
    def base_pw(self):
        return (self.hi - self.lo) / np.asarray(self.base_npx, dtype=np.float64)

    def centres(self):
        """(Ncell,3) leaf-cell centres."""
        return self.cell_lo + 0.5 * self.cell_size

    def volumes(self):
        """(Ncell,) leaf-cell volumes."""
        return self.cell_size[:, 0] * self.cell_size[:, 1] * self.cell_size[:, 2]

    def edges(self):
        """(cell_lo, cell_hi): the lower/upper corners of every leaf."""
        return self.cell_lo, self.cell_lo + self.cell_size

    @classmethod
    def build(cls, particles, base_npx=16, max_level=4, criterion=None,
              lo=None, hi=None, periodic=None):
        """Build an AMR octree over `particles` by refining a base grid.

        particles : Particles (uses pos, mass, h).
        base_npx  : base-grid cells per axis (int or (3,)); the level-0 tiling.
        max_level : maximum refinement depth (leaf size = base / 2**level).
        criterion : callable(count, mass, hmin, size, level) -> bool[Mcells],
                    True where a cell should be refined. Defaults to
                    `refine_smoothing_length()`. See `refine_mass` /
                    `refine_count`. Empty cells (count==0) are never refined.
        lo, hi    : domain corners (default: centered particles.box, else the
                    particle extent).
        periodic  : periodicity (default: particles.periodic).
        """
        if criterion is None:
            criterion = refine_smoothing_length()
        pos = np.ascontiguousarray(particles.pos, dtype=np.float64)
        mass = np.ascontiguousarray(particles.mass, dtype=np.float64)
        h = np.ascontiguousarray(particles.h, dtype=np.float64)
        N = pos.shape[0]

        if lo is None or hi is None:
            if particles.box is not None:
                half = 0.5 * np.asarray(particles.box, dtype=np.float64).reshape(3)
                lo, hi = -half, half
            else:
                lo, hi = pos.min(axis=0), pos.max(axis=0)
        lo = np.asarray(lo, dtype=np.float64).reshape(3)
        hi = np.asarray(hi, dtype=np.float64).reshape(3)
        if periodic is None:
            periodic = particles.periodic

        nb = base_npx
        if np.isscalar(nb):
            nb = (int(nb), int(nb), int(nb))
        nb = tuple(int(x) for x in nb)
        nbx, nby, nbz = nb
        nbase = nbx * nby * nbz
        bpw = (hi - lo) / np.array(nb, dtype=np.float64)

        # particle -> base cell (clip strays into the edge cells)
        ix = np.clip(((pos[:, 0] - lo[0]) / bpw[0]).astype(np.int64), 0, nbx - 1)
        iy = np.clip(((pos[:, 1] - lo[1]) / bpw[1]).astype(np.int64), 0, nby - 1)
        iz = np.clip(((pos[:, 2] - lo[2]) / bpw[2]).astype(np.int64), 0, nbz - 1)
        base_flat = (ix * nby + iy) * nbz + iz

        # active cells start as the full base grid (so the mesh tiles the domain)
        idx = np.arange(nbase)
        ac_lo = np.empty((nbase, 3))
        ac_lo[:, 0] = lo[0] + (idx // (nby * nbz)) * bpw[0]
        ac_lo[:, 1] = lo[1] + ((idx // nbz) % nby) * bpw[1]
        ac_lo[:, 2] = lo[2] + (idx % nbz) * bpw[2]
        ac_size = np.tile(bpw, (nbase, 1))
        ac_base = np.arange(nbase, dtype=np.int64)

        pidx = np.arange(N)              # particles still in play
        pcell = base_flat.copy()         # their active-cell index

        leaf_lo, leaf_size, leaf_level, leaf_base = [], [], [], []
        for level in range(max_level + 1):
            M = ac_lo.shape[0]
            if M == 0:
                break
            count = np.bincount(pcell, minlength=M).astype(np.int64)
            msum = np.bincount(pcell, weights=mass[pidx], minlength=M)
            hmn = np.full(M, np.inf)
            if pidx.size:
                np.minimum.at(hmn, pcell, h[pidx])
            if level == max_level:
                refine = np.zeros(M, dtype=bool)
            else:
                refine = np.asarray(criterion(count, msum, hmn, ac_size, level),
                                    dtype=bool)
                refine &= count > 0      # never refine empty cells
            leafm = ~refine
            if leafm.any():
                leaf_lo.append(ac_lo[leafm])
                leaf_size.append(ac_size[leafm])
                leaf_level.append(np.full(int(leafm.sum()), level, dtype=np.int64))
                leaf_base.append(ac_base[leafm])
            if not refine.any():
                break

            rc = np.where(refine)[0]
            nref = rc.shape[0]
            newsize = ac_size[rc] * 0.5
            child_lo = np.empty((nref * 8, 3))
            child_size = np.empty((nref * 8, 3))
            child_base = np.repeat(ac_base[rc], 8)
            for k in range(8):           # child k = ox + 2*oy + 4*oz
                ox, oy, oz = k & 1, (k >> 1) & 1, (k >> 2) & 1
                child_lo[k::8, 0] = ac_lo[rc, 0] + ox * newsize[:, 0]
                child_lo[k::8, 1] = ac_lo[rc, 1] + oy * newsize[:, 1]
                child_lo[k::8, 2] = ac_lo[rc, 2] + oz * newsize[:, 2]
                child_size[k::8] = newsize

            keepp = refine[pcell]        # particles whose cell was refined
            pidx2 = pidx[keepp]
            parent = pcell[keepp]
            rel = (pos[pidx2] - ac_lo[parent]) / ac_size[parent]
            oct_ = ((rel[:, 0] >= 0.5).astype(np.int64)
                    + 2 * (rel[:, 1] >= 0.5).astype(np.int64)
                    + 4 * (rel[:, 2] >= 0.5).astype(np.int64))
            parent_to_child = np.full(M, -1, dtype=np.int64)
            parent_to_child[rc] = np.arange(nref, dtype=np.int64) * 8
            pcell = parent_to_child[parent] + oct_
            ac_lo, ac_size, ac_base, pidx = child_lo, child_size, child_base, pidx2

        if leaf_lo:
            cell_lo = np.concatenate(leaf_lo)
            cell_size = np.concatenate(leaf_size)
            cell_level = np.concatenate(leaf_level)
            cell_base = np.concatenate(leaf_base)
        else:
            cell_lo = np.zeros((0, 3)); cell_size = np.zeros((0, 3))
            cell_level = np.zeros(0, dtype=np.int64); cell_base = np.zeros(0, dtype=np.int64)

        order = np.argsort(cell_base, kind="stable")   # group leaves by base cell
        cell_lo, cell_size = cell_lo[order], cell_size[order]
        cell_level, cell_base = cell_level[order], cell_base[order]
        base_count = np.bincount(cell_base, minlength=nbase).astype(np.int64)
        base_start = np.zeros(nbase, dtype=np.int64)
        base_start[1:] = np.cumsum(base_count)[:-1]

        return cls(lo=lo, hi=hi, base_npx=nb, cell_lo=cell_lo, cell_size=cell_size,
                   cell_level=cell_level, cell_base=cell_base,
                   base_start=base_start, base_count=base_count, periodic=periodic)


# --------------------------------------------------------------------------- #
#  Refinement criteria for AMRGrid.build — pluggable, add more freely.
#  Each is criterion(count, mass, hmin, size, level) -> bool[Mcells]:
#    count (M,)   particles in the cell        mass  (M,)   summed particle mass
#    hmin  (M,)   min smoothing length in cell size  (M,3) cell edge lengths
#    level int    current refinement level
#  Return True where the cell should be split. Empty cells are filtered by the
#  builder, so a criterion may ignore the count==0 case.
# --------------------------------------------------------------------------- #

def refine_smoothing_length(factor=1.0):
    """Refine until a cell is no larger than `factor` * the local smoothing
    length, i.e. resolution tracks the SPH resolution (fine where h is small /
    particles are dense, coarse in voids). factor=1 -> cell <~ h; factor=2 ->
    cell <~ 2h (the kernel support)."""
    def crit(count, mass, hmin, size, level):
        return size.max(axis=1) > factor * hmin
    return crit


def refine_mass(threshold):
    """Refine while a cell's summed particle mass exceeds `threshold` (with
    equal-mass particles this is equivalent to a particles-per-cell cap of
    threshold / m_particle)."""
    def crit(count, mass, hmin, size, level):
        return mass > threshold
    return crit


def refine_count(n_target):
    """Refine while a cell holds more than `n_target` particles."""
    def crit(count, mass, hmin, size, level):
        return count > n_target
    return crit


@dataclass
class VoronoiGrid:
    """Voronoi polyhedral cells (via pyvoro) — a third target geometry.

    A cell is an arbitrary convex polyhedron about a generator point, exactly
    what the Petkova exact integral was designed for (it integrates the kernel
    over any polyhedron face-by-face; the cube of UniformGrid/AMRGrid is just a
    special case). The SPH method evaluates the kernel at the generator.

    Build one with `VoronoiGrid.from_points(points, bounds)` (arbitrary mesh
    generators) or `VoronoiGrid.from_particles(particles)` (tessellate the
    particles themselves). pyvoro is required only for the build, not the deposit.

    Geometry is stored flat (numba-friendly), two-level CSR so cells/faces can
    have any number of faces/vertices:
      gen        (Ncell,3)  generator point of each cell
      volume     (Ncell,)   cell volume
      cell_rad   (Ncell,)   circumradius (max vertex distance from generator);
                            a particle's kernel can reach a cell only if
                            |x_particle - gen| < 2h + cell_rad (used to prune).
      cf_start   (Ncell+1,) CSR offsets: cell c's faces are [cf_start[c],cf_start[c+1])
      fv_start   (Nface+1,) CSR offsets: face f's vertices are [fv_start[f],fv_start[f+1])
      fvx/fvy/fvz(Nfv,)     ABSOLUTE vertex coords, inlined per face occurrence
                            (shared vertices are duplicated — the integral only
                            ever reads them in face-winding order, so this keeps
                            the kernel index-indirection-free).
    Face winding is pyvoro's own (the Petkova sign logic depends on it, the same
    convention `_CUBE_FACES` reproduces for the cube).

    NOTE: the deposit is NON-PERIODIC (absolute coordinates). For a periodic box,
    `from_particles` builds a periodic tessellation if asked, but boundary cells
    receive only the in-domain part of a kernel; interior mass is conserved.
    """
    gen: np.ndarray
    volume: np.ndarray
    cell_rad: np.ndarray
    cf_start: np.ndarray
    fv_start: np.ndarray
    fvx: np.ndarray
    fvy: np.ndarray
    fvz: np.ndarray
    lo: np.ndarray
    hi: np.ndarray
    periodic: object = (False, False, False)
    cell_parent: object = None      # (Ntotal,) ghost->central index; None = none
    ncentral: object = None         # number of real (central) cells

    def __post_init__(self):
        self.gen = np.ascontiguousarray(self.gen, dtype=np.float64).reshape(-1, 3)
        self.volume = np.ascontiguousarray(self.volume, dtype=np.float64).ravel()
        self.cell_rad = np.ascontiguousarray(self.cell_rad, dtype=np.float64).ravel()
        self.cf_start = np.ascontiguousarray(self.cf_start, dtype=np.int64).ravel()
        self.fv_start = np.ascontiguousarray(self.fv_start, dtype=np.int64).ravel()
        self.fvx = np.ascontiguousarray(self.fvx, dtype=np.float64).ravel()
        self.fvy = np.ascontiguousarray(self.fvy, dtype=np.float64).ravel()
        self.fvz = np.ascontiguousarray(self.fvz, dtype=np.float64).ravel()
        self.lo = np.asarray(self.lo, dtype=np.float64).reshape(3)
        self.hi = np.asarray(self.hi, dtype=np.float64).reshape(3)
        p = self.periodic
        if np.isscalar(p):
            p = (bool(p), bool(p), bool(p))
        self.periodic = tuple(bool(x) for x in p)
        if self.cell_parent is None:                  # plain (non-ghost) mesh
            self.ncentral = self.gen.shape[0]
        else:
            self.cell_parent = np.ascontiguousarray(self.cell_parent,
                                                    dtype=np.int64).ravel()
            self.ncentral = int(self.ncentral)

    @property
    def ncell(self):
        """Number of REAL cells (central; ghosts are deposit bookkeeping only)."""
        return self.ncentral

    @property
    def ntotal(self):
        """Total stored cells including periodic ghost images (for the deposit)."""
        return self.gen.shape[0]

    @property
    def box(self):
        return self.hi - self.lo

    def fold_to_central(self, arr):
        """Sum a per-(total-)cell deposit array down onto the central cells:
        ghost-image contributions are added back into their parent cell. Returns
        the first `ncentral` rows when there are no ghosts."""
        arr = np.asarray(arr)
        if self.cell_parent is None:
            return arr[:self.ncentral]
        out = np.zeros((self.ncentral,) + arr.shape[1:], dtype=arr.dtype)
        np.add.at(out, self.cell_parent, arr)
        return out

    def centres(self):
        """(Ncell,3) generator points of the REAL cells (central)."""
        return self.gen[:self.ncentral]

    def volumes(self):
        """(Ncell,) volumes of the REAL cells (central)."""
        return self.volume[:self.ncentral]

    def edges(self):
        """Voronoi cells are polyhedra, not axis-aligned boxes — there is no
        per-axis edge array. Returns the flat face/vertex CSR instead."""
        return self.cf_start, self.fv_start, self.fvx, self.fvy, self.fvz

    @classmethod
    def from_points(cls, points, bounds=None, periodic=False, dispersion=None,
                    ghost_pad=0.25):
        """Tessellate `points` (Ncell,3) into Voronoi cells via pyvoro.

        bounds     : [[xlo,xhi],[ylo,yhi],[zlo,zhi]] domain (default: padded
                     point extent). pyvoro requires every point strictly inside.
        periodic   : bool or (3,) — a true PERIODIC tessellation via GHOST IMAGES
                     (see `_from_points_periodic`): boundary generators are
                     replicated into their neighbouring periodic images, the lot
                     is tessellated non-periodically, and the central cells are
                     then the exact periodic Voronoi cells (Sum volume = box
                     volume exactly, no boundary padding). The ghost cells are
                     retained for the deposit (folded back into their parent by
                     `fold_to_central`), so a kernel crossing a face wraps.
        dispersion : pyvoro block-sizing hint (~ the typical generator spacing);
                     default = the mean spacing (domain_volume / N)**(1/3). This
                     is what sizes voro++'s internal block grid — setting it to a
                     large value (e.g. the box size) collapses the mesh into one
                     giant block and makes the build O(N^2) (and memory-hungry).
        ghost_pad  : periodic only — the image shell width as a fraction of the
                     box per axis (generators within this distance of a boundary
                     are replicated across it). Acts as a FLOOR: the shell is
                     auto-grown to at least 3x the mean generator spacing, so
                     SPARSE meshes (few generators -> large cells) still close
                     their boundary cells correctly (otherwise Sum cell vol > box
                     vol, e.g. +5% at N=32). 0.25 already covers dense/near-
                     uniform meshes; raise it only for strongly clustered ones.
        """
        import pyvoro
        pts = np.ascontiguousarray(np.asarray(points, dtype=np.float64).reshape(-1, 3))
        N = pts.shape[0]
        if N == 0:
            raise ValueError("need at least one generator point")
        if bounds is None:
            plo = pts.min(axis=0); phi = pts.max(axis=0)
            span = phi - plo
            pad = 0.05 * np.where(span > 0, span, 1.0) + 1e-6
            lo = plo - pad; hi = phi + pad
        else:
            lo = np.array([float(bounds[a][0]) for a in range(3)])
            hi = np.array([float(bounds[a][1]) for a in range(3)])
        per = periodic
        if np.isscalar(per):
            per = [bool(per)] * 3
        else:
            per = [bool(x) for x in per]

        if any(per):
            return cls._from_points_periodic(pts, lo, hi, per, dispersion,
                                             ghost_pad)

        limits = [[lo[a], hi[a]] for a in range(3)]
        if dispersion is None:
            dispersion = float((np.prod(hi - lo) / N) ** (1.0 / 3.0))
        cells = pyvoro.compute_voronoi(pts.tolist(), limits, dispersion,
                                       radii=[], periodic=per)
        return cls._from_pyvoro(cells, lo, hi, per)

    @classmethod
    def _from_points_periodic(cls, pts, lo, hi, per, dispersion, ghost_pad):
        """Periodic Voronoi by GHOST IMAGES: replicate near-boundary generators
        into their periodic images, tessellate the lot non-periodically over a
        padded box, keep the N central cells as the real (periodic) cells and the
        images as ghosts (`cell_parent` -> central index). The central cells are
        complete and exact, and tile EXACTLY one period (Sum volume = box vol)."""
        import pyvoro
        N = pts.shape[0]
        L = hi - lo
        # image shell half-width. The shell MUST span the Voronoi-neighbour reach
        # (~ a few generator spacings) or boundary cells are not cut by their
        # periodic images -> they bloat and Sum(cell vol) > box vol (badly wrong
        # for SPARSE meshes, e.g. +5% at N=32 with the fixed 0.25). So grow the
        # shell to at least 3x the mean generator spacing; the user's ghost_pad
        # fraction acts as a floor (and still covers near-uniform dense meshes).
        spacing = float((np.prod(L) / max(N, 1)) ** (1.0 / 3.0))
        margin = np.maximum(ghost_pad * L, 3.0 * spacing)
        # central generators first, then the periodic images that fall inside the
        # padded box [lo-margin, hi+margin].
        gens = [pts]
        parents = [np.arange(N, dtype=np.int64)]
        offs = [(-1, 0, 1) if per[a] else (0,) for a in range(3)]
        for ox in offs[0]:
            for oy in offs[1]:
                for oz in offs[2]:
                    if ox == 0 and oy == 0 and oz == 0:
                        continue
                    shift = np.array([ox * L[0], oy * L[1], oz * L[2]])
                    q = pts + shift
                    keep = np.ones(N, dtype=bool)
                    for a in range(3):
                        keep &= (q[:, a] >= lo[a] - margin[a]) & \
                                (q[:, a] <= hi[a] + margin[a])
                    if keep.any():
                        gens.append(q[keep])
                        parents.append(np.where(keep)[0].astype(np.int64))
        gen_all = np.ascontiguousarray(np.concatenate(gens, axis=0))
        parent_all = np.concatenate(parents)
        ntot = gen_all.shape[0]
        plo = lo - margin; phi = hi + margin
        limits = [[plo[a], phi[a]] for a in range(3)]
        if dispersion is None:
            dispersion = float((np.prod(phi - plo) / ntot) ** (1.0 / 3.0))
        cells = pyvoro.compute_voronoi(gen_all.tolist(), limits, dispersion,
                                       radii=[], periodic=[False, False, False])
        vg = cls._from_pyvoro(cells, lo, hi, per)        # lo/hi = the PHYSICAL box
        vg.cell_parent = np.ascontiguousarray(parent_all, dtype=np.int64)
        vg.ncentral = N
        return vg

    @classmethod
    def from_particles(cls, particles, bounds=None, periodic=None, dispersion=None,
                       ghost_pad=0.25):
        """Tessellate the particle positions themselves into Voronoi cells.

        bounds default to the particle box (centered) when known, else the padded
        particle extent; periodic defaults to particles.periodic. For a periodic
        box the bounds ARE the period (no edge padding — ghost images handle the
        boundary), so Sum cell volume = box volume exactly.
        """
        pos = np.ascontiguousarray(particles.pos, dtype=np.float64)
        if periodic is None:
            periodic = particles.periodic
        per_any = (periodic if np.isscalar(periodic) else any(periodic))
        if bounds is None and particles.box is not None:
            half = 0.5 * np.asarray(particles.box, dtype=np.float64).reshape(3)
            if per_any:                                  # period == the box, exactly
                bounds = [[-half[a], half[a]] for a in range(3)]
            else:
                plo = pos.min(axis=0); phi = pos.max(axis=0)
                lo = np.minimum(-half, plo) - 1e-6
                hi = np.maximum(half, phi) + 1e-6
                bounds = [[lo[a], hi[a]] for a in range(3)]
        return cls.from_points(pos, bounds=bounds, periodic=periodic,
                               dispersion=dispersion, ghost_pad=ghost_pad)

    @staticmethod
    def _from_pyvoro(cells, lo, hi, periodic):
        """Flatten a pyvoro `compute_voronoi` result into the CSR arrays."""
        N = len(cells)
        gen = np.empty((N, 3)); volume = np.empty(N); cell_rad = np.empty(N)
        cf_start = np.zeros(N + 1, dtype=np.int64)
        fv_start_list = [0]
        fvx, fvy, fvz = [], [], []
        nface = 0
        for c, cell in enumerate(cells):
            verts = np.asarray(cell["vertices"], dtype=np.float64)
            g = np.asarray(cell["original"], dtype=np.float64)
            gen[c] = g
            volume[c] = float(cell["volume"])
            cell_rad[c] = (np.sqrt(((verts - g) ** 2).sum(axis=1)).max()
                           if verts.shape[0] else 0.0)
            for face in cell["faces"]:
                vi = face["vertices"]
                if len(vi) < 3:
                    continue                       # skip degenerate faces
                for idx in vi:
                    fvx.append(verts[idx, 0])
                    fvy.append(verts[idx, 1])
                    fvz.append(verts[idx, 2])
                fv_start_list.append(len(fvx))
                nface += 1
            cf_start[c + 1] = nface
        return VoronoiGrid(
            gen=gen, volume=volume, cell_rad=cell_rad, cf_start=cf_start,
            fv_start=np.asarray(fv_start_list, dtype=np.int64),
            fvx=np.asarray(fvx), fvy=np.asarray(fvy), fvz=np.asarray(fvz),
            lo=lo, hi=hi, periodic=periodic)


def voronoi_candidates(gen, pos, radii):
    """Per-particle candidate cell lists for a Voronoi deposit, as flat CSR.

    For each particle, every generator within `radii[i]` of it is a candidate
    cell (a BallTree range search). Returns (cand_start, cand_count, cand_cell)
    where particle i's candidate cell indices are
    cand_cell[cand_start[i] : cand_start[i] + cand_count[i]].
    """
    from sklearn.neighbors import BallTree
    pos = np.ascontiguousarray(pos, dtype=np.float64)
    N = pos.shape[0]
    tree = BallTree(np.ascontiguousarray(gen, dtype=np.float64))
    ind = tree.query_radius(pos, r=np.ascontiguousarray(radii, dtype=np.float64))
    counts = np.fromiter((a.shape[0] for a in ind), dtype=np.int64, count=N)
    start = np.zeros(N, dtype=np.int64)
    if N:
        start[1:] = np.cumsum(counts)[:-1]
    cand = (np.concatenate(ind).astype(np.int64) if counts.sum()
            else np.zeros(0, dtype=np.int64))
    return start, counts, cand


@dataclass
class Projection2D:
    """A 2D pixel grid in the (axis1, axis2) plane for direct rendering.

    Fast SPLASH-style deposit (scales as npx^2 * N, no 3D cube):
      mode='column'     line-of-sight integral (surface density for 'rho');
      mode='average'    kernel-weighted LOS average (column(A)/column(1));
      mode='rhocolumn'  density-weighted LOS average (mass-weighted);
      mode='slice'      field value on the plane at `zslice` (3D kernel).

    axis1, axis2 : world axes (0,1,2) mapped to image x,y; the remaining axis is
                   the line of sight (`los_axis`).
    lo, hi       : (2,) image extent in the (axis1, axis2) plane.
    npx          : int or (2,) pixels.
    zslice       : LOS coordinate for mode='slice'.
    """
    axis1: int
    axis2: int
    lo: np.ndarray
    hi: np.ndarray
    npx: object = 512
    mode: str = "column"
    zslice: float = 0.0

    def __post_init__(self):
        if self.axis1 == self.axis2 or self.axis1 not in (0, 1, 2) or self.axis2 not in (0, 1, 2):
            raise ValueError("axis1, axis2 must be distinct axes in {0,1,2}")
        if self.mode not in ("column", "average", "rhocolumn", "slice"):
            raise ValueError(f"unknown mode {self.mode!r}")
        self.lo = np.asarray(self.lo, dtype=np.float64).reshape(2)
        self.hi = np.asarray(self.hi, dtype=np.float64).reshape(2)
        npx = self.npx
        if np.isscalar(npx):
            npx = (int(npx), int(npx))
        self.npx = tuple(int(x) for x in npx)
        if any(n <= 0 for n in self.npx):
            raise ValueError("npx must be positive")
        if np.any(self.hi <= self.lo):
            raise ValueError("hi must exceed lo on both axes")

    @classmethod
    def centered(cls, box, axis1, axis2, npx=512, mode="column", zslice=0.0):
        """Plane spanning the (axis1, axis2) extent of a centered box [-box/2,box/2]."""
        box = np.asarray(box, dtype=np.float64).reshape(-1)
        if box.size == 1:
            box = np.repeat(box, 3)
        half = 0.5 * np.array([box[axis1], box[axis2]])
        return cls(axis1=axis1, axis2=axis2, lo=-half, hi=half, npx=npx,
                   mode=mode, zslice=zslice)

    @property
    def los_axis(self):
        return ({0, 1, 2} - {self.axis1, self.axis2}).pop()

    @property
    def box(self):
        return self.hi - self.lo

    @property
    def pixwidth(self):
        return self.box / np.asarray(self.npx, dtype=np.float64)

    def edges(self):
        return [np.linspace(self.lo[a], self.hi[a], self.npx[a] + 1) for a in range(2)]

    def centres(self):
        pw = self.pixwidth
        return [self.lo[a] + (np.arange(self.npx[a]) + 0.5) * pw[a] for a in range(2)]
