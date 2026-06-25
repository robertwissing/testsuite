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
import os
import contextlib
import subprocess
import tempfile
import numpy as np
import numba

try:                                  # optional fast parallel mesher (OpenMP voro++)
    import multivoro as _multivoro
    _HAS_MULTIVORO = True
except Exception:                     # pragma: no cover - optional dependency
    _multivoro = None
    _HAS_MULTIVORO = False

# optional distributed mesher (MadVoro, MPI/serial) — used via a shell-out driver
# compiled in _madvoro/ (no Python binding; see _madvoro/README.md). The driver
# is a subprocess so a mesher crash on pathological input cannot take down the
# interpreter (a deliberate robustness win over an in-process binding).
_MADVORO_DRIVER = os.path.join(os.path.dirname(__file__), "_madvoro", "madvoro_driver")
_HAS_MADVORO = os.path.isfile(_MADVORO_DRIVER) and os.access(_MADVORO_DRIVER, os.X_OK)

# optional AREPO Voronoi mesher (codes/arepo) — also a shell-out (write HDF5 snapshot
# -> run `Arepo param.txt 14` mesh-dump -> parse), in _arepo/arepo_mesh.py. AREPO is
# the reference moving-mesh code; it tessellates collapsed high-dynamic-range cores
# robustly and ~7.5x faster than serial madvoro on them. See `arepo-mesh-backend`.
try:
    from ._arepo import arepo_mesh as _arepo_mesh
    _HAS_AREPO = _arepo_mesh.has_arepo()
except Exception:                     # pragma: no cover - optional dependency
    _arepo_mesh = None
    _HAS_AREPO = False

_DEFAULT_THREADS = os.cpu_count() or 1

# De-sliver floor as a fraction of the median nearest-neighbour spacing, PER
# MESHER. multivoro (parallel voro++) needs ~4x pyvoro's floor on a collapsed,
# high-dynamic-range core: its independent per-cell search loses completeness when
# the local spacing collapses abruptly, so cells are UNDER-CUT and OVERLAP (mass
# error +2260% / 100% core overlap at frac 0.05, vs pyvoro's graceful +0.21%).
# A larger floor pulls the core dynamic range back into the regime multivoro can
# resolve. Serial pyvoro keeps voro++'s completeness guarantee, so 0.05 (just
# enough to clear true slivers) conserves. See `multivoro-dynamic-range-defect`.
# madvoro starts at 0.0 (de-sliver OFF) deliberately: an open question whether its
# distributed Delaunay handles the collapsed high-dynamic-range core natively or
# breaks like multivoro. Measure first, then set the floor from data.
# arepo at 0.0: AREPO's Delaunay/Voronoi handles the collapsed high-dynamic-range
# core natively (its native regime), like madvoro -- no de-sliver needed; the
# ConvexHull reconstruction is exact on the resulting cells. Measured below.
_DESLIVER_FRAC = {"pyvoro": 0.05, "multivoro": 0.20, "madvoro": 0.0, "arepo": 0.0}


def _resolve_backend(backend):

    if backend == "arepo":
        return "arepo"
    if backend == "auto":
        return "multivoro"
    if backend == "multivoro" or (backend == "auto" and _HAS_MULTIVORO):
        return "multivoro"
    return "pyvoro"


@numba.njit(cache=True)
def _mv_flatten(VV, vstart, FF, fstart, genx, geny, genz):
    """Build the Voronoi CSR (cf_start, fv_start, fvx/fvy/fvz, volume, cell_rad)
    from multivoro's concatenated per-cell vertices (VV, offsets vstart) and flat
    face-vertex runs (FF = [k, local_ids..k, ...] per cell, offsets fstart). Face
    vertex indices are LOCAL to the cell. Volume = signed tets fanned from the
    generator; degenerate faces (k<3) skipped. The hot inner loop, in numba."""
    N = vstart.shape[0] - 1
    nface = 0
    nfv = 0
    for c in range(N):                          # pass 1: count faces (k>=3) + verts
        i = fstart[c]; end = fstart[c + 1]
        while i < end:
            k = FF[i]; i += 1 + k
            if k >= 3:
                nface += 1; nfv += k
    cf_start = np.zeros(N + 1, np.int64)
    fv_start = np.zeros(nface + 1, np.int64)
    fvx = np.empty(nfv); fvy = np.empty(nfv); fvz = np.empty(nfv)
    volume = np.empty(N); cell_rad = np.empty(N)
    fidx = 0; vptr = 0
    for c in range(N):                          # pass 2: fill
        gx = genx[c]; gy = geny[c]; gz = genz[c]
        vs = vstart[c]; nv = vstart[c + 1] - vs
        rmax = 0.0
        for j in range(nv):
            dx = VV[vs + j, 0] - gx; dy = VV[vs + j, 1] - gy; dz = VV[vs + j, 2] - gz
            d = dx * dx + dy * dy + dz * dz
            if d > rmax:
                rmax = d
        cell_rad[c] = rmax ** 0.5
        vol6 = 0.0
        i = fstart[c]; end = fstart[c + 1]
        while i < end:
            k = FF[i]; base = i + 1; i += 1 + k
            if k < 3:
                continue
            for t in range(k):
                gi = vs + FF[base + t]
                fvx[vptr] = VV[gi, 0]; fvy[vptr] = VV[gi, 1]; fvz[vptr] = VV[gi, 2]
                vptr += 1
            fidx += 1
            fv_start[fidx] = vptr
            i0 = vs + FF[base]
            a0x = VV[i0, 0] - gx; a0y = VV[i0, 1] - gy; a0z = VV[i0, 2] - gz
            for t in range(1, k - 1):
                i1 = vs + FF[base + t]; i2 = vs + FF[base + t + 1]
                a1x = VV[i1, 0] - gx; a1y = VV[i1, 1] - gy; a1z = VV[i1, 2] - gz
                a2x = VV[i2, 0] - gx; a2y = VV[i2, 1] - gy; a2z = VV[i2, 2] - gz
                vol6 += (a0x * (a1y * a2z - a1z * a2y)
                         - a0y * (a1x * a2z - a1z * a2x)
                         + a0z * (a1x * a2y - a1y * a2x))
        volume[c] = abs(vol6) / 6.0
        cf_start[c + 1] = fidx
    return cf_start, fv_start, fvx, fvy, fvz, volume, cell_rad


@contextlib.contextmanager
def _silence_cstd():
    """Mute C-level stdout+stderr (voro++ prints 'Order 4 vertex memory scaled
    up...' progress noise per cell) for the duration of the block."""
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved = (os.dup(1), os.dup(2))
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        os.close(devnull)
        os.close(saved[0])
        os.close(saved[1])


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

    @staticmethod
    def _median_nn(pts):
        """Median nearest-neighbour distance of the generators — a robust local
        resolution scale for the de-sliver floor. Preferred over the mean spacing
        (domain_vol/N)**(1/3) because that is inflated by empty volume on a
        CONCENTRATED system (a galaxy/merger with a dense disk in a large empty
        halo), where it overshoots the real spacing ~100x and the floor would
        merge away genuine structure. The median NN distance instead tracks the
        typical spacing of where the generators actually are. Ignores exact
        duplicates (NN==0) so a handful of coincident points don't drag it to 0."""
        from sklearn.neighbors import BallTree
        n = pts.shape[0]
        if n < 2:
            return 0.0
        d, _ = BallTree(pts).query(pts, k=2)
        nn = d[:, 1]
        pos = nn[nn > 0.0]
        if pos.size == 0:
            return 0.0
        return float(np.median(pos))

    @staticmethod
    def _merge_close(pts, tol):
        """Collapse near-coincident generators (sub-`tol` clusters) into one
        point at their centroid, so voro++ does not have to cut degenerate sliver
        faces between them. Returns `(newpts, nmerged)`.

        WHY: on a collapsed core (high dynamic range) the densest particles sit
        at separations << the smoothing length. pyvoro/voro++ cannot tessellate
        such near-coincident generators cleanly — the bisecting faces are
        ill-conditioned and the stored polyhedra develop local GAPS and OVERLAPS
        (they no longer partition space). The signed-tet integrator is exact on
        well-formed cells, so a broken partition is the sole remaining source of
        the Voronoi Petkova conservation error (~0.2% on mhdcollapse). Merging the
        coincident generators removes the slivers and restores Sum frac = 1.
        Connected components are merged transitively (union-find over all pairs
        within `tol`), so a tight chain collapses to a single representative."""
        from scipy.spatial import cKDTree
        n = pts.shape[0]
        if tol <= 0.0 or n < 2:
            return pts, 0
        pairs = cKDTree(pts).query_pairs(tol, output_type='ndarray')
        if len(pairs) == 0:
            return pts, 0
        parent = np.arange(n)

        def find(a):
            root = a
            while parent[root] != root:
                root = parent[root]
            while parent[a] != root:           # path compression
                parent[a], a = root, parent[a]
            return root

        for a, b in pairs:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
        roots = np.array([find(i) for i in range(n)], dtype=np.int64)
        _, inv = np.unique(roots, return_inverse=True)
        m = inv.max() + 1
        newpts = np.zeros((m, 3))
        np.add.at(newpts, inv, pts)
        newpts /= np.bincount(inv, minlength=m)[:, None]
        return np.ascontiguousarray(newpts), n - m

    @classmethod
    def from_points(cls, points, bounds=None, periodic=False, dispersion=None,
                    ghost_pad=0.25, merge_tol=None, backend="auto", n_threads=-1):
        """Tessellate `points` (Ncell,3) into Voronoi cells.

        backend    : mesher — 'auto' (multivoro if installed for parallel OpenMP
                     voro++, else pyvoro), 'multivoro', or 'pyvoro'. multivoro is
                     ~2x faster single-threaded and up to ~40x with threads, and is
                     the default when available; identical cells out.
        n_threads  : multivoro thread count; <=0 -> all cores.

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
        merge_tol  : collapse generators closer than this distance into one (see
                     `_merge_close`) BEFORE tessellating, to keep voro++ from
                     emitting degenerate sliver cells (or failing outright on exact
                     duplicates) on a high-dynamic-range core. Default None -> a
                     BACKEND-AWARE fraction of the MEDIAN nearest-neighbour distance
                     (`_median_nn` x `_DESLIVER_FRAC`): 5% for pyvoro, 20% for
                     multivoro (whose parallel solver under-cuts/overlaps core cells
                     under abrupt-spacing dynamic range and so needs the spacing
                     pulled up further). An absolute resolution floor that bites only
                     near-coincident generators and is a no-op on well-resolved
                     meshes. The larger multivoro floor merges ~10% more core cells,
                     which caps the very densest cell's peak density (the volume-
                     weighted field PDF is unchanged); set merge_tol explicitly to
                     override. NOTE: a
                     LOCAL/relative floor cannot work here: the slivers form exactly
                     where the local spacing collapses, so a local threshold shrinks
                     with it and never reaches them (and merging them unevenly makes
                     voro++ fail) — the floor must be absolute. median-NN (not the
                     mean spacing (vol/N)**(1/3)) keeps it robust on CONCENTRATED
                     systems where empty volume inflates the mean ~100x. 0 disables
                     the sliver de-merge but EXACT/near-coincident generators are
                     STILL collapsed (an unconditional ~1e-9*span floor) so voro++
                     never crashes on duplicates. NOTE: merging means Ncell may be
                     < N (the core is represented by fewer cells); intended — those
                     were duplicate sample points.
        """
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

        # EXACT/near-coincident generators are ALWAYS collapsed, independent of
        # merge_tol (incl. merge_tol=0): two generators voro++ cannot tell apart
        # have no valid distinct cells and make it FAIL outright ("number of cells
        # found was not equal to the number of particles"). This floor sits ~1e-9
        # of the domain span -- ~1e3x above voro++'s internal ~1e-12 relative
        # tolerance, far below any real spacing -- so it removes only true
        # duplicates, never structure. (merge_tol=0 thus still BUILDS on data with
        # duplicates; it only switches off the larger sliver de-merge below.)
        tol_dup = 1e-9 * float(np.max(hi - lo))
        # de-sliver: merge near-coincident generators so the mesh stays a clean
        # partition (the collapsed-core conservation fix). Default is a BACKEND-
        # AWARE fraction of median-NN: multivoro needs ~4x pyvoro's floor (its
        # parallel solver under-cuts/overlaps cells on an abrupt-spacing core where
        # serial pyvoro stays valid) -- see `_DESLIVER_FRAC`.
        if merge_tol is None:
            merge_tol = _DESLIVER_FRAC[_resolve_backend(backend)] * cls._median_nn(pts)
        tol_eff = max(float(merge_tol), tol_dup)
        if tol_eff > 0.0:
            pts, _ = cls._merge_close(pts, tol_eff)
            N = pts.shape[0]

        if any(per):
            return cls._from_points_periodic(pts, lo, hi, per, dispersion,
                                             ghost_pad, backend, n_threads)

        return cls._tessellate(pts, lo, hi, lo, hi, dispersion, backend, n_threads)

    @classmethod
    def _from_points_periodic(cls, pts, lo, hi, per, dispersion, ghost_pad,
                              backend="auto", n_threads=-1):
        """Periodic Voronoi by GHOST IMAGES: replicate near-boundary generators
        into their periodic images, tessellate the lot non-periodically over a
        padded box, keep the N central cells as the real (periodic) cells and the
        images as ghosts (`cell_parent` -> central index). The central cells are
        complete and exact, and tile EXACTLY one period (Sum volume = box vol)."""
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
        plo = lo - margin; phi = hi + margin
        # tessellate the padded ghosted set; store the grid on the PHYSICAL box.
        vg = cls._tessellate(gen_all, plo, phi, lo, hi, dispersion,
                             backend, n_threads)
        vg.periodic = tuple(bool(x) for x in per)
        vg.cell_parent = np.ascontiguousarray(parent_all, dtype=np.int64)
        vg.ncentral = N
        return vg

    @classmethod
    def from_particles(cls, particles, bounds=None, periodic=None, dispersion=None,
                       ghost_pad=0.25, merge_tol=None, backend="auto", n_threads=-1):
        """Tessellate the particle positions themselves into Voronoi cells.

        bounds default to the particle box (centered) when known, else the padded
        particle extent; periodic defaults to particles.periodic. For a periodic
        box the bounds ARE the period (no edge padding — ghost images handle the
        boundary), so Sum cell volume = box volume exactly.

        merge_tol forwards to `from_points` (de-sliver coincident particles on a
        collapsed core; default a backend-aware fraction of median-NN, 0 to disable).
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
                               dispersion=dispersion, ghost_pad=ghost_pad,
                               merge_tol=merge_tol, backend=backend,
                               n_threads=n_threads)

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

    @staticmethod
    def _from_multivoro(cells, pts, lo, hi, periodic):
        """Flatten a multivoro result into the CSR arrays. multivoro returns a
        list of `Cell` (get_vertices/get_face_vertices/get_neighbors) in INPUT
        ORDER, but gives neither the generator nor the cell volume — so the
        generator is the input point and the volume is computed here from the
        polyhedron faces (signed tets fanned from the generator). `get_face_
        vertices()` is a flat list [k0, v..(k0), k1, v..(k1), ...] of per-face
        vertex-index runs. The signed-tet Petkova integrator re-derives face
        orientation from the generator, so the face winding need not match
        pyvoro's."""
        N = len(cells)
        gen = np.ascontiguousarray(np.asarray(pts, dtype=np.float64)[:N])
        # Collect the per-cell vertices/faces via the (fast) nanobind accessors,
        # concatenated with offsets, then do the CSR assembly + volume in numba
        # (the per-vertex/per-tet work is the cost; the accessors are ~0.1s/40k).
        vlist = [None] * N
        flist = [None] * N
        vcounts = np.empty(N, dtype=np.int64)
        fcounts = np.empty(N, dtype=np.int64)
        for c in range(N):
            cell = cells[c]
            v = np.ascontiguousarray(cell.get_vertices(), dtype=np.float64)
            f = np.ascontiguousarray(cell.get_face_vertices(), dtype=np.int64)
            vlist[c] = v; flist[c] = f
            vcounts[c] = v.shape[0]; fcounts[c] = f.shape[0]
        VV = (np.concatenate(vlist, axis=0) if N else np.zeros((0, 3)))
        FF = (np.concatenate(flist) if N else np.zeros(0, dtype=np.int64))
        vstart = np.zeros(N + 1, dtype=np.int64); vstart[1:] = np.cumsum(vcounts)
        fstart = np.zeros(N + 1, dtype=np.int64); fstart[1:] = np.cumsum(fcounts)
        cf_start, fv_start, fvx, fvy, fvz, volume, cell_rad = _mv_flatten(
            VV, vstart, FF, fstart,
            np.ascontiguousarray(gen[:, 0]), np.ascontiguousarray(gen[:, 1]),
            np.ascontiguousarray(gen[:, 2]))
        return VoronoiGrid(
            gen=gen, volume=volume, cell_rad=cell_rad, cf_start=cf_start,
            fv_start=fv_start, fvx=fvx, fvy=fvy, fvz=fvz,
            lo=lo, hi=hi, periodic=periodic)

    @classmethod
    def _from_madvoro(cls, pts, dom_lo, dom_hi, grid_lo, grid_hi, periodic,
                      n_threads):
        """Tessellate via the MadVoro shell-out driver (subprocess: a mesher crash
        on pathological input is isolated from the interpreter). Writes points +
        box to a binary temp file, runs the driver, reads back the per-cell CSR.
        The driver outputs ABSOLUTE inlined face vertices + per-cell volume and
        circumradius (so the assembly here is just cumsums -- no per-cell Python
        loop, fast at 1e6 cells). See `_madvoro/madvoro_driver.cpp`.

        n_threads >1 launches `mpirun -np n_threads` (needs the MPI-built driver,
        `MPICXX=mpicxx bash _madvoro/build.sh`); otherwise runs the driver serially."""
        if not _HAS_MADVORO:
            raise RuntimeError(
                "backend='madvoro' needs the compiled driver at %s -- build it "
                "with `bash setupfiles/sph_interp/_madvoro/build.sh`" % _MADVORO_DRIVER)
        pts = np.ascontiguousarray(np.asarray(pts, dtype=np.float64).reshape(-1, 3))
        N = pts.shape[0]
        box = np.array([dom_lo[0], dom_lo[1], dom_lo[2],
                        dom_hi[0], dom_hi[1], dom_hi[2]], dtype=np.float64)
        with tempfile.TemporaryDirectory() as td:
            fin = os.path.join(td, "in.bin")
            fout = os.path.join(td, "out.bin")
            with open(fin, "wb") as f:
                np.array([N], dtype=np.int64).tofile(f)
                box.tofile(f)
                pts.tofile(f)
            cmd = [_MADVORO_DRIVER, fin, fout]
            if n_threads and n_threads > 1:
                cmd = ["mpirun", "-np", str(int(n_threads))] + cmd
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0 or not os.path.isfile(fout):
                raise RuntimeError("madvoro driver failed (rc=%d): %s"
                                   % (r.returncode, r.stderr[-500:]))
            # the manifest holds the rank count; each rank wrote `${out}.${rank}`
            # with its disjoint subset of cells -> read and concatenate them (every
            # cell/face self-describes its counts + inlines absolute coords, so the
            # CSR concatenates with no re-indexing; cumsums are built below).
            with open(fout, "rb") as f:
                nranks = int(np.fromfile(f, dtype=np.int64, count=1)[0])
            gens_l, vol_l, crad_l, nface_l, nvert_l, fvxyz_l = [], [], [], [], [], []
            for rk in range(nranks):
                with open(fout + "." + str(rk), "rb") as f:
                    nc, tf, tv = np.fromfile(f, dtype=np.int64, count=3)
                    gens_l.append(np.fromfile(f, dtype=np.float64, count=nc * 3).reshape(nc, 3))
                    vol_l.append(np.fromfile(f, dtype=np.float64, count=nc))
                    crad_l.append(np.fromfile(f, dtype=np.float64, count=nc))
                    nface_l.append(np.fromfile(f, dtype=np.int64, count=nc))
                    nvert_l.append(np.fromfile(f, dtype=np.int64, count=tf))
                    fvxyz_l.append(np.fromfile(f, dtype=np.float64, count=tv * 3).reshape(tv, 3))
        gen = np.concatenate(gens_l, axis=0) if gens_l else np.zeros((0, 3))
        volume = np.concatenate(vol_l) if vol_l else np.zeros(0)
        cell_rad = np.concatenate(crad_l) if crad_l else np.zeros(0)
        nface = np.concatenate(nface_l) if nface_l else np.zeros(0, np.int64)
        nvert = np.concatenate(nvert_l) if nvert_l else np.zeros(0, np.int64)
        fvxyz = np.concatenate(fvxyz_l, axis=0) if fvxyz_l else np.zeros((0, 3))
        ncell = gen.shape[0]; tface = nvert.shape[0]
        cf_start = np.zeros(ncell + 1, dtype=np.int64); cf_start[1:] = np.cumsum(nface)
        fv_start = np.zeros(tface + 1, dtype=np.int64); fv_start[1:] = np.cumsum(nvert)
        return VoronoiGrid(
            gen=np.ascontiguousarray(gen), volume=volume, cell_rad=cell_rad,
            cf_start=cf_start, fv_start=fv_start,
            fvx=np.ascontiguousarray(fvxyz[:, 0]),
            fvy=np.ascontiguousarray(fvxyz[:, 1]),
            fvz=np.ascontiguousarray(fvxyz[:, 2]),
            lo=grid_lo, hi=grid_hi, periodic=periodic)

    @classmethod
    def _from_arepo(cls, pts, dom_lo, dom_hi, grid_lo, grid_hi, periodic,
                    n_threads):
        """Tessellate via AREPO (codes/arepo) through the `_arepo` shell-out helper
        (write HDF5 snapshot -> run `Arepo param.txt 14` mesh dump -> parse). AREPO
        runs PERIODIC over a cubic box (its native regime); each cell is rebuilt as
        the ConvexHull of its periodic-unwrapped vertices -> exact volume + clean
        faces (AREPO's raw face lists need it; see `_arepo/arepo_mesh.py`).
        n_threads>1 launches `mpirun -np n_threads` (AREPO's MPI works)."""
        if not _HAS_AREPO:
            raise RuntimeError(
                "backend='arepo' needs the built Arepo executable in codes/arepo "
                "(see setupfiles/sph_interp/_arepo/arepo_mesh.py)")
        m = _arepo_mesh.tessellate(pts, dom_lo, dom_hi, n_threads=n_threads)
        return VoronoiGrid(
            gen=np.ascontiguousarray(m["gen"]), volume=m["volume"],
            cell_rad=m["cell_rad"], cf_start=m["cf_start"], fv_start=m["fv_start"],
            fvx=m["fvx"], fvy=m["fvy"], fvz=m["fvz"],
            lo=grid_lo, hi=grid_hi, periodic=periodic)

    @classmethod
    def _tessellate(cls, pts, dom_lo, dom_hi, grid_lo, grid_hi, dispersion,
                    backend, n_threads):
        """Raw NON-periodic tessellation of `pts` over the mesher domain
        [dom_lo, dom_hi] -> CSR VoronoiGrid stored with bounds [grid_lo, grid_hi]
        (the two differ for the periodic ghost path: padded mesher box vs physical
        box). `backend`: 'auto' (madvoro if its driver is built, else multivoro if
        installed, else pyvoro), 'madvoro' (distributed, via the compiled shell-out
        driver), 'multivoro', or 'pyvoro'. `n_threads`: multivoro thread count, or
        madvoro MPI rank count (>1 -> mpirun); <=0 -> all cores (multivoro)."""
        N = pts.shape[0]
        eff = _resolve_backend(backend)          # auto -> madvoro/multivoro/pyvoro
        if eff == "arepo":
            return cls._from_arepo(pts, dom_lo, dom_hi, grid_lo, grid_hi,
                                   (False, False, False), n_threads)
        if eff == "madvoro":
            return cls._from_madvoro(pts, dom_lo, dom_hi, grid_lo, grid_hi,
                                     (False, False, False), n_threads)
        use_mv = eff == "multivoro"
        if backend == "multivoro" and not _HAS_MULTIVORO:
            raise ImportError("backend='multivoro' but the multivoro package is "
                              "not installed")
        if use_mv:
            nt = int(n_threads) if n_threads and n_threads > 0 else _DEFAULT_THREADS
            lim = np.array([dom_lo, dom_hi], dtype=np.float64)   # (2,3): lower,upper
            with _silence_cstd():
                cells = _multivoro.compute_voronoi(
                    np.ascontiguousarray(pts, dtype=np.float64), limits=lim,
                    radii=np.zeros(N), periodic_boundaries=(False, False, False),
                    n_threads=nt)
            return cls._from_multivoro(cells, pts, grid_lo, grid_hi,
                                       (False, False, False))
        import pyvoro
        if dispersion is None:
            dispersion = float((np.prod(np.asarray(dom_hi) - np.asarray(dom_lo))
                                / max(N, 1)) ** (1.0 / 3.0))
        limits = [[dom_lo[a], dom_hi[a]] for a in range(3)]
        cells = pyvoro.compute_voronoi(pts.tolist(), limits, dispersion,
                                       radii=[], periodic=[False, False, False])
        return cls._from_pyvoro(cells, grid_lo, grid_hi, (False, False, False))


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
