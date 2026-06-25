"""Export an AMR interpolation result to other codes / mesh formats.

The library's input side converts any code *into* the canonical `Particles`
container (`from_tipsy`, `from_gadget`, ...). The output side is the mirror: an
AMR `GridResult` is a code-agnostic description of an octree mesh — a flat list
of leaf cells with centre, size, refinement level, base-cell linkage and one
value per field — and per-code **export adapters** translate *out of* it.

`to_cell_table` is the universal intermediate every exporter builds from: a plain
dict of equal-length arrays (one row per leaf). From it a writer for a specific
code is a thin function — group by `level` for patch/grid codes (Enzo/FLASH),
walk `base`+`level` for octree codes (RAMSES), or dump straight to HDF5/npz/VTK
for post-processing/visualization. `to_ramses` writes a self-describing HDF5 of
the RAMSES hydro variables now; `to_vtk`/`to_enzo` and the RAMSES native
amr_/hydro_ binary backend are declared as documented extension points.

Magnetic fields & RAMSES (constrained transport): RAMSES MHD stores B
FACE-CENTERED — six components per cell (left/right faces in x,y,z) — to keep
div B = 0. The AMR interpolation here is cell-CENTERED, so a CT-correct B export
needs a face-centered deposit (sample each B component at the relevant cell-face
centres) feeding six extra hydro variables. The hydro-variable list is
order-driven and extensible, so adding those is additive — no restructuring.
"""

import numpy as np

from .targets import AMRGrid, UniformGrid, VoronoiGrid


def to_cell_table(result, fields=None):
    """Flatten a `GridResult` into a code-agnostic per-cell table.

    result : a `GridResult` whose target is an `AMRGrid` (one row per leaf) or a
             `UniformGrid` (one row per cell, flattened in C order).
    fields : field names to include (default: every data field except the
             internal `_mass`).

    Returns dict[str -> (Ncell,) ndarray] — the canonical intermediate for any
    mesh-code exporter. Geometry keys: 'x','y','z' (cell centres), 'dx','dy','dz'
    (cell sizes), 'volume', 'level' (0=base/uniform), 'base' (flat base-grid
    index — the octree-root linkage; the cell index for a uniform grid). Plus one
    column per requested field. Trivially saved with `np.savez`, written to HDF5,
    or turned into a structured array / DataFrame by the caller.
    """
    grid = result.target
    if isinstance(grid, AMRGrid):
        c = grid.centres()
        s = grid.cell_size
        table = {
            "x": np.ascontiguousarray(c[:, 0]),
            "y": np.ascontiguousarray(c[:, 1]),
            "z": np.ascontiguousarray(c[:, 2]),
            "dx": np.ascontiguousarray(s[:, 0]),
            "dy": np.ascontiguousarray(s[:, 1]),
            "dz": np.ascontiguousarray(s[:, 2]),
            "volume": grid.volumes(),
            "level": grid.cell_level,
            "base": grid.cell_base,
        }
        flat = (lambda a: np.asarray(a))            # already (Ncell,)
    elif isinstance(grid, UniformGrid):
        cx, cy, cz = grid.centres()
        X, Y, Z = np.meshgrid(cx, cy, cz, indexing="ij")
        n = X.size
        pw = grid.pixwidth
        table = {
            "x": X.ravel(), "y": Y.ravel(), "z": Z.ravel(),
            "dx": np.full(n, pw[0]), "dy": np.full(n, pw[1]),
            "dz": np.full(n, pw[2]),
            "volume": np.full(n, grid.cell_volume),
            "level": np.zeros(n, dtype=np.int64),
            "base": np.arange(n, dtype=np.int64),
        }
        flat = (lambda a: np.asarray(a).ravel())    # (nx,ny,nz) -> (n,), C order
    elif isinstance(grid, VoronoiGrid):
        c = grid.centres()                          # generators (== cell seeds)
        v = grid.volumes()
        s = np.cbrt(v)                              # nominal cell size (vol^(1/3))
        n = c.shape[0]
        table = {
            "x": np.ascontiguousarray(c[:, 0]),
            "y": np.ascontiguousarray(c[:, 1]),
            "z": np.ascontiguousarray(c[:, 2]),
            "dx": s, "dy": s.copy(), "dz": s.copy(),
            "volume": v,
            "level": np.zeros(n, dtype=np.int64),
            "base": np.arange(n, dtype=np.int64),
        }
        flat = (lambda a: np.asarray(a))            # already (Ncell,) [central]
    else:
        raise TypeError("to_cell_table expects an AMRGrid, UniformGrid or "
                        f"VoronoiGrid GridResult (got target {type(grid).__name__})")
    if fields is None:
        fields = [k for k in result.data if not k.startswith("_")]
    for k in fields:
        table[k] = flat(result.data[k])
    return table


def to_structured_array(result, fields=None):
    """`to_cell_table` as a single numpy structured array (one record per leaf).

    Convenient for `np.save`, `tofile`, or as a drop-in mesh record array.
    """
    table = to_cell_table(result, fields=fields)
    dt = [(k, v.dtype) for k, v in table.items()]
    out = np.empty(next(iter(table.values())).shape[0], dtype=dt)
    for k, v in table.items():
        out[k] = v
    return out


# --------------------------------------------------------------------------- #
#  Per-code export adapters — extension points (NOT yet implemented).
#  Each consumes `to_cell_table(result)` and writes the target format. The
#  canonical table already carries everything these need (centre, size, level,
#  base linkage, fields); implementing one is a format-serialization task, not a
#  re-derivation of the mesh.
# --------------------------------------------------------------------------- #

def to_vtk(result, path, fields=None):
    """Write the leaves as a VTK unstructured grid of hexahedra (ParaView/VisIt).

    Contract: 8 corner points per leaf from (x,y,z)±(dx,dy,dz)/2, VTK_HEXAHEDRON
    cells, fields as CELL_DATA. Pure ASCII `.vtk` needs no extra deps.
    """
    raise NotImplementedError(
        "to_vtk: build hexahedra from to_cell_table(result) and write a VTK "
        "UNSTRUCTURED_GRID with the fields as CELL_DATA.")


#: Canonical RAMSES hydro-variable order (ideal hydro, ndim=3). For MHD, RAMSES
#: inserts the SIX face-centered magnetic components (Bx_l,By_l,Bz_l,Bx_r,By_r,
#: Bz_r) BETWEEN the velocities and the pressure -> nvar becomes 11; provide them
#: as named fields once a face-centered deposit exists (see module note on the
#: constrained-transport representation). The writer is order-driven, so adding B
#: is just extending `hydro_fields`.
RAMSES_HYDRO_VARS = ("rho", "vx", "vy", "vz", "P")

#: RAMSES MHD variable order: the six FACE-centred B components (left x/y/z, then
#: right x/y/z) flank the pressure, as in RAMSES's hydro output. Produce the B
#: fields with `solenoidal.project_solenoidal` (cell-centred B -> div-free faces)
#: before exporting with this ordering.
RAMSES_MHD_VARS = ("rho", "vx", "vy", "vz",
                   "Bx_l", "By_l", "Bz_l", "P", "Bx_r", "By_r", "Bz_r")


def to_ramses(result, path, hydro_fields=RAMSES_HYDRO_VARS, gamma=5.0 / 3.0,
              boxlen=None, units=None, time=0.0, fmt="hdf5"):
    """Export an AMR interpolation to a RAMSES-flavored mesh file.

    result       : a `GridResult` whose target is an `AMRGrid`, carrying the
                   hydro fields (interpolate them first, e.g.
                   `interpolate(p, amr, values=["rho","vx","vy","vz","P"])`).
    hydro_fields : ordered RAMSES hydro variables to write (default
                   `RAMSES_HYDRO_VARS`). Each must be a field in `result.data`.
                   For MHD, extend this with the magnetic components in RAMSES
                   order (face-centered B between v and P).
    gamma,boxlen,units,time : run metadata. `boxlen` defaults to the max domain
                   side; `units` is an optional dict (e.g. unit_l/unit_d/unit_t).
    fmt          : "hdf5" (implemented; self-describing, reader-agnostic) or
                   "native" (RAMSES amr_/hydro_ Fortran binary — NOT yet
                   implemented; needs a RAMSES reader to validate against).

    HDF5 layout (cell-centered leaf list — the data, decoupled from RAMSES's
    on-disk oct nesting, which a native backend would rebuild from level+base):
      /cells/{x,y,z,dx,level,base,ramses_level}   leaf geometry
      /hydro/data                                 (Ncell, nvar) in var order
      /hydro  attrs: varnames, nvar
      /        attrs: ndim, boxlen, levelmin, levelmax, gamma, time, ncpu,
                      ordering, mhd, + any units
    """
    grid = result.target
    if isinstance(grid, AMRGrid):
        nb = grid.base_npx
    elif isinstance(grid, UniformGrid):
        nb = grid.npx
    else:
        raise TypeError("to_ramses expects an AMRGrid or UniformGrid GridResult "
                        f"(got target {type(grid).__name__})")
    if fmt == "native":
        raise NotImplementedError(
            "to_ramses(fmt='native'): RAMSES amr_/hydro_ Fortran-binary writer "
            "is not implemented (no RAMSES reader installed here to validate "
            "against). Use fmt='hdf5'; the native backend would rebuild the oct "
            "nesting from the leaves' level+base and serialize per-cpu records.")
    if fmt != "hdf5":
        raise ValueError(f"unknown fmt {fmt!r} (use 'hdf5' or 'native')")

    missing = [f for f in hydro_fields if f not in result.data]
    if missing:
        raise KeyError(f"result is missing hydro field(s) {missing}; interpolate "
                       f"them (have {sorted(result.data)})")

    import h5py
    table = to_cell_table(result, fields=list(hydro_fields))
    hydro = np.ascontiguousarray(
        np.stack([table[f] for f in hydro_fields], axis=1))   # (Ncell, nvar)

    # RAMSES levels: base level l0 = log2(base_npx) (RAMSES needs a power-of-two
    # cubic base grid); a leaf at refinement r sits at RAMSES level l0 + r. For a
    # UniformGrid every cell is at the base level (table 'level' is all 0).
    cell_level = table["level"]
    l0f = np.log2(nb[0])
    cubic_pow2 = (nb[0] == nb[1] == nb[2]) and float(l0f).is_integer()
    if cubic_pow2:
        l0 = int(round(l0f))
        ramses_level = l0 + cell_level
        levelmin, levelmax = l0, int(l0 + cell_level.max())
    else:
        ramses_level = cell_level.copy()
        levelmin = levelmax = -1            # not a RAMSES-compatible base grid

    if boxlen is None:
        boxlen = float(np.max(grid.box))

    with h5py.File(path, "w") as f:
        c = f.create_group("cells")
        c.create_dataset("x", data=table["x"])
        c.create_dataset("y", data=table["y"])
        c.create_dataset("z", data=table["z"])
        c.create_dataset("dx", data=table["dx"])
        c.create_dataset("level", data=table["level"])
        c.create_dataset("base", data=table["base"])
        c.create_dataset("ramses_level", data=ramses_level)
        hg = f.create_group("hydro")
        hg.create_dataset("data", data=hydro)
        hg.attrs["varnames"] = np.array(hydro_fields, dtype=h5py.string_dtype())
        hg.attrs["nvar"] = len(hydro_fields)
        f.attrs["ndim"] = 3
        f.attrs["ncell"] = int(hydro.shape[0])
        f.attrs["boxlen"] = boxlen
        f.attrs["levelmin"] = levelmin
        f.attrs["levelmax"] = levelmax
        f.attrs["gamma"] = gamma
        f.attrs["time"] = time
        f.attrs["ncpu"] = 1
        f.attrs["ordering"] = "cell-list"
        f.attrs["mhd"] = any(v in hydro_fields for v in
                             ("Bx", "By", "Bz", "Bx_l", "Bx_r",
                              "By_l", "By_r", "Bz_l", "Bz_r"))
        if not cubic_pow2:
            f.attrs["warning"] = ("base_npx is not a power-of-two cube; "
                                  "levelmin/levelmax unset (not RAMSES-native).")
        for k, v in (units or {}).items():
            f.attrs[k] = v
    return path


# --------------------------------------------------------------------------- #
#  AREPO (moving-mesh) export.
#
#  AREPO is a finite-volume Voronoi code: an IC is just GENERATORS + per-cell
#  CONSERVED quantities; AREPO rebuilds its own Voronoi mesh from the generator
#  positions at run start and derives the primitives (rho = Mass/V_voronoi,
#  v = Momentum/Mass, ...). So ANY of our targets is a valid AREPO IC -- we hand
#  it the cell centres as generators plus the per-cell mass and (intensive)
#  primitive fields. Two paths:
#
#  * DEPOSIT (to_arepo on a GridResult): generators = cell centres of a
#    Uniform/AMR/Voronoi interpolation; Mass = rho_cell * V_cell (the extensive,
#    conserved quantity -> total mass preserved exactly regardless of how AREPO's
#    rebuilt cell volumes differ from ours); velocity / u / B are the intensive
#    cell means the deposit already produced. Use this when the mesh differs from
#    the particles (regridding / subsampling / structured grid).
#
#  * COPY (particles_to_arepo): generators = the SPH particles themselves, 1:1.
#    Each cell inherits its particle's mass (exact) and v/u/B verbatim -- no
#    kernel, no smoothing. The most faithful conversion; AREPO tessellates and
#    sets rho = m_i / V_voronoi_i. (Note: this swaps SPH's overlapping-kernel
#    volume for the space-filling Voronoi volume -- see particles_to_arepo doc.)
#
#  B is CELL-CENTRED in AREPO (Powell/Dedner cleaning), unlike RAMSES's
#  face-centred CT, so a single cell-centred B vector per generator is all that
#  is needed; the SPH field's div(B) is absorbed by AREPO's cleaning.
# --------------------------------------------------------------------------- #

#: Default per-particle field names mapped onto AREPO gas primitives.
AREPO_VELOCITY = ("vx", "vy", "vz")
AREPO_MAGNETIC = ("Bx", "By", "Bz")
#: tipsy stores the gas thermal state as a TEMPERATURE (tgdata col 8), so that is
#: the only thermodynamic source. AREPO wants specific internal energy u; the
#: gasoline conversion is u = dTuFac * T (read_tufac / `--tu-fac`).
AREPO_TEMPERATURE = ("temp", "T", "t")


def _internal_energy_from_temperature(get, temperature, tu_fac):
    """Specific internal energy u = (temperature field) * tu_fac for the AREPO IC.

    tipsy carries only a temperature column, so u comes from it via gasoline's
    dTuFac (`tu_fac`). `get(name)` returns the flattened cell field or None.
    Returns (u_array, description); raises if no temperature field is present."""
    for t in temperature:
        a = get(t)
        if a is not None:
            return a * tu_fac, f"{t!r} * tu_fac({tu_fac:g})  [u = dTuFac*T]"
    raise KeyError(
        f"to_arepo needs a temperature field (one of {temperature}) to set the "
        "specific internal energy u = dTuFac*T; add it to the interpolated "
        "fields (e.g. --fields rho,vx,vy,vz,temp).")


def _frame_coords(coords, box=None, periodic=False):
    """Place particle coordinates into AREPO's box frame [0, L] and pick L.

    box      : periodic period (3,) from the `.param` (authoritative). When
               given, each axis is centred on the data and coords are translated
               into [0, period_axis] (wrapped if `periodic`), so a gasoline box
               [-L/2,L/2] per axis maps to [0,L]. The returned `boxsize` is the
               X period (AREPO's scalar BoxSize = boxSize_X); a non-cubic box is
               flagged with the LONG_Y/LONG_Z ratios AREPO must be compiled with
               (boxSize_Y = BoxSize*LONG_Y, etc.).
    box=None : fall back to the particle extent (shift by min, L = max span).
    Returns (framed_coords, boxsize). boxsize is the X-axis period.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    if box is not None:
        box = np.asarray(box, dtype=np.float64).reshape(3)
        centre = 0.5 * (coords.min(axis=0) + coords.max(axis=0))
        lo = centre - 0.5 * box                         # per-axis box origin
        c = coords - lo
        if periodic:
            c = np.mod(c, box)                           # wrap strays into [0,L)
        if not np.allclose(box, box[0]):
            print(f"  non-cubic period {box}: AREPO needs LONG_Y="
                  f"{box[1] / box[0]:.8g} LONG_Z={box[2] / box[0]:.8g} "
                  f"(BoxSize={box[0]:g}).")
        return c, float(box[0])
    c = coords - coords.min(axis=0)
    return c, float(c.max() * (1.0 + 1e-6))


def _write_arepo_hdf5(path, coords, mass, vel, u, bfield, boxsize, ids=None,
                      time=0.0):
    """Write a GADGET/AREPO-format type-0 (gas) HDF5 IC snapshot.

    `coords` (N,3) must already be framed into [0, boxsize] (see `_frame_coords`);
    mass (N,), vel (N,3), u (N,) specific internal energy, bfield (N,3) or None.
    """
    import h5py
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    n = coords.shape[0]
    boxsize = float(boxsize)
    # keep strictly inside (0, L): AREPO rejects points on/over the box face
    eps = 1e-9 * boxsize
    coords = np.clip(coords, eps, boxsize - eps)
    if ids is None:
        ids = np.arange(1, n + 1, dtype=np.uint32)
    with h5py.File(path, "w") as f:
        h = f.create_group("Header")
        p0 = f.create_group("PartType0")
        npart = np.array([n, 0, 0, 0, 0, 0], dtype=np.int32)
        h.attrs.create("NumPart_ThisFile", npart)
        h.attrs.create("NumPart_Total", npart)
        h.attrs.create("NumPart_Total_HighWord", np.zeros(6, dtype=np.int32))
        h.attrs.create("MassTable", np.zeros(6, dtype=np.float64))
        h.attrs.create("Time", float(time))
        h.attrs.create("Redshift", 0.0)
        h.attrs.create("BoxSize", np.float64(boxsize))
        h.attrs.create("NumFilesPerSnapshot", 1)
        h.attrs.create("Omega0", 0.0)
        h.attrs.create("OmegaB", 0.0)
        h.attrs.create("OmegaLambda", 0.0)
        h.attrs.create("HubbleParam", 1.0)
        for fl in ("Flag_Sfr", "Flag_Cooling", "Flag_StellarAge", "Flag_Metals",
                   "Flag_Feedback"):
            h.attrs.create(fl, 0)
        h.attrs.create("Flag_DoublePrecision", 1)
        p0.create_dataset("ParticleIDs", data=np.asarray(ids, dtype=np.uint32))
        p0.create_dataset("Coordinates", data=coords)
        p0.create_dataset("Masses", data=np.ascontiguousarray(mass, np.float64))
        p0.create_dataset("Velocities",
                          data=np.ascontiguousarray(vel, np.float64))
        p0.create_dataset("InternalEnergy",
                          data=np.ascontiguousarray(u, np.float64))
        if bfield is not None:
            p0.create_dataset("MagneticField",
                              data=np.ascontiguousarray(bfield, np.float64))
    return path, boxsize


def to_arepo(result, path, density="rho", velocity=AREPO_VELOCITY,
             temperature=AREPO_TEMPERATURE, tu_fac=1.0, magnetic=AREPO_MAGNETIC,
             box=None, periodic=False, time=0.0):
    """Export an interpolation `GridResult` as an AREPO moving-mesh IC (HDF5).

    Works for UniformGrid / AMRGrid / VoronoiGrid results: the cell centres
    become the mesh generators; per cell we write the CONSERVED mass
    `Mass = rho * V_cell` and the intensive primitives (velocity, specific
    internal energy, cell-centred B). AREPO rebuilds the Voronoi mesh from the
    generators and recomputes density from its own cell volumes -- total mass is
    conserved exactly in the handoff.

    density       : field giving cell density rho (default 'rho'); required.
    velocity      : 3 field names for v (missing -> zero velocity).
    temperature   : candidate names for the tipsy temperature field; specific
                    internal energy is u = (temperature) * tu_fac (gasoline
                    dTuFac convention). A temperature field is REQUIRED.
    tu_fac        : dTuFac (u = dTuFac * T); read from the run .log by the CLI.
    magnetic      : 3 field names for cell-centred B (all present -> MHD IC).
    box           : periodic period (3,) for the AREPO BoxSize/frame (from the
                    `.param`); None -> use the generator extent.
    periodic      : wrap generators into [0, period) when `box` is given.
    """
    grid = result.target
    if not isinstance(grid, (UniformGrid, AMRGrid, VoronoiGrid)):
        raise TypeError("to_arepo expects a UniformGrid/AMRGrid/VoronoiGrid "
                        f"GridResult (got {type(grid).__name__})")
    # Flatten geometry + every field to per-cell (Ncell,) arrays via the
    # canonical cell table (handles the uniform C-order ravel, AMR leaves and
    # Voronoi central cells uniformly).
    have = [k for k in result.data if not k.startswith("_")]
    table = to_cell_table(result, fields=have)
    coords = np.stack([table["x"], table["y"], table["z"]], axis=1)
    vol = table["volume"]

    def get(name):
        return table[name] if (name is not None and name in table) else None

    if get(density) is None:
        raise KeyError(f"to_arepo needs the density field {density!r}; have "
                       f"{sorted(have)} (interpolate it, e.g. --fields rho,...)")
    rho = get(density)
    mass = rho * vol

    vel = np.zeros((coords.shape[0], 3))
    vmiss = [c for c in velocity if get(c) is None]
    for k, c in enumerate(velocity):
        if get(c) is not None:
            vel[:, k] = get(c)

    u, udesc = _internal_energy_from_temperature(get, temperature, tu_fac)

    bfield = None
    if all(get(c) is not None for c in magnetic):
        bfield = np.stack([get(c) for c in magnetic], axis=1)

    coords, L = _frame_coords(coords, box=box, periodic=periodic)
    path, L = _write_arepo_hdf5(path, coords, mass, vel, u, bfield, L,
                                time=time)
    print(f"  AREPO IC: {coords.shape[0]} generators ({type(grid).__name__}), "
          f"BoxSize={L:.6g}, total mass={mass.sum():.6g}")
    print(f"           u from {udesc}; velocity "
          f"{'zero' if len(vmiss) == 3 else 'from ' + ','.join(velocity)}; "
          f"{'MHD (cell-centred B)' if bfield is not None else 'hydro (no B)'}")
    return path


def particles_to_arepo(particles, path, velocity=AREPO_VELOCITY,
                       temperature=AREPO_TEMPERATURE, tu_fac=1.0,
                       magnetic=AREPO_MAGNETIC, box=None, periodic=None,
                       time=0.0):
    """COPY path: write the SPH particles directly as an AREPO IC (1 generator
    per particle, no interpolation). Each cell inherits its particle's mass
    (exact) and the intensive fields (v/u/B) verbatim from `particles.values`;
    AREPO tessellates and sets rho = m_i / V_voronoi_i.

    This is the most faithful transcription of the per-particle state, BUT note
    it swaps SPH's overlapping-kernel volume (m/rho_SPH) for the space-filling
    Voronoi cell volume, so the resulting cell DENSITY is AREPO's Voronoi
    estimate, not rho_SPH. For a kernel-consistent, mass-conserving density use
    the deposit path (`to_arepo` on a VoronoiGrid petkova result) instead.
    """
    vals = particles.values

    def get(name):
        return vals[name] if (name is not None and name in vals) else None

    n = particles.pos.shape[0]
    mass = np.ascontiguousarray(particles.mass, np.float64)
    vel = np.zeros((n, 3))
    vmiss = [c for c in velocity if get(c) is None]
    for k, c in enumerate(velocity):
        if get(c) is not None:
            vel[:, k] = get(c)
    u, udesc = _internal_energy_from_temperature(get, temperature, tu_fac)
    bfield = None
    if all(get(c) is not None for c in magnetic):
        bfield = np.stack([get(c) for c in magnetic], axis=1)

    # default the box / periodicity from the particle container (which from_tipsy
    # fills from the .param period).
    if box is None:
        box = particles.box
    if periodic is None:
        periodic = bool(np.any(np.atleast_1d(particles.periodic)))
    coords, L = _frame_coords(particles.pos, box=box, periodic=periodic)
    path, L = _write_arepo_hdf5(path, coords, mass, vel, u, bfield, L,
                                time=time)
    print(f"  AREPO IC (copy): {n} particles -> generators, BoxSize={L:.6g}, "
          f"total mass={mass.sum():.6g}")
    print(f"           u from {udesc}; velocity "
          f"{'zero' if len(vmiss) == 3 else 'from ' + ','.join(velocity)}; "
          f"{'MHD (cell-centred B)' if bfield is not None else 'hydro (no B)'}")
    return path


def to_enzo(result, path, **kw):
    """Write Enzo/FLASH-style patch hierarchy (HDF5).

    Contract: these are patch-based, not pure octrees; tile each refinement
    `level` into rectangular grid patches (runs of same-level leaves sharing a
    base cell) and write the grid hierarchy + fields to HDF5.
    """
    raise NotImplementedError(
        "to_enzo: coalesce same-level leaves into rectangular patches and write "
        "the HDF5 grid hierarchy.")
