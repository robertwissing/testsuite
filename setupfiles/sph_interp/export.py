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

from .targets import AMRGrid, UniformGrid


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
    else:
        raise TypeError("to_cell_table expects an AMRGrid or UniformGrid "
                        f"GridResult (got target {type(grid).__name__})")
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


def to_enzo(result, path, **kw):
    """Write Enzo/FLASH-style patch hierarchy (HDF5).

    Contract: these are patch-based, not pure octrees; tile each refinement
    `level` into rectangular grid patches (runs of same-level leaves sharing a
    base cell) and write the grid hierarchy + fields to HDF5.
    """
    raise NotImplementedError(
        "to_enzo: coalesce same-level leaves into rectangular patches and write "
        "the HDF5 grid hierarchy.")
