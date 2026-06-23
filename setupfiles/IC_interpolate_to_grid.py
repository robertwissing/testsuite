#!/usr/bin/env python
"""Interpolate a tipsy SPH snapshot onto a grid and write it out.

A thin command-line front end for the `sph_interp` package: read a gasoline/tipsy
snapshot (-> Particles via from_tipsy), deposit the requested fields onto a chosen
TARGET geometry with a chosen METHOD, and serialize the result.

    interpolate_to_grid.sh INPUT OUTPUT [options]

Run with NO positional arguments to be prompted interactively for everything.

TARGETS (--grid):   uniform | amr | voronoi | projection
METHODS (--method): petkova (exact, mass-conserving; default) | sph (fast, approximate)
OUTPUT (--output-format): npz | hdf5 | ramses | cell-table | native

Examples
  # 256^3 uniform SPH density grid -> npz
  interpolate_to_grid.sh run/snap.00100 rho.npz --grid uniform --npx 256

  # mass-conserving AMR (octree) of rho+|B|, exported to RAMSES HDF5
  interpolate_to_grid.sh run/snap.00100 mesh.h5 --grid amr --method petkova \
      --fields rho,vx,vy,vz --output-format ramses

  # 50000-cell Voronoi mesh subsampled from the particles, exact deposit
  interpolate_to_grid.sh run/snap.00100 voro.npz --grid voronoi \
      --voronoi-source subsample --ncells 50000 --method petkova

  # column-density projection map in the xy plane
  interpolate_to_grid.sh run/snap.00100 col.npz --grid projection \
      --axes xy --project column
"""

import argparse
import os
import sys

import numpy as np

# Make sure the package next to this script is importable when invoked directly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# --------------------------------------------------------------------------- #
#  Interactive prompts (used when no positional args are given)
# --------------------------------------------------------------------------- #

def _ask(prompt, default=None, choices=None):
    suffix = f" [{default}]" if default is not None else ""
    if choices:
        suffix = f" ({'/'.join(choices)}){suffix}"
    while True:
        ans = input(f"{prompt}{suffix}: ").strip()
        if not ans and default is not None:
            return default
        if not ans:
            print("  (required)")
            continue
        if choices and ans not in choices:
            print(f"  choose one of: {', '.join(choices)}")
            continue
        return ans


def interactive(args):
    """Fill `args` by prompting the user. Returns the populated namespace."""
    print("=" * 64)
    print(" sph_interp : interpolate a tipsy snapshot onto a grid")
    print("=" * 64)
    args.input = _ask("Input tipsy snapshot")
    args.grid = _ask("Grid / target geometry", default="uniform",
                     choices=["uniform", "amr", "voronoi", "projection"])
    args.method = _ask("Interpolation method", default="petkova",
                       choices=["sph", "petkova"])
    args.fields = _ask("Fields to interpolate (comma list)", default="rho")

    if args.grid == "uniform":
        args.npx = int(_ask("Cells per axis (npx)", default="128"))
    elif args.grid == "amr":
        args.base_npx = int(_ask("Base grid cells per axis", default="16"))
        args.npx = int(_ask("Finest resolution / cells per axis (npx)",
                            default="128"))
        args.refine = _ask("Refinement criterion", default="smoothing",
                           choices=["smoothing", "mass", "count"])
    elif args.grid == "voronoi":
        args.voronoi_source = _ask(
            "Voronoi generators from", default="particles",
            choices=["particles", "random", "subsample", "grid"])
        if args.voronoi_source != "particles":
            args.npx = int(_ask("Cells per axis (npx; total npx^3 cells)",
                                default="64"))
    elif args.grid == "projection":
        args.axes = _ask("Image plane", default="xy",
                         choices=["xy", "xz", "yz"])
        args.project = _ask("Projection mode", default="column",
                            choices=["column", "average", "rhocolumn", "slice"])
        args.npx = int(_ask("Pixels per axis (npx)", default="512"))

    default_fmt = {"uniform": "npz", "amr": "ramses", "voronoi": "npz",
                   "projection": "npz"}[args.grid]
    args.output_format = _ask("Output format", default=default_fmt,
                              choices=["npz", "hdf5", "ramses", "cell-table",
                                       "native"])
    args.output = _ask("Output file")
    if _ask("Write a SPH-vs-grid diagnostic figure?", default="y",
            choices=["y", "n"]) == "y":
        args.visualise = _ask("  figure path", default=os.path.splitext(
            args.output)[0] + "_visualise.pdf")
        args.visualise_field = _ask("  field to visualise",
                                    default=args.fields.split(",")[0].strip())
        if args.grid != "projection":
            args.visualise_axes = _ask("  slice plane", default="xy",
                                       choices=["xy", "xz", "yz"])
    print("-" * 64)
    return args


# --------------------------------------------------------------------------- #
#  Building the target geometry
# --------------------------------------------------------------------------- #

def amr_max_level(args):
    """AMR refinement depth: use --max-level if given, else derive it from --npx
    (the finest target resolution) and --base-npx, so base_npx*2^level >= npx."""
    if args.max_level is not None:
        return int(args.max_level)
    import math
    return max(0, math.ceil(math.log2(max(args.npx, 1) / args.base_npx)))


def build_target(args, p):
    from sph_interp import (UniformGrid, AMRGrid, VoronoiGrid, Projection2D,
                            refine_smoothing_length, refine_mass, refine_count)
    box = p.box if p.box is not None else (p.pos.max(0) - p.pos.min(0))

    if args.grid == "uniform":
        print(f"  target: UniformGrid {args.npx}^3 over box {np.round(box, 4)}")
        return UniformGrid.centered(box, args.npx)

    if args.grid == "amr":
        crit = {"smoothing": refine_smoothing_length(),
                "mass": refine_mass(args.refine_value
                                    if args.refine_value is not None
                                    else float(p.mass.sum()) / 1e4),
                "count": refine_count(int(args.refine_value)
                                      if args.refine_value is not None
                                      else 8)}[args.refine]
        max_level = amr_max_level(args)
        fine = args.base_npx * (2 ** max_level)
        print(f"  target: AMRGrid base {args.base_npx}^3, max_level "
              f"{max_level}, refine={args.refine}")
        print(f"          adaptive: finest level resolves like {fine}^3 "
              f"(= base_npx x 2^max_level, derived from --npx {args.npx}); "
              "coarse in voids, fine only where the data warrants it")
        g = AMRGrid.build(p, base_npx=args.base_npx, max_level=max_level,
                          criterion=crit)
        dense = fine ** 3
        print(f"          -> {g.ncell} leaf cells "
              f"({dense/max(g.ncell,1):.0f}x fewer than a dense {fine}^3 grid)")
        return g

    if args.grid == "voronoi":
        half = 0.5 * np.asarray(box).reshape(3)
        # Periodic box -> bounds ARE the period (ghost images give exact periodic
        # cells, Sum volume = box volume); non-periodic -> pad past the box.
        vper = p.periodic if any(np.atleast_1d(p.periodic)) else False
        if vper:
            bounds = [[-half[a], half[a]] for a in range(3)]
            ptag = " (periodic, ghost images)"
        else:
            pad = args.voronoi_pad
            bounds = [[-half[a] - pad, half[a] + pad] for a in range(3)]
            ptag = ""
        if args.voronoi_source == "particles":
            print(f"  target: VoronoiGrid from {p.n} particles{ptag}")
            return VoronoiGrid.from_particles(p, bounds=bounds, periodic=vper)
        if args.voronoi_source == "grid":
            nx = args.npx
            axes = [(-half[a] + (np.arange(nx) + 0.5) * (2.0 * half[a] / nx))
                    for a in range(3)]
            X, Y, Z = np.meshgrid(*axes, indexing="ij")
            pts = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
            print(f"  target: VoronoiGrid from a {nx}^3 = {nx**3} regular "
                  f"lattice (cubic cells){ptag}")
            return VoronoiGrid.from_points(pts, bounds=bounds, periodic=vper)
        ncells = args.npx ** 3                     # npx per-axis-equivalent count
        rng = np.random.default_rng(args.seed)
        if args.voronoi_source == "random":
            pts = rng.uniform(-half + 1e-9, half - 1e-9, size=(ncells, 3))
            print(f"  target: VoronoiGrid from {ncells} (={args.npx}^3) random "
                  f"generators{ptag}")
        else:  # subsample
            n = min(ncells, p.n)
            idx = rng.choice(p.n, size=n, replace=False)
            pts = p.pos[idx]
            print(f"  target: VoronoiGrid from {n} subsampled particles "
                  f"(target {args.npx}^3={ncells}){ptag}")
        return VoronoiGrid.from_points(pts, bounds=bounds, periodic=vper)

    if args.grid == "projection":
        a1, a2 = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}[args.axes]
        print(f"  target: Projection2D {args.axes} {args.npx}^2, "
              f"mode={args.project}")
        return Projection2D.centered(box, a1, a2, npx=args.npx,
                                     mode=args.project)

    raise ValueError(f"unknown grid {args.grid!r}")


# --------------------------------------------------------------------------- #
#  Writing the result
# --------------------------------------------------------------------------- #

def _collect(res):
    """A serialization-ready dict for any target geometry."""
    from sph_interp import UniformGrid, AMRGrid, VoronoiGrid, Projection2D
    grid = res.target
    out = {"_grid": type(grid).__name__, "_method": res.method}
    fields = {k: v for k, v in res.data.items()}

    if isinstance(grid, UniformGrid):
        cx, cy, cz = grid.centres()
        out.update(x=cx, y=cy, z=cz, lo=grid.lo, hi=grid.hi,
                   npx=np.array(grid.npx))
        out.update(fields)                       # 3D arrays as-is
    elif isinstance(grid, Projection2D):
        ex, ey = grid.centres()
        out.update(x=ex, y=ey, lo=grid.lo, hi=grid.hi,
                   axis1=grid.axis1, axis2=grid.axis2, mode=grid.mode)
        out.update(fields)                       # 2D arrays as-is
    elif isinstance(grid, AMRGrid):
        c = grid.centres()
        out.update(x=c[:, 0], y=c[:, 1], z=c[:, 2], dx=grid.cell_size[:, 0],
                   dy=grid.cell_size[:, 1], dz=grid.cell_size[:, 2],
                   volume=grid.volumes(), level=grid.cell_level,
                   base=grid.cell_base)
        out.update(fields)                       # (Ncell,) arrays
    elif isinstance(grid, VoronoiGrid):
        c = grid.centres()
        out.update(x=c[:, 0], y=c[:, 1], z=c[:, 2], volume=grid.volumes())
        out.update(fields)                       # (Ncell,) arrays
    return out


def write_output(res, path, fmt, args, p):
    grid_name = type(res.target).__name__

    if fmt in ("ramses", "native"):
        from sph_interp import to_ramses, RAMSES_HYDRO_VARS
        rfmt = "native" if fmt == "native" else "hdf5"
        hydro = (args.hydro_fields.split(",") if args.hydro_fields
                 else [f for f in RAMSES_HYDRO_VARS if f in res.data])
        if not hydro:
            raise SystemExit(
                "ramses: none of the RAMSES hydro variables "
                f"{RAMSES_HYDRO_VARS} are in the interpolated fields "
                f"{sorted(res.data)}; pass --hydro-fields or add them to "
                "--fields (e.g. --fields rho,vx,vy,vz).")
        print(f"  writing RAMSES ({rfmt}) hydro vars {hydro} -> {path}")
        return to_ramses(res, path, hydro_fields=tuple(hydro),
                         gamma=args.gamma, boxlen=args.boxlen, fmt=rfmt)

    if fmt == "cell-table":
        from sph_interp import to_cell_table
        table = to_cell_table(res)
        print(f"  writing cell table ({len(next(iter(table.values())))} rows) "
              f"-> {path}")
        np.savez(path, **table)
        return path

    payload = _collect(res)
    if fmt == "hdf5":
        import h5py
        print(f"  writing HDF5 -> {path}")
        with h5py.File(path, "w") as f:
            for k, v in payload.items():
                if isinstance(v, str):
                    f.attrs[k] = v
                else:
                    f.create_dataset(k, data=np.asarray(v))
        return path

    # default: npz (works for every target)
    print(f"  writing npz -> {path}")
    np.savez(path, **{k: v for k, v in payload.items()
                      if not isinstance(v, str)},
             **{f"attr_{k}": v for k, v in payload.items()
                if isinstance(v, str)})
    return path


# --------------------------------------------------------------------------- #
#  Mass-conservation report (a quick sanity line)
# --------------------------------------------------------------------------- #

def report_mass(res, p):
    from sph_interp import UniformGrid, AMRGrid, VoronoiGrid
    grid = res.target
    pm = float(p.mass.sum())
    if "_mass" in res.data:                       # petkova exposes conserved mass
        gm = float(np.asarray(res.data["_mass"]).sum())
    elif isinstance(grid, UniformGrid) and "rho" in res.data:
        gm = float(res.data["rho"].sum() * grid.cell_volume)
    elif isinstance(grid, (AMRGrid, VoronoiGrid)) and "rho" in res.data:
        gm = float((res.data["rho"] * grid.volumes()).sum())
    else:
        return
    tag = ("" if "_mass" in res.data
           else "  [approx: SPH density estimate, not mass-conserving; "
                "edge/padding cells inflate sum(rho*vol)]")
    print(f"  mass check: grid {gm:.6g} vs particles {pm:.6g}  "
          f"({100*(gm/pm - 1):+.3f}%){tag}")


# --------------------------------------------------------------------------- #
#  Main
# --------------------------------------------------------------------------- #

def build_parser():
    pa = argparse.ArgumentParser(
        prog="interpolate_to_grid",
        description="Interpolate a tipsy SPH snapshot onto a grid (sph_interp). "
                    "Run with no arguments for an interactive prompt.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    pa.add_argument("input", nargs="?", help="input tipsy snapshot")
    pa.add_argument("output", nargs="?", help="output file/prefix")

    pa.add_argument("--grid", choices=["uniform", "amr", "voronoi", "projection"],
                    default="uniform", help="target geometry")
    pa.add_argument("--method", choices=["sph", "petkova"], default="petkova",
                    help="petkova=exact/mass-conserving (default), sph=fast/approx")
    pa.add_argument("--fields", default="rho",
                    help="comma list of fields to interpolate (rho,Bmag,vx,...)")
    pa.add_argument("--output-format",
                    choices=["npz", "hdf5", "ramses", "cell-table", "native"],
                    default=None, help="default: npz (ramses for --grid amr)")

    # snapshot / domain
    pa.add_argument("--nsmooth", type=int, default=64,
                    help="neighbours for h when no smoothlength aux")
    pa.add_argument("--box", type=float, default=None,
                    help="override periodic box size (cube)")
    g = pa.add_mutually_exclusive_group()
    g.add_argument("--periodic", dest="periodic", action="store_true")
    g.add_argument("--no-periodic", dest="periodic", action="store_false")
    pa.set_defaults(periodic=True)

    # resolution (single knob): cells per axis
    pa.add_argument("--npx", type=int, default=128,
                    help="cells per axis (THE resolution knob). uniform/"
                         "projection: grid/image resolution. amr: the FINEST "
                         "target resolution -> max_level is derived so "
                         "base_npx x 2^max_level >= npx (AMR is adaptive). "
                         "voronoi: lattice size for --voronoi-source grid.")
    # projection
    pa.add_argument("--axes", choices=["xy", "xz", "yz"], default="xy",
                    help="projection image plane")
    pa.add_argument("--project",
                    choices=["column", "average", "rhocolumn", "slice"],
                    default="column", help="projection mode")
    # amr
    pa.add_argument("--base-npx", type=int, default=16, help="AMR base grid/axis")
    pa.add_argument("--max-level", type=int, default=None,
                    help="AMR max level (default: derived from --npx and "
                         "--base-npx). Set to override the npx-derived depth.")
    pa.add_argument("--refine", choices=["smoothing", "mass", "count"],
                    default="smoothing", help="AMR refinement criterion")
    pa.add_argument("--refine-value", type=float, default=None,
                    help="threshold for --refine mass/count")
    # voronoi
    pa.add_argument("--voronoi-source",
                    choices=["particles", "random", "subsample", "grid"],
                    default="particles",
                    help="Voronoi generators: particles (one cell per particle, "
                         "resolution fixed by the data), or subsample / random / "
                         "grid which use --npx (npx^3 cells: subsampled particles, "
                         "random points, or a regular lattice)")
    pa.add_argument("--voronoi-pad", type=float, default=0.05,
                    help="pad cell domain past the box (keeps boundary kernels in)")
    pa.add_argument("--seed", type=int, default=0, help="RNG seed (voronoi)")

    # visualisation
    pa.add_argument("--visualise", "--visualize", dest="visualise", nargs="?",
                    const="", default=None,
                    help="write a diagnostic figure comparing the SPH particle "
                         "field to the interpolated grid (geometry-aware slice, "
                         "density PDFs, per-particle grid-vs-SPH scatter, mass "
                         "conservation). Optionally pass a path; default is "
                         "<output>_visualise.pdf.")
    pa.add_argument("--visualise-field", default=None,
                    help="which field to visualise (default: first of --fields)")
    pa.add_argument("--visualise-axes", choices=["xy", "xz", "yz"], default="xy",
                    help="slice/projection plane for --visualise (3D targets)")

    # ramses
    pa.add_argument("--hydro-fields", default=None,
                    help="ordered RAMSES hydro vars (default rho,vx,vy,vz,P "
                         "intersected with --fields)")
    pa.add_argument("--gamma", type=float, default=5.0 / 3.0)
    pa.add_argument("--boxlen", type=float, default=None)
    return pa


def main(argv=None):
    pa = build_parser()
    args = pa.parse_args(argv)

    # No positional args at all -> interactive mode.
    if args.input is None and args.output is None:
        args = interactive(args)
    elif args.input is None or args.output is None:
        pa.error("provide BOTH input and output, or NEITHER (interactive)")

    if args.output_format is None:
        args.output_format = "ramses" if args.grid == "amr" else "npz"

    from sph_interp import from_tipsy, interpolate

    fields = tuple(f.strip() for f in args.fields.split(",") if f.strip())
    print(f"reading {args.input}  fields={fields}")
    box = (np.repeat(args.box, 3) if args.box is not None else None)
    p = from_tipsy(args.input, fields=fields, nsmooth=args.nsmooth,
                   periodic=args.periodic, box=box)
    print(f"  {p.n} gas particles, box {np.round(p.box, 4)}, "
          f"periodic={p.periodic}")

    target = build_target(args, p)
    print(f"  depositing with method={args.method} ...")
    res = interpolate(p, target, method=args.method, values=list(fields))
    report_mass(res, p)

    write_output(res, args.output, args.output_format, args, p)

    if args.visualise is not None:
        from sph_interp import visualise
        vfield = args.visualise_field or fields[0]
        vpath = args.visualise or os.path.splitext(args.output)[0] + "_visualise.pdf"
        plane = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}[args.visualise_axes]
        print(f"  visualising field '{vfield}' ({args.visualise_axes}) ...")
        visualise(res, p, field=vfield, out=vpath, plane=plane)

    print("done.")


if __name__ == "__main__":
    main()
