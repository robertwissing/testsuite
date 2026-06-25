#!/usr/bin/env python
"""Convert a gasoline/tipsy IC into an AREPO moving-mesh IC (HDF5).

The AREPO analog of `ICconvert_tipsytosphexa.py`. It is a thin front end over the
`sph_interp` interpolation software (`IC_interpolate_to_grid.py --output-format
arepo`): read the tipsy snapshot (-> Particles via from_tipsy, periodic box +
dTuFac from the .param/.log), turn it into AREPO generators + per-cell conserved
quantities, and write the GADGET/AREPO HDF5 IC.

    ICconvert_tipsytoarepo.py <input> <output> [mode] [method]

  input    tipsy snapshot (e.g. datafiles/run.00000).
  output   AREPO IC HDF5 path (e.g. datafiles/run_arepo.hdf5).
  mode     particles (default; 1:1 copy of the SPH particles as generators, the
           most faithful conversion) | voronoi | amr | uniform (deposit the SPH
           fields onto that mesh and use its cell centres as generators).
  method   sph | petkova (default; exact, mass-conserving) -- deposit modes only.

Options: --fields (default rho,vx,vy,vz,Bx,By,Bz), --tu-fac (override dTuFac),
--npx, and the other IC_interpolate_to_grid knobs are forwarded.

If --param-out is given, a fresh AREPO simulation `param.txt` is also written,
with run parameters (timestepping, box, output cadence) taken from the gasoline
`.param` beside the input (see sph_interp/_arepo/arepo_param.py).  This replaces
the old approach of sed-editing AREPO's noh_3d example file.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from IC_interpolate_to_grid import main as interp_main
from sph_interp._arepo.arepo_param import write_run_param


def _find_gasoline_param(tipsy_input):
    """The gasoline `.param` beside a tipsy snapshot `<name>.NNNNN`: try
    `<name>.param` (strip the snapshot index)."""
    base = tipsy_input
    if "." in os.path.basename(base):
        base = base.rsplit(".", 1)[0]
    cand = base + ".param"
    return cand if os.path.exists(cand) else None


def build_argv(args):
    argv = [args.input, args.output, "--output-format", "arepo",
            "--fields", args.fields]
    if args.mode == "particles":
        argv += ["--voronoi-particles"]
    else:
        argv += ["--grid", args.mode, "--method", args.method]
        if args.npx is not None:
            argv += ["--npx", str(args.npx)]
    if args.tu_fac is not None:
        argv += ["--tu-fac", str(args.tu_fac)]
    return argv


def parse(argv=None):
    pa = argparse.ArgumentParser(
        prog="ICconvert_tipsytoarepo",
        description="tipsy -> AREPO moving-mesh IC (via sph_interp).")
    pa.add_argument("input", help="tipsy snapshot")
    pa.add_argument("output", help="AREPO IC HDF5 path")
    pa.add_argument("mode", nargs="?", default="particles",
                    choices=["particles", "voronoi", "amr", "uniform"],
                    help="particles=1:1 copy (default), or a deposit mesh")
    pa.add_argument("method", nargs="?", default="petkova",
                    choices=["sph", "petkova"], help="deposit method")
    pa.add_argument("--fields", default="rho,vx,vy,vz,Bx,By,Bz")
    pa.add_argument("--tu-fac", type=float, default=None)
    pa.add_argument("--npx", type=int, default=None)
    pa.add_argument("--param-out", default=None,
                    help="also write a fresh AREPO sim param.txt here, with run "
                         "parameters taken from the gasoline .param beside INPUT")
    pa.add_argument("--output-dir", default=None,
                    help="AREPO OutputDir to write into the generated param.txt "
                         "(default: '<param-out dir>/')")
    pa.add_argument("--boxsize", type=float, default=None,
                    help="override BoxSize in the param.txt (default: dxPeriod). "
                         "For non-cubic boxes pass the X period; y/z come from "
                         "compile-time LONG_Y/LONG_Z.")
    return pa.parse_args(argv)


def write_param(args):
    """Write the AREPO simulation param.txt for this conversion (if requested)."""
    if not args.param_out:
        return
    gas_param = _find_gasoline_param(args.input)
    if gas_param is None:
        print("WARNING: no gasoline .param found beside %r; cannot write "
              "param.txt." % args.input)
        return
    # InitCondFile = the AREPO IC path without its extension (ICFormat 3 = HDF5).
    init_cond = os.path.splitext(args.output)[0]
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = (os.path.dirname(os.path.abspath(args.param_out)) or ".") \
            + os.sep
    write_run_param(gas_param, args.param_out, init_cond, output_dir,
                    boxsize=args.boxsize)
    print("Wrote AREPO param.txt -> %s (run params from %s)"
          % (args.param_out, gas_param))


if __name__ == "__main__":
    args = parse()
    interp_main(build_argv(args))
    write_param(args)
