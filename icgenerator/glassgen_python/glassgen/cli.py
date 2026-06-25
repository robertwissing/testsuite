"""CLI: functional replacement for the gasoline.gdicgen invocation in
runtest.sh (tipsy pre-IC + density table in -> relaxed tipsy glass out).

Example, mirroring runtest.sh ICCONFIG:
    PYTHONPATH=icgenerator/glassgen_python python -m glassgen.cli \
        -o datafiles/myglass -I datafiles/pre.00000 -s 64 -n 4800 -oi 4800 \
        --param datafiles/pre.param --table densitytable_xdr --direction r_sph

Without --table the target density is uniform (total gas mass / box
volume). Output: ${out}_IC (tipsy, gas positions/density/hsmooth updated,
velocities zeroed) and ${out}.NNNNN snapshots every -oi iterations.
"""
import argparse
import os
import sys

import numpy as np

_SETUPFILES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           '..', '..', '..', 'setupfiles')
sys.path.insert(0, _SETUPFILES)

from .density import MAPPINGS, CallableDensity, TableDensity
from .params import box_from_param, density_from_param, read_param
from .relax import relax


def write_denserr_plot(out, pos, mass, h, rho, density, box, nsmooth,
                       flavour='gdforce', rhopow=2.0, nsmoothmin=32):
    """Write ${out}_denserr.png: the density-error |1-rho/rho0| and residual
    descent-force |icvel| distributions of the final glass (the force whose RMS
    the controller drives to equilibrium / stop), each annotated mean/p95/max.
    Self-contained (matplotlib only); the force reuses diagnostics.residual_force."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from .diagnostics import residual_force

    pos = np.asarray(pos, float)
    if box is not None:
        box = np.asarray(box, float)
        pos = pos - box * np.rint(pos / box)
    rho = np.asarray(rho, float)
    rho0 = density.rho0(pos)
    err = np.abs(1.0 - rho / rho0)
    force = residual_force(pos, mass, h, rho, rho0, box, rhopow=rhopow,
                           nsmooth=nsmooth, flavour=flavour,
                           nsmoothmin=nsmoothmin)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, data, lab in (
            (axes[0], err, r'$|1-\rho/\rho_0|$'),
            (axes[1], force, r'residual force $|\mathrm{icvel}|$')):
        ax.hist(data, bins=80, histtype='step', lw=1.8, color='C0')
        ax.set_yscale('log')
        ax.set_xlabel(lab); ax.set_ylabel('count')
        ax.set_title(f'mean {data.mean():.4g}  p95 {np.percentile(data,95):.4g}  '
                     f'max {data.max():.4g}')
    fig.suptitle(f'{os.path.basename(out)}: glass density error + residual force')
    fig.tight_layout()
    path = f'{out}_denserr.png'
    fig.savefig(path, dpi=110, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {path}  '
          f'(denserr mean/p95/max {err.mean():.4g}/{np.percentile(err,95):.4g}/'
          f'{err.max():.4g})', flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='glassgen: glass-relaxation IC generator')
    ap.add_argument('-o', dest='out', required=True,
                    help='output prefix (writes ${out}_IC)')
    ap.add_argument('-I', dest='infile', required=True,
                    help='input tipsy pre-IC')
    ap.add_argument('-s', dest='nsmooth', type=int, default=64)
    ap.add_argument('-n', dest='max_iter', type=int, default=4800)
    ap.add_argument('-oi', dest='oi', type=int, default=0,
                    help='snapshot output interval (0 = none)')
    ap.add_argument('--table', default=None,
                    help='densitytable_xdr file (default: uniform density)')
    ap.add_argument('--direction', default='r_sph', choices=MAPPINGS,
                    help='coordinate a 1D table is indexed by')
    ap.add_argument('--flavour', default='gdforce')
    ap.add_argument('--rhopow', type=float, default=2.0, help='ICrhopow')
    ap.add_argument('--icr0', type=float, default=1.0, help='ICR0')
    ap.add_argument('--icrloop0', type=float, default=0.5, help='ICRLOOP0')
    ap.add_argument('--icr0rate', type=float, default=1.5, help='ICR0Rate')
    ap.add_argument('--margin', type=int, default=8,
                    help='extra kNN candidates kept per particle; larger '
                         'values widen the Verlet skin (fewer rebuilds) at '
                         'slightly higher per-iteration cost')
    ap.add_argument('--one-sided', action='store_true',
                    help='SPH-EXA-style transpose-free force walk: faster, '
                         'pair forces no longer exactly antisymmetric')
    ap.add_argument('--h-mode', dest='hmode', default='knn',
                    choices=('knn', 'measured', 'target'),
                    help='smoothing length for the force kernel + move step: '
                         'knn (default), measured (=(m*hfact3/rho_sph)^1/3, '
                         'the C density-based h), or target '
                         '(=(m*hfact3/rho0)^1/3, noise-free analytic h)')
    ap.add_argument('--nsmoothmin', type=int, default=32,
                    help='floor on density-based h: never gather fewer than '
                         '~nsmoothmin real neighbours (measured/target only)')
    ap.add_argument('--no-progressive', dest='progressive',
                    action='store_false',
                    help='disable the DEFAULT multi-resolution path and run a '
                         'single full-resolution relax. The default coarsens '
                         'the input to the auto density-chosen coarse count, '
                         'relaxes it, then cascades (x2 splits + confined heals) '
                         'back up to full res.')
    ap.add_argument('--coarse-target', dest='coarse_target', default='auto',
                    help='progressive coarse particle count: "auto" picks '
                         'N_coarse_min from the target density (resolves the '
                         'steepest feature); an integer forces a coarse count. '
                         'The merge factor is floor_pow2(N_full/target); '
                         'factor 1 => no merge / direct relax (default auto)')
    ap.add_argument('--alpha', type=float, default=4.0,
                    help='resolution tolerance for the auto coarse count '
                         '(default 4.0)')
    ap.add_argument('--nx-full', dest='nx_full', type=int, default=0,
                    help='full target resolution (default: round(N**(1/3)) '
                         'from the input)')
    ap.add_argument('--inter-iter', dest='inter_iter', type=int, default=10,
                    help='icgen iterations after each intermediate x2 split in '
                         'the progressive cascade (default 10)')
    ap.add_argument('--momentum', type=float, default=0.5,
                    help='heavy-ball coefficient on the post-split confined '
                         'stages (default 0.5; 0 = off)')
    ap.add_argument('--icr0-stop', dest='icr0_stop', type=float, default=3e-3,
                    help='median-smoothed icr0 stop (max step as a fraction of '
                         'h) - the quality dial: 3e-3 fast/good, ~1e-3 tight. '
                         '0 = off (then -n / --polish-iter ends the run)')
    ap.add_argument('--polish-iter', dest='polish_iter', type=int, default=0,
                    help='fixed polish-iteration budget on the fine stage '
                         '(0 = off; composes with --icr0-stop)')
    ap.add_argument('--no-self-bias', dest='self_bias', action='store_false',
                    help='disable the Wendland C2 self-density bias correction '
                         '(Wzero10); density then uses the bare W0 self term')
    ap.add_argument('--plot', action='store_true',
                    help='after relaxation, write ${out}_denserr.png: the '
                         'density-error |1-rho/rho0| and residual-force |icvel| '
                         'distributions of the final glass (mean/p95/max)')
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument('--param', help='gasoline .param file to take the box from')
    grp.add_argument('--box', nargs=3, type=float, metavar=('LX', 'LY', 'LZ'))
    grp.add_argument('--open', action='store_true', help='no periodic box')
    args = ap.parse_args(argv)

    import readtipsy as tip

    gas, dark, star, header, time = tip.readtipsy(args.infile)
    if len(gas) == 0:
        sys.exit('no gas particles in input')
    print(f'read {args.infile}: {len(gas)} gas particles', flush=True)
    pos = gas[:, 1:4].astype(np.float64)
    mass = gas[:, 0].astype(np.float64)

    param_dict = read_param(args.param) if args.param else None
    if args.open:
        box = None
    elif args.box:
        box = np.array(args.box)
    else:
        box = box_from_param(param_dict)

    if args.table:
        density = TableDensity.from_xdr(args.table, mapping=args.direction)
    elif param_dict is not None:
        # read the density profile straight from the .param (profiles 1/2/3),
        # exactly as the C gdicgen does - no density table needed for these.
        density = density_from_param(param_dict)
    else:
        # --box without a .param: uniform target (total gas mass / box volume)
        if box is None:
            sys.exit('uniform density needs a box (--box)')
        rho_mean = mass.sum() / np.prod(box)
        density = CallableDensity(lambda p, r=rho_mean: np.full(len(p), r))

    def write_snapshot(path, pos_now, rho_now, h_now):
        g = gas.copy()
        g[:, 1:4] = pos_now
        g[:, 4:7] = 0.0
        g[:, 7] = rho_now
        g[:, 9] = h_now
        tip.writetipsy(g.astype(gas.dtype), dark, star, path, header, time)

    if args.progressive:
        if box is None:
            sys.exit('--progressive needs a periodic box (--param or --box)')
        from .progressive import relax_progressive
        nx_full = args.nx_full or int(round(len(gas) ** (1.0 / 3.0)))
        # merge the input down to the density-chosen coarse count (auto
        # N_coarse_min, alpha), relax, then cascade (x2 splits + confined heals,
        # momentum on the post-split stages) back to the input count 1:1.
        # --coarse-target overrides 'auto' with an explicit count.
        ct = args.coarse_target
        ct = ct if ct == 'auto' else int(ct)
        res = relax_progressive(
            pos, mass, density, box, nx_full=nx_full, coarse_seed='merge',
            coarse_target=ct, alpha=args.alpha, inter_iter=args.inter_iter,
            nsmooth=args.nsmooth, max_iter=args.max_iter,
            momentum=args.momentum, icr0_stop=args.icr0_stop,
            polish_iter=args.polish_iter, cap=2.0, start_icr0=2.0,
            flavour=args.flavour, rhopow=args.rhopow, icr0rate=args.icr0rate,
            margin=args.margin, one_sided=args.one_sided, hmode=args.hmode,
            nsmoothmin=args.nsmoothmin, self_bias=args.self_bias, verbose=True)
        nf = len(res.pos)
        mf = mass.sum() / nf  # equal-mass glass
        g = np.zeros((nf, 12), dtype=gas.dtype)
        g[:, 0] = mf
        g[:, 1:4] = res.pos
        g[:, 7] = res.rho
        g[:, 9] = res.h
        hdr = header.copy()
        hdr[0] = nf
        hdr[2] = nf
        tip.writetipsy(g, dark, star, f'{args.out}_IC', hdr, time)
        err = np.abs(1.0 - res.rho / density.rho0(res.pos))
        print(f'done (progressive -> nx_full={nx_full}): final N={nf}, '
              f'mean|1-rho/rho0|={err.mean():.3e}, wrote {args.out}_IC')
        for l in res.levels:
            print(f'  stage {l["name"]:<6} N={l["n"]:<8} niter={l["niter"]:<5} '
                  f'{l["wall"]:.0f}s  denserr(mn/p95/mx)='
                  f'{l["denserr_mean"]:.4f}/{l["denserr_p95"]:.4f}/'
                  f'{l["denserr_max"]:.4f}  force(mn/p95/mx)='
                  f'{l["force_mean"]:.2e}/{l["force_p95"]:.2e}/'
                  f'{l["force_max"]:.2e}')
        if args.plot:
            write_denserr_plot(args.out, res.pos, mf, res.h, res.rho, density,
                               box, args.nsmooth, flavour=args.flavour,
                               rhopow=args.rhopow, nsmoothmin=args.nsmoothmin)
        return

    callback = None
    if args.oi > 0:
        def callback(st):
            write_snapshot(f'{args.out}.{st.it:05d}', st.pos, st.rho, st.h)

    res = relax(pos, mass, density, box=box, nsmooth=args.nsmooth,
                flavour=args.flavour, rhopow=args.rhopow,
                icr0=args.icr0, icrloop0=args.icrloop0,
                icr0rate=args.icr0rate, max_iter=args.max_iter,
                margin=args.margin, one_sided=args.one_sided,
                hmode=args.hmode, nsmoothmin=args.nsmoothmin,
                icr0_stop=args.icr0_stop, polish_iter=args.polish_iter,
                self_bias=args.self_bias,
                callback=callback, callback_every=args.oi or args.max_iter,
                verbose=True)

    write_snapshot(f'{args.out}_IC', res.pos, res.rho, res.h)
    err = np.abs(1.0 - res.rho / density.rho0(res.pos))
    print(f'done: {res.niter} iterations, converged={res.converged}, '
          f'mean|1-rho/rho0|={err.mean():.3e}, wrote {args.out}_IC')
    if args.plot:
        write_denserr_plot(args.out, res.pos, mass, res.h, res.rho, density,
                           box, args.nsmooth, flavour=args.flavour,
                           rhopow=args.rhopow, nsmoothmin=args.nsmoothmin)


if __name__ == '__main__':
    main()
