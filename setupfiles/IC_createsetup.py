"""Create the initial-condition tipsy file (+ aux fields + .param) for one test.

CLI (the runtest.sh contract -- both its lattice/rand and glass invocations):

    python IC_createsetup.py <testsim> <nx> <vm> <entry> <output> [extra ...]

    testsim  test name (run with no arguments to list the choices)
    nx       resolution (number of elements along x)
    vm       variable SPH particle masses (0 = equal-mass particles)
    entry    "0" = lattice, "1" = random, anything else = path to a relaxed
             glass IC file to rebuild the final IC from
    output   output name: writes datafiles/<output>.00000 (+ aux, .param)
    extra    per-test parameters, matched POSITIONALLY to the selected
             setup's create() signature after (nx, distri, vm, entry);
             omitted ones keep the create() defaults

Cosmological-unit conventions (for the comoving tests):

    Boxsize = 1 = length unit. Lengths and velocities are comoving;
    velocities are peculiar from the Hubble flow. The critical density of
    the Universe in simulation units is 1, so the total mass in the box in
    simulation units is rho/rho_crit = Omega_matter at z = 0. From the
    Friedmann equation, H_0^2 = 8 pi rho_crit / 3; hence H_0 =
    sqrt(8 pi/3) = 2.894405 in simulation units. The velocity unit is
    (Hubble velocity across the box)/2.894405; the time unit is
    2.894405/(Hubble parameter).

    Example: Omega0 = 1, box size = 100 Mpc at z=0, h=0.5 (sigma8 and P(k)
    are not important): mass unit = (rho_crit in Msun/Mpc^3) * (100 Mpc)^3
    = 6.398e16 Msun; total mass in the box = 1 code mass unit; box size =
    100 Mpc/(1+z); velocity unit = 1727.47/(1+z) = (50*100/2.894405)/(1+z).

Note: closed-packed distributions are not optimal when sides are equal, as
this generates small asymmetries due to uneven spacing.
"""

import argparse
import importlib
import inspect
import sys

import numpy as np

import readtipsy as tip
from IC_createparamfile import createparamfile
from IC_writeascii import convertfromcodeunits, calculatecosmoMsol

# Every name maps to IC_setup_<name>.setup_<name>, imported lazily on use.
SETUPS = [
    'orzag',
    'evrard',
    'polytrope',
    'planet',
    'planetcollision',
    'areablob',
    'mhdrotor',
    'mhddivtest',
    'mhdblast',
    'mhdblast_garcia',
    'mhdpinch_garcia',
    'kh_garcia',
    'sedov',
    'shocktube',
    'mhdloop',
    'mhdloopshear',
    'cube',
    'alfven',
    'cosmowave',
    'taylorgreen',
    'mhdslabcond',
    'mri',
    'mti',
    'kh',
    'square',
    'khslab',
    'rt',
    'zeldovich',
    'mhdcollapse',
    'rotatingcube',
    'gresho',
    'mhdisowave',
    'wave',
    'mhdbalsaravortex',
    'currentsheet',
    'noh',
    'isentropicvortex',
    'coldkeplerian',
    'accdisk',
    'blob',
]


def list_setups():
    print("Available options:")
    for option in SETUPS:
        print(f"- {option}")


def load_setup(setup_case):
    module = importlib.import_module(f"IC_setup_{setup_case}")
    return getattr(module, f"setup_{setup_case}")()


def resolve_distribution(entry):
    """entry "0"/"1" select lattice/random; any other value is a glass IC file."""
    if entry == "0":
        return 0, "lattice"
    if entry == "1":
        return 1, "rand"
    return 2, "glass"


def coerce_extra(value, default):
    """Coerce a CLI extra to its create() parameter type.

    String-valued parameters (planet preset names, planetcollision
    target/impactor) pass through; everything else follows the historical
    all-floats convention. A non-numeric value for a None-default parameter
    also passes through as a string.
    """
    if isinstance(default, str):
        return value
    try:
        return float(value)
    except ValueError:
        return value


def main():
    if len(sys.argv) < 2:
        list_setups()
        sys.exit(0)

    parser = argparse.ArgumentParser(
        description="Create the IC tipsy file (+ aux, .param) for one test case.",
        usage="%(prog)s testsim nx vm entry output [extra ...]",
        epilog=__doc__.split("Cosmological")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("testsim", help="test name (no arguments lists the choices)")
    parser.add_argument("nx", type=float, help="resolution (elements along x)")
    parser.add_argument("vm", type=float, help="variable particle masses (0 = equal mass)")
    parser.add_argument("entry", help='"0" lattice, "1" random, else a glass IC file')
    parser.add_argument("output", help="written as datafiles/<output>.00000")
    # The per-test extras are sliced off BEFORE argparse so values like
    # "-1e-4" (negative scientific notation, which argparse would reject as
    # an unknown option) and arbitrary preset strings survive untouched.
    args = parser.parse_args(sys.argv[1:6])
    extras = sys.argv[6:]

    if args.testsim not in SETUPS:
        list_setups()
        sys.exit(f"Unknown setup_case: {args.testsim}")

    sim = load_setup(args.testsim)
    distri, distribution_name = resolve_distribution(args.entry)
    glasssetup = 1 if distri == 1 else 0

    standard_args = {
        'nx': args.nx,
        'distri': distri,
        'vm': args.vm,
        'entry': args.entry,
    }
    sig = inspect.signature(sim.create)
    extra_args = {param.name: param.default for param in sig.parameters.values()
                  if param.name not in standard_args}
    for name, value in zip(extra_args, extras):
        extra_args[name] = coerce_extra(value, extra_args[name])
    all_args = {**standard_args, **extra_args}
    print(f"IC_createsetup[{args.testsim}]: {distribution_name} IC, parameters {all_args}")

    ## Generate gas IC
    sim.create(**all_args)

    zero = [0.0] * len(sim.x)
    if args.testsim != 'mti' and args.testsim != 'blob':
        sim.metals = zero
    sim.pot = zero
    sim.Bpsi = zero
    ## (data,kpc,msol,ifnottemp,phystocodevel,phystocodeB)
    if args.testsim == 'mhdcollapse':
        sim = convertfromcodeunits(sim, 0.001, 1000, 1, 1, 1)
    elif args.testsim in ('areablob', 'planet', 'planetcollision'):
        sim = convertfromcodeunits(sim, 1, 1, 0, 0, 0)
    else:
        sim = convertfromcodeunits(sim, sim.dkpcunit, sim.dmsolunit, 1, 0, 0)

    # COSMO: comoving run -- recompute the mass unit from the box size and
    # start the snapshot clock at the IC scale factor.
    if sim.cosmo == 1:
        h0 = 1.0
        ascale = sim.a
        sim.dmsolunit = calculatecosmoMsol(h0, sim.dkpcunit)
        sim = convertfromcodeunits(sim, sim.dkpcunit, sim.dmsolunit, 1, 0, 0)

    npart = len(sim.x)
    tgdata = np.column_stack((sim.mass, sim.x, sim.y, sim.z,
                              sim.vx, sim.vy, sim.vz,
                              sim.rho, sim.u, sim.h, sim.metals, sim.pot))
    B = np.column_stack((sim.Bx, sim.By, sim.Bz))

    # Tipsy header (none of the tests include DM or stars at this point)
    data_header = np.array([npart, 3, npart, 0, 0, 0])

    if sim.cosmo == 1:
        a2 = np.geomspace(ascale, 1, num=300)
        z2 = 1 / a2 - 1
        np.savetxt("mhdcosmodivteststd.red", z2, fmt='%.4f')
        time = np.array([ascale])
    else:
        time = np.array([0.0])

    tddata = []
    tsdata = []

    # Generate dark-matter / star (e.g. central sink/BH) particles if the setup
    # defines them. Each hook returns (data_array, nparticles); the per-test logic
    # lives in the IC_setup_*.py files (keeps this driver generic).
    if glasssetup == 0 and hasattr(sim, 'create_dark'):
        tddata, ndark = sim.create_dark(tgdata)
        data_header[3] = data_header[3] + ndark
        data_header[0] = data_header[0] + ndark

    if glasssetup == 0 and hasattr(sim, 'create_star'):
        tsdata, nstar = sim.create_star(tgdata)
        data_header[4] = data_header[4] + nstar
        data_header[0] = data_header[0] + nstar

    myoutput = "datafiles/" + args.output + ".00000"
    tip.writetipsy(tgdata, tddata, tsdata, myoutput, data_header, time)
    tip.writealltipsyaux(B, myoutput)

    if args.testsim == 'planetcollision':
        tip.writealltipsyauxB(sim.Xg, myoutput, sim.labels)

    if args.testsim == 'planet':
        for i in range(1, 6):
            tip.writetipsyaux(getattr(sim, f"Mat{i}"), f"Material{i}", myoutput)

    createparamfile(args.output, sim, glasssetup, args.testsim)
    print('time', time)
    print("created file: ", args.output)


if __name__ == "__main__":
    main()
