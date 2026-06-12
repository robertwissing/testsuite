import readtipsy as tip
import numpy as np
import h5py
import os
import sys

def get_attr_scalar(attr):
        val = attr[()]
        return val[0] if np.ndim(val) > 0 else val

def convert_logfile(constants_file, output_log, box=None):
    """
    Convert SPHExA constants.txt into a GASOLINE-style .log file.
    Missing quantities are filled with zeros. `box` is the SPH-EXA step attrs
    "box" [x0,x1,y0,y1,z0,z1]; when given, a GASOLINE-style "# BOX PARAMETERS"
    header line is written so the analysis framework can read the domain
    (dxPeriod/dyPeriod/dzPeriod) from the .log.
    """

    header = (
        "# [dTime] [z] [Etot] [Ekin] [Epot] [Eth] [Emag] "
        "[totentrop] [totenstro] "
        "[Lx] [Ly] [Lz] "
        "[Llinx] [Liny] [Llinz] "
        "[cmx] [cmy] [cmz] "
        "[MWxy] [MWyz] [MWzx] "
        "[RSxy] [RSyz] [RSzx] "
        "[divBAvg] [divBMax] [divBerrAvg] [divBerrMax] "
        "[alphaAvg] [alphaMax] "
        "[betaMax] [betaAvg] [betaMin] "
        "[etaresAvg] [kinviscAvg] "
        "[rmsmach] [vrms] "
        "[rhogasAvg] [rhogasMax] "
        "[Q1Avg] [Q1Max] [Q2Avg] [Q2Max] "
        "[E0Avg] [E0Max] [Q4Avg] [Q4Max] "
        "[WallTime] [dWMax] [dImax] [dEMax] [dMultiEff]"
    )

    # SPH-EXA emits MHD columns (eMag, meanDivBError, maxDivBError at cols 9-11)
    # only when md.eMag > 0; otherwise constants.txt has 9 columns.
    data = np.loadtxt(constants_file, ndmin=2)
    if data.size == 0:
        print(f"WARNING: '{constants_file}' is EMPTY -- the .log will have no "
              f"energy/time series (analyses fall back to per-snapshot sums).",
              file=sys.stderr)
        data = data.reshape(0, 1)
    ncols = data.shape[1]

    def col(row, i):
        return row[i] if i < ncols else 0.0

    with open(output_log, "w") as f:
        if box is not None:
            dx, dy, dz = box[1] - box[0], box[3] - box[2], box[5] - box[4]
            f.write(f"# BOX PARAMETERS:  bPeriodic: 1 dPeriod: {dx:.8g} "
                    f"dxPeriod: {dx:.8g} dyPeriod: {dy:.8g} "
                    f"dzPeriod: {dz:.8g}\n")
        f.write(header + "\n")

        for row in data:
            # SPHExA constants.txt mapping
            dTime    = col(row, 1)   # d.ttot
            Etot     = col(row, 3)   # d.etot
            Ekin     = col(row, 4)   # d.ecin
            Eth      = col(row, 5)   # d.eint
            Epot     = col(row, 6)   # d.egrav
            Emag     = col(row, 9)   # md.eMag
            divBAvg  = col(row, 10)  # md.meanDivBError
            divBMax  = col(row, 11)  # md.maxDivBError

            out = [
                dTime, 0.0, Etot, Ekin, Epot, Eth, Emag,
                0.0, 0.0,                     # entropy
                0.0, 0.0, 0.0,                # Lx Ly Lz
                0.0, 0.0, 0.0,                # linear momentum
                0.0, 0.0, 0.0,                # center of mass
                0.0, 0.0, 0.0,                # Maxwell stress
                0.0, 0.0, 0.0,                # Reynolds stress
                0.0, 0.0, divBAvg, divBMax,   # divB
                0.0, 0.0,                     # alpha
                0.0, 0.0, 0.0,                # beta
                0.0, 0.0,                     # resistivity / viscosity
                0.0, 0.0,                     # mach / vrms
                0.0, 0.0,                     # rho gas
                0.0, 0.0, 0.0, 0.0,           # Q1/Q2
                0.0, 0.0, 0.0, 0.0,           # E0/Q4
                0.0, 0.0, 0.0, 0.0, 0.0        # walltime etc.
            ]

            f.write(" ".join(f"{x:.8e}" for x in out) + "\n")


def convert_step(file_exa, step, output, idx, idx_interval):
    group = file_exa[step]

    npart = get_attr_scalar(group.attrs["numParticlesGlobal"])
    time  = np.array(get_attr_scalar(group.attrs["time"]))
    ng0  = np.array(get_attr_scalar(group.attrs["ng0"]))
    data_header = np.array([npart, 3, npart, 0, 0, 0])
    tgdata = np.zeros((npart, 12))
    tddata = []
    tsdata = []

    # Mandatory positions (but still checked defensively)
    if "x" in group:
        tgdata[:, 1] = group["x"]
    if "y" in group:
        tgdata[:, 2] = group["y"]
    if "z" in group:
        tgdata[:, 3] = group["z"]

    # Velocities
    if "vx" in group:
        tgdata[:, 4] = group["vx"]
    if "vy" in group:
        tgdata[:, 5] = group["vy"]
    if "vz" in group:
        tgdata[:, 6] = group["vz"]

    # Density
    if "rho" in group:
        tgdata[:, 7] = group["rho"]

    # Internal energy (scaled)
    if "u" in group:
        tgdata[:, 8] = group["u"][:] / (1921.03 * 1.5)

    if "temp" in group:
        tgdata[:,8] = group["temp"][:] * 8.317e7 / 1921.03

    # Smooth-length
    if "h" in group:
        tgdata[:, 9] = group["h"]
    # mass
    if "m" in group:
        tgdata[:, 0] = group["m"]
    elif "h" in group and "rho" in group:
        tgdata[:, 0] = 3.0*ng0/(np.pi*32.0)*group["h"][:]*group["h"][:]*group["h"][:]*group["rho"][:]

    if "nc" in group:
        tgdata[:,11] = group["nc"]



    out_idx = idx_interval + idx * idx_interval
    outname = os.path.join(output, f"{output}.{out_idx:05d}")

    tip.writetipsy(tgdata, tddata, tsdata, outname, data_header, time)

    if "magneto::Bx" in group:
        tip.writetipsyaux(group["magneto::Bx"][:],"BFieldx",outname)

    if "magneto::By" in group:
        tip.writetipsyaux(group["magneto::By"][:],"BFieldy",outname)

    if "magneto::Bz" in group:
        tip.writetipsyaux(group["magneto::Bz"][:],"BFieldz",outname)

    if "magneto::divB" in group:
        tip.writetipsyaux(group["magneto::divB"][:],"DivB",outname)

    print(f"Wrote {outname} from {step}")


def convertfile(entry, output, idx_interval=1):
    # Ensure output directory exists
    os.makedirs(output, exist_ok=True)
    with h5py.File(entry, "r") as file_exa:

        # Collect and sort Step groups numerically
        steps = sorted(
            [k for k in file_exa.keys() if k.startswith("Step#")],
            key=lambda s: int(s.split("#")[1])
        )

        print(f"Found {len(steps)} steps")

        box = file_exa[steps[0]].attrs["box"][:] if steps else None
        try:
            convert_logfile("constants.txt", output + "/" + output + ".log",
                            box=box)
        except FileNotFoundError:
            convert_logfile(output + "/constants.txt",
                            output + "/" + output + ".log", box=box)

        for idx, step in enumerate(steps):
            convert_step(file_exa, step, output, idx, idx_interval)


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        entry = sys.argv[1]
        output = sys.argv[2]
        idx_interval = int(sys.argv[3])
    else:
        entry = input("File name: ")
        output = input("Output name: ")
        idx_interval = int(input("idx_interval: "))
    convertfile(entry, output, idx_interval)
