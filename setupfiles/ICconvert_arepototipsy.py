#!/usr/bin/env python
"""Convert AREPO snapshot HDF5 output back into gasoline/tipsy snapshots.

The AREPO analog of `ICconvert_sphexatotipsy.py`: AREPO writes its moving-mesh
output as `snap_NNN.hdf5` files (PartType0 gas: Coordinates, Velocities, Masses,
Density, InternalEnergy, Volume, and MagneticField when built with MHD). The
test-suite analysis pipeline reads tipsy, so each AREPO snapshot is turned into a
gasoline tipsy snapshot `<output>.<NNNNN>` plus the usual aux files.

    ICconvert_arepototipsy.py <input> <output> <idx_interval>

  input         AREPO run directory containing snap_*.hdf5 (or a single .hdf5).
  output        run name / output directory for the tipsy snapshots.
  idx_interval  snapshot stride (= iOutInterval): AREPO snap_K -> tipsy
                <output>.{K*idx_interval:05d}, so snap_000 (t=0) -> .00000.

Field mapping (AREPO -> tgdata column):
  Masses->0, Coordinates->1:4 (re-centred to the gasoline [-L/2,L/2] box),
  Velocities->4:7, Density->7, InternalEnergy->8 as temperature (T = u/dTuFac),
  reconstructed smoothing length->9 (getsmooth2, also a `smoothlength` aux).
  MagneticField (MHD) -> BFieldx/y/z aux (+ DivB if present).
"""
import os
import sys
import glob

import numpy as np
import h5py

import readtipsy as tip

# Gasoline .log column header (same as the SPH-EXA converter) so the analysis
# framework's log reader is happy; we only fill the BOX PARAMETERS line + header
# (no per-step energy series -> analyses fall back to per-snapshot sums).
ENERGY_HEADER = (
    "# [dTime] [z] [Etot] [Ekin] [Epot] [Eth] [Emag] "
    "[totentrop] [totenstro] "
    "[Lx] [Ly] [Lz] [Llinx] [Liny] [Llinz] "
    "[cmx] [cmy] [cmz] [MWxy] [MWyz] [MWzx] [RSxy] [RSyz] [RSzx] "
    "[divBAvg] [divBMax] [divBerrAvg] [divBerrMax] "
    "[alphaAvg] [alphaMax] [betaMax] [betaAvg] [betaMin] "
    "[etaresAvg] [kinviscAvg] [rmsmach] [vrms] [rhogasAvg] [rhogasMax] "
    "[Q1Avg] [Q1Max] [Q2Avg] [Q2Max] [E0Avg] [E0Max] [Q4Avg] [Q4Max] "
    "[WallTime] [dWMax] [dImax] [dEMax] [dMultiEff]"
)


def _get_tufac(rundir, default=1.0):
    """dTuFac for T = u / dTuFac (inverse of the IC's u = dTuFac*T)."""
    try:
        from IC_analysis_general import read_tufac
        return read_tufac(rundir, default)
    except Exception:
        return default


def _energy_txt_rows(energy_txt):
    """Yield gasoline .log value-rows from an AREPO energy.txt (one per timestep).

    AREPO energy.txt columns: 0=time, 1=E_int (thermal), 2=E_pot, 3=E_kin, then
    per-type int/pot/kin (4..21), per-type mass (22..27), energy injection (28).
    There is NO magnetic energy in energy.txt (E_int = Mass*Utherm only, see
    src/io/global.c), so Emag is left 0 -- the magnetic-energy analyses detect the
    all-zero Emag and fall back to the (sparse) snapshot BField integration.
    Maps time->dTime, E_int->Eth, E_pot->Epot, E_kin->Ekin, Etot=int+pot+kin."""
    ncol = len(ENERGY_HEADER.split()) - 1        # named columns (drop leading '#')
    data = np.loadtxt(energy_txt, ndmin=2)
    for r in data:
        if r.size < 4:
            continue
        out = [0.0] * ncol
        out[0] = r[0]                  # dTime
        out[2] = r[1] + r[2] + r[3]    # Etot = E_int + E_pot + E_kin
        out[3] = r[3]                  # Ekin
        out[4] = r[2]                  # Epot
        out[5] = r[1]                  # Eth   (E_int = thermal only)
        yield out


def write_log(output_log, periods, energy_txt=None):
    """Write a gasoline-style .log: the periodic-box header (so the framework reads
    dxPeriod/dyPeriod/dzPeriod) + the energy column header, followed by one DENSE
    per-timestep row parsed from AREPO's energy.txt when available (far denser than
    the snapshot dumps -- the preferred source for energy-vs-time diagnostics). If
    no energy.txt is found the .log is header-only and analyses fall back to
    per-snapshot sums. Other AREPO *.txt (cpu/balance/domain/timebins/memory/info)
    are runtime diagnostics with no physics time-series, so they are not used."""
    px, py, pz = (float(p) for p in periods)
    with open(output_log, "w") as f:
        f.write(f"# BOX PARAMETERS:  bPeriodic: 1 dPeriod: {px:.8g} "
                f"dxPeriod: {px:.8g} dyPeriod: {py:.8g} dzPeriod: {pz:.8g}\n")
        f.write(ENERGY_HEADER + "\n")
        if energy_txt and os.path.exists(energy_txt):
            n = 0
            for out in _energy_txt_rows(energy_txt):
                f.write(" ".join(f"{x:.8e}" for x in out) + "\n")
                n += 1
            print(f"  wrote {n} dense energy rows from {energy_txt}")
        else:
            print("  (.log header only -> analyses use per-snapshot sums)")


def convert_snapshot(path, output, out_idx, tu_fac, periods=None):
    import IC_smoothlength as smth
    with h5py.File(path, "r") as f:
        hdr = f["Header"].attrs
        time = np.array([float(np.atleast_1d(hdr["Time"])[0])])
        g = f["PartType0"]
        n = g["Coordinates"].shape[0]

        tgdata = np.zeros((n, 12))
        # AREPO box is [0, period_axis] per axis (BoxSize*LONG_*); gasoline tipsy
        # is centred. Centre EACH axis on its own period (from the .param) -- a
        # single 0.5*BoxSize shift mis-places a thin (LONG_Z<1) slab. Fall back to
        # the per-axis data centroid when periods are unknown.
        coords = np.asarray(g["Coordinates"][:], float)
        if periods is not None:
            coords = coords - 0.5 * np.asarray(periods, float).reshape(3)
        else:
            coords = coords - 0.5 * (coords.min(axis=0) + coords.max(axis=0))
        tgdata[:, 1:4] = coords
        if "Velocities" in g:
            tgdata[:, 4:7] = g["Velocities"][:]
        mass = np.asarray(g["Masses"][:], float) if "Masses" in g \
            else np.zeros(n)
        tgdata[:, 0] = mass
        rho = np.asarray(g["Density"][:], float) if "Density" in g \
            else np.zeros(n)
        tgdata[:, 7] = rho
        if "InternalEnergy" in g:                       # temperature = u/dTuFac
            tgdata[:, 8] = np.asarray(g["InternalEnergy"][:], float) / tu_fac
        # reconstruct an SPH-equivalent smoothing length from mass/density
        hsm = (smth.getsmooth2(mass, rho, 64)
               if np.any(rho > 0) else np.zeros(n))
        tgdata[:, 9] = hsm

    outname = os.path.join(output, f"{output}.{out_idx:05d}")
    data_header = np.array([n, 3, n, 0, 0, 0])
    tip.writetipsy(tgdata, [], [], outname, data_header, time)
    tip.writetipsyaux(hsm, "smoothlength", outname)

    # Magnetic field (only present if AREPO was built with MHD).
    with h5py.File(path, "r") as f:
        g = f["PartType0"]
        if "MagneticField" in g:
            B = np.asarray(g["MagneticField"][:], float)
            tip.writetipsyaux(np.ascontiguousarray(B[:, 0]), "BFieldx", outname)
            tip.writetipsyaux(np.ascontiguousarray(B[:, 1]), "BFieldy", outname)
            tip.writetipsyaux(np.ascontiguousarray(B[:, 2]), "BFieldz", outname)
        if "MagneticFieldDivergence" in g:
            tip.writetipsyaux(np.asarray(g["MagneticFieldDivergence"][:], float),
                              "DivB", outname)
    print(f"Wrote {outname} from {os.path.basename(path)}")


def _periods_from_param(paramfile):
    """Per-axis periods (dxPeriod,dyPeriod,dzPeriod) from a gasoline .param, or
    None. AREPO's snapshot only carries the scalar BoxSize (= x-period); the y/z
    periods live in compile-time LONG_Y/LONG_Z, so the true non-cubic box has to
    come from the original IC's .param."""
    if not paramfile or not os.path.exists(paramfile):
        return None
    vals = {}
    with open(paramfile) as f:
        for line in f:
            line = line.split("#")[0]
            if "=" in line:
                k, _, v = line.partition("=")
                vals[k.strip()] = v.strip()

    def g(k):
        try:
            return float(vals[k])
        except (KeyError, ValueError):
            return None
    px, py, pz = g("dxPeriod"), g("dyPeriod"), g("dzPeriod")
    if px and py and pz:
        return (px, py, pz)
    return None


def convertfile(entry, output, idx_interval=1, paramfile=None):
    os.makedirs(output, exist_ok=True)
    if os.path.isdir(entry):
        snaps = sorted(glob.glob(os.path.join(entry, "snap_*.hdf5")))
    elif os.path.isdir(os.path.join(entry, "output")):
        snaps = sorted(glob.glob(os.path.join(entry, "output", "snap_*.hdf5")))
    else:
        snaps = [entry]
    if not snaps:
        print(f"No AREPO snap_*.hdf5 found under '{entry}'", file=sys.stderr)
        sys.exit(1)
    print(f"Found {len(snaps)} AREPO snapshot(s)")

    tu_fac = _get_tufac(output)
    periods = _periods_from_param(paramfile)
    if periods is None:
        with h5py.File(snaps[0], "r") as f:
            box = float(np.atleast_1d(f["Header"].attrs["BoxSize"])[0])
        periods = (box, box, box)
        print(f"  (no .param) assuming cubic box {box:g}; per-axis centring "
              "falls back to the data centroid.")
    else:
        print(f"  box periods from .param: {periods}")
    # AREPO's energy.txt is thermal/kinetic/potential only -- it has NO magnetic
    # energy (E_int = Mass*Utherm; see src/io/global.c). So for an MHD run the
    # dense energy.txt total would be hydro-only and Emag would be missing. Rely on
    # the snapshots (BField aux) for BOTH Etot and Emag there: a header-only .log
    # forces every energy quantity to come from the snapshots. Hydro runs (no
    # MagneticField in the snapshot) get the dense energy.txt series.
    # NB: the MHD-enabled AREPO binary writes a MagneticField array even for a
    # hydro run (all zeros, since the IC has no B), so test for a NON-ZERO field,
    # not mere presence -- else sedov (B==0) would be mis-flagged as MHD.
    with h5py.File(snaps[0], "r") as f:
        g = f["PartType0"]
        is_mhd = ("MagneticField" in g
                  and bool(np.any(np.asarray(g["MagneticField"][:]) != 0.0)))
    if is_mhd:
        print("  MHD run (snapshots carry B) -> header-only .log; Etot and Emag "
              "are computed from snapshots (AREPO logs no magnetic energy).")
        energy_txt = None
    else:
        # AREPO writes energy.txt into the OutputDir (= run dir, beside the snaps).
        snapdir = os.path.dirname(snaps[0]) or "."
        energy_txt = next((p for p in (os.path.join(snapdir, "energy.txt"),
                                       os.path.join(entry, "energy.txt"))
                           if os.path.exists(p)), None)
    write_log(os.path.join(output, output + ".log"), periods, energy_txt)

    # only pass periods for centring when they came from the .param (else let
    # convert_snapshot use the per-axis data centroid).
    cen = periods if paramfile and _periods_from_param(paramfile) else None
    for idx, snap in enumerate(snaps):
        convert_snapshot(snap, output, idx * idx_interval, tu_fac, periods=cen)


if __name__ == "__main__":
    if len(sys.argv) >= 4:
        entry = sys.argv[1]
        output = sys.argv[2]
        idx_interval = int(sys.argv[3])
        paramfile = sys.argv[4] if len(sys.argv) >= 5 else None
    else:
        entry = input("AREPO run dir / snapshot: ")
        output = input("Output name: ")
        idx_interval = int(input("idx_interval (iOutInterval): "))
        paramfile = input("original tipsy .param (blank to skip): ").strip() or None
    convertfile(entry, output, idx_interval, paramfile)
