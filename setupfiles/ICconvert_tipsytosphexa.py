import readtipsy as tip
import numpy as np
import pdb
import h5py
import sys
from scipy.spatial import KDTree

step_string = "/Step#0/"


def read_param_file(paramfile):
    params = {}

    with open(paramfile, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if "=" in line:
                key, val = line.split("=", 1)
                key = key.strip()
                val = val.strip()

                # try numeric conversion
                try:
                    if "." in val or "e" in val or "E" in val:
                        val = float(val)
                    else:
                        val = int(val)
                except ValueError:
                    pass

                params[key] = val

    return params

def convertfile(entry,output):
    base = entry
    if entry.endswith(".00000"):
        base = entry[:-6]
    params = read_param_file(base + ".param")
    tgdata,tddata,tsdata,data_header,time=tip.readtipsy(entry);
    N=data_header[0]
    ngas=data_header[2]
    ndark=data_header[3]
    nstar=data_header[4]

    positions = tgdata[:,1:4]
    mass = tgdata[:,0]
    velocities = tgdata[:,4:7]
    density = tgdata[:,7]
    temp = tgdata[:,8]
    soft = tgdata[:,9]

    num_particles = ngas
    ids = np.arange(num_particles)
    
    du_m1 = np.zeros(num_particles, dtype=np.float32)

    alphamin = 1.0
    alpha = np.ones(num_particles, dtype=np.float32) * alphamin


    file_exa = h5py.File(output+".h5", "w")

    tree = KDTree(positions)
    smoothings = np.zeros(num_particles)
    for i, p in enumerate(positions):
        distances, indices = tree.query(p, k=100)
        smoothings[i] = np.max(distances) / 2.


    minDt = 1e-9

    file_exa["/Step#0/x"] = positions[:, 0]
    file_exa["/Step#0/y"] = positions[:, 1]
    file_exa["/Step#0/z"] = positions[:, 2]
    file_exa["/Step#0/id"] = ids

    #both a temp and u in sphexa
    #Gas constant in SPHEXA set to 8.317e7, in Gasoline msr->param.dGasConst = msr->param.dKpcUnit*KPCCM*KBOLTZ/MHYDR/GCGS/msr->param.dMsolUnit/MSOLG; 1kpc 1msol = 1921.03
    #temperature needs to be scaled down by difference
    file_exa["/Step#0/temp"] = temp*1921.03/8.317e7
    file_exa["/Step#0/u"] = temp*1921.03*1.5/1.0
    file_exa["/Step#0/du_m1"] = du_m1
    file_exa["/Step#0/alpha"] = alpha
    file_exa["/Step#0/x_m1"] = velocities[:, 0] * minDt #diff in position from last dt
    file_exa["/Step#0/y_m1"] = velocities[:, 1] * minDt
    file_exa["/Step#0/z_m1"] = velocities[:, 2] * minDt
    file_exa["/Step#0/vx"] = velocities[:, 0]
    file_exa["/Step#0/vy"] = velocities[:, 1]
    file_exa["/Step#0/vz"] = velocities[:, 2]

    # MHD detection from the BField aux VALUES: the .BFieldx/y/z aux files are
    # written even for hydro runs, so their presence alone isn't enough -- treat
    # the run as MHD only if the field is actually non-zero. Falls back to hydro
    # if the aux files are absent.
    try:
        Bx, _ = tip.readtipsyaux(entry, 'BFieldx')
        By, _ = tip.readtipsyaux(entry, 'BFieldy')
        Bz, _ = tip.readtipsyaux(entry, 'BFieldz')
        mhd = bool(np.any((Bx != 0) | (By != 0) | (Bz != 0)))
    except FileNotFoundError:
        mhd = False
    if mhd:
        file_exa["/Step#0/magneto::Bx"] = Bx
        file_exa["/Step#0/magneto::By"] = By
        file_exa["/Step#0/magneto::Bz"] = Bz
        file_exa["/Step#0/magneto::dBx"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::dBy"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::dBz"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::dBx_m1"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::dBy_m1"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::dBz_m1"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::psi_ch"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::d_psi_ch"] = np.zeros(num_particles, dtype=np.float32)
        file_exa["/Step#0/magneto::d_psi_ch_m1"] = np.zeros(num_particles, dtype=np.float32)
    file_exa["/Step#0/m"] = mass
    file_exa["/Step#0/h"] = smoothings

    bPeriodic = int(params.get("bPeriodic", 0))
    if bPeriodic == 1:
        file_exa["/Step#0"].attrs["boundaryType"] = np.array([1, 1, 1], dtype=np.int8)
        dxP = float(params.get("dxPeriod", 1.0))
        dyP = float(params.get("dyPeriod", 1.0))
        dzP = float(params.get("dzPeriod", 1.0))
        xmin, xmax = -dxP*0.5, dxP*0.5
        ymin, ymax = -dyP*0.5, dyP*0.5
        zmin, zmax = -dzP*0.5, dzP*0.5
    else:
        file_exa["/Step#0"].attrs["boundaryType"] = np.array([0, 0, 0], dtype=np.int8)

        xmin = np.min(positions[:, 0])
        xmax = np.max(positions[:, 0])
        ymin = np.min(positions[:, 1])
        ymax = np.max(positions[:, 1])
        zmin = np.min(positions[:, 2])
        zmax = np.max(positions[:, 2])

    file_exa["/Step#0"].attrs["time"] = 0.0
    file_exa["/Step#0"].attrs["minDt"] = minDt
    file_exa["/Step#0"].attrs["minDt_m1"] = minDt

    bDoGravity = int(params.get("bDoGravity", 0))

    if bDoGravity == 1:
        file_exa["/Step#0"].attrs["gravConstant"] = 1.0
    else:
        file_exa["/Step#0"].attrs["gravConstant"] = 0.0


    gamma = float(params.get("dConstGamma", 5.0 / 3.0))
    mui = float(params.get("dMeanMolWeight", 1.0))

    gasiso = int(params.get("bGasIsothermal",0))
    gasadi = int(params.get("bGasAdiabatic",1))

    dDiskHr = float(params.get("dDiskHr", 0.0))
    dDiskq = float(params.get("dDiskq", 0.0))
    dPhysViscEta = float(params.get("dPhysViscEta", 0.0))
    dPhysViscXi = float(params.get("dPhysViscXi", 0.0))
    if dPhysViscEta > 0.0:
        file_exa["/Step#0"].attrs["eta_visc"] = dPhysViscEta
        file_exa["/Step#0"].attrs["xi_visc"]  = dPhysViscXi
        file_exa["/Step#0"].attrs["Kvisc"] = 0.2
    if gasiso == 1:
        if dDiskHr > 0.0:
            file_exa["/Step#0"].attrs["eosChoice"] = np.int32(3)  # locally isothermal
            file_exa["/Step#0"].attrs["cs0_isor"] = dDiskHr
            file_exa["/Step#0"].attrs["q_isor"]   = dDiskq
            file_exa["/Step#0"].attrs["gravConstant"] = 0.0 #low mass gas disk
        else:
            file_exa["/Step#0"].attrs["eosChoice"] = np.int32(1)  # isothermal
    else:
        file_exa["/Step#0"].attrs["eosChoice"] = np.int32(0)  # adiabatic
    file_exa["/Step#0"].attrs["gamma"] = gamma
    file_exa["/Step#0"].attrs["muiConst"] = mui

    file_exa["/Step#0"].attrs["etaAcc"] = 0.2
    file_exa["/Step#0"].attrs["Kcour"] = 0.2
    file_exa["/Step#0"].attrs["alphamin"] = alphamin
    file_exa["/Step#0"].attrs["alphamax"] = alphamin+0.01
    file_exa["/Step#0"].attrs["kernelChoice"] = np.int32(0) #WENDLAND C2

    file_exa["/Step#0"].attrs["box"] = np.array([xmin, xmax, ymin, ymax, zmin, zmax], dtype=np.float64)

    file_exa["/Step#0"].attrs["iteration"] = np.uint64(0)
    file_exa["/Step#0"].attrs["numParticlesGlobal"] = np.uint64(num_particles)
    file_exa["/Step#0"].attrs["ng0"] = np.int32(100)
    file_exa["/Step#0"].attrs["ngmax"] = np.int32(104)

    if nstar == 1:
        file_exa["/Step#0"].attrs["star::m"]          = tsdata[0,0]      # star mass (in your code units, GM=1 → m=1/G)
        file_exa["/Step#0"].attrs["star::x"]          = tsdata[0,1]      # star position
        file_exa["/Step#0"].attrs["star::y"]          = tsdata[0,2]
        file_exa["/Step#0"].attrs["star::z"]          = tsdata[0,3]
        file_exa["/Step#0"].attrs["star::x_m1"]       = tsdata[0,1]
        file_exa["/Step#0"].attrs["star::y_m1"]       = tsdata[0,2]
        file_exa["/Step#0"].attrs["star::z_m1"]       = tsdata[0,3]
        sinkradius=float(params.get("dSinkRadius", 1.0))
        file_exa["/Step#0"].attrs["star::inner_size"] = sinkradius     # accretion radius — particles inside are removed
        file_exa["/Step#0"].attrs["star::fixed_star"] = np.int32(1)          # 1 = fixed, 0 = let star move
        file_exa["/Step#0"].attrs["star::potentialType"] = np.int32(0)       # 0=Newtonian, 1=Einstein precession
        file_exa["/Step#0"].attrs["removeUnconvergedParticles"] = np.int32(1)
    file_exa.close()

# Check for arguments
if len(sys.argv) >= 3:
    entry = sys.argv[1]
    output = sys.argv[2]
else:
    entry = input("File name: ")
    output = input("Output name: ")
convertfile(entry,output)
