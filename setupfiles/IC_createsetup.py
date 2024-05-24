from IPython import get_ipython
def __reset__(): get_ipython().magic('reset -sf')
from IC_setup_orzag import setup_orzag
from IC_setup_mhdrotor import setup_mhdrotor
from IC_setup_mhddivtest import setup_mhddivtest
from IC_setup_mhdblast import setup_mhdblast
from IC_setup_mhdblast_garcia import setup_mhdblast_garcia
from IC_setup_mhdpinch_garcia import setup_mhdpinch_garcia
from IC_setup_kh_garcia import setup_kh_garcia
#from IC_setup_mhdloop_garcia import setup_mhdloop_garcia
from IC_setup_sedov import setup_sedov
from IC_setup_shocktube import setup_shocktube
from IC_setup_mhdloop import setup_mhdloop
from IC_setup_mhdloopshear import setup_mhdloopshear
from IC_setup_cube import setup_cube
from IC_setup_accdisk import setup_accdisk
from IC_setup_alfven import setup_alfven
from IC_setup_cosmowave import setup_cosmowave
from IC_setup_taylorgreen import setup_taylorgreen
from IC_setup_mhdslabcond import setup_mhdslabcond
from IC_setup_mri import setup_mri
from IC_setup_wave import setup_wave
from IC_setup_kh import setup_kh
from IC_setup_khslab import setup_khslab
from IC_setup_rt import setup_rt
from IC_setup_zeldovich import setup_zeldovich
from IC_setup_mhdblob import setup_mhdblob
from IC_setup_mhdcollapse import setup_mhdcollapse
from IC_setup_rotatingcube import setup_rotatingcube
from IC_setup_gresho import setup_gresho
from IC_setup_mhdisowave import setup_mhdisowave
from IC_setup_mhdbalsaravortex import setup_mhdbalsaravortex
from IC_writeascii import writeascii,convertfromcodeunits,calculatecosmoMsol
from matplotlib import pyplot as plot
from IC_createparamfile import createparamfile
from mpl_toolkits.mplot3d import Axes3D
#from mayavi import mlab
import numpy as np
from scipy import stats
import readtipsy as tip
import timeit
import pickle
import sys
import inspect
"""
Boxsize=1=length unit

Lengths and velocities are comoving.

Velocities are peculiar from the Hubble flow.

The critical density of the Universe in simulation units is 1.


The total mass in the box in simulation units is rho/rho_crit = Omega_{matter} at z = 0.
From the Friedmann equation, H_0^2 = 8 pi rho_crit / 3; hence, H_0 = sqrt(8 pi/3) = 2.894405 in simulation units.
The velocity unit is gotten from (Hubble velocity across the box/2.894405).
The time unit is 2.894405/(Hubble parameter).

Examples

Omega0 = 1 Box size = 100mpc at z=0 h=0.5

(just in case you are wondering: sigma8 and P(k) are not important)

mass unit = (rho_crit in Msun/Mpc^3) times (100 Mpc)^3 = 6.398e16 solar masses

Total mass in the box = 1 code mass unit

box size = 100Mpc/(1+z)

Velocity unit= 1727.47/(1+z) = (50 100 /2.894405)/(1+z)

"""

"""
 Closed packed distributions are not optimal when sides are equal as this will generate
small asymmetries due to not even spacing everywhere.
"""

#COSMO SETTINGS
h0=1.0
z=127
ascale=1./(1.+z)
#ascale = 0.0

glasssetup = 0



# Dictionary mapping simplified names to setup functions
setup_functions = {
    'orzag': setup_orzag,
    'mhdrotor': setup_mhdrotor,
    'mhddivtest': setup_mhddivtest,
    'mhdblast': setup_mhdblast,
    'mhdblast_garcia': setup_mhdblast_garcia,
    'mhdpinch_garcia': setup_mhdpinch_garcia,
    'kh_garcia': setup_kh_garcia,
    'sedov': setup_sedov,
    'shocktube': setup_shocktube,
    'mhdloop': setup_mhdloop,
    'mhdloopshear': setup_mhdloopshear,
    'cube': setup_cube,
    'alfven': setup_alfven,
    'cosmowave': setup_cosmowave,
    'taylorgreen': setup_taylorgreen,
    'mhdslabcond': setup_mhdslabcond,
    'mri': setup_mri,
    'wave': setup_wave,
    'kh': setup_kh,
    'khslab': setup_khslab,
    'rt': setup_rt,
    'zeldovich': setup_zeldovich,
    'mhdblob': setup_mhdblob,
    'mhdcollapse': setup_mhdcollapse,
    'rotatingcube': setup_rotatingcube,
    'gresho': setup_gresho,
    'mhdisowave': setup_mhdisowave,
    'mhdbalsaravortex': setup_mhdbalsaravortex,
    'accdisk': setup_accdisk,
}


# If sys.argv[0] is not given, list all available options
if len(sys.argv) < 2:
    print("Available options:")
    for option in setup_functions:
        print(f"- {option}")
    sys.exit(0)
    
# Determine the setup function based on setup_case
setup_case = sys.argv[1]

# Determine the setup function based on setup_case
if any(name in setup_case for name in setup_functions):
    # Find the key that matches the substring in setup_case
    matching_key = next(name for name in setup_functions if name in setup_case)
    # Call the corresponding setup function
    sim = setup_functions[matching_key]()
else:
    # Default case if no match is found
    raise ValueError(f"Unknown setup_case: {setup_case}")


n=float(sys.argv[2])
fileinput = sys.argv[3] if len(sys.argv) > 3 else ""
distri = 2 if fileinput else input("Enter wanted distribution (0: lattice 1: random): ")
distri = 0 if fileinput == "0" else 1 if fileinput == "1" else distri

print("distri: ",distri,fileinput)

distribution_name = "lattice" if distri == 0 else "rand" if distri == 1 else "glass" if distri == 2 else None

if distribution_name is None:
    print("Invalid distribution choice. Exiting.")
    sys.exit()

default_output = f"{setup_case}_{n}_{distribution_name}"
fileoutput = sys.argv[4] if len(sys.argv) > 4 else input(f"Enter fileoutput (default: {default_output}): ") or default_output

# Retrieve the signature of the create function
sig = inspect.signature(sim.create)
## Standard arguments for all IC
standard_args = {
    'nx': n,
    'distri': int(distri), 
    'vm': 0,
    'entry': fileinput
}
# Determine extra arguments by filtering out standard arguments
extra_args = {param.name: param.default for param in sig.parameters.values() if param.name not in standard_args}
print("Default value of extra arguments in ", setup_case, "  ", extra_args)
# Update extra_args with values from sys.argv, if provided
argv_index = 5  # Starting index in sys.argv for extra arguments
for arg_name in extra_args.keys():
    if argv_index < len(sys.argv):  # Check if there are enough arguments provided
        # Convert the argument to the correct type (int, float, str, etc.)
        # Here we assume all extra arguments are floats; you may need to adjust this based on your needs
        extra_args[arg_name] = float(sys.argv[argv_index])
    argv_index += 1  # Move to the next argument in sys.argv
all_args = {**standard_args, **extra_args}
print("Your setup parameter list: ", setup_case, "  ", all_args)
    
## Generate gas IC
sim.create(**all_args)
## Generate dark matter IC
#orz.create_dark(**all_args)
## Generate star IC
#orz.create_star(**all_args)

## File name
file_string=fileoutput


zero=[0]*len(sim.x);
sim.metals=zero
sim.pot=zero
sim.Bpsi=zero
## (data,kpc,msol,ifnottemp,phystocodevel,phystocodeB)
if setup_case == 'mhdcollapse':
    sim=convertfromcodeunits(sim,0.001,1000,1,1,1)
else:
    sim=convertfromcodeunits(sim,1,1,1,0,0)

#COSMO
if(sim.cosmo==1):
    ascale=sim.a
    sim.dmsolunit=calculatecosmoMsol(h0,sim.dkpcunit)
    sim=convertfromcodeunits(sim,sim.dkpcunit,sim.dmsolunit,1,0,0)


# Data Initialization
npart = len(sim.x)
tgdata = np.column_stack((sim.mass, sim.x, sim.y, sim.z, sim.vx, sim.vy, sim.vz, sim.rho, sim.u, sim.h, sim.metals, sim.pot))
B = np.column_stack((sim.Bx, sim.By, sim.Bz))


# Tipsy Header(None of the tests include DM or stars atm)
data_header=np.array([npart, 3, npart, 0, 0, 0])

if(sim.cosmo==1):
    a2=np.geomspace(ascale,1,num=300)
    z2=1/a2-1
    np.savetxt("mhdcosmodivteststd.red", z2,fmt='%.4f')
    time=np.array([ascale])
else:
    time=np.array([0.0])


tddata=[]
tsdata=[]
if distri == 1:
    glasssetup = 1
    


if setup_case == 'accdisk' and glasssetup == 0:
    if sim.onlyoneBH == 1:
        massprim = 1.0
        softprim = 0.5
        sinkprim = -1.0
        # Set cm velocity to primary
        tsdata=np.array([massprim,0.0,0.0,0.0,0.0,0.0,0.0,0.0,sinkprim,softprim,0.0])
        data_header[4]=data_header[4]+1;
        data_header[0]=data_header[0]+1;
    else:
        massprim = 1.0
        softprim = 0.5
        sinkprim = -1.0
        xsec = -sim.rorb
        vysec = -sim.rorb**(-sim.q)
        softsec = sim.hr*0.6
        masssec = massprim*sim.massrat
        tsdata = np.array([masssec,xsec,0.0,0.0,0.0,vysec,0.0,0.0,0.0,softsec,0.0])
        # Set cm velocity to primary
        totmass=massprim+masssec+np.sum(tgdata[:,0])
        vcmx = (0.0+np.sum(tgdata[:,0]*tgdata[:,4]))/totmass
        vcmy = (masssec*vysec+np.sum(tgdata[:,0]*tgdata[:,5]))/totmass
        print(np.sum(tgdata[:,0]*tgdata[:,5]))
        print(masssec*vysec)
        print(totmass)
        vcmz = (0.0+np.sum(tgdata[:,0]*tgdata[:,6]))/totmass
        print(vcmx,vcmy,vcmz)
        tsdata=np.array([tsdata,np.array([massprim,0.0,0.0,0.0,-vcmx,-vcmy,-vcmz,0.0,sinkprim,softprim,0.0])])
        data_header[4]=data_header[4]+2;
        data_header[0]=data_header[0]+2;
    
tip.writetipsy(tgdata,tddata,tsdata,"datafiles/" + file_string + ".00000",data_header,time)
tip.writealltipsyaux(B,"datafiles/" + file_string + ".00000")
createparamfile(file_string,sim.dxbound,sim.dybound,sim.dzbound,sim.periodic,sim.deltastep,sim.nsteps,sim.dmsolunit,sim.dkpcunit,sim.freqout,sim.adi,sim.molweight,sim.gamma,sim.grav,sim.cosmo,sim.rhoit,sim.ns,glasssetup,sim.dICdensRsmooth,sim.dICdensprofile,sim.dICdensdir,sim.dICdensR,sim.dICdensinner,sim.dICdensouter)
print('time', time)
print("created file: ",file_string)

