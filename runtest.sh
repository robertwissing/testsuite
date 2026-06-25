#!/bin/bash

export SPHEXA_DIR=/mn/stornext/d17/extragalactic/personal/robertwi/Projects/ICGEN/sims/OPTKERNELSIMS/newtests/sphexa/sphexa/build/main/src/sphexa/
export GASOLINE_DIR=/mn/stornext/d17/extragalactic/personal/robertwi/Projects/ICGEN/sims/OPTKERNELSIMS/newtests/gasoline_claude/gasoline_claude
export CHANGA_DIR=/mn/stornext/d17/extragalactic/personal/robertwi/Projects/ICGEN/sims/OPTKERNELSIMS/a/testsuite/ci/changa_ci_test/changa/
export AREPO_DIR=/mn/stornext/d17/extragalactic/personal/robertwi/Projects/ICGEN/sims/OPTKERNELSIMS/a/testsuite/codes/arepo/

#Directory for the setupfiles
current_dir=$(pwd -P)
analysisdir=$current_dir/setupfiles/
#Python (numba) glass-relaxation IC generator package root (see PLAN.md there)
glassgenpydir=$current_dir/icgenerator/glassgen_python

# AREPO (code=3): how ICconvert_tipsytoarepo.py turns the tipsy IC into an AREPO
# moving-mesh IC. "particles" (default) = a 1:1 copy of the SPH particles as mesh
# generators (most faithful, no interpolation; AREPO tessellates). For a
# regridded IC instead set "voronoi", "amr" or "uniform" (deposits the SPH fields
# onto that mesh, petkova/mass-conserving).
AREPOICMODE="particles"

#Code is selected by the 4th positional argument: 0=Gasoline 1=ChaNGa 2=SPH-EXA 3=AREPO
#The executable directory and binary name per code are set in the code-selection block after arg parsing.

#If running serial/mpi or with charmrun
#mpirun="$dir/charmrun.smp +p 63"
#mpirun=""
mpirun="mpirun "

#An extra word added at the end of the simulation output files and directory
EXTRANAME=""

#An extra word added to the glass file names (e.g. to tag the code used to relax)
GLASSNAME=""

# Format a numeric parameter for file/run names: round to 2 decimals and strip
# trailing zeros (2.0 -> 2, 0.48989 -> 0.49, 100.00 -> 100). A nonzero value
# that would round to 0.00 (e.g. ampl=1e-4, pini=1e-6) keeps its %g form so
# distinct small parameters still produce distinct names (cache keys!).
# Non-numeric arguments (preset names etc.) pass through unchanged.
fmt() {
    awk -v x="$1" 'BEGIN {
        if (x !~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) { print x; exit }
        v = x + 0
        s = sprintf("%.2f", v)
        if ((s == "0.00" || s == "-0.00") && v != 0) { printf "%g\n", v; exit }
        sub(/0+$/, "", s); sub(/\.$/, "", s)
        print s
    }'
}

# Glass relaxation: the Python (numba) glassgen.cli relaxes the pre-IC into a
# glass (multi-resolution cascade + heavy-ball momentum + the Wendland C2
# self-density correction; writes ${glassfile}_IC directly). The target density
# is read AUTOMATICALLY from the pre-IC .param (dICdens* profiles 1/2/3 =
# uniform/step/tanh). For a disk/accdisk profile (>=4) add a density table to
# GLASSMODE (e.g. "--table densitytable_xdr --direction r_sph").
#
# Extra glassgen.cli options appended to the default run, e.g.
#   "--no-progressive"  single full-resolution relax (no cascade)
#   "--icr0-stop 1e-3"  tighter glass (default 3e-3);  "--momentum 0.3"
#   "--no-self-bias"    bare W0 density self term
# Empty = the production defaults (progressive cascade, momentum 0.5,
# icr0-stop 3e-3, self-bias on). ICCONFIGPY sets the per-stage iteration cap and
# the ICR0/ICRLOOP0/rhopow/rate knobs (ICR0 = initial average move fraction of
# the box; rhopow sharpens the density match; -n is the max relaxation steps).
GLASSMODE=""
ICCONFIGPY="-n 4800 -oi 0 --icr0 1.0 --icrloop0 0.5 --rhopow 2.0 --icr0rate 1.5"

#Extra condition during runtime
#EXTRACONDALL=" -facQ 1.0 -fch 1.0 -psidecfac 1.0 -thermaldiff 0.1 -alpha 2.0 -beta 8.0 -nInitCorrIter 50 -smin 32 -smax 40000 "
EXTRACONDALL=" "


#Extra runtime conditions and output fields specific to SPH-EXA (code=2)
EXTRACONDSPHEXA=" --prop std "

fields="m,x,y,z,vx,vy,vz,temp,u,rho,h,nc"


usage() {
    echo "Usage: $0 [-i] [-r] [-a value1] ... [-h value8] <nameofsimulation> <Nsmooth> <Nres> <distrib> <code> [vm] [extra]"
    echo ""
    echo "To see the parameter choices (-a .. -h) for a specific simulation, run it by name with no further arguments, e.g.:"
    echo "    $0 accdisk"
    echo ""
    echo "Options:"
    echo "  -a, -b, -c, ...    These options depend on the simulation. To see the specific options, enter the nameofsimulation without arguments."
    echo "  -i                 Interactive."
    echo "  -r                 Do-reference: after the run, save the analysis result"
    echo "                     as the regression baseline '<runname>_reference.json'"
    echo "                     (requires setupfiles/IC_analysis_<simulation>.py)."
    echo ""
    echo "Arguments:"
    echo "  nameofsimulation   Name of the simulation to be run (e.g., shocktube)."
    echo "  Nsmooth            Number of neigbours within smoothing kernel"
    echo "  Nres               Resolution for the simulation (usually specified as number of elements in the x direction)."
    echo "  distrib            Initial conditions type:"
    echo "                     0: lattice"
    echo "                     1: random"
    echo "                     2: glass"
    echo "  code               Simulation code: 0: Gasoline  1: ChaNGa  2: SPH-EXA  3: AREPO"
    echo "  vm                 (optional) enable variable SPH particle masses in the IC generator (default 0 = equal-mass particles)."
    echo "  extra              (optional) extra flags appended to the runtime command."
    echo "If glass data file exist for a specific Nsmooth it will skipp the glass generation and run the simulation. If you choose an nSmooth that does not have a glass datafile it will do the whole glass generation again."
    echo ""
    echo "List of available simulations:"
    echo "  accdisk            Accretion disk simulation."
    echo "  alfven             Alfvén wave simulation."
    echo "  gresho             Gresho vortex simulation."
    echo "  mhdcollapse        MHD collapse simulation."
    echo "  kh                 Kelvin-Helmholtz instability simulation."
    echo "  orzag              Orszag-Tang vortex simulation (-a 1 for true 3D, Tu 2022 / Helzel et al. 2011)."
    echo "  mhdrotor           MHD rotor simulation."
    echo "  mhdloop            MHD loop advection simulation."
    echo "  shocktube          Different kinds of shock tubes."
    echo "                     (Enter 'shocktube' without arguments for more details on available options)"
    echo "  sedov              Sedov-Taylor blast wave simulation."
    echo "  evrard             Evrard collapse simulation."
    echo "  polytrope          Polytrope / Oscillating Polytrope simulation. (Requires additional graviational relaxation)"
    echo "  taylorgreen        Taylor-Green simulation."
    echo "  areablob           Hydrostatic blob(s) simulation."
    echo "  blob               Cloud crushing (wind-cloud) simulation."
    echo "  rt                 Rayleigh-Taylor instability simulation."
    echo "  mti                Magneto-thermal instability simulation."
    echo "  mhdisowave         Isolated MHD wave simulation."
    echo "  mhddivtest         Magnetic divergence test simulation."
    echo "  planet             Planet (single body) simulation."
    echo "  planetcollision    Planet collision simulation."
    echo "  square             Hydrostatic square (2D) / inner-cube (3D) test."
    echo "  cosmowave          Cosmological MHD wave (Berlok 2022) simulation."
    echo "  mhdbalsaravortex   Balsara MHD vortex simulation."
    echo "  zeldovich          Zeldovich pancake (cosmological) simulation."
    echo "  mri                Magnetorotational instability simulation."
    echo "  mhdloopshear       Magnetized advection loop with shear simulation."
    echo "  currentsheet       MHD current sheet (reconnection robustness) simulation."
    echo "  noh                Noh problem (converging strong shock) simulation."
    echo "  isentropicvortex   Isentropic (Yee) vortex (smooth convergence) simulation."
    echo "  wave               Linear (M)HD wave: sound/fast/slow/Alfven (convergence) simulation."
    echo "  coldkeplerian      Cold Keplerian disk/ring (AV / angular-momentum) simulation."
    echo ""
    exit 1
}


# Parse command line options
interactive=0
doreference=0
declare -a args=("" "" "" "" "" "" "" "")
while getopts ":ira:b:c:d:e:f:g:h:" opt; do
  case $opt in
    i)
      interactive=1
      ;;
    r)
      doreference=1
      ;;
    a|b|c|d|e|f|g|h)
      args[$(( $(printf "%d" "'$opt") - 97 ))]=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage $testsim
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage $testsim
      ;;
  esac
done
shift $((OPTIND -1))

# Check for the testsim argument
testsim=${1:-"default"}
if [ "$testsim" = "default" ]; then
    echo "Test simulation type not specified."
    usage
fi
shift

# Function for mhdcollapse setup with generic parameters
setup_mhdcollapse() {
    # Set default values if not specified
    mu=${args[0]:-10}
    rhodiff=${args[1]:-360}
    Erat=${args[2]:-0.045}

    echo "mu is the magnetic flux ratio -a "
    echo "rhodiff is the density difference between central cloud and medium -b "
    echo "Erat is the rotational energy ratio -c "

    echo "Final settings for MHD Collapse:"
    echo "mu: $mu, rhodiff: $rhodiff, Erat: $Erat"

    # Setup based on the potentially updated values
    EXTRACONDIC="$mu $rhodiff $Erat"
    EXTRACONDRUN=" "
    prefile_extra="rd$(fmt "${rhodiff}")"
    glassfile_extra="rd$(fmt "${rhodiff}")${GLASSNAME}"
    runname_extra="mu$(fmt "${mu}")rd$(fmt "${rhodiff}")Erat$(fmt "${Erat}")"
}

# Function for accdisk setup with generic parameters
setup_accdisk() {
    # Set default values if not specified
    hr=${args[0]:-0.1}
    massrat=${args[1]:-1E-3}
    onlyoneBH=${args[2]:-0}
    q=${args[3]:-0.5}
    rorb=${args[4]:-1.0}
    rdisk=${args[5]:-10.0}
    rinner=${args[6]:-0.5}

    echo "Using settings for accdisk:"
    echo "Disk height ratio (hr) -a : $hr"
    echo "Mass ratio (massrat) -b : $massrat"
    echo "Binary black hole presence (onlyoneBH) -c : $onlyoneBH"
    echo "Mass distribution parameter (q) -d : $q"
    echo "Orbital radius (rorb) -e : $rorb"
    echo "Disk radius (rdisk) -f : $rdisk"
    echo "Inner disk radius (rinner) -g : $rinner"

    # Setup based on the potentially updated values
    EXTRACONDIC="$hr $massrat $onlyoneBH $q $rorb $rdisk $rinner "
    EXTRACONDRUN=" "
    prefile_extra="hr$(fmt "${hr}")mr$(fmt "${massrat}")bh$(fmt "${onlyoneBH}")q$(fmt "${q}")ro$(fmt "${rorb}")rd$(fmt "${rdisk}")ri$(fmt "${rinner}")"
    glassfile_extra="hr$(fmt "${hr}")mr$(fmt "${massrat}")bh$(fmt "${onlyoneBH}")q$(fmt "${q}")ro$(fmt "${rorb}")rd$(fmt "${rdisk}")ri$(fmt "${rinner}")${GLASSNAME}"
    runname_extra="hr$(fmt "${hr}")mr$(fmt "${massrat}")bh$(fmt "${onlyoneBH}")q$(fmt "${q}")ro$(fmt "${rorb}")rd$(fmt "${rdisk}")ri$(fmt "${rinner}")"
}

# Function for alfven setup with generic parameters
setup_alfven() {
    # Set default values if not specified
    rotate=${args[0]:-1}

    echo "Using settings for circularized alfvenwave:"
    echo "rotate wave (rotate) -a : $rotate"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rotate "
    EXTRACONDRUN=" "
    prefile_extra=""
    glassfile_extra="${GLASSNAME}"
    runname_extra="rot$(fmt "${rotate}")"
}

# Function for Shock-tube setup with generic parameters
setup_shocktube() {
    shock=${args[0]:-1.0}
    # Array mapping shock indices to shock names (must cover every `choice` in
    # IC_setup_shocktube.choose_shock, else the index check below rejects them)
    shocks=(
        "Sod shock"
        "Ryu 1a"
        "Ryu 1b"
        "Ryu 2a"
        "Ryu 2b"
        "Brio-Wu Ryu 5a"
        "C-shock"
        "Steady shock"
        "Toth"
        "MHD discontinuity"
        "strong blast"
        "fast slow shock"
        "isothermal B98"
        "Mach 25 DW94"
        "rarefaction"
    )

    # Get the shock name based on the shock index
    shock_index=$(echo "$shock" | awk '{printf "%.0f", $0}')  # Convert shock to integer index
    if (( shock_index < 1 || shock_index > ${#shocks[@]} )); then
        echo "Invalid shock index: $shock"
        return 1
    fi
    shockname=${shocks[$((shock_index - 1))]}

    # Display the settings being used
    echo "Using settings for Shocktube:"
    echo "Different shocks (argument -a): $shock"
    echo "
    # shocks(1) = 'Sod shock'
    # shocks(2) = 'Ryu 1a'
    # shocks(3) = 'Ryu 1b'
    # shocks(4) = 'Ryu 2a'
    # shocks(5) = 'Ryu 2b'
    # shocks(6) = 'Brio-Wu (Ryu 5a)'
    # shocks(7) = 'C-shock'
    # shocks(8) = 'Steady shock'
    "
    echo "Selected shock: $shockname"

    # Setup based on the potentially updated values
    EXTRACONDIC="$shock "
    EXTRACONDRUN=""
    prefile_extra="${shockname// /_}"
    glassfile_extra="${shockname// /_}${GLASSNAME}"
    runname_extra="${shockname// /_}"

    echo "prefile_extra: $prefile_extra"
    echo "glassfile_extra: $glassfile_extra"
    echo "runname_extra: $runname_extra"

}

# Function for Orzag-Tang Vortex setup with generic parameters
setup_orzag() {
    dim3=${args[0]:-0}
    # Set default values if not specified
    echo "Using settings for Orzag-Tang Vortex:"
    echo "3D vortex -a (0: 2D slab, 1: true 3D Tu 2022 / Helzel et al. 2011) : $dim3"

    # Setup based on the potentially updated values
    EXTRACONDIC=" $dim3 "
    EXTRACONDRUN=""
    prefile_extra=""
    glassfile_extra="${GLASSNAME}"
    runname_extra=""
    if [ "$dim3" != "0" ]; then
        prefile_extra="_3d"
        glassfile_extra="_3d${GLASSNAME}"
        runname_extra="_3d"
    fi
}

# Function for Magnetized Rotor setup with generic parameters
setup_mhdrotor() {
    rhodisk=${args[0]:-10.0}
    # Set default values if not specified
    echo "Using settings for Magnetized Rotor:"
    echo "density disk -a : $rhodisk"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodisk "
    EXTRACONDRUN=""
    prefile_extra="rhodiff$(fmt "${rhodisk}")"
    glassfile_extra="rhodiff$(fmt "${rhodisk}")${GLASSNAME}"
    runname_extra="rhodiff$(fmt "${rhodisk}")"
}

# Function for Magnetized advection loop setup with generic parameters
setup_mhdloop() {
    rhodisk=${args[0]:-1.0}
    # Set default values if not specified
    echo "Using settings for Magnetized Advection Loop:"
    echo "density disk -a : $rhodisk"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodisk "
    EXTRACONDRUN=""
    prefile_extra="rhodiff$(fmt "${rhodisk}")"
    glassfile_extra="rhodiff$(fmt "${rhodisk}")${GLASSNAME}"
    runname_extra="rhodiff$(fmt "${rhodisk}")"
}

# Function for gresho setup with generic parameters
setup_gresho() {
    # Set default values if not specified
    result=$(echo "scale=5; sqrt(3/25)" | bc -l)
    echo "$result"
    mach=${args[0]:-$result}
    
    echo "Using settings for gresho:"
    echo "mach number of flow -a : $mach"
    
    # Setup based on the potentially updated values
    EXTRACONDIC="$mach "
    EXTRACONDRUN=""
    prefile_extra=""
    glassfile_extra="${GLASSNAME}"
    runname_extra="mach$(fmt "${mach}")"
}

# Function for kelvin-helmholtz setup with generic parameters
setup_kh() {
    # Set default values if not specified
    result=$(echo "scale=5; sqrt(6)/5" | bc -l)
    echo "$result"
    rhodiff=${args[0]:-2.0}
    mach=${args[1]:-$result}
    smooth=${args[2]:-1}
    B0=${args[3]:-0.0}
    Bdir=${args[4]:-1}

    echo "Using settings for kelvin-helmholtz:"
    echo "rho diff -a : $rhodiff"
    echo "mach number of flow -b : $mach"
    echo "smooth density disc -c : $smooth"
    echo "initial magnetic field strength -d : $B0"
    echo "Direction of field x=1 y=2 z=3 -e : $Bdir"
    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodiff $mach $smooth $B0 $Bdir "
    EXTRACONDRUN=""
    prefile_extra="rhodiff$(fmt "${rhodiff}")smooth$(fmt "${smooth}")"
    glassfile_extra="rhodiff$(fmt "${rhodiff}")smooth$(fmt "${smooth}")${GLASSNAME}"
    runname_extra="rhodiff$(fmt "${rhodiff}")mach$(fmt "${mach}")smooth$(fmt "${smooth}")B$(fmt "${B0}")Bdir$(fmt "${Bdir}")"
}


# Function for Sedov blast wave + mhd setup with generic parameters
setup_sedov() {
    select=${args[0]:-1.0}
    betain=${args[1]:-0.0}
    ublast=${args[2]:-10.0}
    Rin=${args[3]:-0.125}
    Pin=${args[4]:-100.0}
    Pout=${args[5]:-1.0}
    # Set default values if not specified
    echo "Using settings for Sedov blast wave + mhd:"
    echo "    # For -a Select 0 we define and inner and outer pressure and an inner radius of blast. And corresponding inner plasmabeta for mhd.
    #Select 1 and 2 we choose an blast energy and give it to a few particles in centre.
    #In select 1 we have different magnetic strength in the blast and outside blast(need divergence cleaning) and select 2 we have a uniform field."
    echo "Using ublast and setting to few particles(=1,2) or using defined inner pressure and radius(=0) -a : $select"
    echo "Plasma beta of initial blast region -b : $betain"
    echo "Energy of blast(select=1) -c : $ublast"
    echo "Radius of inner blast(select=0) -d : $Rin"
    echo "Pressure inside(select=0) -e : $Pin"
    echo "Pressure outside(select=0) -f : $Pout"

    # Setup based on the potentially updated values
    EXTRACONDIC="$select $betain $ublast $Rin $Pin $Pout"
    EXTRACONDRUN=""
    prefile_extra="select$(fmt "${select}")beta$(fmt "${betain}")ublast$(fmt "${ublast}")"
    glassfile_extra="${GLASSNAME}"
    # Rin/Pin/Pout only shape the IC for select=0 (pressure-jump sphere); for
    # select 1/2 the setup overrides Rin (2*mean(h)) and recomputes Pin/Pout
    # from ublast, so they are kept out of the name there.
    runname_extra="select$(fmt "${select}")beta$(fmt "${betain}")ublast$(fmt "${ublast}")"
    if [ "$(fmt "${select}")" == "0" ]; then
        runname_extra="${runname_extra}Rin$(fmt "${Rin}")Pin$(fmt "${Pin}")Pout$(fmt "${Pout}")"
    fi
}

# Function for alfven setup with generic parameters
setup_evrard() {
    # Set default values if not specifie
    echo "Using settings for evrard collapse:"
    # Setup based on the potentially updated values
    EXTRACONDIC=" "
    EXTRACONDRUN=" "
    prefile_extra=""
    glassfile_extra="${GLASSNAME}"
    runname_extra=""
}

# Function for alfven setup with generic parameters
setup_polytrope() {
        # Set default values if not specifie
        echo "Using settings for polytrope:"
        # Setup based on the potentially updated values
        EXTRACONDIC=" "
        EXTRACONDRUN=" "
        prefile_extra=""
        glassfile_extra="${GLASSNAME}"
        runname_extra=""
}

# Function for static blobs (force test) setup with generic parameters
setup_areablob() {
    # Set default values if not specified
    rhodiff=${args[0]:-1000}
    multcloud=${args[1]:-1}

    echo "mu is the magnetic flux ratio -a "
    echo "rhodiff is the density difference between central cloud and medium -b "
    echo "Erat is the rotational energy ratio -c "

    echo "Final settings for MHD Collapse:"
    echo "mu: $multcloud, rhodiff: $rhodiff"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodiff $multcloud"
    EXTRACONDRUN="  "
    prefile_extra="mult$(fmt "${multcloud}")rd$(fmt "${rhodiff}")${GLASSNAME}"
    glassfile_extra="mult$(fmt "${multcloud}")rd$(fmt "${rhodiff}")${GLASSNAME}"
    runname_extra="mult$(fmt "${multcloud}")rd$(fmt "${rhodiff}")"
}

# Function for Taylor-green setup with generic parameters
setup_taylorgreen() {
    etavisc=${args[0]:-0.1}
    xivisc=${args[1]:-0.0}
    # Set default values if not specified
    echo "Using settings for Taylorgreen Vortex:"
    echo "eta visc -a : $etavisc"
    echo "xi visc -b : $xivisc"

    # Setup based on the potentially updated values
    EXTRACONDIC=" "
    EXTRACONDRUN=" -PhysViscEta $etavisc -PhysViscXi $xivisc "
    prefile_extra=""
    glassfile_extra="${GLASSNAME}"
    runname_extra="eta$(fmt "${etavisc}")xi$(fmt "${xivisc}")"
}

# Function for cloud crushing (wind-cloud) setup with generic parameters
setup_blob() {
    # Set default values if not specified
    rhodiff=${args[0]:-10}
    beta=${args[1]:-0}
    mach=${args[2]:-2.7}
    inflow=${args[3]:-0}
    cool=${args[4]:-0}
    windref=${args[5]:-0}
    winddens=${args[6]:-26}
    tcloud=${args[7]:-4000}
    echo "Using settings for cloud crushing:"
    echo "rho diff -a : $rhodiff"
    echo "plasma beta -b : $beta"
    echo "mach number -c : $mach"
    echo "inflow conditions -d : $inflow"
    echo "with cooling -e : $cool"
    echo "reference frame of wind  -f : $windref"
    echo "density of wind g/cm  -g : $winddens"
    echo "temperature cloud  -h : $tcloud"
    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodiff $beta $mach $inflow $cool $windref $winddens $tcloud "
    EXTRACONDRUN=" "
    prefile_extra="rhodiff$(fmt "${rhodiff}")beta$(fmt "${beta}")mach$(fmt "${mach}")inflow$(fmt "${inflow}")cool$(fmt "${cool}")wr$(fmt "${windref}")wd$(fmt "${winddens}")tc$(fmt "${tcloud}")"
    glassfile_extra="rhodiff$(fmt "${rhodiff}")${GLASSNAME}"
    runname_extra="rhodiff$(fmt "${rhodiff}")beta$(fmt "${beta}")mach$(fmt "${mach}")inflow$(fmt "${inflow}")cool$(fmt "${cool}")wr$(fmt "${windref}")wd$(fmt "${winddens}")tc$(fmt "${tcloud}")"
}

# Function for Rayleigh-Taylor setup with generic parameters
setup_rt() {
    # Set default values if not specified
    rhodiff=${args[0]:-2.0}
    case=${args[1]:-1}
    nlowdens=${args[2]:-2}
    smooth=${args[3]:-0}
    beta=${args[4]:-0.0}
    bdir=${args[5]:-1}
    echo "Using settings for rayleigh-taylor:"
    echo "rho diff -a : $rhodiff"
    echo "Case 1 for symmetric and 2 for unsymmetric -b : $case"
    echo "Number of low density segments (same size as high density) -c : $nlowdens"
    echo "Smooth -d : $smooth"
    echo "Plasma beta at the interface (0 = hydro) -e : $beta"
    echo "B direction 1=x horizontal (default) 2=y 3=z vertical -f : $bdir"
    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodiff $case $nlowdens $smooth $beta $bdir"
    EXTRACONDRUN=" "
    prefile_extra="rhodiff$(fmt "${rhodiff}")smooth$(fmt "${smooth}")nlow$(fmt "${nlowdens}")c$(fmt "${case}")beta$(fmt "${beta}")bdir$(fmt "${bdir}")"
    glassfile_extra="rhodiff$(fmt "${rhodiff}")smooth$(fmt "${smooth}")nlow$(fmt "${nlowdens}")c$(fmt "${case}")beta$(fmt "${beta}")bdir$(fmt "${bdir}")${GLASSNAME}"
    runname_extra="rhodiff$(fmt "${rhodiff}")smooth$(fmt "${smooth}")nlow$(fmt "${nlowdens}")c$(fmt "${case}")beta$(fmt "${beta}")bdir$(fmt "${bdir}")"
}

# Function for magneto-thermal instability setup
setup_mti() {
    # Set default values if not specified
    echo "Using settings for MTI:"
    # Setup based on the potentially updated values
    EXTRACONDIC=" "
    EXTRACONDRUN=" "
    prefile_extra=""
    glassfile_extra="${GLASSNAME}"
    runname_extra=""
}

# Function for isolated MHD wave setup with generic parameters
setup_mhdisowave() {
    # Set default values if not specified
    beta=${args[0]:-0.1}
    echo "Using settings for isolated wave:"
    echo "Plasma beta -a : $beta"
    # Setup based on the potentially updated values
    EXTRACONDIC="$beta "
    EXTRACONDRUN=" "
    prefile_extra="$(fmt "${beta}")"
    glassfile_extra="$(fmt "${beta}")${GLASSNAME}"
    runname_extra="beta$(fmt "${beta}")"
}

# Function for magnetic divergence test setup with generic parameters
setup_mhddivtest() {
    # Set default values if not specified
    rhodiff=${args[0]:-1}
    smooth=${args[1]:-1}
    square=${args[2]:-1}
    cosmo=${args[3]:-0}
    echo "Using settings for divergence test:"
    echo "Density contrast -a : $rhodiff"
    echo "Smooth density -b : $smooth"
    echo "Square test -c : $square"
    echo "Cosmo test -d : $cosmo"
    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodiff $smooth $square $cosmo  "
    EXTRACONDRUN=" "
    prefile_extra="rd$(fmt "${rhodiff}")s$(fmt "${smooth}")sq$(fmt "${square}")"
    glassfile_extra="rd$(fmt "${rhodiff}")s$(fmt "${smooth}")sq$(fmt "${square}")${GLASSNAME}"
    runname_extra="rd$(fmt "${rhodiff}")s$(fmt "${smooth}")sq$(fmt "${square}")cos$(fmt "${cosmo}")"
}

# Function for planet (single body) setup with generic parameters
setup_planet() {
    # Set default values if not specified
    planet=${args[0]:-"Earth_PES095_2layer"}
    echo "Using settings for planet:"
    # Setup based on the potentially updated values
    EXTRACONDIC="$planet"
    EXTRACONDRUN=" "
    prefile_extra="$planet"
    glassfile_extra="$planet${GLASSNAME}"
    runname_extra="$planet"
}

# Function for planet collision setup with generic parameters
setup_planetcollision() {
    # Set default values if not specified
    target=${args[0]:-"None"}
    impactor=${args[1]:-"None"}
    impangle=${args[2]:-0.73}
    vinf=${args[3]:-0.0}
    impgeom=${args[4]:-2}
    echo "Using settings for planet collision:"
    # Setup based on the potentially updated values
    EXTRACONDIC="$target $impactor $impangle $vinf $impgeom"
    EXTRACONDRUN=" "
    prefile_extra="t${target}i${impactor}b$(fmt "${impangle}")g$(fmt "${impgeom}")"
    glassfile_extra="t${target}i${impactor}b$(fmt "${impangle}")g$(fmt "${impgeom}")${GLASSNAME}"
    runname_extra="t${target}i${impactor}b$(fmt "${impangle}")v$(fmt "${vinf}")g$(fmt "${impgeom}")"
}

# Function for hydrostatic-square test
# is3d=0: 2D thin-slab variant (Saitoh-Makino 2013; Hopkins 2013, 2015; Hu et al. 2014)
# is3d=1: 3D inner-cube variant (Sandnes et al. 2025)
setup_square() {
    gamma_default=$(echo "scale=12; 5/3" | bc -l)
    rhodiff=${args[0]:-4.0}
    przero=${args[1]:-2.5}
    gamma=${args[2]:-$gamma_default}
    is3d=${args[3]:-0}

    echo "Using settings for hydrostatic-square test:"
    echo "density contrast chi -a : $rhodiff"
    echo "uniform pressure P -b : $przero"
    echo "gamma -c : $gamma"
    echo "3D inner-cube (Sandnes 2025) -d : $is3d"

    # Order matches IC_setup_square.create(... rhodiff, przero, gamma, is3d) signature
    EXTRACONDIC="$rhodiff $przero $gamma $is3d "
    EXTRACONDRUN=""
    prefile_extra="chi$(fmt "${rhodiff}")P$(fmt "${przero}")d3$(fmt "${is3d}")"
    glassfile_extra="chi$(fmt "${rhodiff}")d3$(fmt "${is3d}")${GLASSNAME}"
    runname_extra="chi$(fmt "${rhodiff}")P$(fmt "${przero}")g$(fmt "${gamma}")d3$(fmt "${is3d}")"
}

# Function for cosmological MHD wave (Berlok 2022) setup with generic parameters
setup_cosmowave() {
    # Set default values if not specified
    case=${args[0]:-3}
    factor_A=${args[1]:-1}
    factor_G=${args[2]:-0}
    echo "Using settings for cosmowave (Berlok 2022):"
    echo "test case -a : $case  (0:standing Alfven 1:traveling Alfven 2:standing compr. 3:traveling compr. gamma=4/3 4:div-clean EdS)"
    echo "factor_A (omega_A = factor_A*pi; 0=MHD off, 1=paper default) -b : $factor_A"
    echo "factor_G (omega_G = factor_G*pi/2; 0=self-grav off, 1=paper default) -c : $factor_G"
    # Order matches IC_setup_cosmowave.create(... case, factor_A, factor_G) signature
    EXTRACONDIC="$case $factor_A $factor_G "
    EXTRACONDRUN=" "
    prefile_extra="case$(fmt "${case}")A$(fmt "${factor_A}")G$(fmt "${factor_G}")"
    glassfile_extra="${GLASSNAME}"
    runname_extra="case$(fmt "${case}")A$(fmt "${factor_A}")G$(fmt "${factor_G}")"
}

# Function for Balsara MHD vortex setup with generic parameters
setup_mhdbalsaravortex() {
    # Set default values if not specified
    u=${args[0]:-1.0}
    k=${args[1]:-1.0}
    echo "Using settings for Balsara MHD vortex:"
    echo "magnetic field amplitude -a : $u"
    echo "velocity perturbation amplitude -b : $k"
    # Order matches IC_setup_mhdbalsaravortex.create(... u, k) signature
    EXTRACONDIC="$u $k "
    EXTRACONDRUN=""
    prefile_extra="u$(fmt "${u}")k$(fmt "${k}")"
    glassfile_extra="${GLASSNAME}"
    runname_extra="u$(fmt "${u}")k$(fmt "${k}")"
}

# Function for Zeldovich pancake (cosmological) setup with generic parameters
setup_zeldovich() {
    # Set default values if not specified
    Lboxinkpc=${args[0]:-1000.0}
    zi=${args[1]:-100}
    echo "Using settings for Zeldovich pancake:"
    echo "box size in kpc (length unit) -a : $Lboxinkpc"
    echo "initial redshift zi (must be > collapse redshift zc=2) -b : $zi"
    # Order matches IC_setup_zeldovich.create(... Lboxinkpc, zi) signature
    EXTRACONDIC="$Lboxinkpc $zi "
    EXTRACONDRUN=" "
    prefile_extra="L$(fmt "${Lboxinkpc}")zi$(fmt "${zi}")"
    glassfile_extra="${GLASSNAME}"
    runname_extra="L$(fmt "${Lboxinkpc}")zi$(fmt "${zi}")"
}

# Function for magnetorotational instability (MRI) setup with generic parameters
setup_mri() {
    # Set default values if not specified
    case=${args[0]:-2}
    echo "Using settings for MRI:"
    echo "configuration case -a : $case  (1:thin unstrat nonetflux 2:regular unstrat nonetflux 3:tall unstrat nonetflux 4:stratified netflux 5:unstrat netflux 6:regular unstrat nonetflux 7:- 8:unstrat netflux tall)"
    # Order matches IC_setup_mri.create(... case) signature
    EXTRACONDIC="$case "
    EXTRACONDRUN=""
    prefile_extra="case$(fmt "${case}")"
    glassfile_extra="case$(fmt "${case}")${GLASSNAME}"
    runname_extra="case$(fmt "${case}")"
}

# Function for Magnetized advection loop with shear setup with generic parameters
setup_mhdloopshear() {
    rhodiff=${args[0]:-1.0}
    # Set default values if not specified
    echo "Using settings for Magnetized Advection Loop with shear:"
    echo "density contrast -a : $rhodiff"
    # Order matches IC_setup_mhdloopshear.create(... rhodiff) signature
    EXTRACONDIC="$rhodiff "
    EXTRACONDRUN=""
    prefile_extra="rhodiff$(fmt "${rhodiff}")"
    glassfile_extra="rhodiff$(fmt "${rhodiff}")${GLASSNAME}"
    runname_extra="rhodiff$(fmt "${rhodiff}")"
}

# Function for MHD current sheet (Gardiner & Stone 2005) setup with generic parameters
setup_currentsheet() {
    # Set default values if not specified
    beta=${args[0]:-0.1}
    v0=${args[1]:-0.1}
    echo "Using settings for MHD current sheet:"
    echo "plasma beta (sets gas pressure = beta*B0^2/2) -a : $beta"
    echo "velocity perturbation amplitude v0 -b : $v0"
    # Order matches IC_setup_currentsheet.create(... beta, v0) signature
    EXTRACONDIC="$beta $v0 "
    EXTRACONDRUN=""
    prefile_extra="beta$(fmt "${beta}")v$(fmt "${v0}")"
    glassfile_extra="${GLASSNAME}"
    runname_extra="beta$(fmt "${beta}")v$(fmt "${v0}")"
}

# Function for the Noh problem (converging strong shock) setup with generic parameters
setup_noh() {
    # Set default values if not specified
    is3d=${args[0]:-1}
    pini=${args[1]:-1e-6}
    v0=${args[2]:-1.0}
    echo "Using settings for Noh problem:"
    echo "3D spherical (1) / 2D cylindrical (0) -a : $is3d"
    echo "initial (tiny) pressure -b : $pini"
    echo "radial inflow speed v0 -c : $v0"
    # Order matches IC_setup_noh.create(... is3d, pini, v0) signature
    EXTRACONDIC="$is3d $pini $v0 "
    EXTRACONDRUN=""
    prefile_extra="d3$(fmt "${is3d}")"
    glassfile_extra="d3$(fmt "${is3d}")${GLASSNAME}"
    runname_extra="d3$(fmt "${is3d}")p$(fmt "${pini}")v$(fmt "${v0}")"
}

# Function for the isentropic (Yee) vortex setup with generic parameters
setup_isentropicvortex() {
    # Set default values if not specified
    beta_v=${args[0]:-5.0}
    vadv=${args[1]:-1.0}
    echo "Using settings for isentropic (Yee) vortex:"
    echo "vortex strength beta_v -a : $beta_v"
    echo "mean-flow (advection) velocity vadv -b : $vadv"
    # Order matches IC_setup_isentropicvortex.create(... beta_v, vadv) signature
    EXTRACONDIC="$beta_v $vadv "
    EXTRACONDRUN=""
    prefile_extra="b$(fmt "${beta_v}")"
    glassfile_extra="b$(fmt "${beta_v}")${GLASSNAME}"
    runname_extra="b$(fmt "${beta_v}")vadv$(fmt "${vadv}")"
}

# Function for linear (M)HD wave setup with generic parameters
setup_wave() {
    # Set default values if not specified
    wavetype=${args[0]:-0}
    ampl=${args[1]:-1e-4}
    echo "Using settings for linear (M)HD wave:"
    echo "wave type -a : $wavetype  (0:sound 1:fast magnetosonic 2:slow magnetosonic 3:Alfven)"
    echo "eigenmode amplitude -b : $ampl"
    # Order matches IC_setup_wave.create(... wavetype, ampl) signature
    EXTRACONDIC="$wavetype $ampl "
    EXTRACONDRUN=""
    prefile_extra="wt$(fmt "${wavetype}")a$(fmt "${ampl}")"
    glassfile_extra="wt$(fmt "${wavetype}")a$(fmt "${ampl}")${GLASSNAME}"
    runname_extra="wt$(fmt "${wavetype}")a$(fmt "${ampl}")"
}

# Function for the cold Keplerian disk/ring setup with generic parameters
setup_coldkeplerian() {
    # Set default values if not specified
    hr=${args[0]:-0.05}
    rinner=${args[1]:-1.0}
    rdisk=${args[2]:-2.0}
    q=${args[3]:-0.5}
    echo "Using settings for cold Keplerian disk/ring:"
    echo "aspect ratio hr (=1/Mach; small=cold) -a : $hr"
    echo "inner ring radius -b : $rinner"
    echo "outer ring radius -c : $rdisk"
    echo "rotation power-law q (v=r^-q) -d : $q"
    # Order matches IC_setup_coldkeplerian.create(... hr, rinner, rdisk, q) signature
    EXTRACONDIC="$hr $rinner $rdisk $q "
    EXTRACONDRUN=""
    prefile_extra="hr$(fmt "${hr}")ri$(fmt "${rinner}")ro$(fmt "${rdisk}")q$(fmt "${q}")"
    glassfile_extra="hr$(fmt "${hr}")ri$(fmt "${rinner}")ro$(fmt "${rdisk}")q$(fmt "${q}")${GLASSNAME}"
    runname_extra="hr$(fmt "${hr}")ri$(fmt "${rinner}")ro$(fmt "${rdisk}")q$(fmt "${q}")"
}


# Main execution block
case $testsim in
    accdisk)
        setup_accdisk
        ;;
    alfven)
        setup_alfven
        ;;
    gresho)
        setup_gresho
        ;;
    mhdcollapse)
        setup_mhdcollapse
        ;;
    kh)
        setup_kh
        ;;
    orzag)
        setup_orzag
        ;;
    mhdrotor)
        setup_mhdrotor
        ;;
    mhdloop)
        setup_mhdloop
        ;;
    shocktube)
        setup_shocktube
        ;;
    sedov)
        setup_sedov
	;;
    evrard)
        setup_evrard
        ;;
    polytrope)
        setup_polytrope
        ;;
    taylorgreen)
        setup_taylorgreen
        ;;
    areablob)
        setup_areablob
        ;;
    blob)
        setup_blob
        ;;
    rt)
        setup_rt
        ;;
    mti)
        setup_mti
        ;;
    mhdisowave)
        setup_mhdisowave
        ;;
    mhddivtest)
        setup_mhddivtest
        ;;
    planet)
        setup_planet
        ;;
    planetcollision)
        setup_planetcollision
        ;;
    square)
        setup_square
        ;;
    cosmowave)
        setup_cosmowave
        ;;
    mhdbalsaravortex)
        setup_mhdbalsaravortex
        ;;
    zeldovich)
        setup_zeldovich
        ;;
    mri)
        setup_mri
        ;;
    mhdloopshear)
        setup_mhdloopshear
        ;;
    currentsheet)
        setup_currentsheet
        ;;
    noh)
        setup_noh
        ;;
    isentropicvortex)
        setup_isentropicvortex
        ;;
    wave)
        setup_wave
        ;;
    coldkeplerian)
        setup_coldkeplerian
        ;;
    *)
        echo "Error: Unknown test simulation '$testsim'."
	usage
        exit 1
        ;;
esac

mkdir test_cases
mkdir test_cases/$testsim
cd test_cases/$testsim

CURRENT_DIR=$(pwd -P)

if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <Nsmooth> <Nres> <0: lattice 1: rand 2: glass> <code 0:Gasoline 1:ChaNGa 2:SPH-EXA 3:AREPO> [vm] [extra]"
    exit 1
fi
echo "$0 $1 $2 $3 $4 $5 $6"

# Append the optional extra positional flags ($6) to the runtime conditions
echo "$EXTRACONDALL"
EXTRACONDRUN=" $EXTRACONDRUN $EXTRACONDALL $6 "

# Varied-mass argument ($5): changes the particle distribution itself, so it
# must key the pre-IC file, the glass AND the run name (empty tag when off).
if [ -z "$5" ] || [ "$5" -le 1 ]; then
            vm=0.0
    else
                vm=$5
fi
if [ "$vm" != "0.0" ]; then
    vmtag="vm$(fmt "${vm}")"
else
    vmtag=""
fi
prefile_extra="${prefile_extra}${vmtag}"
glassfile_extra="${glassfile_extra}${vmtag}"
runname_extra="${runname_extra}${vmtag}"

# Set distribution and d1 based on the third command-line argument.
# prefile = the file IC_createsetup.py writes directly into datafiles/:
#  - lattice/rand: that file IS the final IC (velocities/B field included), so it
#    must be keyed by the FULL parameter set (runname_extra) -- otherwise a rerun
#    with e.g. a different B0 silently reuses the cached IC of the old run.
#  - glass: it is only the pre-relaxation density distribution, so it is keyed by
#    the density-structure params (prefile_extra) and shared across runs that
#    differ only in velocities/B; the final IC (built from the relaxed glass) is
#    keyed by runname_extra instead.
if [ "$3" -eq 0 ]; then
echo "RUNNING LATTICE MODE"
    distri="lattice"
    distri2="lattice"
    d1=0
    simfilename="$testsim$2_$distri${runname_extra}"
    prefile="$simfilename"
elif [ "$3" -eq 1 ]; then
echo "RUNNING RANDOM MODE"
    distri="rand"
    distri2="rand"
    d1=1
    simfilename="$testsim$2_$distri${runname_extra}"
    prefile="$simfilename"
else
echo "RUNNING GLASS MODE"
    distri="rand"
    distri2="glass"
    d1=1
    simfilename="$testsim$2_glassN$1${runname_extra}-I"
    prefile="$testsim$2_${distri}${prefile_extra}"
    echo "$simfilename"
fi

# Code name and binary/launcher selection based on the 4th argument
if [ "$4" -eq 0 ]; then
echo "RUNNING GASOLINE"
    codename="GASOLINE"
    dir=$GASOLINE_DIR
    coderun=gasolinemhd
elif [ "$4" -eq 1 ]; then
echo "RUNNING CHANGA"
    codename="CHANGA"
    dir=$CHANGA_DIR
    coderun=ChaNGa.smp
    runlauncher="$CHANGA_DIR/charmrun +p 4"
elif [ "$4" -eq 2 ]; then
echo "RUNNING SPH-EXA"
    codename="SPHEXA"
    dir=$SPHEXA_DIR
    coderun=sphexa
elif [ "$4" -eq 3 ]; then
echo "RUNNING AREPO"
    codename="AREPO"
    dir=$AREPO_DIR
    coderun=Arepo
fi

glassfile_prefix="$testsim$2_glassN$1${glassfile_extra}"
ICfile="${simfilename}"
runname="$testsim$2_N$1${distri2}${runname_extra}${EXTRANAME}${codename}"
directory="initruns"
glassfile="${directory}/${glassfile_prefix}"

mkdir datafiles
mkdir initruns



    # Check if the file already exists
    if [ -f "datafiles/$prefile.00000" ]; then
        echo "File $prefile.00000 already exists. Skipping creation."
    else
    # File doesn't exist, create it
    preIC="python ${analysisdir}/IC_createsetup.py $testsim $2 $vm $d1 $prefile $EXTRACONDIC"
    echo $preIC
    eval $preIC
    # Check if the command was successful
    if [ $? -ne 0 ]; then
    echo "Error: Command failed to execute successfully."
    exit 1
    fi
    
    fi

##Glass relaxation    
    if [ "$3" -ne 0 ] && [ "$3" -ne 1 ]; then
    if [ -f "${glassfile}_IC" -o -f "datafiles/$ICfile.00000" ]; then
        echo "File ${glassfile}_IC already exists. Skipping creatin.."
    else
	# Python (numba) glass relaxation - glassgen.cli writes ${glassfile}_IC
	# directly; target density is read from the .param.
	relaxglass=(env PYTHONPATH="${glassgenpydir}" python -m glassgen.cli -o "$glassfile" -I datafiles/$prefile.00000 -s $1 $ICCONFIGPY --param datafiles/$prefile.param $GLASSMODE)
	echo "${relaxglass[@]}"
	"${relaxglass[@]}"
	if [ $? -ne 0 ]; then
	    echo "Error: Command failed to execute successfully."
	    exit 1
	fi
    fi

    if [ -f "datafiles/$ICfile.00000" ]; then
        echo "File $ICfile.00000 already exists. Skipping creation."
    else
    ICfromglass="python ${analysisdir}/IC_createsetup.py $testsim $2 $vm "${glassfile}_IC" ${ICfile} $EXTRACONDIC "
    echo $ICfromglass
    eval $ICfromglass

    # Check if the command was successful
    if [ $? -ne 0 ]; then
    echo "Error: Command failed to execute successfully."
    exit 1
    fi
    
    fi
    fi
    mkdir $runname
    echo "Code selection: $4"
    if [ "$4" == 0 ] || [ "$4" == 1 ]; then
    if [ "$4" == 1 ]; then sed -i '/^dCloudDensity/d' datafiles/$ICfile.param; fi
        runline=( ${runlauncher:-$mpirun} $dir/${coderun} -s $1 -o $runname/$runname -I datafiles/$ICfile.00000 $EXTRACONDRUN datafiles/$ICfile.param)
    fi

    if [ "$4" == 2 ]; then
	# SPH-EXA
	paramfile="datafiles/$ICfile.param"

	get_param () {
	    awk -F= -v key="$1" '
		$1 ~ "^[[:space:]]*" key "[[:space:]]*$" {
		    gsub(/[[:space:]]*/, "", $2)
		    print $2
		}
	    ' "$paramfile"
	}

	iOutInterval=$(get_param iOutInterval)
	dDelta=$(get_param dDelta)
	nSteps=$(get_param nSteps)

	# Sanity check
	if [[ -z "$iOutInterval" || -z "$dDelta" || -z "$nSteps" ]]; then
	    echo "Error: missing required parameters in $paramfile"
	    exit 1
	fi

	# Compute physical times
	outinterval=$(awk -v io="$iOutInterval" -v dd="$dDelta" 'BEGIN { printf "%.8f", io * dd }')
	endtime=$(awk -v ns="$nSteps" -v dd="$dDelta" 'BEGIN { printf "%.8f", ns * dd }')

	# load special modules for sph-exa
	module purge
	#module load nvhpc/nvhpc-hpcx-cuda12/24.7
	module load gcc/13.3.1
	module load openmpi/gcc13/5.0.9

	# convert the tipsy IC into the SPH-EXA hdf5 format
	echo "python ${analysisdir}/ICconvert_tipsytosphexa.py datafiles/$ICfile.00000 datafiles/$ICfile "
	python ${analysisdir}/ICconvert_tipsytosphexa.py datafiles/$ICfile.00000 datafiles/$ICfile
	runline=( "$dir/$coderun" --init "./datafiles/$ICfile.h5" -o "$runname/$runname.h5" $EXTRACONDSPHEXA -s "$endtime" -w "$outinterval" -f "$fields")
    fi

    if [ "$4" == 3 ]; then
	# AREPO moving-mesh: convert the tipsy IC -> AREPO HDF5 IC (sph_interp),
	# then build an AREPO param.txt whose BoxSize / run times come from the
	# tipsy .param we generated for this lattice (NOT the particle extent).
	paramfile="datafiles/$ICfile.param"

	get_param () {
	    awk -F= -v key="$1" '
		$1 ~ "^[[:space:]]*" key "[[:space:]]*$" {
		    gsub(/[[:space:]]*/, "", $2); print $2
		}' "$paramfile"
	}

	dxP=$(get_param dxPeriod)
	dyP=$(get_param dyPeriod)
	dzP=$(get_param dzPeriod)
	iOutInterval=$(get_param iOutInterval)
	if [[ -z "$dxP" || -z "$iOutInterval" ]]; then
	    echo "Error: missing dxPeriod/iOutInterval in $paramfile"
	    exit 1
	fi
	[ -z "$dyP" ] && dyP=$dxP
	[ -z "$dzP" ] && dzP=$dxP
	# AREPO BoxSize is the X period; the other axes are set by the COMPILE-TIME
	# LONG_Y/LONG_Z = dyPeriod/dxPeriod, dzPeriod/dxPeriod (boxSize_Y=BoxSize*LONG_Y).
	# (The run times / timestepping go into param.txt on the Python side below.)
	boxsize=$dxP
	LONGY=$(awk -v a="$dyP" -v b="$dxP" 'BEGIN { printf "%.10g", a/b }')
	LONGZ=$(awk -v a="$dzP" -v b="$dxP" 'BEGIN { printf "%.10g", a/b }')

	# Non-cubic box -> AREPO needs LONG_Y/LONG_Z compiled in (public AREPO is
	# cubic by default). Cubic (both ~1): use the prebuilt MHD `Arepo`. Else
	# build (once, cached by aspect tag) Arepo_asp<LY>_<LZ> from Config.sh + the
	# LONG_Y/LONG_Z lines.
	iscubic=$(awk -v y="$LONGY" -v z="$LONGZ" 'BEGIN { d=0.0005; print (y>1-d && y<1+d && z>1-d && z<1+d)?1:0 }')
	if [ "$iscubic" != "1" ]; then
	    tag=$(awk -v y="$LONGY" -v z="$LONGZ" 'BEGIN { printf "%.4f_%.4f", y, z }' | tr '.' 'p')
	    coderun="Arepo_asp${tag}"
	    if [ ! -x "$dir/$coderun" ]; then
		echo "Building AREPO with LONG_Y=$LONGY LONG_Z=$LONGZ -> $coderun"
		cfg="Config_asp${tag}.sh"
		cp "$dir/Config.sh" "$dir/$cfg"
		printf "\nLONG_Y=%s\nLONG_Z=%s\n" "$LONGY" "$LONGZ" >> "$dir/$cfg"
		( cd "$dir" && make -j8 CONFIG="$cfg" BUILD_DIR="build_asp${tag}" EXEC="$coderun" )
		if [ ! -x "$dir/$coderun" ]; then echo "Error: AREPO LONG build failed."; exit 1; fi
	    fi
	    echo "Non-cubic box [$dxP,$dyP,$dzP] -> using $coderun (LONG_Y=$LONGY LONG_Z=$LONGZ)"
	fi

	# tipsy IC -> AREPO IC (sph_interp). $AREPOICMODE picks the conversion
	# (default "particles": 1:1 particle->generator copy). Periodic box + dTuFac
	# come from the .param/.log inside the converter; coords are framed per-axis
	# into [0,period_axis] so they fit the LONG_Y/LONG_Z box.
	# --param-out also writes a FRESH AREPO param.txt with run parameters
	# (timestepping, box, output cadence) derived from the tipsy .param --
	# replacing the old noh-template sed (MaxSizeTimestep is now dDelta, not 0.01).
	arepoIC="datafiles/${ICfile}_arepo.hdf5"
	arepoparam="$runname/param.txt"
	convertcmd="python ${analysisdir}/ICconvert_tipsytoarepo.py datafiles/$ICfile.00000 $arepoIC $AREPOICMODE --param-out $arepoparam --output-dir $runname/ --boxsize $boxsize"
	echo "$convertcmd"
	eval "$convertcmd"
	if [ $? -ne 0 ]; then echo "Error: AREPO IC conversion failed."; exit 1; fi

	# RestartFlag 0 = start from the IC (AREPO rebuilds the Voronoi mesh).
	runline=( $mpirun "$dir/${coderun}" "$arepoparam" 0)
    fi
    echo "${runline[@]}"
    "${runline[@]}"

    if [ "$4" == 2 ]; then
	# convert the SPH-EXA hdf5 output back into tipsy snapshots
	python $current_dir/setupfiles/ICconvert_sphexatotipsy.py $runname/$runname.h5 $runname $iOutInterval
    fi

    if [ "$4" == 3 ]; then
	# convert the AREPO snap_*.hdf5 output back into tipsy snapshots so the
	# analysis pipeline (which reads tipsy) can process the AREPO run. Pass the
	# tipsy .param so the per-axis periods (non-cubic box) centre the snapshots
	# correctly and land in the .log.
	python $current_dir/setupfiles/ICconvert_arepototipsy.py $runname $runname $iOutInterval datafiles/$ICfile.param
    fi

## Reference (regression baseline): -r saves the analysis result of this run as
## '<runname>_reference.json' so later runs can be compared against it. We forward
## the same -a..-h values used for the IC (only those the user set; unset ones
## fall back to defaults that the analysis script shares with this setup), so the
## analysis parameters always match the IC parameters.
if [ "$doreference" -eq 1 ]; then
    analysisscript="${analysisdir}/IC_analysis_${testsim}.py"
    if [ -f "$analysisscript" ]; then
        refflags=""
        letters=(a b c d e f g h)
        for i in 0 1 2 3 4 5 6 7; do
            if [ -n "${args[$i]}" ]; then
                refflags="$refflags -${letters[$i]} ${args[$i]}"
            fi
        done
        echo "Saving analysis reference: python $analysisscript $runname $refflags --save-reference"
        python "$analysisscript" "$runname" $refflags --save-reference
    else
        echo "No analysis script $analysisscript for '$testsim'; skipping -r reference save."
    fi
fi
