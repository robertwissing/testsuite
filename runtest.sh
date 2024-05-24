#!/bin/bash


#Directory for the ic generator and the setupfiles (you have to do make for icgenerator/mdl/mpi and make mpi for icgenerator/gasoline/
current_dir=$(pwd -P)
ICdir=$current_dir/icgenerator/gasoline/
gasolineicgen=gasoline.gdicgen
analysisdir=$current_dir/setupfiles/

# Directory where the executable(gasoline/changa) is and the name of it:

#dir=/mn/stornext/u3/robertwi/Projects/RECON/gasolinenew/
#gasolinerun=gasoline.mhd
dir=$current_dir/../../
gasolinerun=ChaNGa.smp

#If running serial/mpi or with charmrun
mpirun="$dir/charmrun.smp +p 63"
#mpirun=""
#mpirun="mpirun "

#An extra word added at the end of the simulation output files and directory
EXTRANAME="CHANGA"

#Specific config for glass generation. Can leave be, should work for all, but if sharper density is needed increase ICrhopow (though at the cost of higher zero order errors usually), ICR0 determines how much of the periodic box particles should move on average in beginning, can be set to 1, ensures that whatever IC you give it, it will go towards the target density (and remove long-range biases). Initially density gradients are bad and it can be hard for the partilces to correctly fit the large-scale density structure in one go, as such we quickly decay the speed and then increase it to some factor of the initial ICR0 this is eveventually reduced as we get closer to the target density. ICR0rate is how fast the global speed should decay. and -n 4800 is just a maximum relaxation number it usually relaxes within 200-1000 steps(depending on the other parameters and IC). and oi is output incase you want to see what is going on during relaxation.
ICCONFIG="-ICR0 1.0 -n 4800 -oi 4800 -ICRLOOP0 0.5 -ICrhopow 2.0 -ICR0Rate 1.5"

#Extra condition during runtime
#EXTRACONDALL="-SpinEta 0.0  "
EXTRACONDALL=""


usage() {
    echo "Usage: $0 [-i] [-a value1] ... [-h value8] <nameofsimulation> <Nsmooth> <Nres> <0: lattice | 1: rand | 2: glass>"
    echo ""
    echo "Options:"
    echo "  -a, -b, -c, ...    These options depend on the simulation. To see the specific options, enter the nameofsimulation without arguments."
    echo "  -i                 Interactive."
    echo ""
    echo "Arguments:"
    echo "  nameofsimulation   Name of the simulation to be run (e.g., shocktube)."
    echo "  Nsmooth            Number of neigbours within smoothing kernel"
    echo "  Nres               Resolution for the simulation (usually specified as number of elements in the x direction)."
    echo "  0 | 1 | 2          Initial conditions type:"
    echo "                     0: lattice"
    echo "                     1: random"
    echo "                     2: glass"
    echo "If glass data file exist for a specific Nsmooth it will skipp the glass generation and run the simulation. If you choose an nSmooth that does not have a glass datafile it will do the whole glass generation again."
    echo ""
    echo "List of available simulations:"
    echo "  accdisk            Accretion disk simulation."
    echo "  alfven             Alfvén wave simulation."
    echo "  gresho             Gresho vortex simulation."
    echo "  mhdcollapse        MHD collapse simulation."
    echo "  kh                 Kelvin-Helmholtz instability simulation."
    echo "  orzag              Orszag-Tang vortex simulation."
    echo "  mhdrotor           MHD rotor simulation."
    echo "  mhdloop            MHD loop advection simulation."
    echo "  shocktube          Different kinds of shock tubes."
    echo "                     (Enter 'shocktube' without arguments for more details on available options)"
    echo "  sedov              Sedov-Taylor blast wave simulation."
    echo ""
    exit 1
}


# Parse command line options
interactive=0
declare -a args=("" "" "" "" "" "" "" "")
while getopts ":ia:b:c:d:e:f:g:h:" opt; do
  case $opt in
    i)
      interactive=1
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
    EXTRACONDRUN="$EXTRACONDALL -n 900 -dt 0.0001 "
    prefile_extra="rd${rhodiff}"
    glassfile_extra="rd${rhodiff}"
    runname_extra="mu${mu}rd${rhodiff}Erat${Erat}"
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
    EXTRACONDRUN="$EXTRACONDALL "
    prefile_extra="hr$hr"
    glassfile_extra="hr$hr"
    runname_extra="hr$hr"
}

# Function for alfven setup with generic parameters
setup_alfven() {
    # Set default values if not specified
    rotate=${args[0]:-1}

    echo "Using settings for circularized alfvenwave:"
    echo "rotate wave (rotate) -a : $rotate"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rotate "
    EXTRACONDRUN="$EXTRACONDALL "
    prefile_extra=""
    glassfile_extra=""
    runname_extra="rot${rotate}"
}

# Function for Shock-tube setup with generic parameters
setup_shocktube() {
    shock=${args[0]:-1.0}
    avisc=${args[1]:-1.0}
    bvisc=${args[2]:-2}



    
    # Array mapping shock indices to shock names
    shocks=(
        "Sod shock"
        "Ryu 1a"
        "Ryu 1b"
        "Ryu 2a"
        "Ryu 2b"
        "Brio-Wu (Ryu 5a)"
        "C-shock"
        "Steady shock"
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
    echo "alpha visc (argument -b): $avisc"
    echo "beta visc (argument -c): $bvisc"

    # Setup based on the potentially updated values
    EXTRACONDIC="$shock "
    EXTRACONDRUN="$EXTRACONDALL -beta $bvisc -alpha $avisc "
    prefile_extra="${shockname// /_}"
    glassfile_extra="${shockname// /_}"
    runname_extra="${shockname// /_}"

    echo "prefile_extra: $prefile_extra"
    echo "glassfile_extra: $glassfile_extra"
    echo "runname_extra: $runname_extra"

}

# Function for Orzag-Tang Vortex setup with generic parameters
setup_orzag() {
    avisc=${args[0]:-1.0}
    bvisc=${args[1]:-2}
    # Set default values if not specified
    echo "Using settings for Orzag-Tang Vortex:"
    echo "alpha visc -a : $avisc"
    echo "beta visc -b : $bvisc"

    # Setup based on the potentially updated values
    EXTRACONDIC=" "
    EXTRACONDRUN="$EXTRACONDALL -beta $bvisc -alpha $avisc "
    prefile_extra=""
    glassfile_extra=""
    runname_extra=""
}

# Function for Magnetized Rotor setup with generic parameters
setup_mhdrotor() {
    rhodisk=${args[0]:-10.0}
    avisc=${args[1]:-1.0}
    bvisc=${args[2]:-2}
    # Set default values if not specified
    echo "Using settings for Magnetized Rotor:"
    echo "density disk -a : $rhodisk"
    echo "alpha visc -b : $avisc"
    echo "beta visc -c : $bvisc"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodisk "
    EXTRACONDRUN="$EXTRACONDALL -beta $bvisc -alpha $avisc "
    prefile_extra="rhodiff${rhodisk}"
    glassfile_extra="rhodiff${rhodisk}"
    runname_extra="rhodiff${rhodisk}"
}

# Function for Magnetized advection loop setup with generic parameters
setup_mhdloop() {
    rhodisk=${args[0]:-1.0}
    avisc=${args[1]:-1.0}
    bvisc=${args[2]:-2}
    # Set default values if not specified
    echo "Using settings for Magnetized Advection Loop:"
    echo "density disk -a : $rhodisk"
    echo "alpha visc -b : $avisc"
    echo "beta visc -c : $bvisc"

    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodisk "
    EXTRACONDRUN="$EXTRACONDALL -beta $bvisc -alpha $avisc "
    prefile_extra="rhodiff${rhodisk}"
    glassfile_extra="rhodiff${rhodisk}"
    runname_extra="rhodiff${rhodisk}"
}

# Function for gresho setup with generic parameters
setup_gresho() {
    # Set default values if not specified
    result=$(echo "scale=5; sqrt(3/25)" | bc -l)
    echo "$result"
    mach=${args[0]:-$result}
    avisc=${args[1]:-0.01}
    bvisc=${args[2]:-2}
    steps=${args[3]:-300}
    
    echo "Using settings for gresho:"
    echo "mach number of flow -a : $mach"
    echo "alpha visc -b : $avisc"
    echo "beta visc -c : $bvisc"
    echo "number of steps -d : $steps"
    
    # Setup based on the potentially updated values
    EXTRACONDIC="$mach "
    EXTRACONDRUN="$EXTRACONDALL -n $steps -beta $bvisc -alpha $avisc"
    prefile_extra=""
    glassfile_extra=""
    runname_extra="mach${mach}"
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
    avisc=${args[5]:-0.01}
    bvisc=${args[6]:-2}
    steps=${args[7]:-200}

    echo "Using settings for kelvin-helmholtz:"
    echo "rho diff -a : $rhodiff"
    echo "mach number of flow -b : $mach"
    echo "smooth density disc -c : $smooth"
    echo "initial magnetic field strength -d : $B0"
    echo "Direction of field x=1 y=2 z=3 -e : $Bdir"
    echo "alpha visc -f : $avisc"
    echo "beta visc -g : $bvisc"
    echo "number of steps -h : $steps"
    # Setup based on the potentially updated values
    EXTRACONDIC="$rhodiff $mach $smooth $B0 $Bdir "
    EXTRACONDRUN="$EXTRACONDALL -n $steps -beta $bvisc -alpha $avisc"
    prefile_extra="rhodiff${rhodiff}smooth${smooth}"
    glassfile_extra="rhodiff${rhodiff}smooth${smooth}"
    runname_extra="rhodiff${rhodiff}mach${mach}smooth${smooth}B${B0}Bdir${Bdir}"
}


# Function for Sedov blast wave + mhd setup with generic parameters
setup_sedov() {
    select=${args[0]:-1.0}
    betain=${args[1]:-2.0}
    ublast=${args[2]:-10.0}
    Rin=${args[3]:-0.125}
    Pin=${args[4]:-100.0}
    Pout=${args[5]:-1.0}
    avisc=${args[6]:-1.0}
    bvisc=${args[7]:-2.0}
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
    echo "alpha visc -g : $avisc"
    echo "beta visc -h : $bvisc"

    # Setup based on the potentially updated values
    EXTRACONDIC="$select $betain $ublast $Rin $Pin $Pout"
    EXTRACONDRUN="$EXTRACONDALL -beta $bvisc -alpha $avisc "
    prefile_extra="select${select}beta${betain}ublast${ublast}"
    glassfile_extra=""
    runname_extra="select${select}beta${betain}ublast${ublast}"
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
    echo "Usage: $0 <Nsmooth> <Nres> <0: lattice 1: rand 2: glass> < 1 isph else no isph > <kernel(see file for list of kernels)>"
    exit 1
fi
echo $0 $1 $2 $3 $4 $5
#Specifically for different SPH versions
if [ "$4" -eq 1 ]; then
is=isph
fi
if [ "$4" -eq 2 ]; then
is=corrdens
fi
if [ "$4" -eq 3 ]; then
is=corrdenssym
fi
if [ "$4" -eq 4 ]; then
is=correntrop
fi
if [ "$4" -eq 5 ]; then
is=corrdenssymisph
fi
if [ "$4" -eq 6 ]; then
is=corrdenssymentrop
fi
if [ "$4" -eq 7 ]; then
is=corrdensisph
fi
if [ "$4" -eq 8 ]; then
is=corrdensentrop
fi

# Specifically done for special version of gasoline, which allow for easy changing of kernel
kern_list=("32D" "64D" "128D" "256D" "32DIS" "64DIS" "128DIS" "256DIS" "32" "64" "128" "256" "CM04" "WuC2" "BUH" "CM05" "C2" "C4" "C6" "C8" "WuC4")
alpha_values=(
    "-alp14 1.0"
    "-alp16 1.0"
    "-alp18 1.0"
    "-alp20 1.0"
    "-alp15 1.0"
    "-alp17 1.0"
    "-alp19 1.0"
    "-alp21 1.0"
    "-alp 1.0"
    "-alp2 1.0"
    "-alp3 1.0"
    "-alp4 1.0"
    "-alp5 1.0"
    "-alp6 1.0"
    "-alp7 1.0"
    "-alp8 1.0"
    "-alp9 1.0"
    "-alp10 1.0"
    "-alp11 1.0"
    "-alp12 1.0"
    "-alp13 1.0"
)

if [[ -z $5 ]]; then
    kern=""
    alpha=""
else
if [[ $5 -ge 0 && $5 -lt ${#kern_list[@]} ]]; then
    kern="${kern_list[$5]}"
    alpha="${alpha_values[$5]}"
    else
    kern=""
    alpha=""
fi
fi
    echo $kern

# Set distribution and d1 based on the third command-line argument
if [ "$3" -eq 0 ]; then
echo "RUNNING LATTICE MODE"
    distri="lattice"
    distri2="lattice"
    d1=0
    simfilename="$testsim$2_$distri${prefile_extra}"
elif [ "$3" -eq 1 ]; then
echo "RUNNING RANDOM MODE"
    distri="rand"
    distri2="rand"
    d1=1
    simfilename="$testsim$2_$distri${prefile_extra}"
else
echo "RUNNING GLASS MODE"
    distri="rand"
    distri2="glass"
    d1=1
    simfilename="$testsim$2_glassN$1KERN$kern${is}${runname_extra}-I"
    echo "$simfilename"
fi

prefile="$testsim$2_${distri}${prefile_extra}"
glassfile_prefix="$testsim$2_glassN$1KERN$kern${is}${glassfile_extra}"
ICfile="$simfilename"
runname="$testsim$2_N$1${distri2}K$kern${is}${runname_extra}${EXTRANAME}"
directory="initruns"
glassfile="${directory}/${glassfile_prefix}"

mkdir datafiles
mkdir initruns

    # Check if the file already exists
    if [ -f "datafiles/$prefile.00000" -o -f "${glassfile}_IC" -o -f "datafiles/$ICfile.00000" ]; then
        echo "File $prefile.00000 already exists. Skipping creation."
    else
    # File doesn't exist, create it
    preIC="python ${analysisdir}/IC_createsetup.py $testsim $2 $d1 $prefile $EXTRACONDIC"
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
	relaxglass=(mpirun -np $(($(nproc) / 2)) ${ICdir}/${gasolineicgen}$is -o $glassfile -I datafiles/$prefile.00000 -s $1 $alpha $ICCONFIG datafiles/$prefile.param)
	echo "${relaxglass[@]}"
	"${relaxglass[@]}"

	# Check if the command was successful
	if [ $? -ne 0 ]; then
	    echo "Error: Command failed to execute successfully."
	    exit 1
	fi

	regex_pattern="${glassfile_prefix}\.(\d+)"

	# Find the highest file number that matches the specific pattern in the directory
	highest_number_file=$(ls -1 $directory | grep -oP "${regex_pattern}" | sort -nr | head -n1)
	echo "HIGHEST NUMBER ${highest_number_file}"
	
	echo "${directory}/${highest_number_file}" "${glassfile}_IC"
	mv "${directory}/${highest_number_file}" "${glassfile}_IC"

    fi

    if [ -f "datafiles/$ICfile.00000" ]; then
        echo "File $ICfile.00000 already exists. Skipping creation."
    else
    ICfromglass="python ${analysisdir}/IC_createsetup.py $testsim $2 "${glassfile}_IC" ${ICfile} $EXTRACONDIC "
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
    runline=( $mpirun $dir/${gasolinerun}$is -s $1 -o $runname/$runname -I datafiles/$ICfile.00000 $alpha $EXTRACONDRUN datafiles/$ICfile.param)
    echo "${runline[@]}"
    "${runline[@]}"
