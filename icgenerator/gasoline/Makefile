 # Makefile for pkdgrav or gasoline (specify below)
# Then type make <platform>
# <platform> should be one of:
# null, sgi, pvm, pthread, pthread_dec, pthread_sgi, ampi, lam_mpi, qmpi, mpi,
# spx, t3d, t3dmpi, t3empi, ksr, charm
#
# You will get a pair of new files: gasoline.NNNNNNNN and buildlog.NNNNNNNN,
# where NNNNNNNN is an 8 character hash.
# The first is the actual gasoline executable, and the second one contains all
# defines that were used to build gasoline.  You should copy that file into
# the directory of any simulation you run so you know exactly which options were
# used.
PNG_INCL =
PNG_LIB =
PNG_OBJ =
PNG_DEF =

# If you do dump frame stuff with png you need these libs ( -DUSE_PNG )
#PNG_INCL = -I/usr/include
#PNG_LIB = -L/usr/lib -lpng -lz
#PNG_LIB = -lpng -lz
#PNG_OBJ = writepng.o
#PNG_DEF = $(PNG_INCL) -DUSE_PNG

# if you use OUTURBDRIVER you need libgsl
#TURB_OBJ = outurb.o
#GSL_LIB = -lgsl -lgslcblas
TURB_OBJ =
GSL_LIB =

COOLING_OBJ =
COOLING_DEF = -DNOCOOLING

BASE_LD_FLAGS = $(PNG_LIB) $(GSL_LIB)
BASE_LD_FLAGS = $(PNG_LIB) $(GSL_LIB) $(COOLING_LIB)

# Which C compiler should I use?
CC = icc
#CC = gcc

CODE_DEF =
# Gasoline uses a number of defines to enable/disable code features,
# debug sections of the code, and switch out different physics.
# Below is a complete list of the available user defines for gasoline.
# Most of them will NEVER be useful to you, as gasoline includes a lot of
# very case-specific code. If you don't know what a define means, you probably
# don't need it. The uncommented ones that come in the repository are probably fine.
# Simply uncomment a line to enable that define, or comment it to disable it.

CODE_DEF = -DDTADJUST -DEPSACCH  -DGASOLINE 
CODE_DEF += -DSETTRAPFPE
CODE_DEF += -DACCZERO
CODE_DEF += -DGDFORCE
CODE_DEF += -DCHANGESOFT
CODE_DEF += -DDODVDS
#CODE_DEF += -DCOMBOKERNEL
CODE_DEF += -DICGEN
CODE_DEF += -DWENDLAND
#CODE_DEF += -DCORRPARTITION
#CODE_DEF += -DISPH
#CODE_DEF += -DCORRENTROPYPART
CODE_DEF += -DCHECKQUALITY
#CODE_DEF += -DWENDLANDCOMBO

BASE_DEF = $(PNG_DEF) $(CODE_DEF) $(COOLING_DEF)

# The filename of the compiled executable
EXE = gasoline.gdicgen
#
#       NULL defines
#
NULL_MDL		= ../mdl/null
#NULL_CFLAGS		= -g -I$(NULL_MDL) $(BASE_DEF)
NULL_CFLAGS		= -O3 -I$(NULL_MDL) $(BASE_DEF)
NULL_LD_FLAGS	= $(BASE_LD_FLAGS) #-L/usr/lib -L/lib
NULL_XOBJ		= erf.o v_sqrt1.o
NULL_LIBMDL		= $(NULL_MDL)/mdl.o -lm

#
#       DBG defines
#
DBG_MDL		= ../mdl/null
DBG_CFLAGS		= -g -I$(DBG_MDL) $(BASE_DEF)
DBG_LD_FLAGS	= $(BASE_LD_FLAGS) -L/usr/lib -L/lib
DBG_XOBJ		= erf.o v_sqrt1.o
DBG_LIBMDL		= $(NULL_MDL)/mdl.o -lm

#
#       SGI defines
#
SGI_MDL			= ../mdl/mpi
SGI_CFLAGS		= -O3 -I$(SGI_MDL) $(BASE_DEF) -mips4 -64 -r10000
SGI_LD_FLAGS	= -mips4 -64 -r10000 $(BASE_LD_FLAGS)
SGI_XOBJ		=
SGI_LIBMDL		= $(SGI_MDL)/mdl.o -lmpi -lm
SGI_MDL_CFLAGS	= -g2 -O2 -mips4 -64 -r10000

#
#       LINUX AMPI (Charm MPI) defines
#
CHARM=../charm/net-linux/bin/charmc
CHARMLINK=
AMPI_MDL			= ../mdl/ampi
AMPI_CFLAGS		= -O3 -malign-double -mpentiumpro -I$(AMPI_MDL) $(BASE_DEF)
AMPI_LD_FLAGS		=  $(BASE_LD_FLAGS)
AMPI_XOBJ                = erf.o v_sqrt1.o
AMPI_LIBMDL              = $(AMPI_MDL)/mdl.o -language ampi $(CHARMLINK) -lm

#
#       LINUX LAM MPI defines
#
LAM_MDL			= ../mdl/mpi
LAM_CFLAGS		= -fast -I$(LAM_MDL) $(BASE_DEF) -DMPI_LINUX
LAM_LD_FLAGS		=  $(BASE_LD_FLAGS)
LAM_XOBJ                = erf.o v_sqrt1.o
LAM_LIBMDL              = $(LAM_MDL)/mdl.o -lm
LAM_MDL_CFLAGS  = -fast -I$(LAM_MDL) $(BASE_DEF)

#
#       Quadrics MPI defines
#
QMPI_MDL                        = ../mdl/mpi
QMPI_CFLAGS             = -O3 -I$(QMPI_MDL) $(BASE_DEF) -DMPI_LINUX -DCCC
QMPI_LD_FLAGS           =  $(BASE_LD_FLAGS)
QMPI_XOBJ                = erf.o
QMPI_LIBMDL              = $(QMPI_MDL)/mdl.o -lmpi -lm
QMPI_MDL_CFLAGS  = -I$(QMPI_MDL) $(BASE_DEF) -DMPI_LINUX


#       SP1/2 defines
#
SPX_MDL			= ../mdl/mpi
SPX_CFLAGS		= -O3 -march=core-avx2 -I$(SPX_MDL) $(BASE_DEF) -I /usr/include/tirpc -std=gnu89
SPX_LD_FLAGS	= $(BASE_LD_FLAGS) -ltirpc
SPX_XOBJ		= v_sqrt1.o
SPX_LIBMDL		= $(SPX_MDL)/mdl.o -lm
SPX_MDL_CFLAGS	= -O3

#
#		PVM defines
#
PVMDIR	= $(PVM_ROOT)
PVMLIB	= $(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
BDIR	= $(HOME)/pvm3/bin
XDIR	= $(BDIR)/$(PVM_ARCH)
PVM_MDL		= ../mdl/pvm
PVM_CFLAGS	= -O3 -I$(PVMDIR)/include -I$(PVM_MDL) $(BASE_DEF)
PVM_XOBJ	= v_sqrt1.o
PVM_LIBMDL	= $(PVM_MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) -lm
PVM_LD_FLAGS	= $(BASE_LD_FLAGS)

#
#       PTHREAD defines
#
PTHREAD_MDL			= ../mdl/pthread
PTHREAD_CFLAGS		= -O3  -D_REENTRANT -I$(PTHREAD_MDL) $(BASE_DEF)
PTHREAD_LD_FLAGS 	=  $(BASE_LD_FLAGS)
PTHREAD_XOBJ		= erf.o v_sqrt1.o
PTHREAD_LIBMDL 		= $(PTHREAD_MDL)/mdl.o -lm -lpthread

#
#       PTHREAD_SGI defines
#
PTHREAD_SGI_MDL			= ../mdl/pthread
PTHREAD_SGI_CFLAGS		= -O2 -D_REENTRANT -I$(PTHREAD_SGI_MDL) $(BASE_DEF) -mips4 -64 -r10000
PTHREAD_SGI_LD_FLAGS 	= -mips4 -64 -r10000 $(BASE_LD_FLAGS)
PTHREAD_SGI_XOBJ		=
PTHREAD_SGI_LIBMDL 		= $(PTHREAD_SGI_MDL)/mdl.o -lm -lpthread
PTHREAD_SGI_MDL_CFLAGS	= -O2 -mips4 -64 -r10000

#
#       T3D MPP defines
#
T3D_MDL	= ../mdl/mpp
V_SQRT	= ../v_sqrt/lib
V_SQRT1	= ../v_sqrt1/lib
RPC		= ../rpc
T3D_CFLAGS		= -O3 -g -DCRAY_T3D -I$(T3D_MDL) -I$(V_SQRT) -I$(RPC) $(BASE_DEF)
T3D_XOBJ		= hyperlib.o
T3D_LIBMDL		= -O3 -g-L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) \
					$(T3D_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3D_LD_FLAGS	= $(BASE_LD_FLAGS)

#
#       T3DMPI MPP defines
#
T3DMPI_MDL	= ../mdl/t3dmpi
V_SQRT		= ../v_sqrt/lib
V_SQRT1		= ../v_sqrt1/lib
RPC			= ../rpc
T3DMPI_CFLAGS	= -O3 -DCRAY_T3D -I$(T3DMPI_MDL) -I$(V_SQRT) -I$(RPC) $(BASE_DEF)
T3DMPI_XOBJ		= hyperlib.o
T3DMPI_LIBMDL	= -O3 -L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) \
					$(T3DMPI_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3DMPI_LD_FLAGS	= $(BASE_LD_FLAGS)

#
#       T3EMPI MPP defines
#
T3EMPI_MDL		= ../mdl/t3empi
T3EMPI_CFLAGS	= -O3 -DCRAY_T3D -I$(T3EMPI_MDL) -I$(RPC) $(BASE_DEF)
T3EMPI_XOBJ		= hyperlib.o v_sqrt1.t3x.o
T3EMPI_LIBMDL	= $(T3EMPI_MDL)/mdl.o -lmpi -lm
T3DMPI_LD_FLAGS	= $(BASE_LD_FLAGS)

#
#       LINUX Charm  defines
#
CHARM=../charm/net-linux/bin/charmc
CHARMLINK=
CHARM_MDL			= ../mdl/charm
CHARM_CFLAGS		= -verbose -g -I$(CHARM_MDL) $(BASE_DEF)
CHARM_MDL_CFLAGS	= -verbose -g -Wall
CHARM_LD_FLAGS		=  $(BASE_LD_FLAGS) -language charm++ -memory os
CHARM_XOBJ                = erf.o v_sqrt1.o
CHARM_LIBMDL              = $(CHARM_MDL)/mdl.o $(CHARMLINK) -lm

#
#       KSR1 defines
#
KSR_MDL			= ../mdl/ksr
KSR_CFLAGS		= -O2 -w2 -I$(KSR_MDL) $(BASE_DEF)
KSR_LD_FLAGS	= -para $(BASE_LD_FLAGS)
KSR_XOBJ		= erf.o v_sqrt1.ksr.o
KSR_LIBMDL		= $(KSR_MDL)/mdl.o -lm -lrpc

OBJ	= main.o master.o param.o outtype.o pkd.o pst.o grav.o \
	  ewald.o walk.o eccanom.o hypanom.o fdl.o htable.o smooth.o \
	  smoothfcn.o collision.o qqsmooth.o $(COOLING_OBJ) cosmo.o romberg.o \
	  starform.o feedback.o millerscalo.o supernova.o supernovaia.o \
	  startime.o stiff.o runge.o dumpframe.o dffuncs.o dumpvoxel.o \
	  rotbar.o special.o ssio.o $(PNG_OBJ) $(TURB_OBJ) \
	  treezip.o log.o

EXTRA_OBJ = erf.o hyperlib.o v_sqrt1.o v_sqrt1.ksr.o v_sqrt1.t3x.o


default:
	@echo "Please tell me what architecture to make."
	@echo "Choices are null, sgi, pvm, pthread, pthread_dec, pthread_sgi, ampi, lam_mpi, qmpi, mpi, spx, t3d, t3dmpi, t3empi, ksr, charm"

install:
	@echo "No installation rules."

clean:
	-rm -f $(OBJ) $(EXTRA_OBJ)

spotless: clean
	@echo rm -rf $(EXE)
	@rm -f $(EXE)
	@rm -f ../mdl/*/mdl.o

depend:
	makedepend -Y -- $(BASE_DEF) -- *.c

$(XDIR):
	-mkdir $(BDIR)
	-mkdir $(XDIR)

null:
	cd $(NULL_MDL); make "CC=$(CC)" "CFLAGS=$(NULL_CFLAGS)"
	make $(EXE) "CFLAGS=$(NULL_CFLAGS)" "LD_FLAGS=$(NULL_LD_FLAGS)"\
		"MDL=$(NULL_MDL)" "XOBJ=$(NULL_XOBJ)" "LIBMDL=$(NULL_LIBMDL)"
dbg:
	cd $(DBG_MDL); make "CC=$(CC)" "CFLAGS=$(DBG_CFLAGS)"
	make $(EXE) "CFLAGS=$(DBG_CFLAGS)" "LD_FLAGS=$(DBG_LD_FLAGS)"\
		"MDL=$(DBG_MDL)" "XOBJ=$(DBG_XOBJ)" "LIBMDL=$(DBG_LIBMDL)"

sgi:
	cd $(SGI_MDL); make CC=cc "CC_FLAGS=$(SGI_MDL_CFLAGS)"
	make $(EXE) CC=cc "CFLAGS=$(SGI_CFLAGS)" "LD_FLAGS=$(SGI_LD_FLAGS)"\
		"MDL=$(SGI_MDL)" "XOBJ=$(SGI_XOBJ)" "LIBMDL=$(SGI_LIBMDL)"

pvm:
	cd $(PVM_MDL); aimk
	make $(EXE) "CFLAGS=$(PVM_CFLAGS)" "LD_FLAGS=$(PVM_LD_FLAGS)"\
		"MDL=$(PVM_MDL)" "XOBJ=$(PVM_XOBJ)" "LIBMDL=$(PVM_LIBMDL)"
	mv -f $(EXE) $(XDIR)

pthread:
	cd $(PTHREAD_MDL); make "CC=$(CC)" "CFLAGS=$(PTHREAD_CFLAGS)"
	make $(EXE) "EXE=$(EXE)" "CFLAGS=$(PTHREAD_CFLAGS)" "LD_FLAGS=$(PTHREAD_LD_FLAGS)"\
		"MDL=$(PTHREAD_MDL)" "XOBJ=$(PTHREAD_XOBJ)" "LIBMDL=$(PTHREAD_LIBMDL)"

pthread_dec:
	cd $(PTHREAD_MDL); make CC=ccc
	make $(EXE) CC=ccc "CFLAGS=$(PTHREAD_CFLAGS) -fast -arch ev6" "LD_FLAGS=$(PTHREAD_LD_FLAGS)"\
		"MDL=$(PTHREAD_MDL)" "XOBJ=$(PTHREAD_XOBJ)" "LIBMDL=$(PTHREAD_LIBMDL)"

pthread_sgi:
	cd $(PTHREAD_MDL); make CC=cc "CC_FLAGS=$(PTHREAD_SGI_MDL_CFLAGS)"
	make $(EXE) "CFLAGS=$(PTHREAD_SGI_CFLAGS)" "LD_FLAGS=$(PTHREAD_SGI_LD_FLAGS)"\
		"MDL=$(PTHREAD_SGI_MDL)" "XOBJ=$(PTHREAD_SGI_XOBJ)" "LIBMDL=$(PTHREAD_SGI_LIBMDL)"

ampi:
	cd $(AMPI_MDL); make CC=$(CHARM) "CFLAGS=$(AMPI_MDL_CFLAGS)"
	make $(EXE) CC=$(CHARM) "CFLAGS=$(AMPI_CFLAGS)" "LD_FLAGS=$(AMPI_LD_FLAGS)"\
		"MDL=$(AMPI_MDL)" "XOBJ=$(AMPI_XOBJ)" "LIBMDL=$(AMPI_LIBMDL)"

lam_mpi:
	cd $(LAM_MDL); make CC=mpicc "CC_FLAGS=$(LAM_MDL_CFLAGS)"
	make $(EXE) CC=mpicc "CFLAGS=$(LAM_CFLAGS)" "LD_FLAGS=$(LAM_LD_FLAGS)"\
		"MDL=$(LAM_MDL)" "XOBJ=$(LAM_XOBJ)" "LIBMDL=$(LAM_LIBMDL)"

qmpi:
	cd $(QMPI_MDL); make "CC_FLAGS=$(QMPI_MDL_CFLAGS)"
	make $(EXE) "CFLAGS=$(QMPI_CFLAGS)" "LD_FLAGS=$(QMPI_LD_FLAGS)"\
		"MDL=$(QMPI_MDL)" "XOBJ=$(QMPI_XOBJ)" "LIBMDL=$(QMPI_LIBMDL)"

mpi: spx

spx:
	cd $(SPX_MDL); make CC=mpicc "CC_FLAGS=$(SPX_MDL_CFLAGS)"
	make $(EXE) CC=mpicc "CFLAGS=$(SPX_CFLAGS)" "LD_FLAGS=$(SPX_LD_FLAGS)"\
		"MDL=$(SPX_MDL)" "XOBJ=$(SPX_XOBJ)" "LIBMDL=$(SPX_LIBMDL)"

t3d:
	cd $(T3D_MDL); make
	make $(EXE) "CFLAGS=$(T3D_CFLAGS)" "LD_FLAGS=$(T3D_LD_FLAGS)"\
		"MDL=$(T3D_MDL)" "XOBJ=$(T3D_XOBJ)" "LIBMDL=$(T3D_LIBMDL)"

t3dmpi:
	cd $(T3DMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3DMPI_CFLAGS)" "LD_FLAGS=$(T3DMPI_LD_FLAGS)"\
		"MDL=$(T3DMPI_MDL)" "XOBJ=$(T3DMPI_XOBJ)" "LIBMDL=$(T3DMPI_LIBMDL)"

t3empi:
	cd $(T3EMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3EMPI_CFLAGS)" "LD_FLAGS=$(T3EMPI_LD_FLAGS)"\
		"MDL=$(T3EMPI_MDL)" "XOBJ=$(T3EMPI_XOBJ)" "LIBMDL=$(T3EMPI_LIBMDL)"

ksr:
	cd $(KSR_MDL); make
	make $(EXE) "CFLAGS=$(KSR_CFLAGS)" "LD_FLAGS=$(KSR_LD_FLAGS)"\
		"MDL=$(KSR_MDL)" "XOBJ=$(KSR_XOBJ)" "LIBMDL=$(KSR_LIBMDL)"

charm:
	cd $(CHARM_MDL); make CC=$(CHARM) "CFLAGS=$(CHARM_MDL_CFLAGS)"
	make $(EXE) CC=$(CHARM) "CFLAGS=$(CHARM_CFLAGS)" "LD_FLAGS=$(CHARM_LD_FLAGS)"\
		"MDL=$(CHARM_MDL)" "XOBJ=$(CHARM_XOBJ)" "LIBMDL=$(CHARM_LIBMDL)"

$(EXE): $(OBJ) $(XOBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o $@ $(OBJ) $(XOBJ) $(LIBMDL)

$(OBJ) $(EXTRA_OBJ): Makefile

# DO NOT DELETE

cosmo.o: runge.h cosmo.h
ewald.o: ewald.h pkd.h floattype.h cooling.h meval.h qeval.h
fdl.o: htable.h fdl.h
grav.o: pkd.h floattype.h cooling.h grav.h meval.h qeval.h
htable.o: htable.h
integration.o: linalg.h floattype.h integration.h
linalg.o: linalg.h floattype.h
main.o: master.h param.h pst.h pkd.h floattype.h cooling.h smoothfcn.h
main.o: parameters.h cosmo.h outtype.h
master.o: master.h param.h pst.h pkd.h floattype.h cooling.h smoothfcn.h
master.o: parameters.h cosmo.h tipsydefs.h opentype.h fdl.h htable.h
master.o: outtype.h
millerscalo.o: millerscalo.h
outtype.o: pkd.h floattype.h cooling.h outtype.h
outurb.o: pkd.h outurb.h
param.o: param.h
pkd.o: pkd.h floattype.h cooling.h ewald.h grav.h walk.h opentype.h
pkd.o: tipsydefs.h dumpframe.h
pst.o: pst.h pkd.h floattype.h cooling.h smoothfcn.h starform.h feedback.h
pst.o: outtype.h smooth.h dumpframe.h
romberg.o: floattype.h
smooth.o: smooth.h pkd.h floattype.h cooling.h smoothfcn.h
smoothfcn.o: smoothfcn.h pkd.h floattype.h cooling.h
startime.o: floattype.h startime.h
supernova.o: pkd.h floattype.h cooling.h feedback.h supernova.h startime.h
supernova.o: millerscalo.h supernovaia.h
supernovaia.o: supernova.h startime.h millerscalo.h feedback.h pkd.h
supernovaia.o: floattype.h cooling.h supernovaia.h
dumpframe.o: dumpframe.c dumpframe.h
dumpvoxel.o: dumpvoxel.c dumpframe.h dumpvoxel.h
dffuncs.o: dffuncs.c dumpframe.h
treezip.o: treezip.c treezip.h treezipkey.h treeziptypes.h
writepng.o: writepng.c writepng.h
walk.o: walk.h pkd.h floattype.h cooling.h
rotbar.o: master.h rotbar.h
cooling_metal.o: cooling_metal.h cooling_metal.c
cooling_grackle.o: cooling_grackle.h cooling_grackle.c
log.o: log.h
