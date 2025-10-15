# might be looking for the compile-time Makefile options of the code...
#
# They have moved to a separate file.
#
# To build the code, do the following:
#
#  (0) Run "make check" for a defaults compilation test.
#      Run "make test"  for a defaults compilation+running test (needs MPI environment..).
#      You will likely want to use "make SYSTYPE=... check" -- scroll down the comments for this.
#
#  (1) Copy the file "Template-Config.sh"  to  "Config.sh"
#
#        cp Template-Config.sh Config.sh
#
#  (2) Edit "Config.sh" as needed for your application
#
#  (3) Run "make"
#
#
#  New compile-time options should be added to the
#  file "Template-Config.sh" only. Usually, the should be added
#  there in the disabled/default version.
#
#  "Config.sh" should *not* be checked in to the repository
#
#  Note: It is possible to override the default name of the
#  Config.sh file, if desired, as well as the name of the
#  executable. For example:
#
#   make  CONFIG=MyNewConf.sh  EXEC=Gadget_new
#
#-----------------------------------------------------------------
#
# You might also be looking for the target system SYSTYPE option
#
# It has also moved to a separate file.
#
# To build the code, do the following:
#
# (A) set the SYSTYPE variable in your .bashrc (or similar file):
#
#        e.g. export SYSTYPE=Magny
# or
#
# (B) set SYSTYPE in Makefile.systype
#     This file has priority over your shell variable.:
#
#    (1) Copy the file "Template-Makefile.systype"  to  "Makefile.systype"
#
#        cp Template-Makefile.systype Makefile.systype
#
#    (2) Uncomment your system in  "Makefile.systype".
#
# If you add an ifeq for a new system below, also add that systype to
# Template-Makefile.systype


# Set default values for compiler options if no systype options are given or found
CC        = mpiCC                           # sets the C-compiler (default)
CXX       = mpiCC                           # sets the C++-compiler (default)
OPTIMIZE  = -std=c++11 -Wall -g -O2
MPICHLIB  = -lmpich
SWITCHES =

ifneq (,$(wildcard ./Makefile.conf))
$(info Including configure-generated Makefile.conf.)
SYSTYPE=AUTOCONF
include Build/Makefile.conf
RLIBS=$(LDFLAGS)
else
# $(warning aiuto)
endif

ifdef SYSTYPE
SYSTYPE := $(SYSTYPE)
#include Build/Makefile.systype
else
include Build/Makefile.systype
endif

ifeq ($(wildcard Makefile.systype), Makefile.systype)
INCL = Build/Makefile.systype
else
INCL =
endif


CONFIG   =  Config.sh
PERL     =  /usr/bin/perl

RESULT     := $(shell CONFIG=$(CONFIG) PERL=$(PERL) make -f config-makefile)
CONFIGVARS := $(shell cat gadgetconfig.h)



ifeq (SOFTDOUBLEDOUBLE,$(findstring SOFTDOUBLEDOUBLE,$(CONFIGVARS)))
CC       =   mpiCC     # default C++ compiler
OPTIMIZE =   -g
OPT     +=  -DX86FIX   # only needed for 32-bit intel/amd systems
endif

ifeq (HAVE_HDF5,$(findstring HAVE_HDF5,$(CONFIGVARS)))
HDF5_LIB_NAMES = -lhdf5
endif

ifeq (LT_METAL_COOLING_WAL,$(findstring LT_METAL_COOLING_WAL,$(CONFIGVARS)))
HDF5_LIB_NAMES = -lhdf5
endif

-include Build/Makefile.local

#----------------------------------------------------------------------------------------------

include Build/Makefile.$(SYSTYPE)

#----------------------------------------------------------------------------------------------

ifeq (PMGRID,$(findstring PMGRID,$(CONFIGVARS)))                              # need fftw ?
  FFTW_LIBNAMES =  -lfftw3_mpi -lfftw3 -lfftw3f_mpi -lfftw3f
  $(info $$OPTIMIZE is [${OPTIMIZE}])
ifeq (openmp,$(findstring openmp,$(OPTIMIZE)))                                # need OpenMP fftw ?
  FFTW_LIBNAMES +=  -lfftw3_omp -lfftw3f_omp
endif
endif

PM_OBJS  =  ImportFromL3_combined/pm_util.o ImportFromL3_combined/pm_powerspec.o ImportFromL3_combined/pm_density.o ImportFromL3_combined/pm_init.o ImportFromL3_combined/pm_readout.o ImportFromL3_combined/pm.o


# forcetree.o appears in two switches!

ifeq (AR_INCLUDE_OPENACC_FOLDER,$(findstring AR_INCLUDE_OPENACC_FOLDER,$(CONFIGVARS))) # OpenACC
	HYDRO_OBJS += OpenACC/gpuallvars_branch.o
endif

ifeq (ACC_NGB,$(findstring ACC_NGB,$(CONFIGVARS))) # OpenACC
 HYDRO_OBJS += OpenACC/ngb.o
else

	HYDRO_OBJS += Gravity/ngb.o
endif

OPTIONS = $(OPTIMIZE) $(OPT)

FOPTIONS = $(OPTIMIZE) $(FOPT)

EXEC   = P-Gadget3

OBJS  =  Hydro/density.o Hydro/smooth_simple.o \
         Gravity/gravtree.o Gravity/forcetree.o Gravity/darkenergy.o Gravity/accel.o Gravity/forcetree_update.o \
         Gravity/potential.o \
         Integrator/driftfac.o Integrator/predict.o Integrator/timestep.o Integrator/kicks.o \
	 CodeBase/peano.o CodeBase/endrun.o CodeBase/global.o CodeBase/begrun.o \
         CodeBase/run.o CodeBase/allvars.o CodeBase/init.o CodeBase/domain.o \
	 IO/read_ic.o IO/read_ic_cluster.o IO/second_order.o IO/restart.o IO/io.o \
         ImportFromL3/parallel_sort_l3.o\
	 System/system.o System/allocate.o System/sizelimited_sendrecv.o System/checksummed_sendrecv.o \
         System/parallel_sort_special.o System/mpi_util.o System/mymalloc.o System/calc_checksum.o \
         System/Md5.o System/pinning.o \
	 compile_time_info.o

INCL    += Makefile  CodeBase/allvars.h CodeBase/proto.h Gravity/forcetree.h CodeBase/domain.h \
	 System/dd.h System/myqsort.h System/tags.h System/vector.h Hydro/kernel.h CodeBase/switches.h \
         Blackholes/blackhole_begrun.h System/communication.h GreenTree_AR/ar_density.c \
	 GreenTree_AR/ar_hydra.c GreenTree_AR/ar_conduction.c GreenTree_AR/recompose_domain.c
# cooling_begrun.h sfr_begrun.h chemistry_begrun.h muppi_begrun.h sfr_lt_begrun.h


ifeq (COOLING,$(findstring COOLING,$(CONFIGVARS)))
OBJS	+= CoolingSfr/cooling.o
INCL    += CoolingSfr/cooling.h
endif

ifeq (SFR,$(findstring SFR,$(CONFIGVARS)))
OBJS	+= CoolingSfr/snIa_heating.o CoolingSfr/sfr_eff.o
#INCL    += CoolingSfr/sfr_eff.h
endif

ifeq (BLACK_HOLES,$(findstring BLACK_HOLES,$(CONFIGVARS)))
OBJS	+= Blackholes/blackhole.o Blackholes/blackhole_feedback_loop.o Blackholes/blackhole_swallow_loop.o
INCL    += Blackholes/blackhole.h
endif

ifeq (SUBFIND,$(findstring SUBFIND,$(CONFIGVARS)))
OBJS	+= FofSubfind/subfind.o FofSubfind/subfind_vars.o FofSubfind/subfind_collective.o FofSubfind/subfind_serial.o FofSubfind/subfind_so.o \
           FofSubfind/subfind_cont.o FofSubfind/subfind_distribute.o FofSubfind/subfind_findlinkngb.o FofSubfind/subfind_nearesttwo.o FofSubfind/subfind_loctree.o \
           FofSubfind/subfind_alternative_collective.o FofSubfind/subfind_reshuffle.o FofSubfind/subfind_potential.o FofSubfind/subfind_density.o
INCL    += FofSubfind/subfind.h
endif

ifeq (FOF,$(findstring FOF,$(CONFIGVARS)))
OBJS	+= FofSubfind/fof.o
INCL    += FofSubfind/fof.h
endif

ifeq (PMGRID,$(findstring PMGRID,$(CONFIGVARS)))
OBJS	+= $(PM_OBJS) Gravity/longrange.o
endif

ifeq (PHIDOT,$(findstring PHIDOT,$(CONFIGVARS)))
OBJS	+= Gravity/phidot.o
endif

ifeq (FORCETEST,$(findstring FORCETEST,$(CONFIGVARS)))
OBJS	+= Gravity/gravtree_forcetest.o
endif

ifeq (CONDUCTION,$(findstring CONDUCTION,$(CONFIGVARS)))
OBJS	+= Hydro/conduction.o
endif

OBJS  += GreenTree_AR/ar_ngb.o

#ifeq (RECOMPOSE_DOMAIN,$(findstring RECOMPOSE_DOMAIN,$(CONFIGVARS)))
#OBJS    += GreenTree_AR/recompose_domain.o
#endif

#ifeq (AR_XMAS_TREE,$(findstring AR_XMAS_TREE,$(CONFIGVARS)))
#OBJS  += GreenTree_AR/ar_forcetree.o
#endif

ifeq (USE_MESHLESS_FINITE_MASS,$(findstring USE_MESHLESS_FINITE_MASS,$(CONFIGVARS)))
OBJS += Mfm/gradients.o Mfm/mfm_fluxes.o Mfm/mfm_limiters.o
OBJS += Mfm/ExactRiemannSolver.o Mfm/FluidVector.o Mfm/FluxSolver.o
OBJS += Mfm/HllRiemannSolver.o Mfm/RiemannSolver.o
INCL += Mfm/Eos.hpp Mfm/FluidVector.hpp Mfm/FluxSolver.hpp
INCL += Mfm/Matrix.hpp Mfm/RiemannSolver.hpp Mfm/Vector.hpp
TEST_OBJS += tests/UnitTests/EosTest.o
TEST_OBJS += tests/UnitTests/FluidVectorTest.o
TEST_OBJS += tests/UnitTests/FluxSolverTest.o
TEST_OBJS += tests/UnitTests/MatrixTest.o
TEST_OBJS += tests/UnitTests/RiemannSolverTest.o
TEST_OBJS += tests/UnitTests/VectorTest.o
else
HYDRO_OBJS += Hydro/hydra.o
endif


ifeq (LT_METAL_COOLING,$(findstring LT_METAL_COOLING,$(CONFIGVARS)))
OBJS  += CoolingSfr/Sfr_LT/lt_wal_cooling.o CoolingSfr/Sfr_LT/lt_io.o CoolingSfr/Sfr_LT/lt_utils.o CoolingSfr/ImfLib/libimf.o CoolingSfr/ImfLib/libimf_vars.o
INCL  += CoolingSfr/Sfr_LT/lt_sn.h CoolingSfr/Sfr_LT/lt.h CoolingSfr/Sfr_LT/lt_error_codes.h CoolingSfr/Sfr_LT/lt_wal_cooling.h
ifeq (LT_SEv_INFO,$(findstring LT_SEv_INFO,$(CONFIGVARS)))
OBJS  += CoolingSfr/Sfr_LT/lt_sfr.o
else
ifeq (LT_ZAGE,$(findstring LT_ZAGE,$(CONFIGVARS)))
OBJS  += CoolingSfr/Sfr_LT/lt_sfr.o
else
OBJS  += CoolingSfr/Sfr_LT/lt_sfr_light.o
endif
endif
ifeq (LT_EJECTA_IN_HOTPHASE,$(findstring LT_EJECTA_IN_HOTPHASE,$(CONFIGVARS)))
OBJS  += CoolingSfr/Sfr_LT/lt_sn.o
else
OBJS  += CoolingSfr/Sfr_LT/lt_sn_light.o
endif
endif

ifeq (GM_MUPPI,$(findstring GM_MUPPI,$(CONFIGVARS))) # add MUPPI
OBJS += CoolingSfr/Muppi/sfr_muppi.o CoolingSfr/Muppi/muppi_communications.o CoolingSfr/Muppi/muppi_normalize.o
endif


ifeq (LMB_SPECTRAL_CRs,$(findstring LMB_SPECTRAL_CRs,$(CONFIGVARS))) # add spectral cr part
OBJS += CosmicRays/lmb_cr_init.o \
		CosmicRays/lmb_cosmic_rays.o \
		CosmicRays/lmb_cr_conduction.o \
		CosmicRays/lmb_cr_tests.o \
		CosmicRays/lmb_cr_debugging_tools.o

INCL += CosmicRays/lmb_cr_datastructs.h \
		CosmicRays/lmb_cr_e_n_integrals.h \
		CosmicRays/lmb_cr_debugging_tools.h \
		CosmicRays/lmb_cr_conduction.h \
		CosmicRays/lmb_cr_slope_solver.h \
		CosmicRays/lmb_cr_update_functions.h \
		CosmicRays/lmb_cr_init.h \
		CosmicRays/lmb_cr_adiabatic_changes.h \
		CosmicRays/lmb_cr_radiative_changes.h \
		CosmicRays/lmb_cr_injection.h \
		CosmicRays/lmb_cr_Dpp.h \
		CosmicRays/lmb_cr_timestep.h \
		CosmicRays/lmb_cr_tests.h
endif


ifeq (JD_DPP,$(findstring JD_DPP,$(CONFIGVARS))) # add JD cr part
OBJS += CosmicRays/cr_electrons.o
endif

ifeq (UM_CHEMISTRY,$(findstring UM_CHEMISTRY,$(CONFIGVARS)))
INCL    += Chemistry/chemistry.h
OBJS += Chemistry/um_chemistry_noneq.o
endif

ifeq (UM_METAL_COOLING,$(findstring UM_METAL_COOLING,$(CONFIGVARS)))
OBJS += Chemistry/um_metal_cooling.o
endif

ifeq (ADAPTGRAVSOFT,$(findstring ADAPTGRAVSOFT,$(CONFIGVARS)))
INCL += Gravity/AdaptGravSoft/ads_ngb.c Gravity/AdaptGravSoft/ads_density.c
endif

CFLAGS = $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) $(HDF5INCL) $(GMP_INCL)


ifeq (ALTERNATIVE_PSORT,$(findstring ALTERNATIVE_PSORT,$(CONFIGVARS)))
OBJS  += System/fof_alt_psort.o System/psort-1.0/error_handling.o
endif

LIBS   = -lm $(HDF5LIB)  $(HDF5_LIB_NAMES) -lz -g $(MPICHLIB) $(GSL_LIBS) -lgsl -lgslcblas $(FFTW_LIBS) $(FFTW_LIBNAMES)

ifeq (NVTX,$(findstring NVTX,$(CONFIGVARS)))
OPT += -DUSE_NVTX
NVTXLIB = -L$(CUDA_ROOT)/lib64 -lnvToolsExt
LIBS += $(NVTXLIB)
endif

ifeq (GTEST,$(findstring GTEST,$(CONFIGVARS)))
LIBS   += $(GTESTLIBS)
CFLAGS += $(GTESTINCL)
OBJS   += $(TEST_OBJS) tests/UnitTests/TestMain.o
EXEC   = TestGadget
else
OBJS   += CodeBase/main.o
endif


OBJS	+= $(HYDRO_OBJS)


CXXFLAGS = $(CFLAGS)


$(EXEC): $(OBJS)
	$(CC) $(SWITCHES) $(OPTIMIZE) $(OBJS) $(FOBJS) $(LIBS) $(RLIBS) -o $(EXEC)


$(OBJS):  $(INCL)  $(CONFIG)  compile_time_info.c


$(FOBJS): $(FINCL)

compile_time_info.c: $(CONFIG)
	$(PERL) prepare-config.perl $(CONFIG)

.PHONY: clean
clean:
	rm -f $(OBJS) $(FOBJS) $(TEST_OBJS) $(EXEC) gadgetconfig.h

.PHONY: check
check:
	OPENGADGET3_TEST_PREAMBLE_COMMAND="set -e" ./tests/run_all_config_tests.sh

.PHONY: test
test:
	OPENGADGET3_TEST_PREAMBLE_COMMAND="set -e" ./run_all_tests.sh

.PHONY: tests
tests: test
