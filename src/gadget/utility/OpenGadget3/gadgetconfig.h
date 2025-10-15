#ifdef GADGET3_IO_LIB

#define ASMTH 1.25
#define AUTO_SWAP_ENDIAN_READIC 
#define DOUBLEPRECISION 3
//#define HAVE_HDF5 

#ifdef GADGET_READ_ID64
#define LONGIDS
#endif

#define MOREPARAMS 
#define MULTIPLEDOMAINS 2
#define MYSORT 
#define NOTEST_FOR_IDUNIQUENESS 
#define NO_ISEND_IRECV_IN_DOMAIN 
#define NO_ISEND_IRECV_IN_PM 
#define PEANOHILBERT 
//#define PEDANTIC_CHECK 
#define PERIODIC 
//#define PMGRID 2048
#define TOPNODEFACTOR 128
#define WALLCLOCK 
#define NOCALLSOFSYSTEM 
//#define CONFIG_FILE_NAME "../scratch/L500N2048/Config.sh"

#define BLACK_HOLES
//#define COOLING
////#define CHEMISTRY //TODO
//#define UM_CHEMISTRY
//#define STELLARAGE
//#define GM_MUPPI
//#define RESCALEVINI
//#define JD_VTURB
//#define JD_DECOMPOSE_VTURB
//#define JD_DPP
//#define MAGNETIC
//#define TRACEDIVB
//#define MAGNETICZERO
//#define DIVBCLEANING_DEDNER
//#define JD_SHOCK
//#define READ_MACH
// TODO LT_*

#ifndef GADGET_MAX_HSML

#define LT_STELLAREVOLUTION
#define LT_NMet 16
#define LT_PM_LIFETIMES
#define LT_STARS_GUESSHSML
#define LT_WIND_VELOCITY 350.0

#define LT_ZAGE
#define LT_ZAGE_LLV
//#define LT_TRACK_CONTRIBUTES
#define LT_METAL_COOLING_WAL

#define WINDS
#endif

//#define GL_CR_DUST

//#define VERBOSE 3
#define DISABLE_MEMORY_MANAGER 1

#else

#define ASMTH 1.25
#define AUTO_SWAP_ENDIAN_READIC 
#define DOUBLEPRECISION 3
#define HAVE_HDF5 
#define LONGIDS 
#define MOREPARAMS 
#define MULTIPLEDOMAINS 2
#define MYSORT 
#define NOTEST_FOR_IDUNIQUENESS 
#define NO_ISEND_IRECV_IN_DOMAIN 
#define NO_ISEND_IRECV_IN_PM 
#define PEANOHILBERT 
#define PEDANTIC_CHECK 
#define PERIODIC 
#define PMGRID 2048
#define TOPNODEFACTOR 128
#define WALLCLOCK 
#define NOCALLSOFSYSTEM 
#define CONFIG_FILE_NAME "../scratch/L500N2048/Config.sh"

#endif