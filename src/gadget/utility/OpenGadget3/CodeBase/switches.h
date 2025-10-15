
#ifndef __SWITCHES_H__
    #define __SWITCHES_H__

    #ifndef AR_GREEN_TREE_HSML_FRACTION
    #define AR_GREEN_TREE_HSML_FRACTION 0.5
    #endif

    #ifndef DENSITY_LESS_NGB_FACTOR
    #define DENSITY_LESS_NGB_FACTOR 1
    #endif

    #ifndef KD_BUFFER_MANAGEMENT
    #define KD_BUFFER_MANAGEMENT 0.1
    #endif

    #ifdef ACC
    #define ACC_RUN
    #define ACC_PROTO
    #define ACC_GRAVITY
    #define ACC_FORCETREE

    #ifndef ACC_ACTIVEPART_CPU
    #define ACC_ACTIVEPART_CPU 1000
    #endif

    #define ACC_ACTIVEPART_II_CPU ACC_ACTIVEPART_CPU



    #define AR_FIX_ASSIGN

    #ifdef LT_STELLAREVOLUTION
    #define AR_FIX_TIME_DEP_ART_COND
    #endif

    #define ACC_NGB
    #define ACC_DENSITY
    #define ACC_HYDRA
    #ifndef ACC_GPU_NGB
    #define ACC_GPU_NGB 32
    #endif

    #ifndef ACC_GPU_NGB_CONDUCTION
    #define ACC_GPU_NGB_CONDUCTION ACC_GPU_NGB
    #endif

    #define ACC_CONDUCTION
    #endif

    #if defined(ONEDIM)
        #define NUMDIMS 1
    #elif defined(TWODIMS)
        #define NUMDIMS 2
    #else
        #define NUMDIMS 3
    #endif

    #define HYDRO_SPH 1
    #define HYDRO_PESPH 2
    #define HYDRO_MFM 3

    #if defined(USE_MESHLESS_FINITE_MASS)
        #define GADGET_HYDRO HYDRO_MFM
    #elif defined(USE_PRESSURE_ENTROPY_SPH)
        #define GADGET_HYDRO HYDRO_PESPH
    #else
        #define GADGET_HYDRO HYDRO_SPH
    #endif

    #define LIMITER_GIZMO 0
    #define LIMITER_TVD_SCALAR 1
    #define LIMITER_SCALAR 2
    #define LIMITER_SPRINGEL_2009 3
    #define LIMITER_NULL 4
    #define LIMITER_ZERO_SLOPES 5

    #if defined(SLOPE_LIMITER_GIZMO)
        #define MFM_SLOPE_LIMITER LIMITER_GIZMO

        #if !defined(SLOPE_LIMITER_TOLERANCE)
            #if defined(AGGRESSIVE_SLOPE_LIMITERS)
                #define SLOPE_LIMITER_TOLERANCE 2
            #else
                #define SLOPE_LIMITER_TOLERANCE 1
            #endif
        #endif /*SLOPE_LIMITER_TOLERANCE */

    #elif defined(SLOPE_LIMITER_TVD_SCALAR)
        #define MFM_SLOPE_LIMITER LIMITER_TVD_SCALAR
    #elif defined(SLOPE_LIMITER_SCALAR)
        #define MFM_SLOPE_LIMITER LIMITER_SCALAR
    #elif defined(SLOPE_LIMITER_SPRINGEL_2009)
        #define MFM_SLOPE_LIMITER LIMITER_SPRINGEL_2009
    #elif defined(SLOPE_LIMITER_NULL)
        #define MFM_SLOPE_LIMITER LIMITER_NULL
    #else
        #define MFM_SLOPE_LIMITER LIMITER_ZERO_SLOPES
    #endif

    #if !defined(HLLC_RIEMANN_VERSION)
        #define HLLC_RIEMANN_VERSION 3
    #endif

    #if (!defined(ROE_AVERAGE) && !defined(PVRS_APPROXIMATION) && !defined(TRRS_APPROXIMATION))
        #define ROE_AVERAGE		// other choice as standard???
    #endif

    #define SFR_NONE 100
    #define SFR_SPRINGEL 101
    #define SFR_TORNATORE 102
    #define SFR_MURANTE 103

    #if defined(GM_MUPPI)
        #define GADGET_SFR SFR_MURANTE
    #elif defined(LT_STELLAREVOLUTION)
        #define GADGET_SFR SFR_TORNATORE
    #elif defined(SFR)
        #define GADGET_SFR SFR_SPRINGEL
    #else
        #define GADGET_SFR SFR_NONE
    #endif

    #define COOLING_KATZ 101
    #define COOLING_SUTHERLAND 102
    #define COOLING_WAL 103
    #define COOLING_LUEDERS 104

    #if defined(SL_COOLING_TABLES)
        #define GADGET_COOLING COOLING_LUEDERS
    #elif defined(LT_METAL_COOLING_WAL)
        #define GADGET_COOLING COOLING_WAL
    #elif defined(LT_METAL_COOLING)
        #define GADGET_COOLING COOLING_SUTHERLAND
    #elif defined(COOLING)
        #define GADGET_COOLING COOLING_KATZ
    #endif

    #if defined(fSIDM) || defined(rSIDM)
        #ifndef NEW_KERNEL
            #define NEW_KERNEL
        #endif
        
        #if defined(mSIDM_VDEP0)
            #define mSIDM_VDEP
        #endif

        #ifdef fSIDM_COMOTEST
            #ifndef fSIDM_NORANDCOMP
                #define fSIDM_NORANDCOMP
            #endif
        #endif
    #endif

    #ifdef NEW_FILTERS
        #define N_BINS_METAL 6
        #define LOG_AGE_BINS 220
    #else
        #define N_BINS_METAL 7
        #define LOG_AGE_BINS 219
    #endif

        #define CB07_TOT_N_FILERS 90
#endif

#ifdef RT_FLUID_TEST
    #define CONSTANT_GRAVITY_Y -0.5
#endif

#ifndef WAKEUP
    #define WAKEUP 3.
#endif
