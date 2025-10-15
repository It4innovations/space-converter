#!/bin/bash            # this line only there to enable syntax highlighting in this file

##################################################
#  Enable/Disable compile-time options as needed #
##################################################

#--------------------------------------- Basic operation mode of code
#PERIODIC
#COOLING
#SFR
#SINKS
#UNEQUALSOFTENINGS
#NUM_THREADS=4              # Now OpenMP works the same, so don't compile with OpenMP *and* PTHREADS !


#--------------------------------------- Kernel Options
#TWODIMS                    # Switch for 2D test problems
#ONEDIM                     # Switch for 1D test problems

#QUINTIC_KERNEL             # Implementation of the Morris 1996 quintic spline kernel, requires (3/2)^3 more neighbours !

#WENDLAND_C4_KERNEL         # Dehnen & Aly 2012, DesNumNgb=200
#WC4_BIAS_CORRECTION        # subtract self contribution, recommended 


#--------------------------------------- TreePM Options
#PMGRID=256
#GRIDBOOST=2
#ASMTH=1.25
#RCUT=4.5
#PLACEHIGHRESREGION=55
#ENLARGEREGION=1.1


#--------------------------------------- Multi-Domain and Top-Level Tree options
MULTIPLEDOMAINS=4
#TOPNODEFACTOR=3.0
#KD_HMAX_ESTIMATE           # Alternative way to update HMAX within Tree nodes

#--------------------------------------- Things that are always recommended
PEANOHILBERT
WALLCLOCK
MYSORT
#AUTO_SWAP_ENDIAN_READIC        # Enables automatic ENDIAN swapping for reading ICs
#WRITE_KEY_FILES                # Enables writing key index files
#WRITE_INFO_BLOCK               # Enables writing the INFO block
#PROCESS_TIMES_OF_OUTPUTLIST    # Chooses the outputtime closest to any global step
#SYNCRONIZ_OUTPUT               # Writes output only at global time steps

#---------------------------------------- Single/Double Precision
#DOUBLEPRECISION
#DOUBLEPRECISION_FFTW
#OUTPUT_IN_DOUBLEPRECISION # snapshot files will be written in double precision
#INPUT_IN_DOUBLEPRECISION


#---------------------------------------- Invariance Test
#INVARIANCETEST
#INVARIANCETEST_SIZE1=2
#INVARIANCETEST_SIZE2=6
#FLTROUNDOFFREDUCTION      # enables (expensive!) `double-double' round-off reduction in particle sums
#SOFTDOUBLEDOUBLE          # needs to be set if a C++ software implementation of 128bit double-double precision should be used


#---------------------------------------- On the fly FOF groupfinder 
#FOF                                # enable FoF output
#FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
#FOF_SECONDARY_LINK_TYPES=1+16+32	 # 2^type for the types linked to nearest primaries
#FOF_GROUP_MIN_LEN=32               # default is 32
#SUBFIND                            # enables substructure finder
#DENSITY_SPLIT_BY_TYPE=1+2+16+32    # 2^type for whch the densities should be calculated seperately  
#MAX_NGB_CHECK=3                    # Max numbers of neighbours for sattlepoint detection (default = 2)
#SAVE_MASS_TAB                      # Saves the an additional array with the masses of the different components
#SUBFIND_SAVE_PARTICLELISTS         # Saves also phase-space and type variables parallel to IDs
#SO_VEL_DISPERSIONS                 # computes velocity dispersions for as part of FOF SO-properties
#ORDER_SNAPSHOTS_BY_ID
#SAVE_HSML_IN_IC_ORDER              # will store the hsml-values in the order of the particles in the IC file
#ONLY_PRODUCE_HSML_FILES            # only carries out density estimate
#KEEP_HSML_AS_GUESS                 # keep using hsml for gas particles in subfind_density
#LINKLENGTH=0.16                    # Linkinglength for FoF (default=0.2)
#NO_GAS_CLOUDS                      # Do not accept pure gaseous substructures
#WRITE_SUB_IN_SNAP_FORMAT           # Save subfind results in snap format
#LT_ADD_GAL_TO_SUB=12                # Adds optical luminosities in 6 bands to subhalos
#DUSTATT=11                         # Includes dust attenuation into the luminosity calculation (using 11 radial bins)
#OBSERVER_FRAME                     # If defined, use CB07 Observer Frame Luminosities, otherwise CB07 Rest Frame Luminosities
#SO_BAR_INFO                        # Adds temperature, Lx, bfrac, etc to Groups
#SUBFIND_COUNT_BIG_HALOS=1e4        # Adds extra blocks for Halos with M_TopHat > SUBFIND_COUNT_BIG_HALOS
#KD_CHOOSE_PSUBFIND_LIMIT           # Increases the limit for the parallel subfind to the maximum possible
#KD_ALTERNATIVE_GROUP_SORT          # Alternative way to sort the Groups/SubGroupe before writing
#SUBFIND_READ_FOF
#SUBFIND_COLLECTIVE_STAGE1
#SUBFIND_COLLECTIVE_STAGE2
#SUBFIND_ALTERNATIVE_COLLECTIVE
#SUBFIND_RESHUFFLE_CATALOGUE
#SUBFIND_RESHUFFLE_CATALOGUE_WITH_VORONOI
#SUBFIND_RESHUFFLE_AND_POTENTIAL    #needs -DSUBFIND_RESHUFFLE_CATALOGUE and COMPUTE_POTENTIAL_ENERGY
#SUBFIND_DENSITY_AND_POTENTIAL       #only calculated density and potential and write them into snapshot 


#--------------------------------------- Adaptive Gravitational Softening (F. Iannuzzi)
#ADAPTGRAVSOFT
#AGS_UPDATEALLPARTICLES
#AGS_NOCORRECTION
#AGS_OUTPUTGRAVSOFT
#AGS_OUTPUTGRAVNUMDENS
#AGS_OUTPUTZETA
#AGS_OUTPUTOMEGA
#AGS_OUTPUTCORR
#AGS_OUTPUTNGBS
#AGS_MAXSOFT=1000
#AGS_MINSOFT=0.001

#-------------------------------------- Mofified Gravity
#MODGRAV               # modified gravity master switch
#MODGRAV_FORCETEST     # for testing the Poisson solver for the scalar field
#MODGRAV_FoR           # for switching to the f(R) implementation 
#MODGRAV_HU_SAWICKI    # for Hu&Sawicki gravity

#--------------------------------------- SFR/feedback model
#GENERATIONS=1           # the number of stars a gas particle may spawn
#SOFTEREQS
#MOREPARAMS
#METALS
#STELLARAGE
#WINDS
#VARIABLE_WINDS
#QUICK_LYALPHA
#ISOTROPICWINDS
#MHM
#MODIFIED_SFR
#ALTERNATIVE_SFR

#-------------------------------------- AGN stuff
#BLACK_HOLES             # enables Black-Holes (master switch)
#BONDI                   # Bondi-Hoyle style accretion model
#BH_VARIABLE_ACCRETION_FACTOR   # variable-alpha model as in Booth&Schaye 2009
#ENFORCE_EDDINGTON_LIMIT # put a hard limit on the maximum accretion rate
#BH_THERMALFEEDBACK      # couple a fraction of the BH luminosity into surrounding gas
#BH_DRAG                 # Drag on black-holes due to accretion
#SWALLOWGAS              # Enables stochastic accretion of gas particles consistent with growth rate of hole
#EVALPOTENTIAL           # computes gravitational potential
#REPOSITION_ON_POTMIN    # repositions hole on potential minimum (requires EVALPOTENTIAL)
#BH_COUNTPROGS	      # carries a counter for each BH that gives the total number of seeds that merged into it
#BH_USE_GASVEL_IN_BONDI  # only when this is enabled, the surrounding gas velocity is used in addition to the sounds speed in the Bondi rate 
#CSND_FRAC_BH_MERGE=0.5  # Relative levocity fraction (in units of soundspeed) for merging blackholes, default=0.5
#DETACH_BLACK_HOLES      # Insert an independent data structure for BHs (currently exlicitly depends on LT_STELLAREVOLUTION)


#-------------------------------------- Some use-full adds
#KD_BHSEED_ON_POTMIN             # Seed in minimal potential instead of max density
#KD_SEED_STAR_MASS_FRACTION=0.02 # minimum star mass fraction for BH seeding 
#KD_FRICTION                     # Turns on dynamical friction for BHs
#KD_FRICTION_DYNAMIC             # Calculates environment dependent terms dynamically
#KD_IGNORE_ACCRETED_GAS_MOMENTUM # Do not transfer momentum for swallowed gas, but follow momentum for swallowed BHs
#KD_TAKE_CENTER_OF_MASS_FOR_BH_MERGER   # Do preserve center of mass when swallowing

#-------------------------------------- AGN-Bubble feedback
#BUBBLES                 # generation of hot bubbles in an isolated halo or the the biggest halo in the run
#MULTI_BUBBLES 	      # hot bubbles in all haloes above certain mass threshold (works only with FOF and without BUBBLES)
#EBUB_PROPTO_BHAR        # Energy content of the bubbles with cosmic time evolves as an integrated BHAR(z) over a Salpeter time (Di Matteo 2003 eq. [11])

#BH_BUBBLES              # calculate bubble energy directly from
                                      # the black hole accretion rate

#UNIFIED_FEEDBACK        # activates BH_THERMALFEEDBACK at high
                                      # Mdot and BH_BUBBLES FEEDBACK al low Mdot

#------------------------------------ # [options from LT and DF]
#LT_BH_ACCRETE_SLICES
#LT_BH_GUESSHSML
#LT_DF_BH
#LT_DF_BH_BHAR_SWITCH=4 	# switch feedback by BHAR

#-------------------------------------- Viscous gas treatment 
#NAVIERSTOKES            # Braginskii-Spitzer parametrization of the shear viscosity: mu = f x T^{5/2}
#NAVIERSTOKES_CONSTANT   # Shear viscosity set constant for all gas particles
#NAVIERSTOKES_BULK       # Bulk viscosity set constant for all gas particles. To run with bulk visocity only one has to set shear viscosity to zero in the parameterfile.
#VISCOSITY_SATURATION    # Both shear and bulk viscosities are saturated, so that unphysical accelerations and entropy increases are avoided. Relevant for the cosmological simulations.
#NS_TIMESTEP             # Enables timestep criterion based on entropy increase due to internal friction forces
#OUTPUTSTRESS            # Outputs diagonal and offdiagonal components of viscous shear stress tensor
#OUTPUTBULKSTRESS        # Outputs viscous bulk stress tensor
#OUTPUTSHEARCOEFF        # Outputs variable shear viscosity coefficient in internal code units


#-------------------------------------------- Things for special behaviour
#WINDTUNNEL
#TWOPOINT_ON_OUTPUT
#POWERSPEC_ON_OUTPUT
#POWERSPEC_ON_OUTPUT_EACH_TYPE
#DO_NOT_CREATE_STAR_PARTICLES
#TRADITIONAL_SPH_FORMULATION
#NOTEST_FOR_IDUNIQUENESS
#SNIA_HEATING
#FIXEDTIMEINFIRSTPHASE=1000.0
NO_ISEND_IRECV_IN_DOMAIN
FIX_PATHSCALE_MPI_STATUS_IGNORE_BUG
#MPISENDRECV_SIZELIMIT=100
#MPISENDRECV_CHECKSUM
#NOGRAVITY
#NOACCEL
#NOISMPRESSURE
#NOVISCOSITYLIMITER
#NOTREERND
#NOTYPEPREFIX_FFTW
#ISOTHERM=200                  # adds potential of an isothermal sphere
#COMPUTE_POTENTIAL_ENERGY
#LONGIDS
#ENLARGE_DYNAMIC_RANGE_IN_TIME  # NOT TESTED !!!
#ASSIGN_NEW_IDS
#INHOMOG_GASDISTR_HINT         # if the gas is distributed very different from collisionless particles, this can helps to avoid problems in the domain decomposition
#LONG_X=140
#LONG_Y=1
#LONG_Z=1
#SPH_BND_PARTICLES
#NEW_RATES                     # switches in updated cooling rates from Naoki
#RADIATIVE_RATES               # used in non-equilibrium chemistry model
#READ_HSML                     # reads hsml from IC file
#ADAPTIVE_GRAVSOFT_FORGAS      # allows variable softening length for gas particles (requires UNEQUALSOFTENINGLENGTH)
#ADAPTIVE_GRAVSOFT_FORGAS_HSML # this sets the gravitational softening for SPH particles equal to the SPH smoothing (requires ADAPTIVE_GRAVSOFT_FORGAS)
#GENERATE_GAS_IN_ICS
#SPLIT_PARTICLE_TYPE=4+8
#NEUTRINOS                     # Option for special integration of light neutrino species 
#OMIT_NEUTRINOS_IN_SNAPS
#KSPACE_NEUTRINOS             #Fourier-space integration of neutrinos
#KSPACE_NEUTRINOS_2           #Improved Fourier-space integration of neutrinos from Ali-Haimoud & Bird 2013.
#INCLUDE_RADIATION
#START_WITH_EXTRA_NGBDEV        # Uses special MaxNumNgbDeviation for starting
#ISOTHERM_EQS                  # isothermal equation of state
#NO_UTHERM_IN_IC_FILE
#SPECIAL_GAS_TREATMENT_IN_HIGHRESREGION
#DONOTUSENODELIST
#MAXHSML
#SUB_TURB_DRIVING
#ADJ_BOX_POWERSPEC         # compiles in a code module that allows via restart-flag 6 the calculation of a gas velocity power spectrum of a snapshot with an adjustable box (user defined center and size)

#---------------------------------------- New shockfinding method
#AB_SHOCK                       # Shock- and Machfinder
#AB_SHOCK_VELDIV                # Additional velocity divergence filter criterion
#AB_SHOCK_MAGNETIC              # Account for Alfvenwaves in signal velocity


#---------------------------------------- Imports from L-Gadget3 (testing phase, do not use)
#SORT_FROM_L3                   # Uses L-Gadget3 version of sort
#PM_FROM_L3                     # Uses L-Gadget3 cersion of sort 

#FFTW3
#FFTW3_THREADS
#FFTW2_THREADS

#OPENMP=4                       # Masterswitch for explicit OpenMP implementation
#OMP_MYSORT                     # Uses a OpenMP version of the self-written mergesort 
#OMP_SORT=2000	                # Replaces standard quicksort with a self-written version with OpenMP support


#--------------------------------------- Time integration options
#ALTERNATIVE_VISCOUS_TIMESTEP
#NOSTOP_WHEN_BELOW_MINTIMESTEP
#NOWINDTIMESTEPPING            # Disable wind reducing timestep (not recommended)
#NOPMSTEPADJUSTMENT
#FORCE_EQUAL_TIMESTEPS

#--------------------------------------- Output/Input options
#OUTPUTPOTENTIAL
#OUTPUT_LONGRANGE_POTENTIAL
#RECOMPUTE_POTENTIAL_ON_OUTPUT # update potential every output even it EVALPOTENTIAL is set
#OUTPUTACCELERATION
#OUTPUTCHANGEOFENTROPY
#OUTPUT_DIV_CURL
#OUTPUT_VORTICITY              # outputs the vorticity vector
#OUTPUTTIMESTEP
#OUTPUTCOOLRATE                # outputs cooling rate, and conduction rate if enabled
#HAVE_HDF5                     # needed when HDF5 I/O support is desired
#OUTPUTBSMOOTH
#OUTPUTDENSNORM
#XXLINFO                       # Enables additional output for viscosityand bfield
#OUTPUTLINEOFSIGHT             # enables on-the-fly output of Ly-alpha absorption spectra
#OUTPUTLINEOFSIGHT_SPECTRUM
#OUTPUTLINEOFSIGHT_PARTICLES

#--------------------------------------- Testing and Debugging options
#FORCETEST=0.1
DEBUG                     # enables core-dumps and FPU exceptions
#PARTICLE_DEBUG            # auxiliary communication of IDs
#VERBOSE
#CHECKSUM_DEBUG

#---------------------------------------- MCs as perturber particles to study spiral structures
#MOL_CLOUDS
#MOL_CLOUDS_LIFETIME=0.01
#--------------------------------------- SCF
#SCFPOTENTIAL               # turn SCF on/off 
#SCF_HYBRID=1               # =1:tree:stars<->stars,DM<->DM,stars->DM/SCF:stars<-DM =2:tree:stars<->stars,stars->DM/SCF:stars<-DM,DM<->DM 
#SCF_HQ_MASS=95.2401        # mass of halo of expansion basis 
#SCF_HQ_A=29.7754           # scale length of halo of expansion basis
#SCF_NEFF=5000000           # effective particle number in halo
#SCF_NMAX=1                 # maximum n for expansion cut-off
#SCF_LMAX=1                 # maximum l for expansion cut-off
#SCF_SCALEFAC               # readin scale factors for coefficients from file scf_scalefac.dat

#--------------------------------------- Static NFW Potential
#STATICNFW
#NFW_C=12
#NFW_M200=100.0
#NFW_Eps=0.01
#NFW_DARKFRACTION=0.87
#NFW_BOXCENTERED

#--------------------------------------- Static Isothermal Potential                      
#STATICISO
#ISO_M200=95.21
#ISO_R200=160.0
#ISO_Eps=0.1
#ISO_FRACTION=0.9

#--------------------------------------- Growing Disk Potential                      
#GROWING_DISK_POTENTIAL

#--------------------------------------- Disk Potential                      
#DISKPOT

#--------------------------------------- Static unit Plummer sphere (G=M=a=1)
#STATICPLUMMER   

#--------------------------------------- Static Hernquist Potential
#STATICHQ
#HQ_M200=95.2401  
#HQ_C=9.0       
#HQ_DARKFRACTION=0.9


#--------------------------------------- Static Brandt Potential
#STATICBRANDT
#BRANDT_OmegaBr=0.5
#BRANDT_R0=2.0


#--------------------------------------- Thermal conduction
#CONDUCTION
#CONDUCTION_CONSTANT
#CONDUCTION_SATURATION


#--------------------------------------- Dark energy
#DARKENERGY # Enables Dark Energy
#TIMEDEPDE  # read w(z) from a DE file
#RESCALEVINI # rescale v_ini in read_ic / read_ic_cluster
#EXTERNALHUBBLE # reads the hubble function from the DE file
#TIMEDEPGRAV # resacles H and G according to DE model
#DARKENERGY_DEBUG # enable writing of drift/kick table


#--------------------------------------- Long-range scalar field
#SCALARFIELD

#--------------------------------------- SPH turbulence on kernel-scale
#JD_VTURB				    # Compute vrms, vblk

#JD_DECOMPOSE_VTURB         # Decompose vrms into vrad and vtan
#JD_TIME_DEP_ART_VISC       # use csound instead of max signal in Dtalpha

#JD_RELAX_CLUSTERS          # keep U=U(t=t0) and set Vel=0 before every snapdump

#VSMOOTH                    # dissipate velocity noise below kernel scale


#--------------------------------------- Real Cosmic Rays
#JD_DPP				        # Compute Reacceleration Coefficient Dpp/p^2
                            #, needs -DMAGNETIC -DJD_VTURB
#JD_DPPONSNAPSHOTONLY	    # Compute Dpp ONLY for the snapshot, not for every timestep 
#KD_RESTRICT_NEIGHBOURS
#BP_REAL_CRs=1				# Number of energy bins for CRs
#BP_SEED_CRs				# enables seed CRs
#BP_REAL_CRs_ARTIFICIAL_CONDUCTIVITY	# enables Price-Monaghan artificial conductivity for CRs

#--------------------------------------- SPH viscosity options
#CONVENTIONAL_VISCOSITY     # enables the old viscosity
#TIME_DEP_ART_VISC          # Enables time dependend viscosity
#NO_SHEAR_VISCOSITY_LIMITER # Turns of the shear viscosity supression
#HIGH_ART_VISC_START        # Start with high rather than low viscosity
#ALTVISCOSITY               # enables alternative viscosity based on div(v)
#ARTIFICIAL_CONDUCTIVITY    # enables Price-Monaghan artificial conductivity

#--------------------------------------- Magnetic Field options
#MAGNETIC                   # Turns on B Field (including induction eqn) 
#HIGH_ORDER_INDUCTION       # Turns on linear error terms in induction equation. Watch energy conservation!!
#EULERPOTENTIALS            # Evolves alpha,betha instead of B
#MAGFORCE                   # Turns on B force
#MAGNETIC_SIGNALVEL         # Extend definition of signal velocity 
                                         # by the magneto sonic waves
#ALFVEN_VEL_LIMITER=10      # Limits the contribution of the Alfven waves
                                         # to the singal velocity 

#....................................... Force equation stuff
#MU0_UNITY                  # Sets \mu_0 to unity
#DIVBFORCE                  # Subtract div(B) force from M-tensor
#DIVBFORCE3=1.0

#....................................... Smoothing Options
#BSMOOTH                    # Turns on B field smoothing
#BSMOOTH_BND

#....................................... Artificial magnetic dissipation options
#MAGNETIC_DISSIPATION       # Turns on artificial magnetic dissipation 
#MAGDISSIPATION_PERPEN      # Uses only the perpendicular magnetic 
                                         # field component for dissipation (REMOVED)
#TIME_DEP_MAGN_DISP         # Enables time dependent coefficients
#HIGH_MAGN_DISP_START       # Starts from high coefficient
#ROT_IN_MAG_DIS             # Adds the RotB term in dissipation 

#....................................... DivB cleaning options
#DIVBCLEANING_DEDNER        # Turns on hyp/par cleaning (Dedner 2002)
#SMOOTH_PHI                 # Turns on smoothing of Phi

# DivBcleaningParabolicSigma = 2
# DivBcleaningHyperbolicSigma = 1
# DivBcleaningQ = 0.5

#....................................... magnetic diffusion options
#MAGNETIC_DIFFUSION         # Turns on magnetic diffusion
#MAGNETIC_DIFFUSION_LIMIT   # Limits transport distance for B
#MAGNETIC_DIFFUSION_HEAT    # Converts diffuesd B field into internal energy
#SMOOTH_ROTB                # Turns on smoothing of rot(B)

#....................................... Debugging stuff
#MAGNETIC_STATISTICS        # Outputs in statistics file for magnetic quantities
#TRACEDIVB                  # Writes div(B) into snapshot
#RELDIVB                    # Use smoothed h*divb/B in snapshot and calculations
#DBOUTPUT                   # Writes dB/dt into snapshot
#OUTPUT_ROTB                # Writes rotB into snapshot
#OUTPUT_SROTB               # Writes smoothed rotB into snapshot
#OUTPUT_XPHI                # Write x cold phase fraction into snapshot

#....................................... Initil condition stuff
#MAGNETIC_SN_SEEDING        # Supernova dipole seeding of magnetic energy. Requires enabled SFR.
#BINISET                    # Allows to set Bx,By,Bz in parameter file
#BFROMROTA                  # Allows to five vector potential instead of
                                         # B within the IC file
#IGNORE_PERIODIC_IN_ROTA    # Don't use the periodic mapping when 
                                         # calculating rot(A)
                                         # Note A might not be periodic even if B is. 
#BRIOWU                     # Extrapolate A outside simulation in case 
                                         # of Brio & Wu test problem
#....................................... Healpix stuff
#HEALPIX=1.05		 # Generates a healpix map, to retain particles in the 
					 # HiRes region. It should be a tolerance value.
#HEALPIX_INNERBOUND=2
#HEALPIX_OUTERBOUND=4

#--------------------------------------- Glass making/ 2nd-order initial conditions
#MAKEGLASS

#--------------------------------------- Raditive Transfer
#RADTRANSFER
#RADTRANSFER_FLUXLIMITER
#RADTRANSFER_MODIFY_EDDINGTON_TENSOR
#RT_COOLING_PHOTOHEATING
#RT_RAD_PRESSURE
#EDDINGTON_TENSOR_GAS
#EDDINGTON_TENSOR_STARS
#EDDINGTON_TENSOR_SFR
#EDDINGTON_TENSOR_BH
#RT_OUTPUT_ET
#RT_OUTPUT_RAD_ACCEL
#HYDROGEN_ONLY
#RT_INCLUDE_HE
#RT_SILENT

#RT_TEST_SST
#RT_MULTI_FREQUENCY
#RT_POPIII

#--------------------------------------- Fine-grained phase space structure analysis 
#DISTORTIONTENSORPS           #main switch: integrate phase-space distortion tensor 
#OUTPUT_DISTORTIONTENSORPS    #write phase-space distortion tensor to snapshot 
#OUTPUT_TIDALTENSORPS         #write configuration-space tidal tensor to snapshot 
#OUTPUT_LAST_CAUSTIC          #write info on last passed caustic to snapshot 
#GDE_TYPES=2+4+8+16+32        #track GDE for these types
#GDE_READIC                   #read initial sheet orientation/initial density/initial caustic count from ICs
#GDE_BIGFLOAT=0               #extended floating point precision -> 0: hpa lib 1: ttmath lib
#GDE_LEAN                     #lean version of GDE
#GDE_LOGOUTPUT                #log output (very slow)

#--------------------------------------- self-interacting DM
#SIDM=2                        #self-interacting particle types
#CONST_CROSS
#SIDM_DO_NOT_SCATTER
#SIDM_DTIME_FAC=0.01
#SIDM_FREEZE
#SIDM_INELASTIC

#---------------------------------------- nonequilibrium proimodal chemisitry
#NONEQUILIBRIUM
#CHEMISTRY
#CMB
#RADIATION


#---------------------------------------- Cosmic Rays (Martin)
#COSMIC_RAYS               # Cosmic Rays Master Switch
#NUMCRPOP=5                # Number of CR populations: Max is 6
#CR_IC                     # IC files contain CR information
#CR_IC_PHYSICAL
#CR_DISSIPATION            # Catastrophic losses
#CR_THERMALIZATION         # Coulomb cooling
#CR_SHOCK=2                # Shock energy is directed into CR
			                # 2 = Mach-Number dependent shocks, Mach-number derived for thermal gas/CR composite
			                # 3 = Mach-Number dependent shocks, Mach-number derived for thermal gas
#CR_DIFFUSION              # Cosmic Ray diffusion
#CR_DIFFUSION_GREEN        # alternative diffusion model
#UPDATE_PARANOIA=1         # 1 = Update on every predict, 2 = Update on every energy injection and on every predict
#CR_OUTPUT_INJECTION       # Output energy injection rate in snapshots
#CR_NO_CHANGE              # Compute changes to CR, but do not execute them, useful for estimating the size of effects
#COSMIC_RAY_TEST           # starts a test routine instead of the simulation
#CR_NOPRESSURE             # computes CRs as usual, but ignores the pressure in the dynamics
#CR_SN_INJECTION           # switches on CRs due to SNe 
#CR_BUBBLES                # CRs in the AGN bubbles
#CR_OUTPUT_TIMESCALES      # returns output from CR_ThermalizationTime and CR_DissipationTime
#CR_OUTPUT_THERMO_VARIABLES # returns output from CR_P0, CR_E0, and CR_n0

#---------------------------------------- Mach number finder (Christoph)
#MACHNUM                   # Mach number Master Switch
#MACHSTATISTIC             # Dissipated thermal energy at shocks
#CR_OUTPUT_JUMP_CONDITIONS # CR: density and thermal energy jump at shocks
#OUTPUT_PRESHOCK_CSND	# Output pre-shock sound speed and pre-shock physical density

#--------------------------------------- Cecilia 's model
#CS_MODEL      
#CS_FEEDBACK
#CS_SNI
#CS_SNII
#CS_ENRICH
#CS_TESTS
#CS_SEED                   # Use seeding for hot phase (Till)
#CS_WSS                    # Use WSS08 cooling (Till)
#CS_WSS_ADDCOMPTON         # Add Compton cooling manually and use corresponding tables (Till)
#CS_ATTEMPTS               # Stop promotion after N attempts (Till)
#CS_TILL_DEBUG             # Print debug info for CS Model (Till)



#--------------------------------------- Voronoi based SPH
#VORONOI
#EXTENDED_GHOST_SEARCH           # This extends the ghost search to the full 3x3 domain instead of the principal domain
#ALTERNATIVE_GHOST_SEARCH        # This switches on the "old" routines that find the ghost neighbours
#VORONOI_MESHOUTPUT
#GAMMA=5.0/3.0
#VORONOI_SHAPESCHEME
#VORONOI_MESHRELAX
#VORONOI_MESHRELAX_KEEPRESSURE
#MESHRELAX_DENSITY_IN_INPUT
#VORONOI_TIME_DEP_ART_VISC
#VORONOI_CFL_COND
#VORONOI_BALSARA
#VORONOI_RELAX_VIA_VISC
#VORONOI_CENTERING
#GRAVITY_CENTROID
#GRAVITY_NOT_PERIODIC
#SHIFT_BY_HALF_BOX

#--------------------------------------- degenerate equation of state
#EOS_DEGENERATE
#EOS_COULOMB_CORRECTIONS
#EOS_NSPECIES=3
#WAKEUP=4.1           /* allows 2 timestep bins within kernel */
#RELAXOBJECT

#--------------------------------------- nuclear network
#NUCLEAR_NETWORK
#FIXED_TEMPERATURE
#NEGLECT_DTDY_TERMS
#NETWORK_OUTPUT
#NETWORK_OUTPUT_BINARY


#-------------------------------------- Test 4 Christian                                                                                   
#GAMMA=1
#ISOTHERM_EQS

#-------------------------------------- Simple turbulence test
#VS_TURB
#POWERSPEC_GRID=128
#GAMMA=1.01

#AB_TURB
#TURB_DRIVING

#--------------------------------------- Luca's model

#LT_METAL_COOLING           # METAL COOLING option
#LT_METAL_COOLING_on_SMOOTH_Z
#LT_METAL_COOLING_WAL

#LT_SMOOTH_Z
#LT_SMOOTH_ALLMETALS
#LT_SMOOTHZ_DETAILS
#LT_SMOOTH_SIZE
#LT_SMOOTH_NGB
#LT_MEAN_Z_ESTIMATE

#LT_SMOOTHZ_IN_IMF_SWITCH

#LT_SMOOTH_XCLD

#....................................

#LT_STARBURSTS

#....................................# STELLAR EVOLUTION options
#LT_STELLAREVOLUTION        # enable stellar evolution
#LT_NMet=9                  # number of species

#LT_SNII                    # enable snII
#LT_SNIa                    # enable SnIa
#LT_AGB                     # enable low mass stars

#LT_PM_LIFETIMES            # use Padovani&matteucci 1993 lifetimes
#LT_MM_LIFETIMES           # use Maeder&Menet 1989 lifetimes

#LT_STARS_GUESSHSML

#LT_STOCHASTIC_IRA

#LT_LOCAL_IRA
#LT_AVOID_SELFENRICHMENT
#LT_AVOID_ENRICH_SFGAS

#LT_TEMP_THRESH_FOR_MULTIPHASE

#LT_WIND_VELOCITY=500.0   # set the winds' velocity in km/s
#LT_HYDROWINDS
#LT_DECOUPLE_POSTWINDS_FROM_SF

#LT_EJECTA_IN_HOTPHASE
#LT_SNegy_IN_HOTPHASE
#LT_HOT_EJECTA
#LT_HOT_WINDS

#LT_CharT                  # extra conditions on SF:
#                                      #  cooling time < sound cross time
#                                      #  ffall time < sound cross time

#.....................................# Other physics options
#LT_WINDS_EXT               # EXTRA WINDS options
#LT_WINDS_EXT_NOSF_BRANCH
#LT_TRACK_WINDS

#LT_USE_TOP_HAT_WEIGHT
#LT_USE_TOP_HAT_WEIGHT_STRICT
#LT_USE_SOLIDANGLE_WEIGHT

#LT_DONTUSE_DENSITY_in_WEIGHT

#.....................................#

#LT_ZAGE
#LT_ZAGE_LLV

#.....................................# logging/debug options
#LT_SEv_INFO                 # infos about enrichment
#LT_SEv_INFO_DETAILS_onSPREAD
#LT_SEv_INFO_DETAILS
#LT_EXTEGY_INFO              # infos about energy
#LT_CharT_INFO               # infos about cooling time,
#                                       #  sound crossing time,
#                                       #  free fall time
#LT_SEvDbg
#LT_TRACK_CONTRIBUTES
#LT_Nbits=8
#LT_power10_Nbits=2
#.....................................#
#LT_COST_SE                 # account for cost of Sn in
#                                       #  domain dec.

#---------------------------------------------

#UM_CHEMISTRY        #<<<<<<< enable met/HD non-eq. chemistry

#UM_METAL_COOLING    #<<<<<<< enable met cooling only
#UM_H_MET_IMPACTS
#UM_e_MET_IMPACTS
#UM_ENABLE_CII_COOLING
#UM_ENABLE_SiII_COOLING
#UM_ENABLE_OI_COOLING
#UM_ENABLE_FeII_COOLING
#UM_WIND_DELAYTIME
#UM_CONTINUE
#UM_CUT_MOLECULES

#NB:
# switch off UM_CHEMISTRY (for non-eq cooling) and use
# UM_METAL_COOLING + UM_MET_IN_LT_COOLING if you want to include low T met. cooling in Luca's  part only and neglect the non-eq options.

#UM_MET_NONEQ_CORRECTIONS
#UM_MET_IN_NONEQ_COOLING
#UM_MET_IN_LT_COOLING
#UM_HD_COOLING
#UM_NEW_RATES
#UM_NEW_RATES_HG97
#UM_CHEMISTRY_INISET  ## Use this with UM_READ_METALS_FROM_IC off

##Some other possibilities to initialize the metals..
#UM_READ_METALS_FROM_IC
#UM_CHECK
#
#	UM_CHEMISTRY: is the main opt; 
#	enables um_chemistry when CHEMISTRY is not defined
#	UM_METAL_COOLING:
#	enables low T metal cooling in cooling.c/um_metal_cooling.c
#	(use with && UM_MET/NONEQ_LT_COOLING and/or UM_MET_IN_NONEQ_COOLING)
#	UM_H_MET_IMPACTS:
#	enables line cooling due to H-impact excitations with metals.
#	UM_e_MET_IMPACTS:
#	enables line cooling due to e-impact excitations with metals.
#	Enable different kind of cooling with:
#            UM_ENABLE_CII_COOLING
#            UM_ENABLE_SiII_COOLING
#            UM_ENABLE_OI_COOLING
#            UM_ENABLE_FeII_COOLING
# 	UM_MET_IN_LT_COOLING:
#	adds low T metal cooling fct in cooling.c
#       UM_MET_IN_NONEQ_COOLING:
# 	adds low T metal cooling in um_chemistry_noneq.c
#       UM_MET_NONEQ_CORRECTIONS:
#	performs corrections in the non-equilibrium stuff
#	when coupled with metals (UM_MET_IN_NONEQ_COOLING + LT_METAL_COOLING)
#	UM_HD_COOLING
#	enable HD cooling
#	UM_H2plus_COOLING
#	enable H2+ cooling
#	UM_NEW_RATES
#	updated rates
#	UM_NEW_RATES_HG97
#	enable Hui&Gnedin case B recombination rates (for high density gas)
#	UM_CHEMISTRY_SETIN
#	set up chemical composition when specified in the parameter file
#	UM_CHECK:
#	just to make some check on the run

#--------------------------------------- TG's options
#CHEMCOOL


#--------------------------------------- Physics implementations of project Eagle

#EAGLE     # master switch

#-------------------------------------- MUPPI SF&FB model by G.Murante & P.Monaco
LT_GM_FBMODEL                    # this switches MUPPI on
LT_GM_FORCE_SINGLE_PHASE         # particle at low T and high pressure are forced into single-phase starburst
LT_GM_DISTRIBUTE_ALL_ENERGY      # forces particles to always distribute thermal energy
EXIT_MP_FOR_WIND_PARTICLES       # particles with abundant hot phase and stellar component (run-up particles) are taken out of multi-phase

LT_GM_EKIN                       # switch for distribution of kinetic energy
LT_GM_KICKS_TO_FOUNTAIN          # fountain type of kinetic feedback
#COUNT_PARTICLES_IN_CONE

#LT_GM_FBMODEL_DEBUG             # outputs feedback energy statistics on file FBEnergy.txt and on stdout
#LT_GM_FBMODEL_RECMAX            # records and writes on auxiliary snap file maximum gas data (e.g., max temperature)



###############################################################################
#
# at compile-time. From the list below, please activate/deactivate the
# options that apply to your run. If you modify any of these options,
# make sure that you recompile the whole code by typing "make clean;
# make".
#
# Main code options:
#
#     These affect the physical model that is simulated.
#
#     - PERIODIC:   Set this if you want to have periodic boundary conditions.
#     - COOLING:    This enables radiative cooling and heating. It also enables
#                   an external UV background which is read from a file.
#     - SFR:        This enables star formation using an effective multiphase
#                   models. This option requires cooling.
#     - METALS:     This model activates the tracking of enrichment in gas and
#                   stars. Note that metal-line cooling is not included yet.
#     - STELLARAGE: This stores the formation redshift of each star particle.
#     - WINDS:      This activates galactic winds. Requires star formation.
#     - ISOTROPICWINDS: This makes the wind isotropic. If not set the wind is
#                       spawned in an axial way. Requires winds to be activated.
#     - NOGRVITY:  This switches off gravity. Makes only sense for pure
#                   SPH simulations in non-expanding space.
#
# Options for SPH:
#
#     - NOFIXEDMASSINKERNEL:  If set, the number of SPH particles in the kernel
#                             is kept constant instead of the mass.
#     - NOGRADHSML:           If actived, an equation of motion without grad(h)
#                             terms is used.
#            Note: To have the default "entropy"-formulation of SPH (Springel &
#                  Hernquist), the switches NOFIXEDMASSINKERNEL and NOGRADHSML
#                  should *not* be set.
#     - NOVISCOSITYLIMITER:   If this is set, there is no explicit upper limit
#                             on the viscosity that tries to prevent particle
#                             'reflection' in case of poor timestepping.
#
# Numerical options:
#
#     - PMGRID:     This enables the TreePM method, i.e. the long-range force
#                   is computed with a PM-algoritthm, and the short range force
#                   with the tree. The parameter has to be set to the size of the
#                   mesh that should be used, (e.g. 64, 96, 128, etc). The mesh
#                   dimensions need not necessarily be a power of two.
#                   Note: If the simulation is not in a periodic box, then a FFT
#                   method for vacuum boundaries is employed, using a mesh with
#                   dimension twice that specified by PMGRID.
#     - PLACEHIGHRESREGION: If this option is set (will only work together
#                   with PMGRID), then the long range force is computed in two
#                   stages: One Fourier-grid is used to cover the whole simulation
#                   volume, allowing the computation of the large-scale force.
#                   A second Fourier mesh is placed on the region occupied by
#                   "high-resolution" particles, allowing the computation of an
#                   intermediate scale force. Finally, the force on very small
#                   scales is supplemented by the tree. This procedure can be useful
#                   for "zoom-simulations", where the majority of particles (the
#                   high-res particles) are occupying only a small fraction of the
#                   volume. To activate this option, the parameter needs to be set
#                   to an integer that encodes the particle types that represent the
#                   high-res particles in the form of a bit mask. For example, if
#                   types 0, 1, and 4 form the high-res particles, set the parameter
#                   to PLACEHIGHRESREGION=1+2+16. The spatial region covered by the
#                   high-res grid is determined automatically from the initial
#                   conditions. Note: If a periodic box is used, the high-res zone
#                   may not intersect the box boundaries.
#     - ENLARGEREGION: The spatial region covered by the high-res zone has a fixed
#                   size during the simulation, which initially is set to the
#                   smallest region that encompasses all high-res particles. Normally, the
#                   simulation will be interrupted, if high-res particles leave this
#                   region in the course of the run. However, by setting this parameter
#                   to a value larger than one, the high-res region can be expanded.
#                   For example, setting it to 1.4 will enlarge its side-length by
#                   40% (it remains centered on the high-res particles). Hence, with
#                   such a setting, the high-res region may expand or move by a
#                   limited amount. If in addition SYNCHRONIZATION is activated, then
#                   the code will be able to continue even if high-res particles
#                   leave the initial high-res grid. In this case, the code will
#                   update the size and position of the grid that is placed onto
#                   the high-resolution region automatically. To prevent that this
#                   potentially happens every single PM step, one should nevertheless
#                   assign a value slightly larger than 1 to ENLARGEREGION.
#     - DOUBLEPRECISION: This makes the code store and compute internal
#                        particle data in double precision. Note that output
#                        files are nevertheless written by converting to single
#                        precision.
#     - NOTREERND:       If this is not set, the tree construction will succeed
#                        even when there are a few particles at identical
#                        locations. This is done by `rerouting' particles once
#                        the node-size has fallen below 1.0e-3 of the softening
#                        length. When this option is activated, this will be
#                        surpressed and the tree construction will always fail
#                        if there are particles at extremely close coordinates.
#     - NOSTOP_WHEN_BELOW_MINTIMESTEP: If this is activated, the code will not
#                        terminate when the timestep falls below the value of
#                        MinSizeTimestep specified in the parameterfile. This
#                        is useful for runs where one wants to enforce a
#                        constant timestep for all particles. This can be done
#                        by activating this option, and by setting Min- and
#                        MaxSizeTimestep to an equal value.
#     - PSEUDOSYMMETRIC: When this option is set, the code will try to "anticipate"
#                        timestep changes by extrapolating the change of the
#                        acceleration into the future. This in general improves the
#                        long-term integration behaviour of periodic orbits.
#     - SYNCHRONIZATION: When this is set, particles may only increase their
#                        timestep if the new timestep will put them into
#                        synchronization with the higher time level. This typically
#                        means that only on half of the timesteps of a particle
#                        an increase of the step may occur.
#     - NOPMSTEPADJUSTMENT: When this is set, the long-range timestep for the
#                        PM force computation is always determined by MaxSizeTimeStep.
#                        Otherwise, it is set to the minimum of MaxSizeTimeStep and
#                        the timestep obtained for the maximum long-range force with
#                        an effective softening scale equal to the PM smoothing-scale.
# - LONG_X/Y/Z:
#     These options can be used together with PERIODIC and NOGRAVITY only.
#     When set, the options define numerical factors that can be used to
#     distorts the periodic simulation cube into a parallelepiped of
#     arbitrary aspect ratio. This can be useful for idealized SPH tests.
#
# - TWODIMS:
#     This effectively switches of one dimension in SPH, i.e. the code
#     follows only 2d hydrodynamics in the xy-, yz-, or xz-plane. This
#     only works with NOGRAVITY, and if all coordinates of the third
#     axis are exactly equal. Can be useful for idealized SPH tests.
#
# - SPH_BND_PARTICLES:
#     If this is set, particles with a particle-ID equal to zero do not
#     receive any SPH acceleration. This can be useful for idealized
#     SPH tests, where these particles represent fixed "walls".
#
#
# Architecture options:
#
#     - T3E:       The code assumes that sizeof(int)=4 holds. A few machines
#                  (like Cray T3E) have sizeof(int)=8. In this case, set the
#                  T3E flag.
#     - NOTYPEPREFIX_FFTW: If this is set, the fftw-header/libraries are accessed
#                  without type prefix (adopting whatever was chosen as default at compile
#                  of fftw). Otherwise, the type prefix 'd' for double is used.
#
# Input options:
#
#     - MOREPARAMS:  Activate this to allow a set of additional parameters in
#                    the parameterfile which control the star formation and
#                    feedback sector. This option must be activated when star
#                    formation is switched on.
#
# Output options:
#
#     - OUTPUTPOTENTIAL: This will force the code to compute gravitational
#                        potentials for all particles each time a snapshot file
#                        is generated. This values are then included in the
#                        snapshot file. Note that the computation of the
#                        values of the potential costs additional time.
#     - OUTPUTACCELERATION: This will include the physical acceleration of
#                        each particle in snapshot files.
#     - OUTPUTCHANGEOFENTROPY: This will include the rate of change of entropy
#                        of gas particles in snapshot files.
#     - OUTPUTTIMESTEP:  This will include an output of the timesteps actually
#                        taken by each particle.
#
# Miscellaneous options:
#
#     - PEANOHILBERT:    This is a tuning option. When set, the code will bring
#                        the particles after each domain decomposition into
#                        Peano-Hilbert order. This improves cache utilization
#                        and performance.
#     - WALLCLOCK:       If set, a wallclock timer is used by the code to
#                        measure internal time consumption (see cpu-log file).
#                        Otherwise a timer that measures consumed processor
#                        ticks is used.
#
# Debugging/testing options:
#
#     - FORCETEST:       This can be set to check the force accuracy of the
#                        code. The option needs to be set to a number between
#                        0 and 1 (e.g. 0.01), which is taken to specify a
#                        random fraction of particles for which at each
#                        timestep forces by direct summation are computed. The
#                        normal tree-forces and the "correct" direct summation
#                        forces are collected in a file. Note that the
#                        simulation itself is unaffected by this option, but it
#                        will of course run much(!) slower
#                        if FORCETEST*NumPart*NumPart >> NumPart. Note: Particle
#                        IDs must be set to numbers >=1 for this to work.
#
###############################################################################




# - FLTROUNDOFFREDUCTION     enables round off reduction in particle sums
#		             if DOUBLEPRECISION is set, these sums are done in 'long double'
#                            if single precision is used, they are done in 'double'
#                            This should in principle allow to make computations
#                            *exactly* invariant to different numbers of CPUs.
#
# - SOFTDOUBLEDOUBLE         when this is set, a software implementation of
#                            128bit double-double addition is used, implemented as a c++ class.
#                            Hence, this option requires compilation with a c++ compiler
#


#     - QUICK_LYALPHA:   This only works for cosmological simulations in periodic boxes
#                        with COOLING & SFR. (WINDS, METALS should be deselected).
#                        It will simply convert all gas particles above overdensity
#                        CritPhysOverdensity and with Temperature below 10^5 K to stars.
#                        This should still leave the Ly-Alpha forest largely unaffected,
#                        but should be faster. It is recommended to set GENERATIONS equal
#                        to 1 for maximum speed-up.
#
#     - TIDALTENSOR:     Calculates the tidal tensor for each particle. 

##################################################################################
#
#  Options for RT:
#
#  - RADTRANSFER        main switch; RT equation solved using the Conjugate Gradient iterative method
#
#  - RADTRANSFER_FLUXLIMITER	intorduces a flux limiter, so that the Ifront speed does not exceed c
#
#  - RADTRANSFER_MODIFY_EDDINGTON_TENSOR	modifies the Eddington tensor to the fully anisotropic version
#
#  - RT_COOLING_PHOTOHEATING	includes photoheating and cooling
#
#  - EDDINGTON_TENSOR_STARS	includes stars as sources of ionising photons
#
#  - EDDINGTON_TENSOR_SFR	includes sf partices
#
#  - EDDINGTON_TENSOR_GAS	includes gas, works only together with RT_MULTI_FREQUENCY
#
#  - EDDINGTON_TENSOR_BH        includes bh, not working atm
#
#  - RT_INCLUDE_HE              includes helium cooling and collisional ionisation ifndef(RT_MULTI_FREQUENCY)
#                               includes helium completely (ionisation, cooling, heating, etc.) ifdef(RT_MULTI_FREQUENCY)
#
#  - HYDROGEN_ONLY              sets hydrogen fraction to 1.0
#
#  - RT_PRINT_ET                outputs the eddington tensor (used for diagnostics)
#
#  - RT_MULTI_FREQUENCY         enables multi-frequency radiation transport. Here the integration
#                               variable is the ionising intensity J_nu
####################################################################################                 
#
#  Overview of options added for fine-grained phase-space analysis: (Mark)
#
#  * Fine-grained phase space structure analysis 
#   - DISTORTIONTENSORPS        
#   - OUTPUT_DISTORTIONTENSORPS 
#   - OUTPUT_TIDALTENSORPS     
#   - CAUSTIC_FINDER=2+4+8+16+32 
#   - ANNIHILATION_RADIATION  
#   - OUTPUT_LAST_CAUSTIC  
####################################################################################
