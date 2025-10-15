/* GENERAL VARIABLES */

#define CM_PER_KPC                ((MyAtLeastDouble)3.085678e21)
#define CM_PER_PC                 ((MyAtLeastDouble)3.085678e18)
#define SOLAR_MASS_IN_CODE_UNITS  (SOLAR_MASS / All.UnitMass_in_g)   
#define SOLAR_MASS_IN_CODE_UNITS_h100  (SOLAR_MASS / All.UnitMass_in_g * All.HubbleParam)   
#define YEAR_IN_CODE_UNITS        (SEC_PER_YEAR / All.UnitTime_in_s) 
/* Constants to transform from Specific Energy (code units) to Temperature (K) and back */
#define SE_TO_T                   (GAMMA_MINUS1*All.UnitEnergy_in_cgs*mu_h*PROTONMASS/BOLTZMANN/All.UnitMass_in_g)
#define T_TO_SE                   (1.0/SE_TO_T)

#define SMALLNEG ((MyAtLeastDouble)-1.e-20)
#define SMALL    ((MyAtLeastDouble)1e-20)


/* MUPPI options */

#define WARNINGS


/* MUPPI PARAMETERS AND VARIABLES */

/* Parameters for the initialization of Multi-Phase particles */
#define FracH ((MyAtLeastDouble)1.0)      /* initial fraction of hot phase */
#define FracC ((MyAtLeastDouble)0.0)      /* initial fraction of cold phase */
#define FracS ((MyAtLeastDouble)0.0)      /* initial fraction of stars */
                                 /* NB: the above three numbers MUST SUM UP AT 1 */
                                 /* PREFERRED CHOICE: FracH=1.0, FracC=0.0, FracS=0.0 */

#define T_c ((MyAtLeastDouble)300.)       /* temperature of the cold phase */

#define CLOCK ((MyAtLeastDouble)1.0)      /* MP particle exits MUPPI after this times tdyn (once it is fixed) */

#define PRESFMOL ((MyAtLeastDouble)20000.) /* pressure at which Fmol=0.5 */

#define f_star ((MyAtLeastDouble)0.02)    /* fraction of cloud that is consumed into stars; 0.01 -> 0.1 */
#define f_evap ((MyAtLeastDouble)0.1)     /* fraction of cloud that is evaporated in a star-formation event */
#define myf_rest ((MyAtLeastDouble)0.2)     /* restored fraction from a generation of young stars in the I.R.A. */

#define NC_COLDMASS_TH 8.5       /* the dynamical time is computed when the cold gas mass is at least equal to 
				    NC_COLDMASS_TH times the hot gas mass*/

#define mybeta_sf ((MyAtLeastDouble)120.0 * SOLAR_MASS_IN_CODE_UNITS)  /* stellar masses of stars formed for each SN */
#define myE_SN_51 ((MyAtLeastDouble)1.e51 / All.UnitEnergy_in_cgs )      /* energy of a single SN, in code units */   

#define TEMPSFTHRES ((MyAtLeastDouble)50000.)                         /* temperature threshold for multi-phase */

#define mu_h  ((MyAtLeastDouble)4./ (5. * HYDROGEN_MASSFRAC + 3.) )   /* molecular weight of the hot phase */
#define mu_c  ((MyAtLeastDouble)4./ (3. * HYDROGEN_MASSFRAC + 1.) )   /* molecular weight of the cold phase */
#define rap_mol_weight   ((MyAtLeastDouble)mu_h/mu_c )                /* ratio of the two */
#define TMAXCOOL ((MyAtLeastDouble)9.e8)                              /* maximal temperature for the cooling time */

#define SCALE_HH ((MyAtLeastDouble)30.)        /* This is used to regulate the first guess for the integration timestep */
                                      /* the higher, the slower and more accurate */


#if defined(MV_GM_STELLAR_KIN_FB2)
#define FOUNTAIN_PROBABILITY ((MyAtLeastDouble)0.03) 
#define FOUNTAIN_DENSITY_THR ((MyAtLeastDouble)0.3)
#define DELAY_IN_TDYN ((MyAtLeastDouble)-15.e6)   /* DELAY_IN_TDYN > 0.0  -->  fountain particle is kicked for N dynamical times of the SF particles 
					  DELAY_IN_TDYN < 0.0  -->  it is kicked for -DELAY_IN_TDYN Myr - the time spent in multi-phase (=CLOCK * tdyn) */
#endif

#if defined (MV_KRUMHOLZ_MOLECULAR_FRACTION) || defined (MV_EARLY_FB_HIGH_SN_ENERGY_FOR_LOW_Z) /* see Krumholz, McKee, and Tumlinson (2009, 699) */
#define SOLAR_METALLICITY ((MyAtLeastDouble)0.01524)  /* Caffau+ 2011 */
#define MIN_METALLICITY ((MyAtLeastDouble)0.05)       /* minimum metallicity for the Krumholz model to be valid, see end of Section 2.1 */
#endif



/* MUPPI_main variables, direct access to memory! */

#define nMuppiInts 6
#define nMuppiDbls 32
#define  CheckMuppiMemory(x1, x2, x3, x4, x5, x6) checkmuppimemory(x1, x2, x3, x4, x5, x6, __FUNCTION__, __FILE__, __LINE__)


#define mnstep            *(intMuppiVars)
#define mpart_ID          *(intMuppiVars + 1)
#define minternal         *(intMuppiVars + 2)
#define mncall            *(intMuppiVars + 3)
#define mmfin             *(intMuppiVars + 4)
#define mmfout            *(intMuppiVars + 5)

#define mMass_tot         *(dblMuppiVars)
#define mEgy_tot          *(dblMuppiVars + 1)
#define mM_h_0            *(dblMuppiVars + 2)
#define mM_c_0            *(dblMuppiVars + 3)
#define mx1               *(dblMuppiVars + 4)
#define mx2               *(dblMuppiVars + 5)
#define mhh               *(dblMuppiVars + 6)
#define mT_h              *(dblMuppiVars + 7)
#define mfill_h           *(dblMuppiVars + 8)
#define mRho_h            *(dblMuppiVars + 9) 
#define mn_h              *(dblMuppiVars + 10)
#define mRho_c            *(dblMuppiVars + 11)
#define mn_c              *(dblMuppiVars + 12)
#define mM_sf0            *(dblMuppiVars + 13)
#define mdE               *(dblMuppiVars + 14)
#define mdM               *(dblMuppiVars + 15)
#define mVol_tot          *(dblMuppiVars + 16)
#define mfill_c           *(dblMuppiVars + 17)
#define mTime             *(dblMuppiVars + 18)
#define mdtime            *(dblMuppiVars + 19)
#define mhh_old           *(dblMuppiVars + 20)
#define mTotMassWtRest    *(dblMuppiVars + 21)

// needed for MV_KRUMHOLZ_MOLECULAR_FRACTION
#define mRho              *(dblMuppiVars + 22)

// needed for  MV_AGNMUPPI_COOLING_OFF
#define mt_startCoolOff   *(dblMuppiVars + 23)

// needed for AGNMUPPI
#define mBHEcold           *(dblMuppiVars + 24)
#define mBHEcoldD          *(dblMuppiVars + 25)
#define mBHEhot            *(dblMuppiVars + 26) //nota: non sarebbe necessario nel buffer mi pare
#define mInitialBHEcold    *(dblMuppiVars + 27)
#define mBHEcoldOrig       *(dblMuppiVars + 28) /* As vars 24,26 but in physical units*/
#define mBHEhotOrig        *(dblMuppiVars + 29)


/* The following variables are used for OUTPUT */
#define mtcool             *(dblMuppiVars + 30) /* INITIAL cooling time */
#define mT_h_0             *(dblMuppiVars + 31)
