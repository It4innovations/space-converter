#ifdef GM_MUPPI	
struct myparams { /* used to pass parameters to derivatives */
  MyAtLeastDouble dE;
  MyAtLeastDouble dM;
  MyAtLeastDouble VolTot;
  int internal;
  int part_ID;
  int filtered;
  int ncall;
  MyAtLeastDouble fill_c;
  MyAtLeastDouble * dblsmuppi;
  int *intsmuppi;
};  



/* sfr_muppi.c */

static inline void starformation(int, int, int *, int *, MyAtLeastDouble, int), cooling(int, MyAtLeastDouble, int), final_statistics(MyAtLeastDouble, int, int);
static inline int check_sf_or_cooling(int);
static inline void make_a_star(int, MyAtLeastDouble, int);
static inline int produce_star_particle(int, MyAtLeastDouble, int);
static inline int convert_gas_into_star_particle(int);
void MUPPI_main(int, MyAtLeastDouble, int*, int*);
static inline int check_exit_muppi(int, MyAtLeastDouble, MyAtLeastDouble *, int *);
static inline int muppi_integration_check_exits(int i, int status, MyAtLeastDouble *dblMuppiVars, int *intMuppiVars, double *y);
static inline void muppi_final_step_operations(int i, MyAtLeastDouble Rho, MyAtLeastDouble *dblMuppiVars, int *intMuppiVars, double *y);
static inline void muppi_initialize_integration(int i, MyAtLeastDouble Rho, MyAtLeastDouble *dblMuppiVars, int *intMuppiVars, double *y, struct myparams *param);
static inline void initialize_muppi_vars(int i, MyAtLeastDouble *Rho, MyAtLeastDouble *dblMuppiVars, int *intMuppiVars);
static inline void start_MUPPI(int i, MyAtLeastDouble *dblMuppiVars, int *intMuppiVars);
static inline void muppi_exec_exits(int i, MyAtLeastDouble *dblMuppiVars, int *intMuppiVar, double *y);
#ifdef MV_GM_AGNMUPPI
static inline int derivatives_initial_check(double x, MyAtLeastDouble E_h, MyAtLeastDouble M_h, MyAtLeastDouble M_c, MyAtLeastDouble dEhydro, MyAtLeastDouble FBVol_tot, MyAtLeastDouble M_s, MyAtLeastDouble E_c_AGN_used, MyAtLeastDouble dBHEcold, int FBpart_ID, int FBinternal, int *FBfiltered, double *dydx, MyAtLeastDouble *dblMuppiVars);
#else
static inline int derivatives_initial_check(double x, MyAtLeastDouble E_h, MyAtLeastDouble M_h, MyAtLeastDouble M_c, MyAtLeastDouble dEhydro, MyAtLeastDouble FBVol_tot, MyAtLeastDouble M_s, int FBpart_ID, int FBinternal, int *FBfiltered, double *dydx, MyAtLeastDouble *dblMuppiVars);
#endif
  static inline void derivatives_maximal_starburst(int dM, MyAtLeastDouble n_h, MyAtLeastDouble T_h, MyAtLeastDouble E_h, MyAtLeastDouble M_h, MyAtLeastDouble Rho_h, MyAtLeastDouble M_c, int FBpart_ID, int FBinternal, const double *y, double *dydx, MyAtLeastDouble *dblMuppiVars);



void checkmuppimemory(int i, MyAtLeastDouble Rho, int *mfin, int *mfout, MyAtLeastDouble *dblMuppiVars, int *intMuppiVars, const char *func, const char *file, int line);

#ifdef LT_STELLAREVOLUTION
/* XXXX warning, these are defined in double (both formal  parameter and output) */

#ifdef LT_METAL_COOLING
double get_cooling_time_fromT(double, double, double, double *, double);
#endif
#ifdef LT_METAL_COOLING_WAL
#ifdef GL_DUST_COOLING
double get_cooling_time_fromT(double, double, double, double *, double, double, double, double, double *);
#else
double get_cooling_time_fromT(double, double, double, double *, double, double, double *);
#endif
void set_cosmo_factors_for_current_time(void);
#endif

#else

double get_cooling_time_fromT(double, double, double, MyAtLeastDouble *);

#endif


//void FB_update_part(int, char);
void check_muppi_nan(int);

//static inline void compute_ism_properties(MyAtLeastDouble M_h, MyAtLeastDouble E_h, MyAtLeastDouble M_c, MyAtLeastDouble Vol_tot,

void compute_ism_properties(MyAtLeastDouble M_h, MyAtLeastDouble E_h, MyAtLeastDouble M_c, MyAtLeastDouble Vol_tot,
					  MyAtLeastDouble *T_h, MyAtLeastDouble *fill_h, MyAtLeastDouble *Rho_h, MyAtLeastDouble *Rho_c, MyAtLeastDouble *n_h, MyAtLeastDouble *n_c);
static inline MyAtLeastDouble cooling_time(MyAtLeastDouble, MyAtLeastDouble, MyAtLeastDouble, MyAtLeastDouble);
static inline MyAtLeastDouble molecular_fraction(MyAtLeastDouble);


#ifdef MV_KRUMHOLZ_MOLECULAR_FRACTION 
static inline MyAtLeastDouble molecular_fraction_Krum(MyAtLeastDouble, MyAtLeastDouble);
#endif


/* muppi_kicks.h */
static inline void muppi_kinetic(int);


#endif

/* warning, some functions used in MUPPI are defined in other 
   proto files, e.g. cooling.h, proto.h, forcetree.h */

