/*
* @file
* This file is part of the developer version of GADGET3 and contains
* the license conditions for its usage.
*
* @author GADGET-development team, led by Volker Springel and Klaus Dolag.
*
* @section LICENSE
* Copyright (c) 2016, Volker Springel, Klaus Dolag, and all contributing authors
* (see change logs). All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*
* 1. Received source code may be modified and used as convenient.
*
* 2. Redistributions of source code or in binary form is only possible with
*    explicit agreement of the copyright holders.
*
* 3. Redistributions of source code must retain the above copyright notice,
*    this list of conditions and the following disclaimer.
*
* 4. Redistributions in binary form must reproduce the above copyright notice,
*    this list of conditions and the following disclaimer in the documentation
*    and/or other materials provided with the distribution.
*
* 5. Neither the name of the copyright holder nor the names of its
*    contributors may be used to endorse or promote products derived from this
*    software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
* IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
* ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
* LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
* CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
* SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
* INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
* CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
* ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
* POSSIBILITY OF SUCH DAMAGE.
*
*/
#ifndef INLINE_FUNC
#ifdef INLINE
#define INLINE_FUNC inline
#else
#define INLINE_FUNC
#endif
#endif


                                                                         /* :: -------------------------------------------- ::  */
                                                                         /*     variables related to reading the directory      */


/* :: -------------------------------------------- ::  */
/*    SEGMENT .DATA                                    */

                                                                         /* :: -------------------------------------------- ::  */
                                                                         /*     variables related to reading the directory      */
struct my_direntry
{
  char   redshift_str[10];
  double redshift;
};


static struct my_direntry *WalCool_CoolTables_redshifts;
static int    WalCool_CoolTables_num, WalCool_redshift_index;

#ifndef ACC_IO_C_INCLUDES_LT
#ifdef ACC
#else
static
#endif
int    myHyd, myHel, HydPos, HelPos;
#endif


                                                                         /* :: -------------------------------------------- ::  */
                                                                         /*     define tables                                   */
static int HighZTable, LowZTable;

#ifndef ACC_IO_C_INCLUDES_LT
#ifdef ACC
#else 
static 
#endif
int WalCool_n[6];
#endif

#define WalCool_n_Rho WalCool_n[0]                                       /*!< number of density bins      */
#define WalCool_n_T   WalCool_n[1]                                       /*!< number of temperature bins  */
#define WalCool_n_Hef WalCool_n[2]                                       /*!< number of helium fractions  */
#define WalCool_n_El  WalCool_n[3]                                       /*!< number of elements          */
#define WalCool_n_Ab  WalCool_n[4]                                       /*!< number of solar abundances  */

static float *WalCool_indexes_and_arrays;
static float *WalCool_Rho;                                              /*!< stores the density bins        */
static float *WalCool_T;                                                /*!< stores the temperature bins    */
static float *WalCool_U;                                                /*!< stores the energy density bins */
static float *WalCool_Hef;                                              /*!< stores the Helium fractions    */
static float *WalCool_Ab;                                               /*!< stores the solar abundances    */
static float *WalCool_redshifts;                                        /*!< stores the redshifts           */

static double Trange, Tmin_wal;
static double Urange, Umin;
static double Rrange, Rmin;
static double Herange, Hemin;

static float *WalCool_enH;                                              /*!< stores the electron abundances over n_H     */
static float *WalCool_enHS;                                             /*!< stores the electron abundances over n_H for solar composition  */
static float *WalCool_mu;                                               /*   mean molecular w and temperature conversion */
static float *WalCool_UtoT;                                             /*   IN PHOTOIONIZATION EQUILIBRIUM              */
/* static float *WalCool_enH_collis;                                       /\*!< stores the electron abundances over n_H     *\/ */
/* static float *WalCool_mu_collis;                                        /\*   mean molecular w and temperature conversion *\/ */
/* static float *WalCool_UtoT_collis;                                      /\*   with NO IONIZING BACKGROUND                 *\/ */
                                                                        




static char **WalCool_names_tab;
static char *WalCool_names, **WalCool_El_names, **WalCool_Ab_names, **WalCool_Ab_symbols;


static float *WALCOOLTABLES, *WalMfreeCoolTables, *WalCoolTables;
static int  WalCoolMfreeSize, WalCoolSize, WalCool_arrays_size, WalCool_arraysS_size;
static int  UtoT_Size, enH_mu_Size;

static double max_cool_redshift;

//#define ENH_MU_IDX(T, nH, Hef) ((nH) * (WalCool_n_T * WalCool_n_Hef) + (Hef) * WalCool_n_T + (T))
//#define ENH_MU_COLLIS_IDX(T, Hef) ((Hef) * WalCool_n_T + (T))


                                                                         /*  --- ACCESSING TABLES ---  */

                                                                         /* [ COOLING RATES with UV BCKGRND ]  */

                                                                         /* ... starting point of a whole table at a given redshift */
#define MFREE_COOLRATE_z_IDX(Redsh_idx) ((Redsh_idx) * WalCoolMfreeSize)
#define COOLRATE_z_IDX(Redsh_idx) ((Redsh_idx) * WalCoolSize)

                                                                         /* ... starting point of an element's / Hef table at a given redshift */
#define COOLRATE_el_IDX(Redsh_idx, el_idx) ((Redsh_idx) * WalCoolSize + (el_idx) * WalCool_n_T * WalCool_n_Rho)
#define MFREE_COOLRATE_Hef_IDX(Redsh_idx, He_idx) ((Redsh_idx) * WalCoolMfreeSize + (He_idx) * WalCool_n_T * WalCool_n_Rho)

                                                                         /* ... access both element and Hef table at a given redshift */
                                                                         /* ... note: must use the previous index as 0 point          */
#define COOLRATE_IIDX(T_idx, R_idx) ((T_idx) * WalCool_n_Rho + (R_idx))

                                                                         /* ... complete index */
#define COOLRATE_IDX(Redsh_idx, el_idx, T_idx, R_idx) ((Redsh_idx) * WalCoolSize + (el_idx) * WalCool_n_T * WalCool_n_Rho + (T_idx) * WalCool_n_Rho + (R_idx))
#define MFREE_COOLRATE_IDX(Redsh_idx, He_idx, T_idx, R_idx) ((Redsh_idx) * WalCoolMfreeSize + (He_idx) * WalCool_n_T * WalCool_n_Rho + (T_idx) * WalCool_n_Rho + (R_idx))


                                                                         /* [ COOLING RATES with NO UV BCKGRND ]  */

#define COOLRATE_collis_IDX(T_idx, el_idx) ((T_idx) * WalCool_n_El + el_idx)
#define MFREE_COOLRATE_collis_IDX(He_idx, T_idx) ((He_idx) * WalCool_n_T + (T_idx))


                                                                         /* [ U to T ]  */

#define UtoT_z_IDX(Redsh_idx) ((Redsh_idx) * UtoT_Size)
#define UtoT_IDX(Redsh_idx, He_idx, U_idx, R_idx) ((Redsh_idx) * UtoT_Size + (He_idx) * WalCool_n_Rho * WalCool_n_T + (U_idx) * WalCool_n_Rho + (R_idx))

#define UtoT_collis_IDX(He_idx, U_idx) ((He_idx) * WalCool_n_T + (U_idx))


                                                                         /* [ ne over nH ]  */

#define ENH_z_IDX(Redsh_idx) ((Redsh_idx) * enH_mu_Size)
#define ENH_IDX(Redsh_idx, He_idx, T_idx, nH) ((Redsh_idx) * enH_mu_Size + (He_idx) * WalCool_n_Rho * WalCool_n_T + (T_idx) * WalCool_n_Rho + (nH))

#define ENH_collis_IDX(He_idx, T_idx) (((He_idx) * WalCool_n_T) + (T_idx))

#define ENHS_z_IDX(Redsh_idx) ((Redsh_idx) * (WalCool_n_Rho * WalCool_n_T))
#define ENHS_IDX(Redsh_idx, T_idx, nH) ((Redsh_idx) * (WalCool_n_Rho * WalCool_n_T) + (T_idx) * WalCool_n_Rho + (nH))

#define ENHS_collis_IDX(T_idx) ((T_idx))

                                                                         /* [ mu ]  */

#define MU_z_IDX(Redsh_idx) ((Redsh_idx) * enH_mu_Size)
#define MU_IDX(Redsh_idx, He_idx, T_idx, nH) ((Redsh_idx) * enH_mu_Size + (He_idx) * WalCool_n_Rho * WalCool_n_T + (T_idx) * WalCool_n_Rho + (nH))



                                                                         /*  --- ---------------- ---  */
#define LIN 0
#define LOG 1

#define z_IDX  0
#define He_IDX 1
#define U_IDX  2
#define H_IDX  3

#define iz  IDX[z_IDX]
#define iHe IDX[He_IDX]
#define iU  IDX[U_IDX]
#define inH IDX[H_IDX]

#define dz  dX[z_IDX]
#define dHe dX[He_IDX]
#define dU  dX[U_IDX]
#define dnH dX[H_IDX]

#define CHECK_ABND_and_IDX {if(!abundances_ad_indexes_are_set) endrun(LT_ERR_WALCOOL_AB_IDX_NOT_SET);}

#define T_CMB 2.75

#define TCMB(z) ( T_CMB * (1 + (z)) )

#ifndef ACC_IO_C_INCLUDES_LT
#ifdef ACC
int    *SpeciesIdx;
int    *SpeciesPos;
#else
static int    *SpeciesIdx;
static int    *SpeciesPos;
#endif
#endif


static int    *Specie_is_coolant;

static int    collisional_table_loaded = 0;

#ifndef ACC_IO_C_INCLUDES_LT
#ifdef ACC

#else
static int    UseHeNumberRatio = 0;
#endif
#endif




void WalCool_Initialize(void);

int WalCool_get_redshift_table(void);
int is_it_a_tablefile(char *, int);

int compare_dir_redshifts(const void *, const void *);
int compare_redshifts(const void *, const void *);

//int get_cool_redshift(double, double*);

void WalCool_get_table(int, int);

void WalCool_Initialize_get_header(void);

float get_Ndim_interp(int , double *, double *);

int set_indexes(double, double, double *, double *, int *);
int set_indexes_T(double, double, double *, double *, int *);
int INLINE_FUNC find_index(int, float*, double, double, int, float, double*);

double InnerInterpolation(float*, float*, double*, int*);

#define LT_ERR_WALCOOL_Z_FOR_NOTGAS 982276
#define LT_ERR_WALCOOL_IMPOSSIBLE_TO_READ_ABUNDANCE 982277
#define LT_ERR_WALCOOL_EL_MISSED 982278
