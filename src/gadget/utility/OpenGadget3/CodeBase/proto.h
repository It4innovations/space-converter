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
#ifndef PROTO_H
#define PROTO_H

#ifndef ALLVARS_H
#include "../CodeBase/allvars.h"
#endif
#include "../Gravity/forcetree.h"
#include "../CodeBase/domain.h"
#ifdef COOLING
#include "../CoolingSfr/cooling.h"
#endif
#ifdef LT_STELLAREVOLUTION
#include "../CoolingSfr/Sfr_LT/lt.h"
#endif

#define MAXTAGS 300
#define MAXTAGLEN 50

#ifdef HAVE_HDF5
#include <hdf5.h>
void write_header_attributes_in_hdf5(hid_t handle);
void read_header_attributes_in_hdf5(char *fname);
void write_parameters_attributes_in_hdf5(hid_t handle);
void write_units_attributes_in_hdf5(hid_t handle);
void write_constants_attributes_in_hdf5(hid_t handle);
#endif
void report_VmRSS(void);
void write_ps_files(void);
void output_compile_time_options(void);

void report_pinning(void);
void pin_to_core_set(void);
void get_core_set(void);

void init_turb(void);
void set_turb_ampl(void);
void add_turb_accel(void);
void reset_turb_temp(void);
void log_turb_temp(void);
void init_static_nfw();

void create_snapshot_if_desired(void);

void mpi_report_comittable_memory(void);
void report_task_list(void);
long long report_comittable_memory(long long *MemTotal,
				   long long *Committed_AS, long long *SwapTotal, long long *SwapFree);

void do_first_halfstep_kick(void);
void do_second_halfstep_kick(void);
void find_timesteps(void);
void compute_densities(void);
void compute_grav_accelerations(void);
void calc_memory_checksum(void *base, size_t bytes);

void get_disk_forces(double RR, double zz, double *f_R, double *f_z);
double get_disk_mass(double time);
void growing_disk_init(void);

double get_turb_pot(double x, double y, double z);

void sub_turb_move_perturbers(double t0, double t1);
void sub_turb_add_forces(void);
void sub_turb_read_table(void);
void sub_turb_parent_halo_accel(double dx, double dy, double dz, double *acc);
double sub_turb_enclosed_mass(double r, double msub, double vmax, double radvmax, double c);


int powerspec_turb_treefind(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			    int *nexport, int *nsend_local);
#ifdef AR_GREEN_TREE_FOF
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *ngblist, int extended_numngb_inbox,
					 int *i_startnode, int node, int i_node);
#else
int powerspec_turb_find_nearest_evaluate(int target, int mode, int *nexport, int *nsend_local);
#endif
void powerspec_turb_calc_dispersion(void);
double powerspec_turb_obtain_fields(void);
void powerspec_turb_save(char *fname, double *disp);
void powerspec_turb_collect(void);
void powerspec_turb(int filenr);


void set_cosmo_factors_for_current_time(void);
void drift_sph_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void set_non_standard_physics_for_current_time(void);
void calculate_non_standard_physics(void);
void compute_statistics(void);
void execute_resubmit_command(void);
void make_list_of_active_particles(void);
void output_extra_log_messages(void);

int check_stop_condition(void);


#ifdef ACC_DMAX			//_PROTO
#define DMAX(a,b) (a > b) ? a : b
#define DMIN(a,b) -DMAX(-a,-b)
#define IMAX(a,b) DMAX(a,b)
#define IMIN(a,b) DMIN(a,b)
#else

static inline double DMAX(double a, double b)
{
  return (a > b) ? a : b;
}

#ifndef GADGET3_IO_LIB
static inline long double DMAX(long double a, long double b)
{
  return (a > b) ? a : b;
}

static inline long double DMAX(double a, long double b)
{
  return (a > b) ? a : b;
}

static inline long double DMAX(long double a, double b)
{
  return (a > b) ? a : b;
}

#endif

static inline double DMIN(double a, double b)
{
  return (a < b) ? a : b;
}

static inline int IMAX(int a, int b)
{
  return (a > b) ? a : b;
}

static inline int IMIN(int a, int b)
{
  return (a < b) ? a : b;
}

#endif //else of ifdef ACC
void do_distortion_tensor_kick(int i, double dt_gravkick);
void set_predicted_sph_quantities_for_extra_physics(int i);
void do_sph_kick_for_extra_physics(int i, integertime tstart, integertime tend, double dt_entr);

void check_particle_for_temperature_minimum(int i);

double get_pressure(int i);
double get_gamma(int i);
double get_gamma_minus1(int i);

void read_fof(int num);
int fof_compare_ID_list_ID(const void *a, const void *b);

void myfree_msg(void *p, char *msg);
void kspace_neutrinos_init(void);

#ifndef TWOPOINT_ON_OUTPUT
void twopoint(void);
void twopoint_save(void);
#else
void twopoint(int);
void twopoint_save(int);
#endif
int twopoint_ngb_treefind_variable(MyLongDouble searchcenter[3], MyFloat rsearch, int target, int *startnode,
				   int mode, int *nexport, int *nsend_local);
int twopoint_count_local(int target, int mode, int *nexport, int *nsend_local);

void powerspec(int flag, int *typeflag);
double PowerSpec_Efstathiou(double k);
void powerspec_save(void);
void foldonitself(int *typelist);
void dump_potential(void);
void read_and_dump_potential(void);

int snIaheating_evaluate(int target, int mode, int *nexport, int *nSend_local);
void snIa_heating(void);
void voronoi_setup_exchange(void);

#ifdef KSPACE_NEUTRINOS
double get_neutrino_powerspec(double k, double ascale);
double get_powerspec(double k, double ascale);
void init_transfer_functions(void);
#endif

#ifdef KSPACE_NEUTRINOS_2
void get_delta_nu_update(double a, int nk_in, double wavenum[], double P_cdm_curr[], double delta_nu_curr[]);
double get_neutrino_powerspec(double kk_in, double SmoothK[], double SmoothPowerNu[],
			      gsl_interp * SplinePowNu, gsl_interp_accel * acc_nu, double SmoothPower[],
			      gsl_interp * SplinePow, gsl_interp_accel * acc, int nbins);
void transfer_init_tabulate(int nk_in);
void rebin_power(double SmoothK[], double SmoothPower[], int nbins, double Kbin[], double SumPower[],
		 long long CountModes[], int bins_ps);
/* Return the total matter density in all neutrino species.*/
double OmegaNu(double a);
#endif

#if !defined(GADGET3_IO_LIB) || defined(GADGET3_IO_LIB_MPI)
int MPI_Check_Sendrecv(void *sendbuf, int sendcount, MPI_Datatype sendtype,
		       int dest, int sendtag, void *recvbufreal, int recvcount,
		       MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm, MPI_Status * status);

int MPI_Sizelimited_Sendrecv(void *sendbuf, size_t sendcount, MPI_Datatype sendtype,
			     int dest, int sendtag, void *recvbuf, size_t recvcount,
			     MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
			     MPI_Status * status);
#endif

int mpi_calculate_offsets(int *send_count, int *send_offset, int *recv_count, int *recv_offset,
			  int send_identical);
void sort_based_on_field(void *data, int field_offset, int n_items, int item_size, void **data2ptr);
void mpi_distribute_items_to_tasks(void *data, int task_offset, int *n_items, int *max_n, int item_size);

void parallel_sort_special_P_GrNr_ID(void);
void calculate_power_spectra(int num, long long *ntot_type_all);
#ifdef ACC_PROTO
#pragma acc routine
#endif
int pmforce_is_particle_high_res(int type, MyLongDouble * pos);

void compare_partitions(void);
void assign_unique_ids(void);
int permut_data_compare(const void *a, const void *b);
void generate_permutation_in_active_list(void);
void get_particle_numbers(char *fname, int num_files);

void conduction(void);
void conduction_matrix_multiply(double *in, double *out);
double conduction_vector_multiply(double *a, double *b);
int conduction_evaluate(int target, int mode, int *ngblist, int extended_numngb_inbox, int *i_startnode,
			int node, int i_node, void *data);



void *conduction_evaluate_primary(void *p);
void *conduction_evaluate_secondary(void *p);

void fof_get_group_center(double *cm, int gr);
void fof_get_group_velocity(double *cmvel, int gr);
int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);
void fof_compute_group_properties(int gr, int start, int len);

void parallel_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
#if !defined(GADGET3_IO_LIB) || defined(GADGET3_IO_LIB_MPI)
void parallel_sort_comm(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *),
			MPI_Comm comm);
#endif
int compare_IDs(const void *a, const void *b);
void test_id_uniqueness(void);

void check_particles_info(const char *func, const char *file, int linenr);
void check_node_info(const char *func, const char *file, int linenr);

int io_compare_P_ID(const void *a, const void *b);
int io_compare_P_GrNr_SubNr(const void *a, const void *b);


void drift_particle(int i, integertime time1);
int ShouldWeDoDynamicUpdate(void);

void put_symbol(double t0, double t1, char c);
void write_cpu_log(void);

int get_timestep_bin(integertime ti_step);

const char *svn_version(void);

void find_particles_and_save_them(int num);
void lineofsight_output(void);
void sum_over_processors_and_normalize(void);
void absorb_along_lines_of_sight(void);
void output_lines_of_sight(int num);
integertime find_next_lineofsighttime(integertime time0);
integertime find_next_gridoutputtime(integertime ti_curr);
void add_along_lines_of_sight(void);
double los_periodic(double x);
void do_the_kick(int i, integertime tstart, integertime tend, integertime tcurrent);

void x86_fix(void);

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file,
				int line);

void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);

void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);

void mymalloc_init(void);
void dump_memory_table(void);
void report_detailed_memory_usage_of_largest_task(size_t *OldHighMarkBytes, const char *label,
						  const char *func, const char *file, int line);

double get_shear_viscosity(int i);

void kinetic_feedback_mhm(void);
int kin_compare_key(const void *a, const void *b);
void kinetic_evaluate(int target, int mode);

void bubble(void);
void multi_bubbles(void);
void bh_bubble(double bh_dmass, MyFloat center[3], MyIDType BH_id);
void find_CM_of_biggest_group(void);
int compare_length_values(const void *a, const void *b);
double rho_dot(double z, void *params);
double bhgrowth(double z1, double z2);

void second_order_ics(void);
double F1_Omega(double);
double F2_Omega(double);

//int smoothed_evaluate(int target, int mode, int *nexport, int *nsend_local);
int smoothed_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		      int *ngblist);
void *smoothed_evaluate_primary(void *p);
void *smoothed_evaluate_secondary(void *p);
void smoothed_values(void);

int fof_find_dmparticles_evaluate(int target, int mode, int *nexport, int *nsend_local);

double INLINE_FUNC hubble_function(double a);
#ifdef DARKENERGY
double DarkEnergy_a(double);
double DarkEnergy_t(double);
#ifdef TIMEDEPDE
void fwa_init(void);
double INLINE_FUNC fwa(double);
double INLINE_FUNC get_wa(double);
#ifdef TIMEDEPGRAV
double INLINE_FUNC dHfak(double a);
double INLINE_FUNC dGfak(double a);
#endif
#endif
#endif

#ifdef EXTERNALHUBBLE
double hubble_function_external(double a);
#endif

void blackhole(void);
void report_internal_model_for_blackhole();

int ngb_treefind_fof_nearest(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode, int mode,
			     int *nexport, int *nsend_local, int MyFOF_PRIMARY_LINK_TYPES);

void ngb_init(void);

void fof_fof(int num);
void fof_import_ghosts(void);
void fof_course_binning(void);
void fof_find_groups(void);
void fof_check_cell(int p, int i, int j, int k);
void fof_find_minids(void);
int fof_link_accross(void);
void fof_exchange_id_lists(void);
int fof_grid_compare(const void *a, const void *b);
void fof_compile_catalogue(void);
void fof_save_groups(int num);
void fof_save_local_catalogue(int num);
double fof_periodic(double x);
double fof_periodic_wrap(double x);
void fof_find_nearest_dmparticle(void);
int fof_find_nearest_dmparticle_evaluate(int target, int mode, int *nexport, int *nsend_local);

int fof_compare_key(const void *a, const void *b);
void fof_link_special(void);
void fof_link_specialpair(int p, int s);
void fof_make_black_holes(void);

int io_compare_P_GrNr_ID(const void *a, const void *b);

void write_file(char *fname, int readTask, int lastTask);

void distribute_file(int nfiles, int firstfile, int firsttask, int lasttask, int *filenr, int *master,
		     int *last);

int get_values_per_blockelement(enum iofields blocknr);

int get_datatype_in_block(enum iofields blocknr);
void get_dataset_name(enum iofields blocknr, char *buf);


int blockpresent(enum iofields blocknr);
void fill_write_buffer(enum iofields blocknr, int *pindex, int pc, int type);
#ifdef TEST_COPY_BOX_PERIODICALLY
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type, int ncopy);
#else
//void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type);
void empty_read_buffer(enum iofields blocknr, int offset, int pc, int type, int offset_gas, int offset_stars,
		       int offset_bhs);
#endif
int get_particles_in_block(enum iofields blocknr, int *typelist);

int get_bytes_per_blockelement(enum iofields blocknr, int mode);

#if defined(TEST_COPY_BOX_PERIODICALLY) || defined(GADGET3_IO_LIB)
void read_file(char *fname, int readTask, int lastTask, int ncopy);
#else
void read_file(char *fname, int readTask, int lastTask);
#endif

void get_Tab_IO_Label(enum iofields blocknr, char *label);


void long_range_init_regionsize(void);

int find_files(char *fname);

int metals_compare_key(const void *a, const void *b);
void enrichment_evaluate(int target, int mode);



#if GADGET_HYDRO == HYDRO_SPH || GADGET_HYDRO == HYDRO_PESPH
void compute_hydro_accelerations(void);
int hydro_evaluate(int target, int mode, int *ngblist, int extended_numngb_inbox, int *i_startnode, int node,
		   int i_node, void *data);

#ifdef ACC_HYDRA
void *hydro_evaluate_primary(void *p, int OnGPU);
#else
void *hydro_evaluate_primary(void *p);
#endif
void *hydro_evaluate_secondary(void *p);

#ifdef AXION_DM
void ax_hydro_force(void);
int ax_hydro_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		      int *ngblist);
void *ax_hydro_evaluate_primary(void *p);
void *ax_hydro_evaluate_secondary(void *p);
#endif

#elif GADGET_HYDRO == HYDRO_MFM
void compute_fluxes(void);
void compute_gradients(void);
void compute_limiters(void);

void gradients(void);
void mfm_fluxes(void);
void mfm_limiters(void);
void mfm_fluxes_to_code(int i, integertime tstart_step, integertime tend_step);

int flux_evaluate(int, int, int *, int *, int *, int *);
void *flux_evaluate_primary(void *p);
void *flux_evaluate_secondary(void *p);

int gradient_evaluate(int, int, int *, int *, int *, int *);
void *gradient_evaluate_primary(void *p);
void *gradient_evaluate_secondary(void *p);

int slope_evaluate(int, int, int *, int *, int *, int *);
void *slope_evaluate_primary(void *p);
void *slope_evaluate_secondary(void *p);

void mfm_synchronize_entropy(int i);
void mfm_update_mass(int i, MyFloat old_mass, MyFloat new_mass);
void mfm_add_mass(int i, MyFloat mass_change);
#endif



void pm_init_nonperiodic_allocate(void);

void pm_init_nonperiodic_free(void);

double get_random_number(MyIDType id);
void set_random_numbers(void);

int grav_tree_compare_key(const void *a, const void *b);
int dens_compare_key(const void *a, const void *b);
int hydro_compare_key(const void *a, const void *b);

int data_index_compare(const void *a, const void *b);
int peano_compare_key(const void *a, const void *b);

#ifdef PMGRID
void create_plans(void);
#endif

void quicksort(void *a, size_t n, size_t es, int (*cmp)(const void *, const void *));
void serial_sort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
#ifdef _OPENMP
void omp_qsort(void *a, size_t n, size_t es, int (*cmp)(const void *, const void *));
void omp_mysort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
void serial_sort_omp(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
int floorLog2(unsigned int n);
#endif


void mysort_dataindex(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_domain(void *b, size_t n, size_t s);
void mysort_idlist(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_pmperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_pmnonperiodic(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));
void mysort_peano(void *b, size_t n, size_t s, int (*cmp)(const void *, const void *));

void check_wind_creation(void);
void treat_outflowing_particles(void);
void set_injection_accel(void);

#ifdef AR_GREEN_TREE_DENSITY
int density_evaluate(int target, int mode, int *ngblist, int extended_numngb_inbox, int *i_startnode,
		     int node, int i_node, void *additional_data);
#else
int density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		     int *ngblist);
#endif

#if (defined(fSIDM) || defined(rSIDM)) && !defined(ADAPTGRAVSOFT)
int pt1_density_evaluate(int target, int mode, int *ngblist, int extended_numngb_inbox, int *i_startnode,
			 int node, int i_node, void *additional_data);
void pt1_force_update_hmax(void);
#endif

void on_beginning_of_timestep();
void on_end_of_timestep();




#ifdef ACC_DENSITY
void *density_evaluate_primary(void *p, int OnGPU);
#else
void *density_evaluate_primary(void *p);
#endif
void *density_evaluate_secondary(void *p);
int density_isactive(int n);

#if (defined(fSIDM) || defined(rSIDM)) && !defined(ADAPTGRAVSOFT)
int pt1_density_isactive(int n);
#endif

#ifdef AXION_DM
int ax_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			int *ngblist, _Bool cycle);
void *ax_density_evaluate_primary(void *p, _Bool cycle);
void *ax_density_evaluate_secondary(void *p, _Bool cycle);
int axion_isactive(int n);
#endif

#ifdef SUBFIND
int ngb_treefind_fof_nearest_openmp(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				    int mode, int *exportflag, int *exportnodecount, int *exportindex,
				    int *ngblist, int MyFOF_PRIMARY_LINK_TYPES);
#ifdef AR_GREEN_TREE_FOF
int ar_green_ngb_treefind_fof_variable_threads(MyLongDouble searchcenter[3], MyFloat * phsml, int *startnode,
					       int mode, int *exportflag, int *exportnodecount,
					       int *exportindex, int *ngblist, int *targets, int n_targets,
					       int should_drift, int);
#endif
#ifdef DENSITY_SPLIT_BY_TYPE
void *subfind_density_evaluate_primary(int thread_id, char *todo, int j_in, int j_target);
#else
void *subfind_density_evaluate_primary(int thread_id, char *todo);
#endif
#ifdef DENSITY_SPLIT_BY_TYPE
void *subfind_density_evaluate_secondary(int thread_id, int j_in);
#else
void *subfind_density_evaluate_secondary(int thread_id);
#endif
#endif

void GetMachNumberCR(struct sph_particle_data *Particle);
void GetMachNumber(struct sph_particle_data *Particle);
void GetShock_DtEnergy(struct sph_particle_data *Particle);

size_t sizemax(size_t a, size_t b);


void reconstruct_timebins(void);

void init_peano_map(void);
peanokey peano_hilbert_key(int x, int y, int z, int bits);
peanokey peano_and_morton_key(int x, int y, int z, int bits, peanokey * morton);
peanokey morton_key(int x, int y, int z, int bits);

void catch_abort(int sig);
void catch_fatal(int sig);
void terminate_processes(void);
void enable_core_dumps_and_fpu_exceptions(void);
void write_pid_file(void);

#ifdef PAUSE_RUN_TO_ATTACH_DEBUGGER
void pause_run_to_attach_debugger();
#endif

void pm_init_periodic_allocate(void);

void pm_init_periodic_free(void);

void move_particles(integertime time1);


void find_next_sync_point_and_drift(void);
void find_dt_displacement_constraint(double hfac);
#ifdef MAKEGLASS
void do_glass_making_step(void);
#endif
#ifdef RELAXOBJECT
void determine_relaxfac(void);
#endif
#ifdef WAKEUP
void process_wake_ups(void);
#endif

void set_units_sfr(void);

void gravity_forcetest(void);

void allocate_commbuffers(void);
void allocate_memory(void);
void begrun(void);
void check_omega(void);
void close_outputfiles(void);
void compute_accelerations(void);
void compute_global_quantities_of_system(void);
void compute_potential(void);
void construct_timetree(void);
void cooling_and_starformation(void);
void cl_CoolingOnly(void);
void count_hot_phase(void);
void delete_node(int i);
void density(void);
#if (defined(fSIDM) || defined(rSIDM)) && !defined(ADAPTGRAVSOFT)
void pt1_density();
#endif
#ifdef AXION_DM
void ax_density(_Bool cycle);
#endif
void density_decouple(void);
void determine_interior(void);
int dissolvegas(void);
void do_box_wrapping(void);
double enclosed_mass(double R);
#ifndef LT_STELLAREVOLUTION
void endrun(int);
#else
void EndRun(int, const char *, const char *, const int);
#endif


void EndRunString(int do_output, const char *, const char *, const int, const char *, ...);

void energy_statistics(void);
void ensure_neighbours(void);

void output_log_messages(void);
void ewald_corr(double dx, double dy, double dz, double *fper);

void ewald_force(int ii, int jj, int kk, double x[3], double force[3]);
void ewald_force_ni(int iii, int jjj, int kkk, double x[3], double force[3]);

void ewald_init(void);
double ewald_psi(double x[3]);
#ifdef ACC_PROTO
#pragma acc routine
#endif
double ewald_pot_corr(double dx, double dy, double dz);
int find_ancestor(int i);
integertime find_next_outputtime(integertime time);
void find_next_time(void);
integertime find_next_time_walk(int node);
void free_memory(void);
void advance_and_find_timesteps(void);
integertime get_timestep(int p, double *a, int flag);

void determine_PMinterior(void);
void gravity_tree(void);
void hydro_force(void);
void init(void);

#if !defined(LT_STELLAREVOLUTION) || defined(GM_MUPPI)
double get_starformation_rate(int i);
#else
double get_starformation_rate(int i, float *Temperature, float *xclouds);
#endif


#ifdef GM_MUPPI
/* muppi_communications.c */
void thermalenergy(void);
int thermalenergy_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			   int *ngblist);
void *thermalenergy_evaluate_primary(void *p);
void *thermalenergy_evaluate_secondary(void *p);
/* muppi_normalize.c */
void muppinorm(void);
int muppinorm_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		       int *ngblist);
void *muppinorm_evaluate_primary(void *p);
void *muppinorm_evaluate_secondary(void *p);
#endif


#ifndef LT_STELLAREVOLUTION
void init_clouds(void);
void integrate_sfr(void);
#else
void init_clouds(int, double, double, double, double, double *, double *);
void integrate_sfr(double, double, double, double, double, double);
#endif
void insert_node(int i);
int mark_targets(void);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE * stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE * stream);
void mpi_printf(const char *fmt, ...);
void open_outputfiles(void);
void write_outputfiles_header(void);
void peano_hilbert_order(void);
double pot_integrand(double xx);
void predict(double time);
void predict_collisionless_only(double time);
void predict_sph_particles(double time);
void prepare_decouple(void);
void read_ic(char *fname);
void read_ic_cluster(char *fname);
void read_ic_cluster_gas(char *fname);
void read_ic_cluster_wimp(char *fname);
int read_outputlist(char *fname);
void read_parameter_file(char *fname, char *tags[], void **addr, int *id, int nt);
void rearrange_particle_sequence(void);
void reorder_gas(void);
void reorder_particles(void);
#ifdef _OPENMP
void kd_reorder_particles2(void);
void kd_reorder_gas2(void);
#endif
void restart(int mod);
void run(void);
void savepositions(int num);
void savepositions_ioformat1(int num);
#ifdef OUTPUT_LIGHTCONES
void lightcone_output(int num);
void read_step_file(void);
void free_step_data(void);
#endif
double second(void);
void set_softenings(void);
void set_sph_kernel(void);
void set_units(void);
void setup_smoothinglengths(void);

void minimum_large_ints(int n, long long *src, long long *res);
void sumup_large_ints(int n, int *src, long long *res);
void sumup_longs(int n, long long *src, long long *res);

void statistics(void);
double timediff(double t0, double t1);
void veldisp(void);
void veldisp_ensure_neighbours(int mode);

void gravity_tree_shortrange(void);

#ifdef AXION_DM
double get_ax_hydrokick_factor(integertime time0, integertime time1);
static double ax_hydrokick_integ(double a, void *param);
#endif

double get_hydrokick_factor(integertime time0, integertime time1);
double get_gravkick_factor(integertime time0, integertime time1);
static double drift_integ(double a, void *param);
static double gravkick_integ(double a, void *param);
// double growthfactor_integ(double a, void *param);
static double hydrokick_integ(double a, void *param);
void init_drift_table(void);
double get_drift_factor(integertime time0, integertime time1);
#ifdef MAGNETIC
double get_magkick_factor(integertime time0, integertime time1);
#endif
#if defined(fSIDM) || defined(rSIDM)
static double mSIDM_integ(double a, void *param);
double get_mSIDMkick_factor(integertime time0, integertime time1);
#endif
double measure_time(void);
double report_time(void);

/* on some DEC Alphas, the correct prototype for pow() is missing,
   even when math.h is included ! */
#ifndef ACC_PROTO
double pow(double, double);
#endif

void long_range_init(void);
void long_range_force(void);

void pm_init(void);
void pm_init_regionsize(void);
int pmforce(int mode_periodic);	//, int mode=0, int *typelist=NULL);

int pmpotential_nonperiodic(int grnr);
void pmpotential_periodic(void);

void readjust_timebase(double TimeMax_old, double TimeMax_new);

double enclosed_mass(double R);
void pm_setup_nonperiodic_kernel(void);

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
double dmax(double, double);
double dmin(double, double);
#endif

#ifdef LT_STELLAREVOLUTION
void fsolver_error_handler(const char *, const char *, int, int);
int get_Yset(int);

void ReadYields(int, int, int);
void reading_thresholds_for_thermal_instability(void);
double GetMetalLambda(double, double);

void init_SN(void);
void TestStellarEvolution();
void TestStarsEvolution(int, int);
void get_metals_checksum(int, long double *);
unsigned int evolve_SN(void);
double INLINE_FUNC get_NextChemTime(double, int, int *);
double INLINE_FUNC get_cost_SE(int);
int INLINE_FUNC get_chemstep(int, int, double *, double);
int INLINE_FUNC get_current_chem_bin(double, int);
int get_chemstep_bin(double, double, int *, int);
int compare_steps(const void *, const void *);

int append_chemicallyactive_particles(unsigned int *);
void drop_chemicallyactive_particles(void);


void read_metals(void);

void init_clouds_cm(int, double *, double *, double, double, double, int, double *);

int load_SFs_IMFs(void);
#ifdef ACC_PROTO
#pragma acc routine
#endif
int get_SF_index(int, int *, int *);

int write_eff_model(int, int);
int read_eff_model(int, int);

void read_SnIa_yields(void);
void read_SnII_yields(void);
void read_AGB_yields(void);

void read_metalcool_table(void);

void initialize_star_lifetimes(void);
double INLINE_FUNC get_age(double);

void recalculate_stellarevolution_vars(void *, int);
void recalc_eff_model(void);

int calculate_effective_yields(double, double, int, int);
double calculate_FactorSN(double, double, void *);

double get_imf_params(double);
double INLINE_FUNC normalize_imf(double, double, void *);

double *get_meanZ(void);
void write_metallicity_stat(void);
double INLINE_FUNC get_metalmass(float *);
double INLINE_FUNC get_metallicity(int, int);
#ifdef ACC_PROTO
#pragma acc routine
#endif
double INLINE_FUNC get_PhysDensThresh(int);
#ifdef SUBFIND
double INLINE_FUNC get_metallicity_subfind(int);
double get_starformation_rate_subfind(int i, float *Temperature, float *xclouds);
int get_SF_index_subfind(int, int *, int *);
#endif
double INLINE_FUNC get_metallicity_solarunits(MyFloat);

double INLINE_FUNC get_entropy(int);
#ifdef ACC_PROTO
#pragma acc routine
#endif

int INLINE_FUNC getindex(double *, int, int, double *, int *);

int INLINE_FUNC perturb(double *, double *);
int is_chemically_active(int);
#endif

#if defined (UM_METAL_COOLING)	/* not necessarily UM_CHEMISTRY */

#ifdef UM_e_MET_IMPACTS
double um_metal_line_cooling(double, double, int, double);
#else
double um_metal_line_cooling(double, double, int);
#endif

#endif

#if defined(CHEMISTRY)
int compute_abundances(int mode, int ithis, double a_start, double a_end);
#endif

#if defined(CHEMISTRY) || defined(UM_CHEMISTRY)
int InitChem(void);
int init_rad(double);
#endif


#ifdef UM_CHEMISTRY

/* main call of chemistry treatment */
int compute_abundances(int mode, int ithis, double a_start, double a_end, double *um_energy);

/* main prototypes of chemistry treatment */
double INLINE_FUNC Um_Compute_MeanMolecularWeight(int);
double Um_AbundanceRatios(double, double, double *, double *, double *, double, double, int);


#if (defined(SFR) || defined(COOLING)) && !(defined(LT_METAL_COOLING) || defined(LT_METAL_COOLING_WAL))
double Um_DoCooling(double, double, double, double *, int, int);	// std SFR with um_chemistry
double Um_GetCoolingTime(double, double, double *, int);	// std SFR + um_chem
#else
#ifndef LT_METAL_COOLING_WAL
double Um_DoCooling(double, double, double, double *, double, int, int);
double Um_GetCoolingTime(double, double, double *, double, int);
#else
double Um_DoCooling(double, double, double, double *, double *, double, double, int, int);
double Um_GetCoolingTime(double, double, double *, double *, double, double, int);
#endif
double Um_AbundanceRatios(double, double, double *, double *, double *, double, double, int);

#endif // close std SFR with um_chemistry


/* additional prototypes of chemistry treatment */
double INLINE_FUNC Um_Compute_basenumden(double, int);
double INLINE_FUNC Um_Compute_n_den(double, int);
double INLINE_FUNC Um_Compute_sum_nuclei(int);
double INLINE_FUNC Um_Compute_mean_nden(int);

#if defined (UM_METAL_COOLING) && (defined (LT_METAL_COOLING) || defined(LT_METAL_COOLING_WAL)) && defined (UM_MET_IN_NONEQ_COOLING) && defined(UM_MET_NONEQ_CORRECTIONS)
void Um_met_corrections(int);
#else
void Std_corrections(int);
#endif

#if defined (UM_METAL_COOLING) && defined (LT_METAL_COOLING) && defined (UM_MET_IN_NONEQ_COOLING)
double INLINE_FUNC Um_Compute_invden(int i);
#endif

/* these fcts are not passed anywhere outside chemistry (see um_chemistry_noneq.c) */
//int energy_solver(int, int, double, double double, double, double *, double *, double *);
//int energy_solver(int, int, double, double, double, double, double *, double *)
//int species_solver(int, int, double, double, double, double, double *);
//int radiative_rates(int, double, double);

#endif

#ifdef UM_CHECK
void Um_cooling_check(void);
#endif





#ifdef AUTO_SWAP_ENDIAN_READIC
void swap_Nbyte(char *data, int n, int m);
void swap_header(void);
#endif
void find_block(char *label, FILE * fd);

#ifdef DISTORTIONTENSORPS
void get_half_kick_distortion(int pindex, MyBigFloat half_kick_add[6][6]);
void analyse_phase_space(int pindex,
			 MyBigFloat * s_1, MyBigFloat * s_2, MyBigFloat * s_3,
			 MyBigFloat * smear_x, MyBigFloat * smear_y, MyBigFloat * smear_z,
			 MyBigFloat * second_deriv, MyBigFloat * sigma);
MyBigFloat get_analytic_annihilation(MyBigFloat s_1, MyBigFloat s_2, MyBigFloat s_3,
				     MyBigFloat s_1_prime, MyBigFloat s_2_prime, MyBigFloat s_3_prime,
				     MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
				     MyBigFloat dt, MyBigFloat sigma);
MyBigFloat get_max_caustic_density(MyBigFloat s_1, MyBigFloat s_2,
				   MyBigFloat s_1_prime, MyBigFloat s_2_prime,
				   MyBigFloat second_deriv, MyBigFloat second_deriv_prime,
				   MyBigFloat sigma, int pindex);
void get_current_ps_info(int pindex, MyBigFloat * flde, MyBigFloat * psde);
void do_phase_space_drift(int i, double dt_drift);
void do_the_phase_space_kick(int i, double dt_gravkick);
void do_long_range_phase_space_kick(int i, double dt_gravkick);
/* some math functions we need from phasespace_math.c */
void ludcmp(MyBigFloat ** a, int n, int *indx, MyBigFloat * d);
MyBigFloat **matrix(long nrl, long nrh, long ncl, long nch);
MyBigFloat *vector(long nl, long nh);
void free_matrix(MyBigFloat ** m, long nrl, long nrh, long ncl, long nch);
void free_vector(MyBigFloat * v, long nl, long nh);
void mult_matrix(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension, MyBigFloat ** matrix_result);
void mult_matrix_transpose_A(MyBigFloat ** matrix_a, MyBigFloat ** matrix_b, int dimension,
			     MyBigFloat ** matrix_result);
void luinvert(MyBigFloat ** input_matrix, int n, MyBigFloat ** inverse_matrix);
void eigsrt(MyBigFloat d[], MyBigFloat ** v, int n);
void jacobi(MyBigFloat ** a, int n, MyBigFloat d[], MyBigFloat ** v, int *nrot);
#ifdef PMGRID
void pmtidaltensor_periodic_diff(void);
void pmtidaltensor_periodic_fourier(int component);
int pmtidaltensor_nonperiodic_diff(int grnr);
int pmtidaltensor_nonperiodic_fourier(int component, int grnr);
void check_tidaltensor_periodic(int particle_ID);
void check_tidaltensor_nonperiodic(int particle_ID);
#endif
#endif


#ifdef JD_DPP
void compute_Dpp(int ipart);
#endif


#ifdef SCFPOTENTIAL
void SCF_do_center_of_mass_correction(double fac_rad, double start_rad, double fac_part, int max_iter);
void SCF_write(int task);
void SCF_calc_from_random(long *seed);
void SCF_calc_from_particles(void);
void SCF_init(void);
void SCF_reset(void);
void SCF_free(void);
void SCF_evaluate(MyDouble x, MyDouble y, MyDouble z, MyDouble * potential, MyDouble * ax, MyDouble * ay,
		  MyDouble * az);
void SCF_collect_update(void);

void sphere_acc(double x, double y, double z, double *xa, double *ya, double *za);
void to_unit(double x, double y, double z, double *xs, double *ys, double *zs);
double ran1(long *idum);
double gasdev(long *idum);
double factrl(int n);
int nlm_all(int num, int n, int l, int m);
int nlm(int n, int l, int m);
int nl(int n, int l);
int lm(int l, int m);
int kdelta(int a, int b);
double gnlm_var(int n, int l, int m);
double hnlm_var(int n, int l, int m);
#endif


#ifdef MOL_CLOUDS
int do_mol_clouds();
void init_mol_clouds();
#endif

#if !defined(LT_STELLAREVOLUTION) && defined(GM_MUPPI)
double get_cooling_time_fromT(double, double, double, double *);
#endif


#ifdef LMB_SPECTRAL_CRs

// init
void init_spectral_crs_bounds();
void init_spectral_crs_arrays(int i);
void init_spectral_crs_pressure(int i, MyFloat Pth);
// read from ic
void init_spectral_crs_from_file(int i);

// injection
void lmb_cr_sn_injection(int i);
void spectral_crs_shock_injection(int i);


// main loop
void evolve_spectral_crs();
void inject_crs_postprocess();

// adiabatic changes
void spectral_crs_adiabatic_changes(int i, int mode);

// conduction
void spectral_crs_art_cond_update(int i, double dt);

// timestep constraint
//MyAtLeastDouble compute_bp_cr_timestep_from_cooling_times(int n);

// debugging
void write_cr_file_headers();
void report_cr_properties(int i);

// tests
void report_internal_model_for_spectral_crs();

// Spectral CR diffusion
void spectral_cr_diffusion();
void cr_diffusion_matrix_multiply(double *in_e, double *out_e);
double cr_diffusion_vector_multiply(double *a, double *b);
int cr_diffusion_evaluate(int target, int mode, int *ngblist, int extended_numngb_inbox, int *i_startnode,
			  int node, int i_node, void *data);

void *ar_cr_diffusion_evaluate_primary(void *p);
void *ar_cr_diffusion_evaluate_secondary(void *p);

#ifdef LMB_CR_OUTPUT_SYNCHROTRON
// synchrotron emission
void synchrotron_emissivity();

// test
void report_synchrotron_kernel();
#endif

#endif // LMB_SPECTRAL_CRs

#ifdef SIDM
void sidm_Init_Particles(void);
void sidm_Init_CrossSection(void);
void sidm_DoScatter(void);
inline MyDouble sidm_cross_sigma(MyDouble rel_vel);
#endif


int ngb_treefind_variable_threads_check(MyLongDouble searchcenter[3], MyFloat hsml, int target,
					int *startnode, int mode);
int ngb_treefind_variable_threads_noexport(MyLongDouble searchcenter[3], MyFloat hsml, int target,
					   int *startnode, int mode, int *ngblist, int debug);
int test_ngb_inconsitency(int *ngblist, int numngb, int *check_ngblist, MyLongDouble * local_pos,
			  MyFloat hsml, int startnode, int mode, int target);


#ifdef ADAPTGRAVSOFT
void ags_setup_smoothinglengths(void);
int ags_ngb_treefind_variable_threads(MyLongDouble searchcenter[3], MyFloat hsml, int target, int *startnode,
				      int mode, int *exportflag, int *exportnodecount, int *exportindex,
				      int *ngblist);
void ags_density(void);
int ags_density_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			 int *ngblist);
void *ags_density_evaluate_primary(void *p);
void *ags_density_evaluate_secondary(void *p);
int ags_density_isactive(int n);
void ags_force_update_hmax(void);
#endif

#ifdef ALTERNATIVE_PSORT
void init_sort_ID(MyIDType * data, int ndata);
#endif

#ifdef STATICBRANDT
inline double OmegaR(int i, int mode);
#endif

#ifdef DISKPOT
void gravity_tree_subfunc_diskpot(void);
#endif

#ifdef _OPENMP
void build_active_particle_list(void);
#endif

#ifdef WRITE_KEY_FILES
int key_compare(const void *a, const void *b);
int key_list_compare(const void *a, const void *b);
void write_key_index_file(int num, int filenr);
#endif

#ifdef PHIDOT
void phidot(void);
int phidot_treeevaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex);
int phidot_treeevaluate_ewald_correction(int target, int mode, int *exportflag, int *exportnodecount,
					 int *exportindex);
void *phidot_primary_loop(void *p);
void *phidot_secondary_loop(void *p);
double phidot_ewald_corr(double dx, double dy, double dz, double vx, double vy, double vz);
#endif

#ifdef ACC_PROTO
void on_beginning_of_timestep();
void on_end_of_timestep();

#endif

void omp_mysort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
void serial_sort_omp(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));
void omp_qsort(void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));

void domain_diffuse();

#ifdef GM_STARDENSITY
void density_stars(void);
int density_evaluate_stars(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
			   int *ngblist);
void *density_evaluate_primary_stars(void *p);
void *density_evaluate_secondary_stars(void *p);
int density_isactive_stars(int n);
int ngb_treefind_variable_threads_stars(MyLongDouble searchcenter[3], MyFloat hsml, int target,
					int *startnode, int mode, int *exportflag, int *exportnodecount,
					int *exportindex, int *ngblist);
#ifdef GM_REPOSITION_ON_STARDENSITY_MAX
void blackhole_stardens_reposition(int);
#endif
#endif


#ifdef GL_CR_DUST
void evolve_dust(void);
void compute_ism_properties_effective_model(double, double, double,
			    double, double *, double *, double *, double *, double *);
#endif


#endif //PROTO_H
