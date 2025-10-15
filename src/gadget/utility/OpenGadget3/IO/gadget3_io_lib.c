#include "gadget3_io_lib.h"

#include <mpi.h>

#include "../CodeBase/allvars.h"
#include "../CodeBase/proto.h"

#include <string.h>

void init_lib(int world_rank, int world_size)
{
	//MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
	//MPI_Comm_size(MPI_COMM_WORLD, &NTask);
	ThisTask = world_rank;
	NTask = world_size;

	MPI_Comm_dup(MPI_COMM_WORLD, &MYMPI_COMM_WORLD);
}

void finish_lib()
{
	MPI_Comm_free(&MYMPI_COMM_WORLD);
}

void set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor, int TotBHs)
{
	memset(&All, 0, sizeof(struct global_data_all_processes));
	All.ICFormat = ICFormat;
	All.SnapFormat = SnapFormat;
	All.NumFilesWrittenInParallel = NumFilesWrittenInParallel;
	All.MaxMemSize = MaxMemSize;
	All.BufferSize = BufferSize;
	All.PartAllocFactor = PartAllocFactor;
#ifdef LT_STELLAREVOLUTION
	All.SFfactor = PartAllocFactor;
#endif

	// BH
	if (TotBHs > 0) {
		All.TotBHs = TotBHs;
		All.BHfactor = PartAllocFactor;
	}	

	/* System of units */
	All.UnitLength_in_cm = 3.0856780e+24;//  3.085678e21 for kpc / h 3.085678e24 for Mpc / h
	All.UnitMass_in_g = 1.989e43;  //1.0e10 Msun / h
	All.UnitVelocity_in_cm_per_s = 1e5;  //1 km / sec
	All.GravityConstantInternal = 0;

	/* Cosmological parameters */
	All.Omega0 = 0.3158;
	All.OmegaLambda = 0.6842;
	All.OmegaBaryon = 0.0494;
	All.HubbleParam = 0.6732;
	//All.BoxSize  500.0

	/* initialize OpenMP thread pool and bind (implicitly though OpenMP runtime) */
#ifdef _OPENMP
#pragma omp parallel
	{
		ThisThread = omp_get_thread_num();
#pragma omp master
		{
			maxThreads = omp_get_num_threads();
		}
}
#else
	ThisThread = 0;
	maxThreads = 1;
#endif

	RestartFlag = 2; //Restart from specified snapshot dump and continue simulation
	RestartSnapNum = 0;

	/* initialize CPU-time/Wallclock-time measurement */
	for (int i = 0; i < CPU_PARTS; i++)
		All.CPU_Sum[i] = CPU_Step[i] = 0;

	CPUThisRun = 0;
	WallclockTime = second();
}