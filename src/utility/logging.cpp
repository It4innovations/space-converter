/*
 * Copyright(C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
 *
 * This program is free software : you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#ifdef WITH_MERIC
#	include <meric.h>
#endif

#include <map>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <cstdio>
#include <string>
#include <mpi.h>

namespace space_converter {

    std::map<std::string, double> logging_measurements;
    int logging_world_rank;

    void LOG_Init(int world_rank) {
        logging_world_rank = world_rank;
    }

    void LOG_MeasureStart(const char* name) {
    #ifdef WITH_MERIC
        MERIC_MeasureStart(name);
    #endif

        MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes start measurement at the same time
        
    #ifdef WITH_OPENMP
        logging_measurements[name] = omp_get_wtime();
    #endif
    }

    void LOG_MeasureStop(const char* name) {
    #ifdef WITH_MERIC
        MERIC_MeasureStop(name);
    #endif

        if (logging_measurements.find(name) == logging_measurements.end()) {
            printf("Measurement '%s' was not started.\n", name);
            return;
        }

        double start_time = logging_measurements[name];

        MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes stop measurement at the same time
    
        double end_time = 0.0;
    #ifdef WITH_OPENMP
        end_time = omp_get_wtime();
    #endif

        if (logging_world_rank == 0) {
            printf("rank #%d: '%s' elapsed time: %f [s]\n", logging_world_rank, name, end_time - start_time);
        }
    }
}