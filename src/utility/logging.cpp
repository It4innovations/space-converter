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
#include <cstdlib>
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

        // NOTE: no MPI_Barrier here (same reason as in LOG_MeasureStop below):
        // measured sections entered by a subset of ranks (save_raw_volume, ...)
        // would desynchronize the barrier counts and deadlock the other ranks.

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

        // NOTE: no MPI_Barrier here. A logging call must not synchronize: several
        // measured sections (save_raw_volume, save_raw_particles_to_vdb, ...) are
        // entered by rank 0 only, so a barrier inside Stop mismatched the barrier
        // counts across ranks and deadlocked multi-rank dense/--dense-file runs.
        // The reported time is therefore rank 0's local elapsed time.

        double end_time = 0.0;
    #ifdef WITH_OPENMP
        end_time = omp_get_wtime();
    #endif

        if (logging_world_rank == 0) {
            printf("rank #%d: '%s' elapsed time: %f [s]\n", logging_world_rank, name, end_time - start_time);
        }
    }

    bool LOG_Verbose() {
        static const bool verbose = [] {
            const char* v = getenv("SPACE_CONVERTER_VERBOSE");
            return v != nullptr && v[0] != '\0' && !(v[0] == '0' && v[1] == '\0');
        }();
        return verbose;
    }
}