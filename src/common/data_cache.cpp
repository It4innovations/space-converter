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

#include "data_cache.h"
#include <algorithm>
#include <numeric>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace common {
    namespace cache {    

        CacheManager::CacheManager() {
            // Initialize GPU pointers to nullptr
            d_particles_ptype_offset = nullptr;
            d_values_particles = nullptr;
            
            // Initialize GPU vector pointers (will be empty vectors)
            // d_pos_particles_per_ptype, d_particles_id_ordered_per_ptype, etc. are std::vectors
            // and will be automatically initialized as empty
            
            // Initialize scalar values
            radius_particle_const = 0.0;
            use_gpu_cuda = false;
        }

        CacheManager::~CacheManager() {
            // Free GPU memory before destruction
            free_gpu_memory();
        }

        void CacheManager::initialize() {
            // Initialize cache manager
        }

        void CacheManager::clear_cache() {
            // Clear all cached data
            free_gpu_memory();
            
            // Clear CPU data
            pos_particles_per_ptype.clear();
            particles_id_ordered_per_ptype.clear();
            radius_particles_per_ptype.clear();
            particles_ptype_offset.clear();
            rho_particles_per_ptype.clear();
            mass_particles_per_ptype.clear();
            values_particles.clear();
            
            // Clear GPU vector pointers
            d_pos_particles_per_ptype.clear();
            d_particles_id_ordered_per_ptype.clear();
            d_radius_particles_per_ptype.clear();
            d_rho_particles_per_ptype.clear();
            d_mass_particles_per_ptype.clear();
            
            // Reset scalar values
            radius_particle_const = 0.0;
        }

        void CacheManager::sort_particles_by_radius_cpu() {

            // Sort particle IDs by radius for each particle type (OpenMP parallel version)
#ifdef WITH_OPENMP
            #pragma omp parallel for schedule(dynamic)
#endif
            for (int ptype = 0; ptype < static_cast<int>(radius_particles_per_ptype.size()); ++ptype) {
                if (ptype >= static_cast<int>(particles_id_ordered_per_ptype.size())) {
                    continue;
                }

                auto& radii = radius_particles_per_ptype[ptype];
                auto& ids = particles_id_ordered_per_ptype[ptype];

                if (radii.size() != ids.size() || radii.empty()) {
                    continue;
                }

                // Create index array
                std::vector<size_t> indices(ids.size());
                std::iota(indices.begin(), indices.end(), 0);

                // Sort indices by radius values (ascending)
                std::sort(indices.begin(), indices.end(),
                    [&radii](size_t i1, size_t i2) {
                        return radii[i1] < radii[i2];
                    });

                // Reorder particle IDs according to sorted indices
                std::vector<size_t> sorted_ids(ids.size());
                for (size_t i = 0; i < indices.size(); ++i) {
                    sorted_ids[i] = ids[indices[i]];
                }
                ids = std::move(sorted_ids);
            }
        }

#ifndef WITH_GPU_CUDA
        // Stub implementation for CPU-only builds
        void CacheManager::free_gpu_memory() {
            // No GPU memory to free in CPU-only mode
        }

        void CacheManager::copy_particles_to_gpu() {
            // Empty stub for CPU-only builds
        }

        void CacheManager::copy_values_to_gpu() {
            // Empty stub for CPU-only builds
        }

        //void CacheManager::sort_particles_by_radius_gpu() {
        //    // Fallback to CPU sorting
        //    sort_particles_by_radius_cpu();
        //}

        void CacheManager::sort_particles_by_radius_gpu_inplace() {
            // Fallback to CPU sorting
            sort_particles_by_radius_cpu();
        }
#endif

    } // namespace cache
} // namespace common