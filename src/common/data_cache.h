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

#pragma once

#include <vector>

namespace common {
    namespace cache {
        // CacheManager is responsible for managing cached data, including GPU caches and other temporary storage.
        class CacheManager {
        public:
            // Constructor and destructor
            CacheManager();
            ~CacheManager();

            // Initialize the cache manager, setting up necessary resources.
            void initialize();

            // Clear all cached data, freeing up resources.
            void clear_cache();

        public:
            // CPU data
            std::vector< std::vector<float> > pos_particles_per_ptype;     ///< Particle positions per type
            std::vector< std::vector<size_t> > particles_id_ordered_per_ptype;     ///< Particle positions ordered by radius

			std::vector< std::vector<float> > radius_particles_per_ptype;  ///< Particle radii per type
			std::vector<size_t> particles_ptype_offset;                    ///< Particle count offsets per type
			std::vector< std::vector<float> > rho_particles_per_ptype;     ///< Density values per type
            std::vector< std::vector<float> > mass_particles_per_ptype;     ///< Mass values per type

            std::vector<float> values_particles; ///< Particle values for extracted type (e.g., density, temperature)
			int cached_value_ptype = -1; ///< Cached particle type for values_particles (to track which type's values are currently cached)
			int cached_value_block_name_id = -1; ///< Cached block name ID for values_particles (to track which block's values are currently cached)

            // GPU data pointers
            std::vector<float*> d_pos_particles_per_ptype;              ///< GPU particle positions per type
            std::vector<size_t*> d_particles_id_ordered_per_ptype;      ///< GPU particle IDs ordered by radius per type
            std::vector<float*> d_radius_particles_per_ptype;           ///< GPU particle radii per type
            size_t* d_particles_ptype_offset;                           ///< GPU particle count offsets per type
            std::vector<float*> d_rho_particles_per_ptype;              ///< GPU density values per type
            std::vector<float*> d_mass_particles_per_ptype;             ///< GPU mass values per type

            float* d_values_particles;                    ///< GPU particle values for extracted type (e.g., density, temperature)

            // Others
            double radius_particle_const = 0.0; ///< Constant particle radius value

            /// Flag to enable GPU acceleration for CUDA-based computations
            bool use_gpu_cuda = false;

        public:

            // Methods for GPU data transfer
            void copy_particles_to_gpu();                               ///< Copy all particle data from CPU to GPU
            void copy_values_to_gpu();                                  ///< Copy particle values from CPU to GPU
            void free_gpu_memory();                                     ///< Free all GPU memory allocations

            // Sorting methods
            void sort_particles_by_radius_cpu();                        ///< Sort particle IDs by radius on CPU (ascending)
            //void sort_particles_by_radius_gpu();                        ///< Sort particle IDs by radius on GPU (ascending)
            void sort_particles_by_radius_gpu_inplace();                ///< Sort particle IDs by radius on GPU using device pointers (no CPU<->GPU copy)
        };
    }

} // namespace common