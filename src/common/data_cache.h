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
            std::vector< std::vector<float> > pos_particles_per_ptype;     ///< Particle positions per type
            std::vector< std::vector<size_t> > particles_id_ordered_per_ptype;     ///< Particle positions ordered by radius

			std::vector< std::vector<float> > radius_particles_per_ptype;  ///< Particle radii per type
			std::vector<size_t> particles_ptype_offset;                    ///< Particle count offsets per type
			std::vector< std::vector<float> > rho_particles_per_ptype;     ///< Density values per type
            std::vector< std::vector<float> > mass_particles_per_ptype;     ///< Mass values per type
        };
    }

} // namespace common