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
#include <limits>
#include <cstdint>

#ifdef WITH_OPENMP
#include <omp.h>
#endif

namespace space_converter {
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
                particle_radius_const = 0.0;
                use_gpu = false;
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
                particles_reader_id_per_ptype.clear();
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
                particle_radius_const = 0.0;
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

                    // ids are compact per-type indices — sort them by the radius of the
                    // particle they refer to (ascending). Composes with other sorts.
                    std::sort(ids.begin(), ids.end(),
                        [&radii](size_t a, size_t b) {
                            return radii[a] < radii[b];
                        });
                }
            }

            // Helper function to compute Morton code (Z-order curve) for 3D position
            // This interleaves the bits of x, y, z coordinates to create a spatial ordering
            static inline uint64_t compute_morton_code_3d(float x, float y, float z, float min_coord, float max_coord) {
                // Normalize coordinates to [0, 1] range
                float range = max_coord - min_coord;
                if (range <= 0.0f) range = 1.0f;

                float nx = (x - min_coord) / range;
                float ny = (y - min_coord) / range;
                float nz = (z - min_coord) / range;

                // Clamp to [0, 1]
                nx = std::max(0.0f, std::min(1.0f, nx));
                ny = std::max(0.0f, std::min(1.0f, ny));
                nz = std::max(0.0f, std::min(1.0f, nz));

                // Convert to 21-bit integers (total 63 bits for Morton code)
                uint64_t ix = static_cast<uint64_t>(nx * ((1 << 21) - 1));
                uint64_t iy = static_cast<uint64_t>(ny * ((1 << 21) - 1));
                uint64_t iz = static_cast<uint64_t>(nz * ((1 << 21) - 1));

                // Interleave bits
                uint64_t morton = 0;
                for (int i = 0; i < 21; ++i) {
                    morton |= ((ix & (1ULL << i)) << (2 * i));
                    morton |= ((iy & (1ULL << i)) << (2 * i + 1));
                    morton |= ((iz & (1ULL << i)) << (2 * i + 2));
                }

                return morton;
            }

            void CacheManager::sort_particles_by_nonoverlap_cpu() {
                // Sort particles using Morton codes (Z-order curve) to ensure spatial coherence
                // This minimizes overlapping voxel regions, reducing/eliminating atomic operations

#ifdef WITH_OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
                for (int ptype = 0; ptype < static_cast<int>(pos_particles_per_ptype.size()); ++ptype) {
                    if (ptype >= static_cast<int>(particles_id_ordered_per_ptype.size())) {
                        continue;
                    }

                    auto& positions = pos_particles_per_ptype[ptype];
                    auto& ids = particles_id_ordered_per_ptype[ptype];

                    if (positions.empty() || ids.empty()) {
                        continue;
                    }

                    size_t n_particles = ids.size();
                    if (positions.size() < n_particles * 3) {
                        continue; // Invalid data
                    }

                    // Find bounding box for position normalization (ids are compact
                    // per-type indices into the positions array)
                    float min_coord = std::numeric_limits<float>::max();
                    float max_coord = std::numeric_limits<float>::lowest();

                    for (size_t i = 0; i < n_particles; ++i) {
                        size_t idx = i * 3;
                        for (int dim = 0; dim < 3; ++dim) {
                            float coord = positions[idx + dim];
                            min_coord = std::min(min_coord, coord);
                            max_coord = std::max(max_coord, coord);
                        }
                    }

                    // Morton code per compact particle slot
                    std::vector<uint64_t> morton_codes(n_particles);
                    for (size_t i = 0; i < n_particles; ++i) {
                        size_t idx = i * 3;
                        morton_codes[i] = compute_morton_code_3d(
                            positions[idx + 0], positions[idx + 1], positions[idx + 2],
                            min_coord, max_coord);
                    }

                    // Sort the iteration order by the Morton code of the referenced particle
                    std::sort(ids.begin(), ids.end(),
                        [&morton_codes](size_t a, size_t b) {
                            return morton_codes[a] < morton_codes[b];
                        });
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

            void CacheManager::sort_particles_by_nonoverlap_gpu_inplace() {
                // Fallback to CPU sorting
                sort_particles_by_nonoverlap_cpu();
            }
#endif

        } // namespace cache
    } // namespace common
} // namespace space_converter