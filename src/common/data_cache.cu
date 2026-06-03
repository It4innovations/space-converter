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
 #include "../utility/gpu_utility.h"
 #include <cuda_runtime.h>
 #include <thrust/device_vector.h>
 #include <thrust/sort.h>
 #include <thrust/copy.h>
 #include <thrust/execution_policy.h>
 #include <iostream>

namespace common {
    namespace cache {

        void CacheManager::copy_values_to_gpu() {
            // Free values_particles
            if (d_values_particles != nullptr) {
                CUDA_CHECK_ERROR(cudaFree(d_values_particles));
                d_values_particles = nullptr;
            }

			// Allocate and copy values_particles
            if (!values_particles.empty()) {
                size_t size_bytes = values_particles.size() * sizeof(float);
                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_values_particles, size_bytes));
                CUDA_CHECK_ERROR(cudaMemcpy(d_values_particles, values_particles.data(), 
                         size_bytes, cudaMemcpyHostToDevice));
            } else {
                d_values_particles = nullptr;
            }
            // Check for CUDA errors
			CUDA_CHECK_LAST_ERROR();
        }
        
        void CacheManager::copy_particles_to_gpu() {
            // Free existing GPU memory if any
            free_gpu_memory();

            // Allocate and copy pos_particles_per_ptype
            d_pos_particles_per_ptype.resize(pos_particles_per_ptype.size());
            for (size_t i = 0; i < pos_particles_per_ptype.size(); ++i) {
                size_t size_bytes = pos_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(CUDA_MALLOC(&d_pos_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_pos_particles_per_ptype[i], pos_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_pos_particles_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy particles_id_ordered_per_ptype
            d_particles_id_ordered_per_ptype.resize(particles_id_ordered_per_ptype.size());
            for (size_t i = 0; i < particles_id_ordered_per_ptype.size(); ++i) {
                size_t size_bytes = particles_id_ordered_per_ptype[i].size() * sizeof(size_t);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(CUDA_MALLOC(&d_particles_id_ordered_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_particles_id_ordered_per_ptype[i], particles_id_ordered_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_particles_id_ordered_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy radius_particles_per_ptype
            d_radius_particles_per_ptype.resize(radius_particles_per_ptype.size());
            for (size_t i = 0; i < radius_particles_per_ptype.size(); ++i) {
                size_t size_bytes = radius_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(CUDA_MALLOC(&d_radius_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_radius_particles_per_ptype[i], radius_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_radius_particles_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy particles_ptype_offset
            d_particles_ptype_offset = nullptr;
            if (!particles_ptype_offset.empty()) {
                size_t size_bytes = particles_ptype_offset.size() * sizeof(size_t);
                CUDA_CHECK_ERROR(CUDA_MALLOC(&d_particles_ptype_offset, size_bytes));
                CUDA_CHECK_ERROR(cudaMemcpy(d_particles_ptype_offset, particles_ptype_offset.data(), 
                         size_bytes, cudaMemcpyHostToDevice));
            }

            // Allocate and copy rho_particles_per_ptype
            d_rho_particles_per_ptype.resize(rho_particles_per_ptype.size());
            for (size_t i = 0; i < rho_particles_per_ptype.size(); ++i) {
                size_t size_bytes = rho_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(CUDA_MALLOC(&d_rho_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_rho_particles_per_ptype[i], rho_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_rho_particles_per_ptype[i] = nullptr;
                }
            }

            // Allocate and copy mass_particles_per_ptype
            d_mass_particles_per_ptype.resize(mass_particles_per_ptype.size());
            for (size_t i = 0; i < mass_particles_per_ptype.size(); ++i) {
                size_t size_bytes = mass_particles_per_ptype[i].size() * sizeof(float);
                if (size_bytes > 0) {
                    CUDA_CHECK_ERROR(CUDA_MALLOC(&d_mass_particles_per_ptype[i], size_bytes));
                    CUDA_CHECK_ERROR(cudaMemcpy(d_mass_particles_per_ptype[i], mass_particles_per_ptype[i].data(), 
                             size_bytes, cudaMemcpyHostToDevice));
                } else {
                    d_mass_particles_per_ptype[i] = nullptr;
                }
            }

            // Check for CUDA errors
            CUDA_CHECK_LAST_ERROR();
        }

        void CacheManager::free_gpu_memory() {
            // Free pos_particles_per_ptype
            for (auto ptr : d_pos_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_pos_particles_per_ptype.clear();

            // Free particles_id_ordered_per_ptype
            for (auto ptr : d_particles_id_ordered_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_particles_id_ordered_per_ptype.clear();

            // Free radius_particles_per_ptype
            for (auto ptr : d_radius_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_radius_particles_per_ptype.clear();

            // Free particles_ptype_offset
            if (d_particles_ptype_offset != nullptr) {
                CUDA_CHECK_ERROR(cudaFree(d_particles_ptype_offset));
                d_particles_ptype_offset = nullptr;
            }

            // Free rho_particles_per_ptype
            for (auto ptr : d_rho_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_rho_particles_per_ptype.clear();

            // Free mass_particles_per_ptype
            for (auto ptr : d_mass_particles_per_ptype) {
                if (ptr != nullptr) {
                    CUDA_CHECK_ERROR(cudaFree(ptr));
                }
            }
            d_mass_particles_per_ptype.clear();

            // Free values_particles
            if (d_values_particles != nullptr) {
                CUDA_CHECK_ERROR(cudaFree(d_values_particles));
                d_values_particles = nullptr;
            }

            // Check for CUDA errors
            CUDA_CHECK_LAST_ERROR();
        }

        void CacheManager::sort_particles_by_radius_gpu_inplace() {
            // Sort particle IDs in-place by radius without modifying radii array
            // Zero additional memory allocation - optimal performance
            
            for (size_t ptype = 0; ptype < d_radius_particles_per_ptype.size(); ++ptype) {
                if (ptype >= d_particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                float* d_radii_ptr = d_radius_particles_per_ptype[ptype];
                size_t* d_ids_ptr = d_particles_id_ordered_per_ptype[ptype];

                if (d_radii_ptr == nullptr || d_ids_ptr == nullptr) {
                    continue;
                }

                // Get the size from the CPU vectors
                if (ptype >= radius_particles_per_ptype.size() || 
                    ptype >= particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                size_t n = radius_particles_per_ptype[ptype].size();
                if (n != particles_id_ordered_per_ptype[ptype].size() || n == 0) {
                    continue;
                }

                // Wrap device pointers
                thrust::device_ptr<float> d_radii_thrust(d_radii_ptr);
                thrust::device_ptr<size_t> d_ids_thrust(d_ids_ptr);

                // Sort ids in-place based on radii values (radii array remains unchanged)
                // d_ids_ptr already contains [0, 1, 2, ..., n-1]
                thrust::sort(d_ids_thrust, d_ids_thrust + n,
                    [d_radii_thrust] __device__ (size_t i, size_t j) {
                        return d_radii_thrust[i] < d_radii_thrust[j];
                    });
            }

            // Check for CUDA errors
            //CUDA_CHECK_LAST_ERROR();
            CUDA_SYNC_CHECK();
        }

        // Device function to compute Morton code (Z-order curve) for 3D position
        __device__ inline uint64_t compute_morton_code_3d_gpu(float x, float y, float z, float min_coord, float max_coord) {
            // Normalize coordinates to [0, 1] range
            float range = max_coord - min_coord;
            if (range <= 0.0f) range = 1.0f;
            
            float nx = (x - min_coord) / range;
            float ny = (y - min_coord) / range;
            float nz = (z - min_coord) / range;
            
            // Clamp to [0, 1]
            nx = fmaxf(0.0f, fminf(1.0f, nx));
            ny = fmaxf(0.0f, fminf(1.0f, ny));
            nz = fmaxf(0.0f, fminf(1.0f, nz));
            
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

        void CacheManager::sort_particles_by_nonoverlap_gpu_inplace() {
            // Sort particles using Morton codes (Z-order curve) to ensure spatial coherence
            // This minimizes overlapping voxel regions, reducing/eliminating atomic operations
            
            for (size_t ptype = 0; ptype < d_pos_particles_per_ptype.size(); ++ptype) {
                if (ptype >= d_particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                float* d_pos_ptr = d_pos_particles_per_ptype[ptype];
                size_t* d_ids_ptr = d_particles_id_ordered_per_ptype[ptype];

                if (d_pos_ptr == nullptr || d_ids_ptr == nullptr) {
                    continue;
                }

                // Get the size from the CPU vectors
                if (ptype >= pos_particles_per_ptype.size() || 
                    ptype >= particles_id_ordered_per_ptype.size()) {
                    continue;
                }

                size_t n_particles = particles_id_ordered_per_ptype[ptype].size();
                if (n_particles == 0 || pos_particles_per_ptype[ptype].size() < n_particles * 3) {
                    continue;
                }

                // Find bounding box on CPU (or use pre-computed values)
                const auto& positions = pos_particles_per_ptype[ptype];
                const auto& ids = particles_id_ordered_per_ptype[ptype];
                
                float min_coord = std::numeric_limits<float>::max();
                float max_coord = std::numeric_limits<float>::lowest();
                
                for (size_t i = 0; i < n_particles; ++i) {
                    size_t idx = ids[i] * 3;
                    for (int dim = 0; dim < 3; ++dim) {
                        float coord = positions[idx + dim];
                        min_coord = std::min(min_coord, coord);
                        max_coord = std::max(max_coord, coord);
                    }
                }

                // Wrap device pointers
                thrust::device_ptr<float> d_pos_thrust(d_pos_ptr);
                thrust::device_ptr<size_t> d_ids_thrust(d_ids_ptr);

                // Sort ids by Morton codes computed from positions
                thrust::sort(d_ids_thrust, d_ids_thrust + n_particles,
                    [d_pos_thrust, min_coord, max_coord] __device__ (size_t i, size_t j) {
                        // Get positions for particle i
                        float xi = d_pos_thrust[i * 3 + 0];
                        float yi = d_pos_thrust[i * 3 + 1];
                        float zi = d_pos_thrust[i * 3 + 2];
                        
                        // Get positions for particle j
                        float xj = d_pos_thrust[j * 3 + 0];
                        float yj = d_pos_thrust[j * 3 + 1];
                        float zj = d_pos_thrust[j * 3 + 2];
                        
                        // Compute Morton codes
                        uint64_t morton_i = compute_morton_code_3d_gpu(xi, yi, zi, min_coord, max_coord);
                        uint64_t morton_j = compute_morton_code_3d_gpu(xj, yj, zj, min_coord, max_coord);
                        
                        return morton_i < morton_j;
                    });
            }

            // Check for CUDA errors
            CUDA_SYNC_CHECK();
        }
    }
}