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

#include "convert_vdb_kernel.h"
#include <float.h>

#ifdef WITH_OPENMP
# include <omp.h>
#endif

namespace common {
    namespace vdb {
        namespace kernel {
            
            /**
             * @brief CPU implementation of bounding box calculation using OpenMP
             * 
             * Parallel reduction to compute min/max coordinates across all particles.
             */
            void find_bbox_cpu(
                const float* pos_particles,
                size_t num_particles,
                const float* offset_position,
                float* bbox_min,
                float* bbox_max
            ) {
                if (num_particles == 0) {
                    bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0f;
                    bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0f;
                    return;
                }
                
                float min_x = FLT_MAX, min_y = FLT_MAX, min_z = FLT_MAX;
                float max_x = -FLT_MAX, max_y = -FLT_MAX, max_z = -FLT_MAX;
                
                // Parallel reduction to find min/max coordinates
#ifdef WITH_OPENMP
#pragma omp parallel for reduction(min : min_x, min_y, min_z) reduction(max : max_x, max_y, max_z)
#endif
                for (size_t i = 0; i < num_particles; ++i) {
                    // Load position (x, y, z are interleaved: [x0, y0, z0, x1, y1, z1, ...])
                    float px = pos_particles[i * 3 + 0] - offset_position[0];
                    float py = pos_particles[i * 3 + 1] - offset_position[1];
                    float pz = pos_particles[i * 3 + 2] - offset_position[2];
                    
                    if (px < min_x) min_x = px;
                    if (py < min_y) min_y = py;
                    if (pz < min_z) min_z = pz;
                    
                    if (px > max_x) max_x = px;
                    if (py > max_y) max_y = py;
                    if (pz > max_z) max_z = pz;
                }
                
                bbox_min[0] = min_x;
                bbox_min[1] = min_y;
                bbox_min[2] = min_z;
                
                bbox_max[0] = max_x;
                bbox_max[1] = max_y;
                bbox_max[2] = max_z;
            }
            
        } // namespace kernel
    } // namespace vdb
} // namespace common