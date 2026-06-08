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

#include <cmath>
#include <vector>
// float4 is a CUDA built-in type; provide a lightweight fallback for
// non-CUDA (e.g. IDE IntelliSense) translation units so the header parses.
#ifdef __CUDACC__
#  include <cuda_runtime.h>
#else
#  ifndef __CUDA_RUNTIME_H__
// Plain-C struct matching CUDA's layout; the real definition wins if
// cuda_runtime.h is included before this header.
struct float4 { float x, y, z, w; };
#  endif
#endif
#include "data_common.h"

namespace utility {
	namespace cudakdtree {

        void run_knn(float* points, size_t N, int k, std::vector<float>& radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_gpu, bool use_cycling, float max_radius, common::SpaceData::DenseType& rho_kernel);

        /**
         * @brief Build a float4 KD-tree from particle positions and store into out_tree.
         *
         * Each output node stores (offset-adjusted world position, w=bit-cast particle index).
         * Calls cukd::buildTree_host so the resulting array is a left-balanced KD-tree.
         *
         * @param positions     Flat array of raw positions [x0,y0,z0, x1,y1,z1, ...].
         * @param ids           Per-entry original particle indices.
         * @param N             Number of particles.
         * @param offset        3-element offset subtracted from positions.
         * @param out_tree      Output tree (resized and filled by this function).
         */
        void build_float4_kdtree(
            const float*  positions,
            const size_t* ids,
            size_t        N,
            const float*  offset,
            std::vector<float4>& out_tree
        );

	}// cudakdtree
} //utility