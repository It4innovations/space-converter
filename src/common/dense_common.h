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

 // data_common.h defines VoxelDenseManager (base class) and raw_common.h.
 // Do NOT create a reverse include of dense_common.h inside data_common.h.
#include "data_common.h"

namespace common {
	namespace vdb {
		namespace dense {

			/**
			 * @brief CPU-only concrete implementation of VoxelDenseManager.
			 *
			 * Overrides clear() and create() to confirm the host-side (CPU) path
			 * with explicit, non-virtual-dispatched implementations.  No device
			 * allocations are performed.
			 */
			class VoxelCPUDenseManager : public VoxelDenseManager {
			public:
				/**
				 * @brief Clear all host grid data and reset dimensions.
				 */
				void clear() override;

				/**
				 * @brief Allocate and zero-initialise host grid buffers.
				 * @param x Width   @param y Height   @param z Depth
				 * @param allocate_data_temp  If true, allocate data_temp buffer for normalization.
				 */
				void create(size_t x, size_t y, size_t z, bool allocate_data_temp = false) override;
			};

#ifdef WITH_GPU_CUDA
			/**
			 * @brief CPU + GPU concrete implementation of VoxelDenseManager.
			 *
			 * Extends the host buffers in VoxelDenseManager with matching device
			 * (CUDA) allocations.  Provides to_device() / from_device() / free_device()
			 * to manage GPU memory explicitly.
			 */
			class VoxelGPUDenseManager : public VoxelDenseManager {
			public:
				// ── Device (GPU) data ─────────────────────────────────────────────
				float* d_data_density = nullptr;  ///< Device pointer to density data
				float* d_data_temp = nullptr;  ///< Device pointer to temp accumulation buffer (only allocated when needed)
				size_t* d_dims = nullptr;  ///< Device pointer to grid dimensions
				size_t* d_offset = nullptr;  ///< Device pointer to grid offset

				/// Atomic counter: number of particles that passed should_process
				uint64_t* d_particle_count = nullptr;

				// ── Lifecycle ─────────────────────────────────────────────────────

				/**
				 * @brief Free all device memory and clear all host grid data.
				 */
				void clear() override;

				/**
				 * @brief Allocate/zero-initialise host grid buffers.
				 * @param x Width   @param y Height   @param z Depth
				 * @param allocate_data_temp  If true, allocate data_temp buffer for normalization.
				 */
				void create(size_t x, size_t y, size_t z, bool allocate_data_temp = false) override;

				// ── GPU memory management ─────────────────────────────────────────

				/**
				 * @brief Allocate device memory and upload all host data to the GPU.
				 *
				 * Allocates d_data_density, d_data_temp, d_dims, d_offset on the device
				 * (reusing existing allocations when the size matches) and copies the
				 * current host data to them via cudaMemcpy.
				 */
				void to_device();

				/**
				 * @brief Download GPU results back to host vectors.
				 *
				 * Copies d_data_density (and d_data_temp if enabled) from device back
				 * to the corresponding host std::vector.
				 */
				void from_device();

				/**
				 * @brief Free all device memory allocations.
				 *
				 * Resets all device pointers to nullptr after calling cudaFree.
				 */
				void free_device();

				/**
				 * @brief Compute the minimum and maximum values of the density data on the GPU.
				 *
				 * Uses a warp-shuffle + shared-memory parallel reduction followed by global
				 * atomic updates.  Requires that to_device() has already been called so that
				 * d_data_density is valid.  Throws std::runtime_error on CUDA errors or if
				 * the device buffer has not been allocated.
				 *
				 * @param min_value  Output: minimum voxel density value.
				 * @param max_value  Output: maximum voxel density value.
				 */
				void find_min_max(float& min_value, float& max_value);
			};
#endif
		} // dense
	} // vdb
} // common
