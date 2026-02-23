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

#include "dense_common.h"
#include "gpu_utility.h"
#include <algorithm>
#include <cfloat>
#include <cstring>
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
// File-local kernel helpers (anonymous namespace → no external linkage)
// ─────────────────────────────────────────────────────────────────────────────
namespace {

	/**
	 * Atomic float min via CAS – works for all finite float values.
	 */
	__device__ __forceinline__ void atomicMinFloat(float* addr, float val) {
		int* addr_i = reinterpret_cast<int*>(addr);
		int  old    = *addr_i, assumed;
		do {
			assumed = old;
			old = atomicCAS(addr_i, assumed,
				__float_as_int(fminf(val, __int_as_float(assumed))));
		} while (assumed != old);
	}

	/**
	 * Atomic float max via CAS – works for all finite float values.
	 */
	__device__ __forceinline__ void atomicMaxFloat(float* addr, float val) {
		int* addr_i = reinterpret_cast<int*>(addr);
		int  old    = *addr_i, assumed;
		do {
			assumed = old;
			old = atomicCAS(addr_i, assumed,
				__float_as_int(fmaxf(val, __int_as_float(assumed))));
		} while (assumed != old);
	}

	/**
	 * Single-pass min/max reduction kernel.
	 *
	 * Each thread accumulates a private min/max over a grid-stride loop, then
	 * the block reduces with warp-shuffle intrinsics and a shared-memory
	 * inter-warp stage before writing a single atomicCAS per block to global
	 * output scalars.
	 *
	 * Assumptions:
	 *   - blockDim.x is a multiple of 32 (warp size).
	 *   - blockDim.x <= 1024  (shared memory arrays sized for max 32 warps).
	 */
	__global__ void find_min_max_kernel(
		const float* __restrict__ data,
		size_t        n,
		float*        d_min,
		float*        d_max)
	{
		float local_min =  FLT_MAX;
		float local_max = -FLT_MAX;

		// Grid-stride loop: every thread covers multiple elements
		for (size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
			 i < n;
			 i += static_cast<size_t>(gridDim.x) * blockDim.x)
		{
			const float v = data[i];
			local_min = fminf(local_min, v);
			local_max = fmaxf(local_max, v);
		}

		// ── Warp-level reduction (shuffle) ────────────────────────────────────
		constexpr unsigned int FULL_MASK = 0xFFFFFFFFu;
#pragma unroll
		for (int delta = 16; delta > 0; delta >>= 1) {
			local_min = fminf(local_min, __shfl_xor_sync(FULL_MASK, local_min, delta));
			local_max = fmaxf(local_max, __shfl_xor_sync(FULL_MASK, local_max, delta));
		}

		// ── Inter-warp reduction via shared memory ────────────────────────────
		__shared__ float s_min[32];  // one slot per warp (max 1024/32 = 32 warps)
		__shared__ float s_max[32];

		const int lane = static_cast<int>(threadIdx.x) & 31;
		const int wid  = static_cast<int>(threadIdx.x) >> 5;

		if (lane == 0) {
			s_min[wid] = local_min;
			s_max[wid] = local_max;
		}
		__syncthreads();

		// Number of active warps in this block
		const int n_warps = (static_cast<int>(blockDim.x) + 31) >> 5;

		// Only the first warp reduces the per-warp results
		if (wid == 0) {
			local_min = (lane < n_warps) ? s_min[lane] :  FLT_MAX;
			local_max = (lane < n_warps) ? s_max[lane] : -FLT_MAX;
#pragma unroll
			for (int delta = 16; delta > 0; delta >>= 1) {
				local_min = fminf(local_min, __shfl_xor_sync(FULL_MASK, local_min, delta));
				local_max = fmaxf(local_max, __shfl_xor_sync(FULL_MASK, local_max, delta));
			}
			if (lane == 0) {
				atomicMinFloat(d_min, local_min);
				atomicMaxFloat(d_max, local_max);
			}
		}
	}

} // anonymous namespace

namespace common {
	namespace vdb {
		namespace dense {

			// 
			// Internal helper
			// 

			//static void CUDA_CHECK_ERROR(cudaError_t err, const char* context) {
			//	if (err != cudaSuccess)
			//		throw std::runtime_error(std::string(context) + ": " + cudaGetErrorString(err));
			//}

			// 
			// VoxelGPUDenseManager  host lifecycle methods
			// 

			__host__ void VoxelGPUDenseManager::clear() {
				free_device();
				data_density.clear();
#ifndef WITH_NO_DATA_TEMP
				data_temp.clear();
#endif
				memset(dims, 0, 3 * sizeof(size_t));
				memset(offset, 0, 3 * sizeof(size_t));
			}

			__host__ void VoxelGPUDenseManager::create(size_t x, size_t y, size_t z) {
				dims[0] = x;  dims[1] = y;  dims[2] = z;

				data_density.resize(size());
				memset(data_density.data(), 0, memsize());
#ifndef WITH_NO_DATA_TEMP
				data_temp.resize(size());
				memset(data_temp.data(), 0, memsize());
#endif
			}

			// 
			// VoxelGPUDenseManager  GPU memory management
			// 

			__host__ void VoxelGPUDenseManager::to_device() {
				const size_t n = size();
				if (n == 0)
					return;

				// d_data_density
				if (d_data_density == nullptr)
					CUDA_CHECK_ERROR(cudaMalloc(&d_data_density, memsize()));
				CUDA_CHECK_ERROR(cudaMemcpy(d_data_density, data_density.data(), memsize(),
					cudaMemcpyHostToDevice));

#ifndef WITH_NO_DATA_TEMP
				// d_data_temp
				if (d_data_temp == nullptr)
					CUDA_CHECK_ERROR(cudaMalloc(&d_data_temp, memsize()));
				CUDA_CHECK_ERROR(cudaMemcpy(d_data_temp, data_temp.data(), memsize(),
					cudaMemcpyHostToDevice));
#endif

				// d_dims
				if (d_dims == nullptr)
					CUDA_CHECK_ERROR(cudaMalloc(&d_dims, 3 * sizeof(size_t)));
				CUDA_CHECK_ERROR(cudaMemcpy(d_dims, dims, 3 * sizeof(size_t),
					cudaMemcpyHostToDevice));

				// d_offset
				if (d_offset == nullptr)
					CUDA_CHECK_ERROR(cudaMalloc(&d_offset, 3 * sizeof(size_t)));
				CUDA_CHECK_ERROR(cudaMemcpy(d_offset, offset, 3 * sizeof(size_t),
					cudaMemcpyHostToDevice));

				// d_particle_count
				if (d_particle_count == nullptr)
					CUDA_CHECK_ERROR(cudaMalloc(&d_particle_count, sizeof(uint64_t)));
				CUDA_CHECK_ERROR(cudaMemset(d_particle_count, 0, sizeof(uint64_t)));
			}

			__host__ void VoxelGPUDenseManager::from_device() {
				const size_t n = size();
				if (n == 0)
					return;

				if (d_data_density != nullptr) {
					if (data_density.size() != n) data_density.resize(n);
					CUDA_CHECK_ERROR(cudaMemcpy(data_density.data(), d_data_density, memsize(),
						cudaMemcpyDeviceToHost));
				}

#ifndef WITH_NO_DATA_TEMP
				if (d_data_temp != nullptr) {
					if (data_temp.size() != n) data_temp.resize(n);
					CUDA_CHECK_ERROR(cudaMemcpy(data_temp.data(), d_data_temp, memsize(),
						cudaMemcpyDeviceToHost));
				}
#endif
			}

			__host__ void VoxelGPUDenseManager::free_device() {
				if (d_data_density) { cudaFree(d_data_density); d_data_density = nullptr; }
#ifndef WITH_NO_DATA_TEMP
				if (d_data_temp) { cudaFree(d_data_temp);    d_data_temp = nullptr; }
#endif
				if (d_dims) { cudaFree(d_dims);         d_dims = nullptr; }
				if (d_offset) { cudaFree(d_offset);       d_offset = nullptr; }
				if (d_particle_count) { cudaFree(d_particle_count); d_particle_count = nullptr; }
			}

			// 
			// VoxelGPUDenseManager  find_min_max
			// 

			__host__ void VoxelGPUDenseManager::find_min_max(float& min_value, float& max_value) {
				const size_t n = size();
				if (n == 0) {
					min_value = 0.0f;
					max_value = 0.0f;
					return;
				}

				if (d_data_density == nullptr)
					throw std::runtime_error("find_min_max: d_data_density is null; call to_device() first");

				// Allocate and initialise the two scalar result buffers on the device
				float* d_min = nullptr;
				float* d_max = nullptr;
				const float init_min =  FLT_MAX;
				const float init_max = -FLT_MAX;

				CUDA_CHECK_ERROR(cudaMalloc(&d_min, sizeof(float)));
				CUDA_CHECK_ERROR(cudaMalloc(&d_max, sizeof(float)));
				CUDA_CHECK_ERROR(cudaMemcpy(d_min, &init_min, sizeof(float),
					cudaMemcpyHostToDevice));
				CUDA_CHECK_ERROR(cudaMemcpy(d_max, &init_max, sizeof(float),
					cudaMemcpyHostToDevice));

				// Launch: 256 threads/block, up to 1024 blocks
				constexpr int THREADS = 256;
				const int blocks = static_cast<int>(
					std::min<size_t>((n + THREADS - 1) / THREADS, 1024));

				find_min_max_kernel<<<blocks, THREADS>>>(d_data_density, n, d_min, d_max);
				CUDA_CHECK_ERROR(cudaGetLastError());
				CUDA_CHECK_ERROR(cudaDeviceSynchronize());

				// Copy results back to host
				CUDA_CHECK_ERROR(cudaMemcpy(&min_value, d_min, sizeof(float),
					cudaMemcpyDeviceToHost));
				CUDA_CHECK_ERROR(cudaMemcpy(&max_value, d_max, sizeof(float),
					cudaMemcpyDeviceToHost));

				CUDA_CHECK_ERROR(cudaFree(d_min));
				CUDA_CHECK_ERROR(cudaFree(d_max));
			}

		} // dense
	} // vdb
} // common
