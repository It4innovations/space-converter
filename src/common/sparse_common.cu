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

#include "sparse_common.h"
#include <cuda_runtime.h>
#include <cub/cub.cuh>

#include "data_common.h"

namespace common {
	namespace vdb {
		namespace sparse {

			// ---------------------------------------------
			// Kernel: Voxels -> (keys, values)
			// ---------------------------------------------
			__global__ void voxelsToKeyValue(const Voxel* __restrict__ in,
				int n,
				uint64_t* __restrict__ keys,
				float* __restrict__ vals)
			{
				int tid = blockIdx.x * blockDim.x + threadIdx.x;
				if (tid >= n) return;

				Voxel v = in[tid];
				keys[tid] = packCoord3(v.i, v.j, v.k);
				vals[tid] = v.value;
			}

			// ---------------------------------------------
			// Kernel: Pixels -> (keys, values)
			// ---------------------------------------------			
			__global__ void pixelsToKeyValue(const Voxel* __restrict__ in,
				int n,
				uint64_t* __restrict__ keys,
				float* __restrict__ vals)
			{
				int tid = blockIdx.x * blockDim.x + threadIdx.x;
				if (tid >= n) return;

				Voxel v = in[tid];
				keys[tid] = packCoord3(v.i, v.j, v.k);
				vals[tid] = v.value;
			}

			// ---------------------------------------------
			// Kernel: (keys, values) -> Voxels (unique output)
			// ---------------------------------------------
			__global__ void keyValueToVoxels(const uint64_t* __restrict__ keys,
				const float* __restrict__ vals,
				int n,
				Voxel* __restrict__ out)
			{
				int tid = blockIdx.x * blockDim.x + threadIdx.x;
				if (tid >= n) return;

				int x, y, z;
				unpackCoord3(keys[tid], x, y, z);
				out[tid] = Voxel(x, y, z, vals[tid]);
			}

			// CustomSum functor
			struct CustomSum
			{
				template <typename T>
				__device__ __forceinline__
				T operator()(const T &a, const T &b) const {
					return a + b;
				}
			};

			// ---------------------------------------------
			// Fast manager: per-batch accumulate via sort+reduce
			// ---------------------------------------------
			
			VoxelGPUManagerSortReduce::VoxelGPUManagerSortReduce(size_t max_particles)
				: m_max(max_particles)
			{
					// cudaMalloc(&d_inVoxels,   m_max * sizeof(Voxel));

					cudaMalloc(&d_keys, m_max * sizeof(uint64_t));
					cudaMalloc(&d_vals, m_max * sizeof(float));

					cudaMalloc(&d_keys_alt, m_max * sizeof(uint64_t));
					cudaMalloc(&d_vals_alt, m_max * sizeof(float));

					// Worst-case: all unique => output size == input size
					cudaMalloc(&d_keys_out, m_max * sizeof(uint64_t));
					cudaMalloc(&d_vals_out, m_max * sizeof(float));
					cudaMalloc(&d_num_out, sizeof(int));

					// Precompute temp storage sizes for CUB primitives (max)
					size_t sort_bytes = 0;
					cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes,
						d_keys, d_keys_alt,
						d_vals, d_vals_alt,
						(int)m_max);
					m_sort_temp_bytes = sort_bytes;
					cudaMalloc(&d_sort_temp, m_sort_temp_bytes);

					size_t reduce_bytes = 0;
					CustomSum op_sum;
					cub::DeviceReduce::ReduceByKey(nullptr, reduce_bytes,
						d_keys_alt, d_keys_out,
						d_vals_alt, d_vals_out,
						d_num_out,
						op_sum,
						(int)m_max);
				m_reduce_temp_bytes = reduce_bytes;
				cudaMalloc(&d_reduce_temp, m_reduce_temp_bytes);
			}

			VoxelGPUManagerSortReduce::~VoxelGPUManagerSortReduce()
			{
					// cudaFree(d_inVoxels);

					cudaFree(d_keys);
					cudaFree(d_vals);

					cudaFree(d_keys_alt);
					cudaFree(d_vals_alt);

					cudaFree(d_keys_out);
					cudaFree(d_vals_out);

					cudaFree(d_num_out);

				cudaFree(d_sort_temp);
				cudaFree(d_reduce_temp);
			}

			// Accumulate a batch of voxels: output becomes unique (i,j,k) with summed values
			// Returns number of unique voxels in the batch.
			int VoxelGPUManagerSortReduce::insertOrUpdate(const Voxel* h_voxels, int num_voxels, bool print_timing)
			{
					if (num_voxels <= 0) {
						m_last_count = 0;
						return 0;
					}
					if ((size_t)num_voxels > m_max) {
						printf("[VoxelGPUManagerSortReduce] ERROR: num_voxels=%d exceeds max=%zu\n",
							num_voxels, m_max);
						return 0;
					}

					cudaEvent_t e0, e1, e2, e3, e4;
					cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventCreate(&e2);
					cudaEventCreate(&e3); cudaEventCreate(&e4);

					cudaEventRecord(e0);

					// H2D
					cudaMemcpy(d_inVoxels, h_voxels, num_voxels * sizeof(Voxel), cudaMemcpyHostToDevice);
					cudaEventRecord(e1);

					// map -> key/value
					{
						int block = 256;
						int grid = (num_voxels + block - 1) / block;
						voxelsToKeyValue << <grid, block >> > (d_inVoxels, num_voxels, d_keys, d_vals);
					}
					cudaEventRecord(e2);

					// sort by key (pairs)
					cub::DeviceRadixSort::SortPairs(d_sort_temp, m_sort_temp_bytes,
						d_keys, d_keys_alt,
						d_vals, d_vals_alt,
						num_voxels);
					cudaEventRecord(e3);

					// reduce-by-key (sum values for identical keys)
					CustomSum op_sum;
					cub::DeviceReduce::ReduceByKey(d_reduce_temp, m_reduce_temp_bytes,
						d_keys_alt, d_keys_out,
						d_vals_alt, d_vals_out,
						d_num_out,
						op_sum,
						num_voxels);

					int h_num_out = 0;
					cudaMemcpy(&h_num_out, d_num_out, sizeof(int), cudaMemcpyDeviceToHost);

					if (print_timing) {
						float t01 = 0, t12 = 0, t23 = 0, t34 = 0, t04 = 0;
						cudaEventElapsedTime(&t01, e0, e1); // H2D
						cudaEventElapsedTime(&t12, e1, e2); // map
						cudaEventElapsedTime(&t23, e2, e3); // sort
						cudaEventElapsedTime(&t34, e3, e4); // reduce
						cudaEventElapsedTime(&t04, e0, e4); // total
						printf("[SortReduce] H2D: %.3f ms | map: %.3f ms | sort: %.3f ms | reduce: %.3f ms | total: %.3f ms | unique=%d\n",
							t01, t12, t23, t34, t04, h_num_out);
					}

					cudaEventDestroy(e0); cudaEventDestroy(e1); cudaEventDestroy(e2);
					cudaEventDestroy(e3); cudaEventDestroy(e4);

				m_last_count = h_num_out;
				return h_num_out;
			}

			// Extract the last accumulated unique voxels back to host
			int VoxelGPUManagerSortReduce::extractAll(Voxel** h_output_voxels, bool print_timing)
			{
					const int n = m_last_count;
					if (n <= 0) {
						*h_output_voxels = nullptr;
						return 0;
					}

					cudaEvent_t e0, e1, e2;
					cudaEventCreate(&e0); cudaEventCreate(&e1); cudaEventCreate(&e2);

					cudaEventRecord(e0);

					Voxel* d_out_voxels = nullptr;
					cudaMalloc(&d_out_voxels, n * sizeof(Voxel));

					{
						int block = 256;
						int grid = (n + block - 1) / block;
						keyValueToVoxels << <grid, block >> > (d_keys_out, d_vals_out, n, d_out_voxels);
					}

					cudaEventRecord(e1);

					*h_output_voxels = new Voxel[n];
					cudaMemcpy(*h_output_voxels, d_out_voxels, n * sizeof(Voxel), cudaMemcpyDeviceToHost);

					cudaEventRecord(e2);
					cudaDeviceSynchronize();

					if (print_timing) {
						float t01 = 0, t12 = 0, t02 = 0;
						cudaEventElapsedTime(&t01, e0, e1); // pack->voxel kernel
						cudaEventElapsedTime(&t12, e1, e2); // D2H copy
						cudaEventElapsedTime(&t02, e0, e2); // total
						printf("[SortReduce extractAll] kernel: %.3f ms | D2H: %.3f ms | total: %.3f ms\n",
							t01, t12, t02);
					}

					cudaFree(d_out_voxels);

				cudaEventDestroy(e0); cudaEventDestroy(e1); cudaEventDestroy(e2);

				return n;
			}
		}
	}
}