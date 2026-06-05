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
#include <cstdint>
#include <unordered_map>
#include "convert_vdb.h"

#include <omp.h>

namespace common {
	namespace vdb {
		namespace sparse {

#ifdef WITH_NANOVDB
			// =====================================================================
			// VoxelNanoVDBManager Implementation
			// =====================================================================

			VoxelNanoVDBManager::VoxelNanoVDBManager() : insert_count(0), common::vdb::VoxelSparseManager() {
				nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
			}

			VoxelNanoVDBManager::~VoxelNanoVDBManager() {
				// NanoVDB grid is managed by shared_ptr, no manual cleanup needed
			}

			void VoxelNanoVDBManager::init(unsigned int expected_voxels) {
				// Clear and re-initialize the NanoVDB grid
				nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
				insert_count = 0;

				// Set transform if available
				if (transform_scale > 0.0) {
					nano_grid->setTransform(transform_scale, nanovdb::Vec3d(0, 0, 0));
				}
			}

			void VoxelNanoVDBManager::insertOrUpdatePackedSequential(uint64_t key, float value) {
				int x, y, z;
				unpackCoord3(key, x, y, z);

				auto accessor = nano_grid->getAccessor();
				nanovdb::Coord coord(x, y, z);

				// NanoVDB GridBuilder supports accumulation: getValue + setValue
				float current_value = accessor.getValue(coord);
				accessor.setValue(coord, current_value + value);
				insert_count++;
			}

			void VoxelNanoVDBManager::insertOrUpdate(Voxel* h_voxels, int num_voxels) {
				if (num_voxels <= 0) return;

				auto start_total = omp_get_wtime();

				// OPTIMIZATION: Use thread-local hash maps to pre-aggregate voxels
				// This drastically reduces NanoVDB accessor contention
				auto start_preagg = omp_get_wtime();
				const int thread_count = omp_get_max_threads();
				std::vector<std::unordered_map<uint64_t, float>> thread_maps(thread_count);

#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					auto& local_map = thread_maps[tid];
					local_map.reserve((size_t)num_voxels / (size_t)thread_count + 64);

#pragma omp for schedule(static)
					for (int idx = 0; idx < num_voxels; idx++) {
						const Voxel& v = h_voxels[idx];
						const uint64_t key = packCoord3(v.i, v.j, v.k);
						local_map[key] += v.value;
					}
				}
				auto end_preagg = omp_get_wtime();

				// Merge thread-local maps into single map
				auto start_merge = omp_get_wtime();
				size_t total_local_unique = 0;
				for (const auto& m : thread_maps) {
					total_local_unique += m.size();
				}
				std::unordered_map<uint64_t, float> merged;
				merged.reserve(total_local_unique + total_local_unique / 4 + 1);

				for (const auto& m : thread_maps) {
					for (const auto& kv : m) {
						merged[kv.first] += kv.second;
					}
				}
				auto end_merge = omp_get_wtime();

				// Sequential insertion into NanoVDB grid (no contention)
				auto start_grid_insert = omp_get_wtime();
				auto accessor = nano_grid->getAccessor();
				for (const auto& kv : merged) {
					int x, y, z;
					unpackCoord3(kv.first, x, y, z);
					nanovdb::Coord coord(x, y, z);

					// Accumulate values
					float current_value = accessor.getValue(coord);
					accessor.setValue(coord, current_value + kv.second);
					insert_count++;
				}
				auto end_grid_insert = omp_get_wtime();
				auto end_total = omp_get_wtime();

				auto preagg_duration = (end_preagg - start_preagg) * 1000.0;
				auto merge_duration = (end_merge - start_merge) * 1000.0;
				auto grid_insert_duration = (end_grid_insert - start_grid_insert) * 1000.0;
				auto total_duration = (end_total - start_total) * 1000.0;

				printf("[VoxelNanoVDBManager::insertOrUpdate] Pre-agg: %.3f ms, Merge: %.3f ms, Grid insert: %.3f ms, Total: %.3f ms\n",
					preagg_duration, merge_duration, grid_insert_duration, total_duration);
			}

			int VoxelNanoVDBManager::extractAll(Voxel** h_output_voxels) {
				auto start_total = omp_get_wtime();

				// Collect all active voxels from NanoVDB grid
				std::vector<Voxel> voxels;

				// NanoVDB accessor to iterate through active voxels
				auto accessor = nano_grid->getAccessor();

				// Note: GridBuilder doesn't provide a convenient iterator over active values
				// We need to iterate through the bounding box
				// This is less efficient than hash table extraction but necessary with NanoVDB

				// Get bounding box
				auto bbox = nano_grid->getBBox();
				if (bbox.empty()) {
					*h_output_voxels = nullptr;
					auto end_total = omp_get_wtime();
					printf("[VoxelNanoVDBManager::extractAll] Total: %.3f ms (empty grid)\n", (end_total - start_total) * 1000.0);
					return 0;
				}

				// Iterate through bounding box and collect non-zero values
				for (auto iter = bbox.begin(); iter; ++iter) {
					nanovdb::Coord coord = *iter;
					float value = accessor.getValue(coord);
					if (value != 0.0f) {
						voxels.push_back(Voxel(coord[0], coord[1], coord[2], value));
					}
				}

				// Allocate output array
				int total_count = static_cast<int>(voxels.size());
				*h_output_voxels = new Voxel[total_count];
				std::copy(voxels.begin(), voxels.end(), *h_output_voxels);

				auto end_total = omp_get_wtime();
				auto total_duration = (end_total - start_total) * 1000.0;
				printf("[VoxelNanoVDBManager::extractAll] Total: %.3f ms\n", total_duration);

				return total_count;
			}

			void VoxelNanoVDBManager::clear() {
				// Create a new empty grid
				nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
				if (transform_scale > 0.0) {
					nano_grid->setTransform(transform_scale, nanovdb::Vec3d(0, 0, 0));
				}
				insert_count = 0;
			}

			void VoxelNanoVDBManager::serialize(uint8_t* bin_data) {
				// Serialize: insert_count + voxel data
				// We'll extract all voxels and serialize them as (count, [i,j,k,value]...)

				uint8_t* ptr = bin_data;

				// Write insert_count
				memcpy(ptr, &insert_count, sizeof(int));
				ptr += sizeof(int);

				// Extract all voxels
				Voxel* voxels = nullptr;
				int count = extractAll(&voxels);

				// Write actual count
				memcpy(ptr, &count, sizeof(int));
				ptr += sizeof(int);

				// Write voxel data
				if (count > 0 && voxels != nullptr) {
					memcpy(ptr, voxels, sizeof(Voxel) * count);
					delete[] voxels;
				}
			}

			void VoxelNanoVDBManager::deserialize(uint8_t* bin_data) {
				const uint8_t* ptr = bin_data;

				// Read insert_count
				memcpy(&insert_count, ptr, sizeof(int));
				ptr += sizeof(int);

				// Read voxel count
				int count;
				memcpy(&count, ptr, sizeof(int));
				ptr += sizeof(int);

				// Clear current grid
				clear();

				// Read voxel data and insert into grid
				if (count > 0) {
					auto accessor = nano_grid->getAccessor();
					const Voxel* voxels = reinterpret_cast<const Voxel*>(ptr);

					for (int i = 0; i < count; i++) {
						const Voxel& v = voxels[i];
						nanovdb::Coord coord(v.i, v.j, v.k);
						accessor.setValue(coord, v.value);
					}
				}
			}

			void VoxelNanoVDBManager::merge(common::vdb::VoxelSparseManager* _other) {
				// Attempt to dynamic_cast to VoxelNanoVDBManager
				VoxelNanoVDBManager* other = dynamic_cast<VoxelNanoVDBManager*>(_other);
				if (!other) {
					printf("[VoxelNanoVDBManager::merge] Error: Incompatible manager type\n");
					return;
				}

			if (!other->nano_grid) {
				return;
			}

			// Get accessor for destination grid
			auto acc_dst = nano_grid->getAccessor();

			// Iterate directly over active voxels in source grid using tree structure
			// This avoids the overhead of extractAll() and intermediate allocations
			auto bbox = other->nano_grid->getBBox();
			if (bbox.empty()) {
				return;
			}

			auto acc_src = other->nano_grid->getAccessor();
			int merged_count = 0;

			// Iterate through bounding box and merge non-zero values
			for (auto iter = bbox.begin(); iter; ++iter) {
				nanovdb::Coord xyz = *iter;
				float v = acc_src.getValue(xyz);

				if (v != 0.0f) {
					// Accumulate: add source value to destination value
					float current_value = acc_dst.getValue(xyz);
					acc_dst.setValue(xyz, current_value + v);
					merged_count++;
				}
			}

			insert_count += merged_count;
		}

		size_t VoxelNanoVDBManager::mem_size() const {
			// Rough estimate: insert_count + count + voxel data
			// For serialization we need: sizeof(int) + sizeof(int) + sizeof(Voxel) * actual_count
			// But we don't know actual_count without extracting
			// As an estimate, use insert_count (may be over-estimate due to overwrites)
			return sizeof(int) + sizeof(int) + sizeof(Voxel) * insert_count;
		}
#endif // WITH_NANOVDB

		} // sparse
	} // vdb
} // common
