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

#ifdef WITH_OPENVDB
#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#endif

namespace common {
	namespace vdb {
		namespace sparse {

#ifdef WITH_OPENVDB
			// =====================================================================
			// VoxelOpenVDBManager Implementation
			// =====================================================================

			VoxelOpenVDBManager::VoxelOpenVDBManager() : insert_count(0), common::vdb::VoxelSparseManager() {
				// Initialize OpenVDB library (thread-safe, only initializes once)
				openvdb::initialize();
				
				// Create a new OpenVDB FloatGrid
				openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
				grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				grid->setName("density");
				
				// Store as void* to avoid header dependency
				openvdb_grid = new openvdb::FloatGrid::Ptr(grid);
			}

			VoxelOpenVDBManager::~VoxelOpenVDBManager() {
				// Clean up the stored grid pointer
				if (openvdb_grid) {
					delete static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
					openvdb_grid = nullptr;
				}
			}

			void VoxelOpenVDBManager::init(unsigned int expected_voxels) {
				// Clear and re-initialize the OpenVDB grid
				if (openvdb_grid) {
					delete static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
				}
				
				openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
				grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				grid->setName("density");
				
				// Set transform if available
				if (transform_scale > 0.0) {
					openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
					grid->setTransform(transform);
				}
				
				openvdb_grid = new openvdb::FloatGrid::Ptr(grid);
				insert_count = 0;
			}

			void VoxelOpenVDBManager::insertOrUpdatePackedSequential(uint64_t key, float value) {
				int x, y, z;
				unpackCoord3(key, x, y, z);

				openvdb::FloatGrid::Ptr grid = *static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
				auto accessor = grid->getAccessor();
				openvdb::Coord coord(x, y, z);

				// OpenVDB supports accumulation: getValue + setValue
				float current_value = accessor.getValue(coord);
				accessor.setValue(coord, current_value + value);
				insert_count++;
			}

			void VoxelOpenVDBManager::insertOrUpdate(Voxel* h_voxels, int num_voxels) {
				if (num_voxels <= 0) return;

				auto start_total = omp_get_wtime();

				// OPTIMIZATION: Use thread-local hash maps to pre-aggregate voxels
				// This drastically reduces OpenVDB accessor contention
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

				// Sequential insertion into OpenVDB grid (no contention)
				auto start_grid_insert = omp_get_wtime();
				openvdb::FloatGrid::Ptr grid = *static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
				auto accessor = grid->getAccessor();
				
				for (const auto& kv : merged) {
					int x, y, z;
					unpackCoord3(kv.first, x, y, z);
					openvdb::Coord coord(x, y, z);

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

				printf("[VoxelOpenVDBManager::insertOrUpdate] Pre-agg: %.3f ms, Merge: %.3f ms, Grid insert: %.3f ms, Total: %.3f ms\n",
					preagg_duration, merge_duration, grid_insert_duration, total_duration);
			}

			int VoxelOpenVDBManager::extractAll(Voxel** h_output_voxels) {
				auto start_total = omp_get_wtime();

				// Collect all active voxels from OpenVDB grid
				std::vector<Voxel> voxels;

				openvdb::FloatGrid::Ptr grid = *static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
				
				// Use OpenVDB's value iterator to efficiently iterate over active voxels
				for (openvdb::FloatGrid::ValueOnCIter iter = grid->cbeginValueOn(); iter; ++iter) {
					openvdb::Coord coord = iter.getCoord();
					float value = *iter;
					if (value != 0.0f) {
						voxels.push_back(Voxel(coord[0], coord[1], coord[2], value));
					}
				}

				// Allocate output array
				int total_count = static_cast<int>(voxels.size());
				if (total_count > 0) {
					*h_output_voxels = new Voxel[total_count];
					std::copy(voxels.begin(), voxels.end(), *h_output_voxels);
				} else {
					*h_output_voxels = nullptr;
				}

				auto end_total = omp_get_wtime();
				auto total_duration = (end_total - start_total) * 1000.0;
				printf("[VoxelOpenVDBManager::extractAll] Total: %.3f ms\n", total_duration);

				return total_count;
			}

			void VoxelOpenVDBManager::clear() {
				// Create a new empty grid
				if (openvdb_grid) {
					delete static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
				}
				
				openvdb::FloatGrid::Ptr grid = openvdb::FloatGrid::create(0.0f);
				grid->setGridClass(openvdb::GRID_FOG_VOLUME);
				grid->setName("density");
				
				if (transform_scale > 0.0) {
					openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
					grid->setTransform(transform);
				}
				
				openvdb_grid = new openvdb::FloatGrid::Ptr(grid);
				insert_count = 0;
			}

			void VoxelOpenVDBManager::serialize(uint8_t* bin_data) {
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

			void VoxelOpenVDBManager::deserialize(uint8_t* bin_data) {
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
					openvdb::FloatGrid::Ptr grid = *static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
					auto accessor = grid->getAccessor();
					const Voxel* voxels = reinterpret_cast<const Voxel*>(ptr);

					for (int i = 0; i < count; i++) {
						const Voxel& v = voxels[i];
						openvdb::Coord coord(v.i, v.j, v.k);
						accessor.setValue(coord, v.value);
					}
				}
			}

			void VoxelOpenVDBManager::merge(common::vdb::VoxelSparseManager* _other) {
				// Attempt to dynamic_cast to VoxelOpenVDBManager
				VoxelOpenVDBManager* other = dynamic_cast<VoxelOpenVDBManager*>(_other);
				if (!other) {
					printf("[VoxelOpenVDBManager::merge] Error: Incompatible manager type\n");
					return;
				}

			if (!other->openvdb_grid) {
				return;
			}

			openvdb::FloatGrid::Ptr grid_dst = *static_cast<openvdb::FloatGrid::Ptr*>(openvdb_grid);
			openvdb::FloatGrid::Ptr grid_src = *static_cast<openvdb::FloatGrid::Ptr*>(other->openvdb_grid);

			// Use OpenVDB's optimized composite operation for efficient parallel merging
			// This is significantly faster than manual iteration
			openvdb::tools::compSum(*grid_dst, *grid_src);

			insert_count += other->insert_count;
			}
#endif // WITH_OPENVDB

		} // sparse
	} // vdb
} // common
