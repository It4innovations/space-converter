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

// Constants
#define EMPTY_KEY -999999
#define HASH_TABLE_LOAD_FACTOR 0.5f

namespace common {
	namespace vdb {
		namespace sparse {

			// Hash function for 3D coordinates (CPU version)
			inline unsigned int VoxelOpenMPManager::hash3D_cpu(int i, int j, int k) const {
				unsigned int hash = 73856093u * i ^ 19349663u * j ^ 83492791u * k;
				return hash % table_size;
			}

		VoxelOpenMPManager::VoxelOpenMPManager(): VoxelCPUManager(), table_size(0), insert_count(0), hash_table(nullptr) {
		}
		
		VoxelOpenMPManager::~VoxelOpenMPManager() {
			if (hash_table != nullptr) {
				delete[] hash_table;
			}
		}
		
		void VoxelOpenMPManager::init(unsigned int expected_voxels) {
			// Size hash table with load factor consideration
			table_size = (unsigned int)(expected_voxels / HASH_TABLE_LOAD_FACTOR);
			insert_count = 0;
			
			// Allocate hash table on CPU
			if (hash_table != nullptr) {
				delete[] hash_table;
			}
			hash_table = new VoxelHashEntry[table_size];
			
			// Initialize hash table using OpenMP
			#pragma omp parallel for
			for (unsigned int i = 0; i < table_size; i++) {
				hash_table[i].i = EMPTY_KEY;
				hash_table[i].j = EMPTY_KEY;
				hash_table[i].k = EMPTY_KEY;
				hash_table[i].value = 0.0f;
				hash_table[i].occupied = 0;
			}
		}
		
		// Helper for sequential insertion
		void VoxelOpenMPManager::insertOrUpdatePackedSequential(uint64_t key, float value) {
					int x, y, z;
					unpackCoord3(key, x, y, z);
					unsigned int slot = hash3D_cpu(x, y, z);

					for (unsigned int probe = 0; probe < table_size; probe++) {
						if (hash_table[slot].occupied == 0) {
							hash_table[slot].i = x;
							hash_table[slot].j = y;
							hash_table[slot].k = z;
							hash_table[slot].value = value;
							hash_table[slot].occupied = 1;
							insert_count++;
							return;
						}

						if (hash_table[slot].i == x &&
							hash_table[slot].j == y &&
							hash_table[slot].k == z) {
							hash_table[slot].value += value;
							return;
						}

				}
			}
			
			// Insert or update voxels using OpenMP parallelization with thread-local pre-aggregation
			void VoxelOpenMPManager::insertOrUpdate(Voxel* h_voxels, int num_voxels) {
					if (num_voxels <= 0) return;

					auto start_total = omp_get_wtime();

					// OPTIMIZATION: Use thread-local hash maps to pre-aggregate voxels
					// This drastically reduces hash table contention
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

					// Sequential insertion into hash table (no contention)
					auto start_hash_insert = omp_get_wtime();
					for (const auto& kv : merged) {
						insertOrUpdatePackedSequential(kv.first, kv.second);
					}
					auto end_hash_insert = omp_get_wtime();
					auto end_total = omp_get_wtime();
					auto preagg_duration = (end_preagg - start_preagg) * 1000.0;
					auto merge_duration = (end_merge - start_merge) * 1000.0;
					auto hash_insert_duration = (end_hash_insert - start_hash_insert) * 1000.0;
					auto total_duration = (end_total - start_total) * 1000.0;

					printf("[insertOrUpdate] Pre-agg: %.3f ms, Merge: %.3f ms, Hash insert: %.3f ms, Total: %.3f ms\n",
						preagg_duration, merge_duration,
						hash_insert_duration, total_duration);
			}
			
			// Extract all voxels from hash table using OpenMP
			int VoxelOpenMPManager::extractAll(Voxel** h_output_voxels) {
					auto start_total = omp_get_wtime();
					
					// Thread-local storage for voxel collection
					std::vector<std::vector<Voxel>> thread_voxels(omp_get_max_threads());
					
					// Parallel extraction - each thread collects its voxels
					#pragma omp parallel
					{
						int thread_id = omp_get_thread_num();
						
						#pragma omp for schedule(static)
						for (unsigned int idx = 0; idx < table_size; idx++) {
							if (hash_table[idx].occupied == 1) {
								thread_voxels[thread_id].push_back(Voxel(
									hash_table[idx].i,
									hash_table[idx].j,
									hash_table[idx].k,
									hash_table[idx].value
								));
							}
						}
					}
					
					// Combine thread-local results
					int total_count = 0;
					for (const auto& vec : thread_voxels) {
						total_count += vec.size();
					}
					
					*h_output_voxels = new Voxel[total_count];
					int offset = 0;
					for (const auto& vec : thread_voxels) {
						std::copy(vec.begin(), vec.end(), *h_output_voxels + offset);
						offset += vec.size();
					}
					
					auto end_total = omp_get_wtime();
					auto total_duration = (end_total - start_total) * 1000.0;
					printf("[extractAll] Total: %.3f ms\n", total_duration);
					
					return total_count;
			}
			
			// Clear hash table
			void VoxelOpenMPManager::clear() {
				#pragma omp parallel for
				for (unsigned int i = 0; i < table_size; i++) {
					hash_table[i].i = EMPTY_KEY;
					hash_table[i].j = EMPTY_KEY;
					hash_table[i].k = EMPTY_KEY;
					hash_table[i].value = 0.0f;
					hash_table[i].occupied = 0;
				}
				insert_count = 0;
			}

			// Serialization: write current voxel data to binary buffer
			void VoxelOpenMPManager::serialize(uint8_t* bin_data) {
				// int insert_count;
				// unsigned int table_size;
				// VoxelHashEntry* hash_table;

				//size_t expected_size = sizeof(int) + sizeof(unsigned int) + sizeof(VoxelHashEntry) * table_size;
				uint8_t* ptr = bin_data;
				
				memcpy(ptr, &insert_count, sizeof(int));
				ptr += sizeof(int);
				memcpy(ptr, &table_size, sizeof(unsigned int));
				ptr += sizeof(unsigned int);
				memcpy(ptr, hash_table, sizeof(VoxelHashEntry) * table_size);
			}
			
			// Deserialization: read voxel data from binary buffer
			void VoxelOpenMPManager::deserialize(uint8_t* bin_data) {
				// int insert_count;
				// unsigned int table_size;
				// VoxelHashEntry* hash_table;

				// size_t expected_size = sizeof(int) + sizeof(unsigned int) + sizeof(VoxelHashEntry) * table_size;
				// if (bin_data.size() != expected_size) {
				// 	printf("[VoxelOpenMPManager::deserialize] Error: Expected size %zu, got %zu\n", expected_size, bin_data.size());
				// 	return; // Invalid data
				// }
				
				const uint8_t* ptr = bin_data;
				memcpy(&insert_count, ptr, sizeof(int));
				ptr += sizeof(int);
				memcpy(&table_size, ptr, sizeof(unsigned int));
				ptr += sizeof(unsigned int);

				if (hash_table) {
					delete[] hash_table;
				}
				hash_table = new VoxelHashEntry[table_size];
				memcpy(hash_table, ptr, sizeof(VoxelHashEntry) * table_size);
			}
			
			// Merge: combine voxels from another manager (accumulate values)
			void VoxelOpenMPManager::merge(common::vdb::VoxelSparseManager* _other) {

				// Attempt to dynamic_cast to VoxelOpenMPManager
				VoxelOpenMPManager* other = dynamic_cast<VoxelOpenMPManager*>(_other);
				if (!other) {
					printf("[VoxelOpenMPManager::merge] Error: Incompatible manager type\n");
					return;
				}

				// Parallel insertion directly from other's hash table using thread-local maps
				const int thread_count = omp_get_max_threads();
				std::vector<std::unordered_map<uint64_t, float>> thread_maps(thread_count);
				
				#pragma omp parallel
				{
					int tid = omp_get_thread_num();
					auto& local_map = thread_maps[tid];
					
					#pragma omp for schedule(static)
					for (unsigned int i = 0; i < other->table_size; i++) {
						if (other->hash_table[i].occupied == 1) {
							uint64_t key = packCoord3(other->hash_table[i].i, other->hash_table[i].j, other->hash_table[i].k);
							local_map[key] += other->hash_table[i].value;
						}
					}
				}
				
				// Merge thread-local maps
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
				
				// Insert into this manager's hash table
				for (const auto& kv : merged) {
					insertOrUpdatePackedSequential(kv.first, kv.second);
				}
			}

			// ---------------------------------------------
			// VoxelNanoVDBManager implementations
			// ---------------------------------------------
			
			void VoxelNanoVDBManager::serialize(uint8_t* bin_data) {
				// Serialize the NanoVDB grid handle
				const uint8_t* grid_data = (const uint8_t*)gridHandle.data();
				size_t grid_size = gridHandle.size();
				
				// Copy size first, then grid data
				memcpy(bin_data, &grid_size, sizeof(size_t));
				memcpy(bin_data + sizeof(size_t), grid_data, grid_size);
			}
			
			void VoxelNanoVDBManager::deserialize(uint8_t* bin_data) {
				// Deserialize the NanoVDB grid handle
				size_t grid_size;
				memcpy(&grid_size, bin_data, sizeof(size_t));
				
				// Create a new buffer with allocated memory
				auto buffer = nanovdb::HostBuffer::createFull(grid_size, nullptr);
				
				// Copy the serialized data into the buffer
				memcpy(buffer.data(), bin_data + sizeof(size_t), grid_size);
				
				// Move the buffer into the GridHandle
				gridHandle = nanovdb::GridHandle<>(std::move(buffer));
			}
			
			void VoxelNanoVDBManager::merge(common::vdb::VoxelSparseManager* _other) {
				VoxelNanoVDBManager* other = dynamic_cast<VoxelNanoVDBManager*>(_other);
				if (!other) {
					printf("[VoxelNanoVDBManager::merge] Error: Incompatible manager type\n");
					return;
				}
				
				// Extract voxels from other grid and add them to this grid
				Voxel* other_voxels = nullptr;
				int other_count = other->extractAll(&other_voxels);
				
				if (other_count > 0) {
					buildFromVoxels(other_voxels, other_count);
					delete[] other_voxels;
				}
			}
			
			void VoxelNanoVDBManager::buildFromVoxels(Voxel* h_voxels, int num_voxels) {
				if (num_voxels <= 0) return;
				
				/**
				 * @brief Create temporary source grid from input voxels
				 * 
				 * This temporary grid will be merged into the persistent destination grid
				 */
#if OPENVDB_VERSION == 11
				std::shared_ptr<nanovdb::build::FloatGrid> grid_src = std::make_shared<nanovdb::build::FloatGrid>(0.0f, "temp", nanovdb::GridClass::FogVolume);
#else
				std::shared_ptr<nanovdb::tools::build::FloatGrid> grid_src = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "temp", nanovdb::GridClass::FogVolume);
#endif
				
				/**
				 * @brief Build source grid from input voxels
				 * 
				 * First pass: accumulate all input voxels into a temporary source grid
				 */
				std::unordered_map<uint64_t, float> merged;
				merged.reserve((size_t)num_voxels + (size_t)num_voxels / 4 + 1);
				for (int idx = 0; idx < num_voxels; idx++) {
					const Voxel& v = h_voxels[idx];
					merged[packCoord3(v.i, v.j, v.k)] += v.value;
				}

				auto acc_src = grid_src->getAccessor();
				for (const auto& kv : merged) {
					int x, y, z;
					unpackCoord3(kv.first, x, y, z);
					acc_src.setValue(nanovdb::Coord(x, y, z), kv.second);
				}
				
				// Build the final grid handle from source grid
#if OPENVDB_VERSION == 11
				gridHandle = nanovdb::createNanoGrid(*grid_src);
#else
				gridHandle = nanovdb::tools::createNanoGrid(*grid_src);
#endif
			}
			
			int VoxelNanoVDBManager::extractAll(Voxel** h_output_voxels) {
				std::vector<Voxel> voxels;
				
				auto* grid_src_float = (nanovdb::NanoGrid<float>*)getGridHandle().data();
				
				// Iterate over tree hierarchy: root -> upper internal -> lower internal -> leaf nodes
				for (auto it2 = grid_src_float->tree().root().cbeginChild(); it2; ++it2) {
					// Iterate over upper internal nodes (level 2)
					for (auto it1 = it2->cbeginChild(); it1; ++it1) {
						// Iterate over lower internal nodes (level 1)
						for (auto it0 = it1->cbeginChild(); it0; ++it0) {
							// Iterate over active voxels in leaf nodes
							for (auto it = it0->cbeginValueOn(); it; ++it) {
								float v = *it;
								nanovdb::Coord xyz = it.getCoord();
								voxels.push_back(Voxel(xyz[0], xyz[1], xyz[2], v));
							}
						}
					}
				}
				
				int count = voxels.size();
				*h_output_voxels = new Voxel[count];
				std::copy(voxels.begin(), voxels.end(), *h_output_voxels);
				
				return count;
			}

#ifdef WITH_OPENVDB
		// ---------------------------------------------
		// VoxelOpenVDBManager implementations
		// ---------------------------------------------
		
		void VoxelOpenVDBManager::serialize(uint8_t* bin_data) {
				// For OpenVDB serialization, we need to extract voxels and serialize them
				// This is a simplified approach - proper OpenVDB serialization would use iostreams
				Voxel* voxels = nullptr;
				int count = extractAll(&voxels);
				
				uint8_t* ptr = bin_data;
				memcpy(ptr, &count, sizeof(int));
				ptr += sizeof(int);
				memcpy(ptr, voxels, count * sizeof(Voxel));
				
				delete[] voxels;
			}
			
			void VoxelOpenVDBManager::deserialize(uint8_t* bin_data) {
				// Deserialize voxels and rebuild the grid
				const uint8_t* ptr = bin_data;
				int count;
				memcpy(&count, ptr, sizeof(int));
				ptr += sizeof(int);
				
				Voxel* voxels = new Voxel[count];
				memcpy(voxels, ptr, count * sizeof(Voxel));
				
				buildFromVoxels(voxels, count);
				delete[] voxels;
			}
			
			void VoxelOpenVDBManager::merge(common::vdb::VoxelSparseManager* _other) {
				VoxelOpenVDBManager* other = dynamic_cast<VoxelOpenVDBManager*>(_other);
				if (!other) {
					printf("[VoxelOpenVDBManager::merge] Error: Incompatible manager type\n");
					return;
				}
				
				// Extract voxels from other grid and add them to this grid
				Voxel* other_voxels = nullptr;
				int other_count = other->extractAll(&other_voxels);
				
				if (other_count > 0) {
					buildFromVoxels(other_voxels, other_count);
					delete[] other_voxels;
				}
			}
			
			void VoxelOpenVDBManager::buildFromVoxels(Voxel* h_voxels, int num_voxels) {
				if (num_voxels <= 0) return;
				
				// Pre-aggregate voxels
				std::unordered_map<uint64_t, float> merged;
				merged.reserve((size_t)num_voxels + (size_t)num_voxels / 4 + 1);
				for (int idx = 0; idx < num_voxels; idx++) {
					const Voxel& v = h_voxels[idx];
					merged[packCoord3(v.i, v.j, v.k)] += v.value;
				}

				// Insert into OpenVDB grid
				openvdb::FloatGrid::Accessor acc = persistentGrid->getAccessor();
				for (const auto& kv : merged) {
					int x, y, z;
					unpackCoord3(kv.first, x, y, z);
					openvdb::Coord xyz(x, y, z);
					
					// Accumulate values if voxel exists, otherwise set new value
					if (acc.isValueOn(xyz)) {
						float existing_value = acc.getValue(xyz);
						acc.setValue(xyz, existing_value + kv.second);
					} else {
						acc.setValue(xyz, kv.second);
					}
				}
			}
			
			int VoxelOpenVDBManager::extractAll(Voxel** h_output_voxels) {
				std::vector<Voxel> voxels;
				
				// Iterate over all active voxels in the grid
				for (openvdb::FloatGrid::ValueOnCIter iter = persistentGrid->cbeginValueOn(); iter; ++iter) {
					openvdb::Coord xyz = iter.getCoord();
					float value = *iter;
					voxels.push_back(Voxel(xyz[0], xyz[1], xyz[2], value));
				}
				
				int count = voxels.size();
				*h_output_voxels = new Voxel[count];
				std::copy(voxels.begin(), voxels.end(), *h_output_voxels);
				
				return count;
			}
#endif // WITH_OPENVDB
	}
}
}