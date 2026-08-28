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

#include "convert_vdb.h"
#include "convert_vdb_kernel.h"

#include <iostream>
#include <string>
#include <vector>
#include <float.h>
#include <iomanip>

#ifdef WITH_OPENMP
# include <omp.h>
#endif

#ifdef WITH_OPENVDB
#	include <openvdb/openvdb.h>
#	include <openvdb/points/PointConversion.h>
#	include <openvdb/points/PointCount.h>
#	include <openvdb/io/Stream.h>

#	include <openvdb/points/PointDataGrid.h>
#	include <openvdb/tools/ParticlesToLevelSet.h>
#	include <openvdb/tools/PointScatter.h>
#	include <openvdb/tools/Dense.h>

#	include <nanovdb/tools/NanoToOpenVDB.h>
#endif

#include "dense_common.h"

#include <mpi.h>
#include <algorithm>

#ifdef WITH_CUDAKDTREE
#	include "cudakdtree_tool.h"
#endif

#ifdef WITH_NANOFLANN
#	include "nanoflann_tool.h"
#	include <nanoflann.hpp>
#endif

#include "utility/dense_utility.h"

#ifdef WITH_GPU_CUDA
#	include "utility/gpu_utility.h"
#endif

namespace space_converter {
	namespace common {
		namespace vdb {

			/**
			 * @brief Read precomputed particle radii and densities from binary files
			 *
			 * Loads particle neighbor search results from disk for each particle type.
			 * This avoids recomputing expensive neighbor searches.
			 */
			void ConvertVDBBase::read_radius_from_file(std::string& calc_radius_neigh_file)
			{
				int ptype_count = get_num_types(); // Ensure that radius_particles is resized correctly based on the number of particle types
				cache_manager.particles_ptype_offset.resize(ptype_count + 1, 0);

				for (int ptype = 0; ptype < ptype_count; ptype++) {
					std::string calc_radius_neigh_file_ptype = calc_radius_neigh_file + "." + std::to_string(ptype) + ".bin";

					// Check if file exists
					std::ifstream file(calc_radius_neigh_file_ptype, std::ios::binary);
					if (!file.good()) {
						printf("File %s does not exist or cannot be opened\n", calc_radius_neigh_file_ptype.c_str());
						return;
					}

					// Read the size of the vector
					size_t size;
					if (!file.read(reinterpret_cast<char*>(&size), sizeof(size_t))) {
						printf("Error reading size from file %s\n", calc_radius_neigh_file_ptype.c_str());
						return;
					}

					// Resize the vector and read the data
					std::vector<float> radius_particles(size);
					if (!file.read(reinterpret_cast<char*>(radius_particles.data()), size * sizeof(float))) {
						printf("Error reading data from file %s\n", calc_radius_neigh_file_ptype.c_str());
						return;
					}

					printf("Read %zu radius particles for particle type %d from file %s\n", size, ptype, calc_radius_neigh_file_ptype.c_str());

					// Resize the vector and read the data
					std::vector<float> rho_particles(size);
					if (!file.read(reinterpret_cast<char*>(rho_particles.data()), size * sizeof(float))) {
						printf("Error reading data from file %s\n", calc_radius_neigh_file_ptype.c_str());
						return;
					}

					printf("Read %zu rho particles for particle type %d from file %s\n", size, ptype, calc_radius_neigh_file_ptype.c_str());

					// Overwrite the per-type slots that find_particle_positions() already
					// created (an earlier version push_back'ed here, appending the loaded
					// arrays PAST index ptype_count so the kernels kept using the stale
					// hsml-based radii).
					if ((int)cache_manager.radius_particles_per_ptype.size() <= ptype)
						cache_manager.radius_particles_per_ptype.resize(ptype + 1);
					if ((int)cache_manager.rho_particles_per_ptype.size() <= ptype)
						cache_manager.rho_particles_per_ptype.resize(ptype + 1);

					cache_manager.radius_particles_per_ptype[ptype] = std::move(radius_particles);
					cache_manager.rho_particles_per_ptype[ptype] = std::move(rho_particles);

					cache_manager.particles_ptype_offset[ptype + 1] = cache_manager.particles_ptype_offset[ptype] + cache_manager.radius_particles_per_ptype[ptype].size();

					file.close();
				}
			}

			/**
			 * @brief Write precomputed particle radii and densities to binary files
			 *
			 * Saves particle neighbor search results to disk for later reuse.
			 */
			void ConvertVDBBase::write_radius_from_file(std::string& calc_radius_neigh_file)
			{
				int ptype_count = get_num_types(); // Ensure that radius_particles is resized correctly based on the number of particle types

				for (int ptype = 0; ptype < ptype_count; ptype++) {
					std::string calc_radius_neigh_file_ptype = calc_radius_neigh_file + "." + std::to_string(ptype) + ".bin";

					std::ofstream file(calc_radius_neigh_file_ptype, std::ios::binary);
					if (!file.is_open()) {
						printf("Cannot create file %s for writing\n", calc_radius_neigh_file_ptype.c_str());
						return;
					}

					// Write the size of the vector
					size_t size = cache_manager.radius_particles_per_ptype[ptype].size();
					if (!file.write(reinterpret_cast<const char*>(&size), sizeof(size_t))) {
						printf("Error writing size to file %s\n", calc_radius_neigh_file_ptype.c_str());
						file.close();
						return;
					}

					// Write the vector data
					if (!file.write(reinterpret_cast<const char*>(cache_manager.radius_particles_per_ptype[ptype].data()), size * sizeof(float))) {
						printf("Error writing data to file %s\n", calc_radius_neigh_file_ptype.c_str());
						file.close();
						return;
					}

					// Write the vector data
					if (!file.write(reinterpret_cast<const char*>(cache_manager.rho_particles_per_ptype[ptype].data()), size * sizeof(float))) {
						printf("Error writing data to file %s\n", calc_radius_neigh_file_ptype.c_str());
						file.close();
						return;
					}

					file.close();
				}
			}

			void ConvertVDBBase::find_particle_positions()
			{
				// double t_start = omp_get_wtime();

				size_t no_points = get_local_num_particles();
				int ptype_count = get_num_types();

				cache_manager.pos_particles_per_ptype.resize(ptype_count);
				cache_manager.radius_particles_per_ptype.resize(ptype_count);
				cache_manager.rho_particles_per_ptype.resize(ptype_count);
				cache_manager.mass_particles_per_ptype.resize(ptype_count);
				cache_manager.particles_id_ordered_per_ptype.resize(ptype_count);
				cache_manager.particles_reader_id_per_ptype.resize(ptype_count);
				cache_manager.particles_ptype_offset.resize(ptype_count + 1, 0);

				int num_threads = 1;
#ifdef WITH_OPENMP
				num_threads = omp_get_max_threads();
#endif

				cache_manager.num_particles_per_ptype.assign(ptype_count, 0);

				if (cache_manager.skip_cache) {
					// Streaming mode: only the per-type counts are needed here, the particle
					// data itself is pulled from the reader again during conversion. This is
					// what keeps the ~32 bytes per particle of cache from being allocated.
					for (size_t i = 0; i < no_points; i++) {
						int ptype = get_particle_type(i);
						if (ptype >= 0 && ptype < ptype_count)
							cache_manager.num_particles_per_ptype[ptype]++;
					}

					for (int ptype = 0; ptype < ptype_count; ptype++)
						cache_manager.particles_ptype_offset[ptype + 1] =
							cache_manager.particles_ptype_offset[ptype] + cache_manager.num_particles_per_ptype[ptype];

					return;
				}

				for (int ptype = 0; ptype < ptype_count; ptype++) {

					std::vector<float> points_master;
					std::vector < std::vector<float> > points_thread(num_threads);

					std::vector<float> pmass_master;
					std::vector < std::vector<float> > pmass_thread(num_threads);

					std::vector<float> rho_master;
					std::vector < std::vector<float> > rho_thread(num_threads);

					std::vector<float> radius_master;
					std::vector < std::vector<float> > radius_thread(num_threads);

					std::vector<size_t> particles_id_master;
					std::vector < std::vector<size_t> > particles_id_thread(num_threads);

#pragma omp parallel num_threads(num_threads)
					{
#ifdef WITH_OPENMP
						int tid = omp_get_thread_num();
#else
						int tid = 0;
#endif

#pragma omp for
						for (size_t i = 0; i < no_points; i++) {

							if (get_particle_type(i) != ptype)
								continue;

							double Pos[3];
							get_particle_position(i, Pos);

							// Collect particle positions per thread
							points_thread[tid].push_back(Pos[0]);
							points_thread[tid].push_back(Pos[1]);
							points_thread[tid].push_back(Pos[2]);

							double mass = get_particle_mass(i);
							pmass_thread[tid].push_back(mass);

							double rho = get_particle_rho(i);
							rho_thread[tid].push_back(rho);

							double radius = get_particle_hsml(i);
							radius_thread[tid].push_back(radius);

							// Global reader id of this compact slot (needed later to
							// query per-particle values from the reader)
							particles_id_thread[tid].push_back(i);
						}
					}

					for (int t = 0; t < num_threads; ++t) {
						points_master.insert(points_master.end(), points_thread[t].begin(), points_thread[t].end());
						pmass_master.insert(pmass_master.end(), pmass_thread[t].begin(), pmass_thread[t].end());
						rho_master.insert(rho_master.end(), rho_thread[t].begin(), rho_thread[t].end());
						radius_master.insert(radius_master.end(), radius_thread[t].begin(), radius_thread[t].end());
						particles_id_master.insert(particles_id_master.end(), particles_id_thread[t].begin(), particles_id_thread[t].end());
					}

					cache_manager.pos_particles_per_ptype[ptype] = std::move(points_master);
					cache_manager.radius_particles_per_ptype[ptype] = std::move(radius_master);
					cache_manager.mass_particles_per_ptype[ptype] = std::move(pmass_master);
					cache_manager.rho_particles_per_ptype[ptype] = std::move(rho_master);
					// Reader ids: compact slot k of this type <-> global reader id.
					cache_manager.particles_reader_id_per_ptype[ptype] = std::move(particles_id_master);
					// Iteration order: COMPACT indices 0..N_t-1 (later permuted by the
					// sorts). These index the compact per-type arrays above; they must
					// NEVER be passed to reader accessors (use the reader ids for that).
					// A previous version stored the global reader ids here, which made
					// the conversion kernels index the compact arrays out of bounds for
					// every particle type > 0.
					{
						size_t n_t = cache_manager.pos_particles_per_ptype[ptype].size() / 3;
						std::vector<size_t> iter_order(n_t);
						for (size_t k = 0; k < n_t; k++)
							iter_order[k] = k;
						cache_manager.particles_id_ordered_per_ptype[ptype] = std::move(iter_order);
					}
					cache_manager.particles_ptype_offset[ptype + 1] = cache_manager.particles_ptype_offset[ptype] + cache_manager.pos_particles_per_ptype[ptype].size() / 3;
					cache_manager.num_particles_per_ptype[ptype] = cache_manager.pos_particles_per_ptype[ptype].size() / 3;
				}

				// if (cache_manager.world_rank == 0) {
				// 	printf("find_particle_positions: Find positions: %f\n", omp_get_wtime() - t_start);
				// }
			}

			bool ConvertVDBBase::find_particle_values(int ptype, int block_name_id)
			{
				// Streaming mode reads the values per chunk during conversion
				if (cache_manager.skip_cache)
					return false;

				//int cached_value_ptype = -1; ///< Cached particle type for values_particles (to track which type's values are currently cached)
				//int cached_value_block_name_id = -1; ///< Cached block name ID for values_particles (to track which block's values are currently cached)
				if (cache_manager.cached_value_ptype == ptype && cache_manager.cached_value_block_name_id == block_name_id) {
					// Values are already cached for this particle type and block, no need to recompute
					return false;
				}

				// double t_start = omp_get_wtime();

				// Values are cached in COMPACT order (slot k of values_particles is
				// the value of the k-th cached particle of this type). The reader is
				// queried with the GLOBAL reader id of that slot — passing the compact
				// index here (as an earlier version did) fetched the values of the
				// wrong particles for every type > 0.
				const auto& reader_ids = cache_manager.particles_reader_id_per_ptype[ptype];
				size_t num_particles = reader_ids.size();

				if (num_particles == 0) {
					cache_manager.values_particles.clear();
					cache_manager.cached_value_ptype = ptype;
					cache_manager.cached_value_block_name_id = block_name_id;
					return true;
				}

				int num_threads = 1;
#ifdef WITH_OPENMP
				num_threads = omp_get_max_threads();
#endif

				// Pre-allocate values_master with known size
				std::vector<float> values_master(num_particles);

#pragma omp parallel for num_threads(num_threads)
				for (size_t idx = 0; idx < num_particles; idx++) {
					float value = get_particle_norm_value(block_name_id, reader_ids[idx]);
					values_master[idx] = value;
				}

				cache_manager.values_particles = std::move(values_master);
				cache_manager.cached_value_ptype = ptype;
				cache_manager.cached_value_block_name_id = block_name_id;

				return true;
			}

#ifdef WITH_CUDAKDTREE
			/**
			 * @brief Calculate particle radii using CUDA-accelerated k-NN search
			 *
			 * Uses GPU-accelerated KD-tree to find k-nearest neighbors for each particle.
			 * Computes smoothing radius and density based on neighbor distances.
			 */
			void ConvertVDBBase::calculate_radius_by_cudakdtree(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, bool use_cudakdtree_cpu, float maxRadius, common::SpaceData::DenseType& rho_kernel)
			{
				// double t_start = omp_get_wtime();

				if (calc_radius_neigh == -1) {
					read_radius_from_file(calc_radius_neigh_file);
				}
				else {
					size_t no_points = get_local_num_particles();

					int ptype_count = get_num_types();

					int num_threads = omp_get_max_threads();

					for (int ptype = 0; ptype < ptype_count; ptype++) {


											// Run GPU/CPU-based k-Nearest Neighbor search
						utility::cudakdtree::run_knn(cache_manager.pos_particles_per_ptype[ptype].data(), cache_manager.pos_particles_per_ptype[ptype].size() / 3, calc_radius_neigh + 1, cache_manager.radius_particles_per_ptype[ptype], cache_manager.rho_particles_per_ptype[ptype], cache_manager.mass_particles_per_ptype[ptype], !use_cudakdtree_cpu, use_cycling, maxRadius, rho_kernel);
					}

					if (calc_radius_neigh_file.length() > 0) {
						write_radius_from_file(calc_radius_neigh_file);
					}
				}

				// printf("cudakdtree: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
			}

			/**
			 * @brief Build per-type float4 KD-trees and cache them in cache_manager.
			 *
			 * Called once after find_particle_positions() in init_converter().
			 * Tree nodes use offset-adjusted world positions (x,y,z) with the original
			 * particle index bit-cast into w, so they can be used directly in
			 * convert_to_dense_grid_loop_over_voxels_cpu/gpu without rebuilding.
			 */
			void ConvertVDBBase::build_voxel_kdtree(float* offset_position)
			{
				const int ptype_count = (int)cache_manager.pos_particles_per_ptype.size();
				cache_manager.voxel_kdtree_per_ptype.resize(ptype_count);

#ifdef WITH_GPU_CUDA
				if (cache_manager.use_gpu) {
					// GPU path: build directly on GPU, store in d_voxel_kdtree_per_ptype
					cache_manager.d_voxel_kdtree_per_ptype.resize(ptype_count, nullptr);

					for (int ptype = 0; ptype < ptype_count; ++ptype) {
						const auto& pos = cache_manager.pos_particles_per_ptype[ptype];
						//const auto& ids = cache_manager.particles_id_ordered_per_ptype[ptype];
						size_t N = pos.size() / 3;

						if (N == 0) {
							cache_manager.d_voxel_kdtree_per_ptype[ptype] = nullptr;
							continue;
						}

						printf("[build_voxel_kdtree] GPU: ptype=%d  N=%zu  building float4 KD-tree...\n", ptype, N);
						utility::cudakdtree::build_float4_kdtree(
							pos.data(),
							//ids.data(),
							N,
							offset_position,
							true,  // use_gpu
							cache_manager.voxel_kdtree_per_ptype[ptype],  // not used for GPU
							&cache_manager.d_voxel_kdtree_per_ptype[ptype]
						);
						printf("[build_voxel_kdtree] GPU: ptype=%d  tree built on device (%zu nodes)\n", ptype, N);
					}
				}
				else
#endif
				{
					// CPU path: build on CPU, store in voxel_kdtree_per_ptype
					for (int ptype = 0; ptype < ptype_count; ++ptype) {
						const auto& pos = cache_manager.pos_particles_per_ptype[ptype];
						//const auto& ids = cache_manager.particles_id_ordered_per_ptype[ptype];
						size_t N = pos.size() / 3;

						if (N == 0) {
							cache_manager.voxel_kdtree_per_ptype[ptype].clear();
							continue;
						}

						printf("[build_voxel_kdtree] CPU: ptype=%d  N=%zu  building float4 KD-tree...\n", ptype, N);
						utility::cudakdtree::build_float4_kdtree(
							pos.data(),
							//ids.data(),
							N,
							offset_position,
							false,  // use_gpu
							cache_manager.voxel_kdtree_per_ptype[ptype],
							nullptr  // d_out_tree not used for CPU
						);
						printf("[build_voxel_kdtree] CPU: ptype=%d  tree built (%zu nodes)\n",
							ptype, cache_manager.voxel_kdtree_per_ptype[ptype].size());
					}
				}
			}

#endif // WITH_CUDAKDTREE


#ifdef WITH_NANOFLANN
			/**
			 * @brief Calculate particle radii using CPU-based nanoflann k-NN search
			 *
			 * Uses CPU-based KD-tree (nanoflann library) to find k-nearest neighbors.
			 * Computes smoothing radius and density based on neighbor distances.
			 */
			void ConvertVDBBase::calculate_radius_by_nanoflann(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, common::SpaceData::DenseType& rho_kernel)
			{
				// double t_start = omp_get_wtime();

				if (calc_radius_neigh == -1) {
					read_radius_from_file(calc_radius_neigh_file);
				}
				else {
					size_t no_points = get_local_num_particles();

					int ptype_count = get_num_types();

#ifdef WITH_OPENMP
					int num_threads = omp_get_max_threads();
#else
					int num_threads = 1;
#endif

					for (int ptype = 0; ptype < ptype_count; ptype++) {


						if (cache_manager.pos_particles_per_ptype[ptype].size() > 0) {
							// Run k-NN search (calc_radius_neigh + 1 includes the query point itself)
							utility::nanoflann_tool::run_knn(cache_manager.pos_particles_per_ptype[ptype].data(), cache_manager.pos_particles_per_ptype[ptype].size() / 3, calc_radius_neigh + 1, cache_manager.radius_particles_per_ptype[ptype], cache_manager.rho_particles_per_ptype[ptype], cache_manager.mass_particles_per_ptype[ptype], use_cycling, rho_kernel);
						}
					}

					if (calc_radius_neigh_file.length() > 0) {
						write_radius_from_file(calc_radius_neigh_file);
					}
				}

				// printf("nanoflann: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
			}
#endif

			/**
			 * @brief Main conversion function: converts particle data to VDB grids
			 *
			 * This is the core function that:
			 * 1. Iterates through cached particles of a specific type
			 * 2. Filters by type, bounding box, and value ranges
			 * 3. Rasterizes particles into voxel grid (dense or sparse)
			 * 4. Applies smoothing kernels for SPH-like density estimation
			 * 5. Tracks min/max values for normalization
			 */
			void ConvertVDBBase::convert_iolib_to_grid(
				int particle_type,
				float particle_radius_multiplier,
				float* bbox_min,
				float* bbox_max,
				int bbox_dim,
				int* bbox_min_orig,
				double bbox_size_orig,
				common::SpaceData::DenseType dense_type,
				common::SpaceData::DenseNorm dense_norm,
				int block_name_id,
				float object_size,
				float& min_value,
				float& max_value,
				size_t& particles_count,
				VDBParticles& grid,
				double& transform_scale,
				float filter_min,
				float filter_max,
				common::SpaceData::AnimType anim_type,
				int frame_req,
				int frame,
				float* bbox_sphere_pos,
				float bbox_sphere_r,
				bool use_simple_density,
				bool use_norm_value,
				float* offset_position
			)
			{
				//printf("===convert_iolib_to_grid: START\n");
#ifdef WITH_OPENMP
				double t = omp_get_wtime();
				// 			double t_step = omp_get_wtime();
#endif   

							// Calculate spatial scaling factors
							// Compare original cubic bounding box diagonal to actual data bounding box diagonal
				double orig_box_space_diagonal = object_size * sqrt(3.0);  // Diagonal of cube with side length object_size
				double new_box_space_diagonal = sqrt(pow(bbox_max[0] - bbox_min[0], 2) + pow(bbox_max[1] - bbox_min[1], 2) + pow(bbox_max[2] - bbox_min[2], 2));

				// Scale factor to fit data into desired coordinate system
				double scale_space_diagonal = new_box_space_diagonal / orig_box_space_diagonal;

				// Voxel size in world space units
				transform_scale = scale_space_diagonal * object_size / (double)bbox_dim;

				size_t no_points = get_local_num_particles();

				// Prepare raw particle data containers if using raw particle type
				// Raw particles store individual particle data without rasterization
				RawParticles::ParticleData raw_positions;
				RawParticles::ParticleData raw_values;
				RawParticles::ParticleData raw_radius;
				RawParticles::ParticleData raw_frame;

				// Set grid transform (voxel size in world space units)
				if (grid.type == VDBParticleType::eNanoVDB || grid.type == VDBParticleType::eOpenVDB || grid.type == VDBParticleType::eCUB) {
					//grid.nano_grid->setTransform(transform_scale);
					grid.sparse_grid->set_transform(transform_scale);
				}
				// #ifdef WITH_OPENVDB
				// 			else if (grid.type == VDBParticleType::eOpenVDB) {
				// 				// Create linear transform for OpenVDB grid
				// 				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
				// 				grid.vdb_grid->setTransform(transform);
				// 			}
				// #endif
				else if (grid.type == VDBParticleType::eRawParticles) {
					// Reserve space for raw particle arrays
					raw_positions.values.reserve(no_points);
					raw_positions.num_comp = 3;
					raw_positions.name = "position";

					raw_values.values.reserve(no_points);
					raw_values.num_comp = 1;
					raw_values.name = get_dataset_name(block_name_id);

					raw_radius.values.reserve(no_points);
					raw_radius.num_comp = 1;
					raw_radius.name = "radius";

					raw_frame.values.reserve(no_points);
					raw_frame.num_comp = 1;
					raw_frame.name = "frame";
				}

				// Normalize bounding box coordinates
				float bbox_x_min_norm = bbox_min[0] / (float)object_size;
				float bbox_x_max_norm = bbox_max[0] / (float)object_size;
				float bbox_y_min_norm = bbox_min[1] / (float)object_size;
				float bbox_y_max_norm = bbox_max[1] / (float)object_size;
				float bbox_z_min_norm = bbox_min[2] / (float)object_size;
				float bbox_z_max_norm = bbox_max[2] / (float)object_size;

				if (grid.type == VDBParticleType::eDense) {
					// Set grid offset for dense grids (position in global coordinate system)
					// This allows tiled dense grids to know their position in the full domain
					grid.dense_grid->offset[0] = (size_t)(bbox_x_min_norm * (double)bbox_dim / scale_space_diagonal);
					grid.dense_grid->offset[1] = (size_t)(bbox_y_min_norm * (double)bbox_dim / scale_space_diagonal);
					grid.dense_grid->offset[2] = (size_t)(bbox_z_min_norm * (double)bbox_dim / scale_space_diagonal);
				}

				// Get cached particle data for this particle type
				if (particle_type < 0 || particle_type >= (int)cache_manager.pos_particles_per_ptype.size()) {
					printf("Warning: Invalid particle type %d in convert_iolib_to_grid\n", particle_type);
					particles_count = 0;
					return;
				}

				size_t num_particles = cache_manager.skip_cache
					? cache_manager.num_particles_per_ptype[particle_type]
					: cache_manager.pos_particles_per_ptype[particle_type].size() / 3;
				if (num_particles == 0) {
					particles_count = 0;
					return;
				}

				//printf("===convert_iolib_to_grid: INIT: %f [s]\n", omp_get_wtime() - t);

				// Cache value particles if not already cached for this particle type and block
				bool new_values = find_particle_values(particle_type, block_name_id);

				//printf("===convert_iolib_to_grid: FIND: %f [s]\n", omp_get_wtime() - t);

#ifdef WITH_GPU_CUDA
				if (cache_manager.use_gpu) {

					if (new_values) {
						// Upload values_particles to GPU if they were newly computed
						cache_manager.copy_values_to_gpu();
					}

					//printf("===convert_iolib_to_grid: COPY TO GPU: %f [s]\n", omp_get_wtime() - t);

					// Use GPU kernel with cached GPU particle positions
					const float* d_pos_particles = cache_manager.d_pos_particles_per_ptype[particle_type];
					const size_t* d_particles_id_ordered = cache_manager.d_particles_id_ordered_per_ptype[particle_type];
					const float* d_radius_particles = cache_manager.d_radius_particles_per_ptype[particle_type];
					const float* d_value_particles = cache_manager.d_values_particles;

					kernel::convert_iolib_to_grid_gpu(
						d_pos_particles,
						d_particles_id_ordered,
						d_radius_particles,
						d_value_particles,
						num_particles,
						particle_radius_multiplier,
						bbox_dim,
						bbox_min_orig,
						bbox_size_orig,
						dense_type,
						dense_norm,
						block_name_id,
						object_size,
						min_value,
						max_value,
						particles_count,
						grid,
						transform_scale,
						filter_min,
						filter_max,
						anim_type,
						frame_req,
						frame,
						bbox_sphere_pos,
						bbox_sphere_r,
						use_simple_density,
						cache_manager.particle_radius_const,

						bbox_x_min_norm,
						bbox_x_max_norm,
						bbox_y_min_norm,
						bbox_y_max_norm,
						bbox_z_min_norm,
						bbox_z_max_norm,
						scale_space_diagonal,
						offset_position,
#ifdef WITH_CUDAKDTREE
						cache_manager.use_dense_loop_over_voxels,
						cache_manager.calc_radius_neigh,
						// Pre-built device KD-tree (nullptr if not built or wrong ptype)
						(particle_type < (int)cache_manager.d_voxel_kdtree_per_ptype.size())
						? cache_manager.d_voxel_kdtree_per_ptype[particle_type] : nullptr,
						(particle_type < (int)cache_manager.d_voxel_kdtree_per_ptype.size())
						? (int)cache_manager.particles_id_ordered_per_ptype[particle_type].size() : 0
#else
						false,
						16,
						nullptr,
						0
#endif
					);

					//printf("===convert_iolib_to_grid: CONVERT GPU: %f [s]\n", omp_get_wtime() - t);
				}
				else
#endif
				{
				if (cache_manager.skip_cache) {
					// Streaming mode: pull the particles of this type from the reader in
					// chunks and splat each chunk with the very same kernel the cached path
					// uses, so the result is identical. The grid accumulates across calls,
					// but min/max and the particle count are per-call outputs and are
					// therefore combined here.
					const size_t CHUNK = size_t(1) << 20;
					const size_t no_points = get_local_num_particles();

					std::vector<float>  pos_chunk(3 * CHUNK);
					std::vector<float>  radius_chunk(CHUNK);
					std::vector<float>  value_chunk(CHUNK);
					std::vector<size_t> ids_chunk(CHUNK);
					std::vector<size_t> reader_idx(CHUNK);

					// The kernel indexes the arrays through particle_ids; here they are
					// already laid out in chunk order, so the mapping is the identity.
					for (size_t j = 0; j < CHUNK; j++)
						ids_chunk[j] = j;

					float  min_total = FLT_MAX, max_total = -FLT_MAX;
					size_t count_total = 0;
					size_t n = 0;

					auto flush_chunk = [&]() {
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static)
#endif
						for (long long k = 0; k < (long long)n; k++) {
							const size_t src = reader_idx[k];
							double Pos[3];
							get_particle_position(src, Pos);
							pos_chunk[3 * k + 0] = (float)Pos[0];
							pos_chunk[3 * k + 1] = (float)Pos[1];
							pos_chunk[3 * k + 2] = (float)Pos[2];
							radius_chunk[k] = (float)get_particle_hsml(src);
							// Reader accessors take the global reader id
							value_chunk[k] = get_particle_norm_value(block_name_id, src);
						}

						float  min_c = FLT_MAX, max_c = -FLT_MAX;
						size_t count_c = 0;

						kernel::convert_iolib_to_grid_cpu(
							pos_chunk.data(), ids_chunk.data(), radius_chunk.data(), value_chunk.data(),
							n,
							particle_radius_multiplier,
							bbox_dim, bbox_min_orig, bbox_size_orig,
							dense_type, dense_norm,
							block_name_id, object_size,
							min_c, max_c, count_c,
							grid, transform_scale,
							filter_min, filter_max,
							anim_type, frame_req, frame,
							bbox_sphere_pos, bbox_sphere_r, use_simple_density,
							cache_manager.particle_radius_const,
							bbox_x_min_norm, bbox_x_max_norm,
							bbox_y_min_norm, bbox_y_max_norm,
							bbox_z_min_norm, bbox_z_max_norm,
							scale_space_diagonal, offset_position,
							false, 16, nullptr, 0);

						if (min_c < min_total) min_total = min_c;
						if (max_c > max_total) max_total = max_c;
						count_total += count_c;
						n = 0;
					};

					for (size_t i = 0; i < no_points; i++) {
						if (get_particle_type(i) != particle_type)
							continue;

						reader_idx[n++] = i;
						if (n == CHUNK)
							flush_chunk();
					}
					if (n > 0)
						flush_chunk();

					min_value = min_total;
					max_value = max_total;
					particles_count = count_total;
					return;
				}

					// Use CPU kernel with cached CPU particle positions
					const float* pos_particles = cache_manager.pos_particles_per_ptype[particle_type].data();
					const size_t* particle_ids = cache_manager.particles_id_ordered_per_ptype[particle_type].data();
					const float* radius_particles = cache_manager.radius_particles_per_ptype[particle_type].data();
					const float* value_particles = cache_manager.values_particles.data();
					kernel::convert_iolib_to_grid_cpu(
						pos_particles,
						particle_ids,
						radius_particles,
						value_particles,
						num_particles,
						particle_radius_multiplier,
						bbox_dim,
						bbox_min_orig,
						bbox_size_orig,
						dense_type,
						dense_norm,
						block_name_id,
						object_size,
						min_value,
						max_value,
						particles_count,
						grid,
						transform_scale,
						filter_min,
						filter_max,
						anim_type,
						frame_req,
						frame,
						bbox_sphere_pos,
						bbox_sphere_r,
						use_simple_density,
						cache_manager.particle_radius_const,

						bbox_x_min_norm,
						bbox_x_max_norm,
						bbox_y_min_norm,
						bbox_y_max_norm,
						bbox_z_min_norm,
						bbox_z_max_norm,
						scale_space_diagonal,

						offset_position,
#ifdef WITH_CUDAKDTREE
						cache_manager.use_dense_loop_over_voxels,
						cache_manager.calc_radius_neigh,
						// Pre-built CPU KD-tree (nullptr if not built or wrong ptype)
						(particle_type < (int)cache_manager.voxel_kdtree_per_ptype.size()
							&& !cache_manager.voxel_kdtree_per_ptype[particle_type].empty())
						? cache_manager.voxel_kdtree_per_ptype[particle_type].data() : nullptr,
						(particle_type < (int)cache_manager.voxel_kdtree_per_ptype.size())
						? (int)cache_manager.voxel_kdtree_per_ptype[particle_type].size() : 0
#else
						false,
						16,
						nullptr,
						0
#endif
					);

					//printf("===convert_iolib_to_grid: CONVERT CPU: %f [s]\n", omp_get_wtime() - t);
				}


				// #ifdef WITH_OPENMP
				// 			t_step = omp_get_wtime();
				// #endif
				//printf("===convert_iolib_to_grid: FINISHED: %f [s]\n", omp_get_wtime() - t);
			}


			/**
			 * @brief Merge two VDB grids by summing their values
			 *
			 * Handles merging of different grid types (dense, sparse, serialized).
			 * Used in MPI distributed processing to combine partial results.
			 */
			void ConvertVDBBase::merge_grid(
				VDBParticles& grid_dst,
				VDBParticles& grid_recv
			)
			{
				if (grid_dst.type == VDBParticleType::eDense && grid_recv.type == VDBParticleType::eDense) {
					// Merge dense grids by summing voxel values (element-wise addition)
#pragma omp parallel for
					for (int z = 0; z < grid_dst.dense_grid->dims[2]; z++) {
						for (int y = 0; y < grid_dst.dense_grid->dims[1]; y++) {
							for (int x = 0; x < grid_dst.dense_grid->dims[0]; x++) {
								size_t index = grid_dst.dense_grid->get_index(x, y, z);
								grid_dst.dense_grid->data_density[index] += grid_recv.dense_grid->data_density[index];
								// Merge temp buffer if both grids have it allocated
								if (!grid_dst.dense_grid->data_temp.empty() && !grid_recv.dense_grid->data_temp.empty()) {
									grid_dst.dense_grid->data_temp[index] += grid_recv.dense_grid->data_temp[index];
								}
							}
						}
					}
				}
				else if ((grid_dst.type == VDBParticleType::eNanoVDB || grid_dst.type == VDBParticleType::eOpenVDB) && grid_recv.type == VDBParticleType::eVector) {
					grid_dst.sparse_grid->merge(grid_recv.vector_grid);
				}

				else if (grid_dst.type == VDBParticleType::eRawParticles && grid_recv.type == VDBParticleType::eVector) {
					RawParticles temp;
					temp.deserialize(grid_recv.vector_grid);
					grid_dst.raw_particles.merge(temp);
				}
			}

			/**
			 * @brief Find axis-aligned bounding box containing all particles
			 *
			 * Computes min/max coordinates across all particles of a given type.
			 * Used to determine grid extents before conversion.
			 * Uses cached particle positions and dispatches to CPU or GPU kernel.
			 */
			void ConvertVDBBase::iolib_find_bbox(
				int particle_type,
				float* bbox_min,
				float* bbox_max,
				float* offset_position
			)
			{
				// Initialize with extreme values in case no particles match
				bbox_min[0] = bbox_min[1] = bbox_min[2] = FLT_MAX;
				bbox_max[0] = bbox_max[1] = bbox_max[2] = -FLT_MAX;

				// Validate particle type
				if (particle_type >= (int)cache_manager.pos_particles_per_ptype.size()) {
					printf("Warning: Invalid particle type %d in iolib_find_bbox\n", particle_type);
					bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0f;
					bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0f;
					return;
				}

				if (cache_manager.skip_cache) {
					// Streaming mode: positions are not cached, so walk the reader once. Same
					// reduction as find_bbox_cpu, just fed straight from the snapshot.
					const size_t no_points = get_local_num_particles();
					float mnx = FLT_MAX, mny = FLT_MAX, mnz = FLT_MAX;
					float mxx = -FLT_MAX, mxy = -FLT_MAX, mxz = -FLT_MAX;

#ifdef WITH_OPENMP
#pragma omp parallel for reduction(min : mnx, mny, mnz) reduction(max : mxx, mxy, mxz) schedule(dynamic, 4096)
#endif
					for (size_t i = 0; i < no_points; i++) {
						if (particle_type != -1 && get_particle_type(i) != particle_type)
							continue;

						double Pos[3];
						get_particle_position(i, Pos);

						const float x = (float)Pos[0] - offset_position[0];
						const float y = (float)Pos[1] - offset_position[1];
						const float z = (float)Pos[2] - offset_position[2];

						if (x < mnx) mnx = x;
						if (y < mny) mny = y;
						if (z < mnz) mnz = z;
						if (x > mxx) mxx = x;
						if (y > mxy) mxy = y;
						if (z > mxz) mxz = z;
					}

					bbox_min[0] = mnx; bbox_min[1] = mny; bbox_min[2] = mnz;
					bbox_max[0] = mxx; bbox_max[1] = mxy; bbox_max[2] = mxz;
					return;
				}

				for (int ptype = 0; ptype < (int)cache_manager.pos_particles_per_ptype.size(); ++ptype) {
					if (particle_type != -1 && ptype != particle_type)
						continue;

					// Get particle positions for this type
					size_t num_particles = cache_manager.pos_particles_per_ptype[ptype].size() / 3;

					if (num_particles == 0) {
						//// No particles of this type
						//bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0f;
						//bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0f;
						continue;
					}

#ifdef WITH_GPU_CUDA
					if (cache_manager.use_gpu) {
						// Use GPU kernel with cached GPU particle positions
						const float* d_pos_particles = cache_manager.d_pos_particles_per_ptype[ptype];
						kernel::find_bbox_gpu(d_pos_particles, num_particles, offset_position, bbox_min, bbox_max);
					}
					else
#endif
					{
						// Use CPU kernel with cached CPU particle positions
						const float* pos_particles = cache_manager.pos_particles_per_ptype[ptype].data();
						kernel::find_bbox_cpu(pos_particles, num_particles, offset_position, bbox_min, bbox_max);
					}
				}
			}

#ifdef WITH_NANOVDB

			std::shared_ptr<nanovdb::tools::build::FloatGrid> ConvertVDBBase::dense_to_nanovdb(VoxelDenseManager* dense_manager, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm)
			{
				std::shared_ptr<nanovdb::tools::build::FloatGrid> nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);

				// Create NanoVDB transform with scale factor and translation offset

				nano_grid->setTransform(
					transform_scale,
					nanovdb::Vec3d(
						dense_manager->offset[0] * transform_scale,
						dense_manager->offset[1] * transform_scale,
						dense_manager->offset[2] * transform_scale
					)
				);

#ifdef WITH_GPU_CUDA
				common::vdb::dense::VoxelGPUDenseManager* dense_manager_gpu = dynamic_cast<common::vdb::dense::VoxelGPUDenseManager*>(dense_manager);
				if (dense_manager_gpu)
				{
					// Use GPU-specific operations if needed
					dense_manager_gpu->from_device();
				}
#endif
				// eVoxelVolume: convert per-voxel accumulated mass to physical density (value / voxel world volume).
				// Result is independent of bbox size and grid resolution -> consistent intensity across zoom levels.
				if (dense_norm == common::SpaceData::DenseNorm::eVoxelVolume && transform_scale > 0.0) {
					const float inv_voxel_volume = (float)(1.0 / (transform_scale * transform_scale * transform_scale));
#pragma omp parallel for
					for (long long i = 0; i < (long long)dense_manager->data_density.size(); i++) {
						float density = dense_manager->data_density[i] * inv_voxel_volume;
						dense_manager->data_density[i] = std::isfinite(density) ? density : 0.0f;
					}
				}
				// Normalize density values using temp buffer if available
				else if (!dense_manager->data_temp.empty()) {
#pragma omp parallel for
					for (int z = 0; z < dense_manager->z(); z++) {
						for (int y = 0; y < dense_manager->y(); y++) {
							for (int x = 0; x < dense_manager->x(); x++) {

								// Get raw density and temp values from dense grid
								size_t index = dense_manager->get_index(x, y, z);
								float density = dense_manager->data_density[index];

								float temp = 0.0f;
								temp = dense_manager->data_temp[index];

								// Apply normalization: divide accumulated density by accumulated weights (temp buffer)
								// This computes the weighted average for SPH-like density estimation
								if (dense_norm != common::SpaceData::DenseNorm::eNone) {
									density = density / temp;
								}

								// If the value is non-zero, set it in the grid
								// isfinite also rejects inf from division by a zero weight (voxel with density but no accumulated norm)
								if (std::isfinite(density)) {
									//accessor.setValue(openvdb::Coord(x + dense_manager->offset[0], y + dense_manager->offset[1], z + dense_manager->offset[2]), density);							
									//if (dense_type == common::SpaceData::DenseType::eType2)
									//	dense_manager->data_density[index] = std::log10(density);
									//else
									dense_manager->data_density[index] = density;
								}
								else {
									dense_manager->data_density[index] = 0.0f;
								}
							}
						}
					}
				}

				/**
				 * @brief Copy normalized dense data to sparse NanoVDB grid
				 *
				 * Iterates through dense array and copies non-zero values to sparse structure.
				 * This maintains memory efficiency by only storing occupied voxels.
				 */
				auto acc_dst = nano_grid->getAccessor();
				for (int z = 0; z < dense_manager->z(); z++) {
					for (int y = 0; y < dense_manager->y(); y++) {
						for (int x = 0; x < dense_manager->x(); x++) {
							nanovdb::Coord xyz(x, y, z);
							float value = dense_manager->data_density[(size_t)x + (size_t)y * (size_t)dense_manager->x() + (size_t)z * (size_t)dense_manager->x() * (size_t)dense_manager->y()];
							// Only store non-zero values to maintain sparse storage efficiency
							// Background value (0.0f) is implicit in NanoVDB
							if (value != 0.0f) {
								acc_dst.setValue(xyz, value);
							}
						}
					}
				}

				return nano_grid;
			}

			std::shared_ptr<nanovdb::tools::build::FloatGrid> ConvertVDBBase::sparse_to_nanovdb(VoxelSparseManager* voxel_manager)
			{
				std::shared_ptr<nanovdb::tools::build::FloatGrid> nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);

				// Create NanoVDB transform with scale factor and translation offset

				nano_grid->setTransform(
					voxel_manager->transform_scale,
					nanovdb::Vec3d(0, 0, 0)
				);

				// Attempt to dynamic_cast to VoxelOpenMPManager
				common::vdb::sparse::VoxelOpenMPManager* voxel_omp_manager = dynamic_cast<common::vdb::sparse::VoxelOpenMPManager*>(voxel_manager);
				if (voxel_omp_manager) {
					auto acc_dst = nano_grid->getAccessor();

					size_t total_occupied = 0;
					for (unsigned int i = 0; i < voxel_omp_manager->table_size; i++) {
						if (voxel_omp_manager->hash_table[i].occupied != 1)
							continue;

						nanovdb::Coord xyz(voxel_omp_manager->hash_table[i].i, voxel_omp_manager->hash_table[i].j, voxel_omp_manager->hash_table[i].k);
						float value = voxel_omp_manager->hash_table[i].value;
						// Only store non-zero values to maintain sparse storage efficiency
						// Background value (0.0f) is implicit in NanoVDB
						if (value != 0.0f) {
							acc_dst.setValue(xyz, value);
							total_occupied++;
						}
					}

					if (cache_manager.world_rank == 0) {
						printf("rank #%d: VoxelOpenMPManager: Total occupied voxels: %f %%\n", cache_manager.world_rank, 100.0f * (float)total_occupied / (float)voxel_omp_manager->insert_count);
					}
				}


#ifdef WITH_GPU_CUDA
				// Attempt to dynamic_cast to VoxelGPUManagerSortReduce
				common::vdb::sparse::VoxelGPUManagerSortReduce* voxel_gpu_manager = dynamic_cast<common::vdb::sparse::VoxelGPUManagerSortReduce*>(voxel_manager);
				if (voxel_gpu_manager) {
					uint64_t* h_keys = new uint64_t[voxel_gpu_manager->m_last_count];
					float* h_vals = new float[voxel_gpu_manager->m_last_count];
					voxel_gpu_manager->get_keys_values_from_device(h_keys, h_vals);

					size_t total_occupied = 0;
					auto acc_dst = nano_grid->getAccessor();
					for (unsigned int i = 0; i < voxel_gpu_manager->m_last_count; i++) {
						int x, y, z;
						common::vdb::sparse::unpackCoord3(h_keys[i], x, y, z);
						nanovdb::Coord xyz(x, y, z);
						float value = h_vals[i];
						// Only store non-zero values to maintain sparse storage efficiency
						// Background value (0.0f) is implicit in NanoVDB
						if (value != 0.0f) {
							acc_dst.setValue(xyz, value);
							total_occupied++;
						}
					}
					delete[] h_keys;
					delete[] h_vals;

					if (cache_manager.world_rank == 0) {
						printf("rank #%d: VoxelGPUManagerSortReduce: Total occupied voxels: %f %%\n", cache_manager.world_rank, 100.0f * (float)total_occupied / (float)voxel_gpu_manager->m_last_count);
					}
				}
#endif

				return nano_grid;
			}

#endif // WITH_NANOVDB

#if defined(WITH_OPENVDB)
			/**
			 * @brief Convert dense grid to sparse OpenVDB grid
			 *
			 * Converts regular dense grid to memory-efficient sparse OpenVDB format.
			 * Applies normalization and uses OpenVDB's optimized dense-to-sparse conversion.
			 */
			openvdb::FloatGrid::Ptr ConvertVDBBase::dense_to_openvdb(VoxelDenseManager* dense_manager, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm)
			{
				openvdb::FloatGrid::Ptr floatgrid = openvdb::FloatGrid::create(0.0f);

				// Configure grid metadata
				floatgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
				std::string grid_name("density");
				floatgrid->setName(grid_name);

				// Create transform with scale and translation offset
				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
				transform->postTranslate(openvdb::Vec3d(
					dense_manager->offset[0] * transform_scale,
					dense_manager->offset[1] * transform_scale,
					dense_manager->offset[2] * transform_scale));

				floatgrid->setTransform(transform);


#ifdef WITH_GPU_CUDA
				common::vdb::dense::VoxelGPUDenseManager* dense_manager_gpu = dynamic_cast<common::vdb::dense::VoxelGPUDenseManager*>(dense_manager);
				if (dense_manager_gpu)
				{
					// Use GPU-specific operations if needed
					dense_manager_gpu->from_device();
				}
#endif

				// eVoxelVolume: convert per-voxel accumulated mass to physical density (value / voxel world volume).
				// Result is independent of bbox size and grid resolution -> consistent intensity across zoom levels.
				if (dense_norm == common::SpaceData::DenseNorm::eVoxelVolume && transform_scale > 0.0) {
					const float inv_voxel_volume = (float)(1.0 / (transform_scale * transform_scale * transform_scale));
#pragma omp parallel for
					for (long long i = 0; i < (long long)dense_manager->data_density.size(); i++) {
						float density = dense_manager->data_density[i] * inv_voxel_volume;
						dense_manager->data_density[i] = std::isfinite(density) ? density : 0.0f;
					}
				}
				// Normalize density values using temp buffer if available
				else if (!dense_manager->data_temp.empty()) {
#pragma omp parallel for
					for (int z = 0; z < dense_manager->z(); z++) {
						for (int y = 0; y < dense_manager->y(); y++) {
							for (int x = 0; x < dense_manager->x(); x++) {

								// Get raw density and temp values from dense grid
								size_t index = dense_manager->get_index(x, y, z);
								float density = dense_manager->data_density[index];

								float temp = 0.0f;
								temp = dense_manager->data_temp[index];

								// Apply normalization: divide accumulated density by accumulated weights (temp buffer)
								// This computes the weighted average for SPH-like density estimation
								if (dense_norm != common::SpaceData::DenseNorm::eNone) {
									density = density / temp;
								}

								// If the value is non-zero, set it in the grid
								// isfinite also rejects inf from division by a zero weight (voxel with density but no accumulated norm)
								if (std::isfinite(density)) {
									//accessor.setValue(openvdb::Coord(x + dense_manager->offset[0], y + dense_manager->offset[1], z + dense_manager->offset[2]), density);							
									//if (dense_type == common::SpaceData::DenseType::eType2)
									//	dense_manager->data_density[index] = std::log10(density);
									//else
									dense_manager->data_density[index] = density;
								}
								else {
									dense_manager->data_density[index] = 0.0f;
								}
							}
						}
					}
				}

				// Use OpenVDB's optimized dense-to-sparse conversion
				// This efficiently identifies non-zero regions and builds the sparse tree
				openvdb::math::CoordBBox bbox(openvdb::Coord(0, 0, 0), openvdb::Coord(dense_manager->x() - 1, dense_manager->y() - 1, dense_manager->z() - 1));
				openvdb::tools::Dense<const float, openvdb::tools::LayoutXYZ> dense(bbox, dense_manager->data_density.data());
				openvdb::tools::copyFromDense(dense, floatgrid->tree(), 0.0f);

				return floatgrid;
			}

			openvdb::FloatGrid::Ptr ConvertVDBBase::sparse_to_openvdb(VoxelSparseManager* voxel_manager) {
				openvdb::FloatGrid::Ptr floatgrid = openvdb::FloatGrid::create(0.0f);

				// Configure grid metadata
				floatgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
				std::string grid_name("density");
				floatgrid->setName(grid_name);

				// Create transform with scale and translation offset
				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(voxel_manager->transform_scale);
				floatgrid->setTransform(transform);

				//// Attempt to dynamic_cast to VoxelCPUManager (works for OpenMP, NanoVDB, OpenVDB managers)
				//common::vdb::sparse::VoxelOpenMPManager* voxel_omp_manager = dynamic_cast<common::vdb::sparse::VoxelOpenMPManager*>(voxel_manager);
				//if (voxel_omp_manager) {
				//	auto acc_dst = floatgrid->getAccessor();

				//	// Use common interface method to extract all voxels
				//	common::vdb::sparse::Voxel* voxels = nullptr;
				//	int voxel_count = voxel_omp_manager->extractAll(&voxels);


				//	// Clean up extracted voxel array
				//	if (voxels != nullptr) {
				//		delete[] voxels;
				//	}

				//	return floatgrid;
				//}

				// Attempt to dynamic_cast to VoxelOpenMPManager
				common::vdb::sparse::VoxelOpenMPManager* voxel_omp_manager = dynamic_cast<common::vdb::sparse::VoxelOpenMPManager*>(voxel_manager);
				if (voxel_omp_manager) {
					auto acc_dst = floatgrid->getAccessor();

					size_t total_occupied = 0;
					for (unsigned int i = 0; i < voxel_omp_manager->table_size; i++) {
						if (voxel_omp_manager->hash_table[i].occupied != 1)
							continue;

						openvdb::Coord xyz(voxel_omp_manager->hash_table[i].i, voxel_omp_manager->hash_table[i].j, voxel_omp_manager->hash_table[i].k);
						float value = voxel_omp_manager->hash_table[i].value;
						// Only store non-zero values to maintain sparse storage efficiency
						// Background value (0.0f) is implicit in OpenVDB
						if (value != 0.0f) {
							acc_dst.setValue(xyz, value);
							total_occupied++;
						}
					}

					printf("VoxelOpenMPManager: Total occupied voxels: %f %%\n", 100.0f * (float)total_occupied / (float)voxel_omp_manager->insert_count);
				}


#ifdef WITH_GPU_CUDA
			// Attempt to dynamic_cast to VoxelGPUManagerSortReduce
				common::vdb::sparse::VoxelGPUManagerSortReduce* voxel_gpu_manager = dynamic_cast<common::vdb::sparse::VoxelGPUManagerSortReduce*>(voxel_manager);
				if (voxel_gpu_manager) {
					uint64_t* h_keys = new uint64_t[voxel_gpu_manager->m_last_count];
					float* h_vals = new float[voxel_gpu_manager->m_last_count];
					voxel_gpu_manager->get_keys_values_from_device(h_keys, h_vals);

					size_t total_occupied = 0;
					auto acc_dst = floatgrid->getAccessor();
					for (unsigned int i = 0; i < voxel_gpu_manager->m_last_count; i++) {
						int x, y, z;
						common::vdb::sparse::unpackCoord3(h_keys[i], x, y, z);
						openvdb::Coord xyz(x, y, z);
						float value = h_vals[i];
						// Only store non-zero values to maintain sparse storage efficiency
						// Background value (0.0f) is implicit in OpenVDB
						if (value != 0.0f) {
							acc_dst.setValue(xyz, value);
							total_occupied++;
						}
					}
					delete[] h_keys;
					delete[] h_vals;

					printf("VoxelGPUManagerSortReduce: Total occupied voxels: %f %%\n", 100.0f * (float)total_occupied / (float)voxel_gpu_manager->m_last_count);
				}
#endif

				return floatgrid;
			}

			/**
			 * @brief Serialize OpenVDB grid to binary buffer
			 *
			 * Converts OpenVDB grid to binary format for MPI transfer or file I/O.
			 * Uses OpenVDB's streaming API for efficient serialization.
			 */
			void ConvertVDBBase::openvdb_to_vector(openvdb::FloatGrid::Ptr grid, std::vector<uint8_t>& file_content)
			{
				// Serialize grid to string stream in binary mode
				std::ostringstream stream(std::ios_base::binary);
				openvdb::io::Stream(stream).write({ grid });
				stream.flush();

				// Convert string stream to byte vector
				const std::string& str = stream.str();
				file_content.assign(str.begin(), str.end());
			}

			/**
			 * @brief Serialize two OpenVDB grids to binary buffer
			 *
			 * Similar to openvdb_to_vector but handles two grids in one stream.
			 */
			void ConvertVDBBase::openvdb_to_vector2(openvdb::FloatGrid::Ptr grid1, openvdb::FloatGrid::Ptr grid2, std::vector<uint8_t>& file_content)
			{
				// Serialize both grids to single stream
				std::ostringstream stream(std::ios_base::binary);
				openvdb::io::Stream(stream).write({ grid1, grid2 });
				stream.flush();

				// Convert to byte vector
				const std::string& str = stream.str();
				file_content.assign(str.begin(), str.end());
			}

			/**
			 * @brief Deserialize OpenVDB grid from binary buffer
			 *
			 * Reconstructs OpenVDB grid from binary data received via MPI or loaded from file.
			 */
			openvdb::FloatGrid::Ptr ConvertVDBBase::vector_to_openvdb(std::vector<uint8_t>& file_content)
			{
				// Convert byte vector back to string stream
				std::istringstream stream(std::string(file_content.begin(), file_content.end()), std::ios_base::binary);

				// Create VDB input stream and read grids
				openvdb::io::Stream vdbStream(stream);
				openvdb::GridPtrVecPtr grids = vdbStream.getGrids();

				// Find and return the first FloatGrid in the collection
				for (auto& grid : *grids) {
					if (grid->isType<openvdb::FloatGrid>()) {
						return openvdb::gridPtrCast<openvdb::FloatGrid>(grid);
					}
				}

				// No FloatGrid found in stream
				return nullptr;
			}

#endif

			/**
			 * @brief Get particle density (rho) value
			 *
			 * Returns density from precomputed k-NN results if available,
			 * otherwise falls back to internal particle data.
			 */
			double ConvertVDBBase::get_particle_rho(uint64_t id) {
				double mass = get_particle_mass(id);
				if (mass != 0.0) {
					// Check if we have precomputed density for this particle
					int particle_type = get_particle_type(id);
					if (cache_manager.rho_particles_per_ptype.size() > 0 && cache_manager.rho_particles_per_ptype[particle_type].size() > 0) {
						return cache_manager.rho_particles_per_ptype[particle_type][id - cache_manager.particles_ptype_offset[particle_type]];
					}
				}

				return get_particle_rho_internal(id);
			}

			/**
			 * @brief Get types and available data blocks
			 *
			 * Returns availability of data blocks for each particle type.
			 * Includes precomputed density from k-NN search if available.
			 */
			void ConvertVDBBase::get_types_and_blocks(std::vector<int>& types_and_blocks) {
				get_types_and_blocks_internal(types_and_blocks);

				int num_types = get_num_types();
				int rho_blocknr = get_particle_rho_blocknr();

				// A reader without a rho block (e.g. HACC GenericIO without --rho-name)
				// reports -1 — nothing to advertise then
				if (rho_blocknr < 0) {
					return;
				}

				// Mark density block as available if we have precomputed rho values
				for (int type = 0; type < num_types; type++) {
					if (cache_manager.rho_particles_per_ptype.size() > 0 && cache_manager.rho_particles_per_ptype[type].size() > 0) {
						types_and_blocks[num_types * rho_blocknr + type]++;
					}
				}
			}

			/**
			 * @brief Get normalized particle value for a specific data block
			 *
			 * Returns particle property value (density, temperature, etc.).
			 * Handles special case for density block with k-NN precomputed values.
			 */
			float ConvertVDBBase::get_particle_norm_value(int blocknr, uint64_t id) {
				if (blocknr == get_particle_rho_blocknr()) {
					return get_particle_rho(id);
				}

				return get_particle_norm_value_internal(blocknr, id);
			}

			int ConvertVDBBase::get_particle_value(int blocknr, uint64_t id, float* out_value) {
				if (blocknr == get_particle_rho_blocknr()) {
					// Return precomputed or internal density value
					*out_value = get_particle_rho(id);
					return 1; // Density is a scalar (1 component)
				}

				return get_particle_value_internal(blocknr, id, out_value);
			}

			/**
			 * @brief Get number of components for a particle value
			 *
			 * Returns 1 for scalars, 3 for vectors, etc.
			 */
			int ConvertVDBBase::get_particle_value_comp(int blocknr, uint64_t id) {
				return get_particle_value_comp_internal(blocknr, id);
			}

			/**
			 * @brief Format filename with zero-padded number
			 *
			 * Replaces "{}" placeholders in pattern with formatted number.
			 * Numbers < 1000 are zero-padded to 3 digits (e.g., "001", "042").
			 */
			std::string ConvertVDBBase::format_filename(const std::string& pattern, int number, int zero_pad) {
				// Create the formatted number as a string
				std::ostringstream formattedNumber;
				if (zero_pad > 0) {
					// Zero-pad numbers to specified width for consistent file sorting
					formattedNumber << std::setw(zero_pad) << std::setfill('0') << number;
				}
				else {
					// Use longer format for numbers >= 1000
					formattedNumber << number;
				}

				// Replace all "{}" placeholders with the formatted number
				std::string result = pattern;
				std::string placeholder = "{}";
				size_t pos = result.find(placeholder);

				while (pos != std::string::npos) {
					// Replace the first occurrence of the placeholder with the formatted number
					result.replace(pos, placeholder.length(), formattedNumber.str());

					// Find the next occurrence of the placeholder
					pos = result.find(placeholder, pos + formattedNumber.str().length());
				}

				return result;
			}

		}//vdb
	}  //common
}  // namespace space_converter