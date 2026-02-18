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

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/NanoToOpenVDB.h>
#else
#	include <nanovdb/tools/NanoToOpenVDB.h>
#endif
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
					printf("File %s does exist or cannot be opened\n", calc_radius_neigh_file_ptype.c_str());
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

				cache_manager.radius_particles_per_ptype.push_back(radius_particles); // Store the radius particles for this type
				cache_manager.rho_particles_per_ptype.push_back(rho_particles); // Store the rho particles for this type

				cache_manager.particles_ptype_offset[ptype + 1] = cache_manager.particles_ptype_offset[ptype] + radius_particles.size();

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
			double t_start = omp_get_wtime();

			size_t no_points = get_local_num_particles();
			int ptype_count = get_num_types();
			
			cache_manager.pos_particles_per_ptype.resize(ptype_count);
			cache_manager.radius_particles_per_ptype.resize(ptype_count);
			cache_manager.rho_particles_per_ptype.resize(ptype_count);
			cache_manager.mass_particles_per_ptype.resize(ptype_count);
			cache_manager.particles_ptype_offset.resize(ptype_count + 1, 0);

			int num_threads = omp_get_max_threads();

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
					int tid = omp_get_thread_num();

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
				cache_manager.particles_id_ordered_per_ptype[ptype] = std::move(particles_id_master);
				cache_manager.particles_ptype_offset[ptype + 1] = cache_manager.particles_ptype_offset[ptype] + cache_manager.pos_particles_per_ptype[ptype].size() / 3;
			}

			printf("find_particle_positions: Find positions: %f\n", omp_get_wtime() - t_start);
		}

#ifdef WITH_CUDAKDTREE
		/**
		 * @brief Calculate particle radii using CUDA-accelerated k-NN search
		 * 
		 * Uses GPU-accelerated KD-tree to find k-nearest neighbors for each particle.
		 * Computes smoothing radius and density based on neighbor distances.
		 */
		void ConvertVDBBase::calculate_radius_by_cudakdtree(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, bool use_cudakdtree_cpu, float maxRadius, common::SpaceData::DenseType &rho_kernel)
		{
			double t_start = omp_get_wtime();

			if (calc_radius_neigh == -1) {
				read_radius_from_file(calc_radius_neigh_file);
			}
			else {
				size_t no_points = get_local_num_particles();

				int ptype_count = get_num_types();
				// cache_manager.radius_particles_per_ptype.resize(ptype_count);
				// cache_manager.rho_particles_per_ptype.resize(ptype_count);
				// cache_manager.particles_ptype_offset.resize(ptype_count + 1, 0);

				int num_threads = omp_get_max_threads();

				for (int ptype = 0; ptype < ptype_count; ptype++) {

// 					std::vector<float> points;
// 					std::vector < std::vector<float> > points_thread(num_threads);

// 					std::vector<float> pmass;
// 					std::vector < std::vector<float> > pmass_thread(num_threads);

// #pragma omp parallel num_threads(num_threads) 
// 					{
// 						int tid = omp_get_thread_num();

// #pragma omp for
// 						for (size_t i = 0; i < no_points; i++) {

// 							if (get_particle_type(i) != ptype)
// 								continue;

// 							double Pos[3];
// 							get_particle_position(i, Pos);

// 							// Collect particle positions per thread
// 							points_thread[tid].push_back(Pos[0]);
// 							points_thread[tid].push_back(Pos[1]);
// 							points_thread[tid].push_back(Pos[2]);

// 							double mass = get_particle_mass(i);
// 							pmass_thread[tid].push_back(mass);
// 						}
// 					}

// 					// Merge thread-local results in order to maintain particle index consistency
// 					// Thread order must be preserved to match particle IDs with computed radii
// 					std::vector<float> result;
// 					for (int t = 0; t < num_threads; ++t) {
// 						points.insert(points.end(), points_thread[t].begin(), points_thread[t].end());
// 						pmass.insert(pmass.end(), pmass_thread[t].begin(), pmass_thread[t].end());
// 					}

// 					cache_manager.particles_ptype_offset[ptype + 1] = cache_manager.particles_ptype_offset[ptype] + points.size() / 3;

// 					printf("cudakdtree: init: %f\n", omp_get_wtime() - t_start);

					// Run GPU/CPU-based k-Nearest Neighbor search
					utility::cudakdtree::run_knn(cache_manager.pos_particles_per_ptype[ptype].data(), cache_manager.pos_particles_per_ptype[ptype].size() / 3, calc_radius_neigh + 1, cache_manager.radius_particles_per_ptype[ptype], cache_manager.rho_particles_per_ptype[ptype], cache_manager.mass_particles_per_ptype[ptype], !use_cudakdtree_cpu, use_cycling, maxRadius, rho_kernel);
				}

				if (calc_radius_neigh_file.length() > 0) {
					write_radius_from_file(calc_radius_neigh_file);
				}
			}

			printf("cudakdtree: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
		}

#endif		


#ifdef WITH_NANOFLANN
		/**
		 * @brief Calculate particle radii using CPU-based nanoflann k-NN search
		 * 
		 * Uses CPU-based KD-tree (nanoflann library) to find k-nearest neighbors.
		 * Computes smoothing radius and density based on neighbor distances.
		 */
		void ConvertVDBBase::calculate_radius_by_nanoflann(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, common::SpaceData::DenseType& rho_kernel)
		{
			double t_start = omp_get_wtime();

			if (calc_radius_neigh == -1) {
				read_radius_from_file(calc_radius_neigh_file);
			}
			else {
				size_t no_points = get_local_num_particles();

				int ptype_count = get_num_types();
				// cache_manager.radius_particles_per_ptype.resize(ptype_count);
				// cache_manager.rho_particles_per_ptype.resize(ptype_count);
				// cache_manager.particles_ptype_offset.resize(ptype_count + 1, 0);

				int num_threads = omp_get_max_threads();

				for (int ptype = 0; ptype < ptype_count; ptype++) {

// 					std::vector<float> points;
// 					std::vector < std::vector<float> > points_thread(num_threads);

// 					std::vector<float> pmass;
// 					std::vector < std::vector<float> > pmass_thread(num_threads);

// #pragma omp parallel num_threads(num_threads) 
// 					{
// 						int tid = omp_get_thread_num();

// #pragma omp for
// 						for (size_t i = 0; i < no_points; i++) {

// 							if (get_particle_type(i) != ptype)
// 								continue;

// 							double Pos[3];
// 							get_particle_position(i, Pos);

// 							// Collect particle positions per thread (x, y, z interleaved)
// 							points_thread[tid].push_back(Pos[0]);
// 							points_thread[tid].push_back(Pos[1]);
// 							points_thread[tid].push_back(Pos[2]);

// 							// Store mass for density computation
// 							double mass = get_particle_mass(i);
// 							pmass_thread[tid].push_back(mass);
// 						}
// 					}

// 					// Merge thread-local results in thread order to maintain particle index consistency
// 					// This is critical for correct radius lookup later
// 					std::vector<float> result;
// 					for (int t = 0; t < num_threads; ++t) {
// 						points.insert(points.end(), points_thread[t].begin(), points_thread[t].end());
// 						pmass.insert(pmass.end(), pmass_thread[t].begin(), pmass_thread[t].end());
// 					}

// 					// Store particle count offset for this type (cumulative sum)
// 					cache_manager.particles_ptype_offset[ptype + 1] = cache_manager.particles_ptype_offset[ptype] + points.size() / 3;

// 					printf("nanoflann: init: %f\n", omp_get_wtime() - t_start);

					if (cache_manager.pos_particles_per_ptype[ptype].size() > 0) {
						// Run k-NN search (calc_radius_neigh + 1 includes the query point itself)
						utility::nanoflann_tool::run_knn(cache_manager.pos_particles_per_ptype[ptype].data(), cache_manager.pos_particles_per_ptype[ptype].size() / 3, calc_radius_neigh + 1, cache_manager.radius_particles_per_ptype[ptype], cache_manager.rho_particles_per_ptype[ptype], cache_manager.mass_particles_per_ptype[ptype], use_cycling, rho_kernel);
					}
				}

				if (calc_radius_neigh_file.length() > 0) {
					write_radius_from_file(calc_radius_neigh_file);
				}
			}

			printf("nanoflann: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
		}
#endif

		/**
		 * @brief Main conversion function: converts particle data to VDB grids
		 * 
		 * This is the core function that:
		 * 1. Iterates through all particles
		 * 2. Filters by type, bounding box, and value ranges
		 * 3. Rasterizes particles into voxel grid (dense or sparse)
		 * 4. Applies smoothing kernels for SPH-like density estimation
		 * 5. Tracks min/max values for normalization
		 */
		void ConvertVDBBase::convert_iolib_to_grid(
			int particle_type,
			float particle_fix_size,
			std::string grid_name,
			//int* grid_dims,
			float grid_transform,
			float* bbox_min,
			float* bbox_max,
			int bbox_dim,
			int* bbox_min_orig,
			double bbox_size_orig,
			common::SpaceData::ExtractedType extracted_type,
			common::SpaceData::DenseType dense_type,
			common::SpaceData::DenseNorm dense_norm,
			int block_name_id,
			float object_size,
			float& min_value,
			float& max_value,
			float min_value_global,
			float max_value_global,
			size_t& particles_count,
			VDBParticles& grid,
			double& transform_scale,
			float filter_min,
			float filter_max,
			float min_rho,
			float max_rho,
			common::SpaceData::AnimType anim_type,
			int frame_req,
			int frame,
			float *bbox_sphere_pos,
			float bbox_sphere_r,
			bool use_simple_density,
			float *offset_position
		)
		{
#ifdef WITH_OPENMP
			double t = omp_get_wtime();
			double t_step = omp_get_wtime();
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
			if (grid.type == VDBParticles::VDBParticleType::eNanoVDB) {
				grid.nano_grid->setTransform(transform_scale);
			}
#ifdef WITH_OPENVDB
			else if (grid.type == VDBParticles::VDBParticleType::eOpenVDB) {
				// Create linear transform for OpenVDB grid
				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
				grid.vdb_grid->setTransform(transform);
			}
#endif
			else if (grid.type == VDBParticles::VDBParticleType::eRawParticles) {
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

			if (grid.type == VDBParticles::VDBParticleType::eDense) {
				// Set grid offset for dense grids (position in global coordinate system)
				// This allows tiled dense grids to know their position in the full domain
				grid.dense_grid.offset[0] = (size_t)(bbox_x_min_norm * (double)bbox_dim / scale_space_diagonal);
				grid.dense_grid.offset[1] = (size_t)(bbox_y_min_norm * (double)bbox_dim / scale_space_diagonal);
				grid.dense_grid.offset[2] = (size_t)(bbox_z_min_norm * (double)bbox_dim / scale_space_diagonal);
			}

			// Track value range and particle count (for parallel reduction)
			float min = FLT_MAX;
			float max = -FLT_MAX;

			size_t particles_count_temp = 0;

#ifdef WITH_OPENMP			
			int nthreads = 1;
			if (grid.type == VDBParticles::VDBParticleType::eDense) {
				// Dense grids support parallel writes (no race conditions)
				// Sparse grids (NanoVDB, OpenVDB) require serial access
				nthreads = omp_get_max_threads();
			}

#pragma omp parallel for reduction(min : min) reduction(max : max) reduction(+ : particles_count_temp) num_threads(nthreads) schedule(dynamic, 256)
#endif
			for (size_t i = 0; i < no_points; i++) {
				// Filter by particle type
				if (get_particle_type(i) != particle_type)
					continue;

				// Filter by animation frame if extracting specific frame
				if (anim_type == common::SpaceData::AnimType::eFrameExtract) {
					if (frame_req != frame)
						continue;
				}

				double Pos[3];
				get_particle_position(i, Pos);

				// Apply position offset if specified
				if (offset_position[0] != 0.0f || offset_position[1] != 0.0f || offset_position[2] != 0.0f) {
					Pos[0] -= offset_position[0];
					Pos[1] -= offset_position[1];
					Pos[2] -= offset_position[2];
				}

				// Normalize position to [0,1] range based on bounding box
				double px_norm = ((double)Pos[0] - (double)bbox_min_orig[0]) / bbox_size_orig;
				double py_norm = ((double)Pos[1] - (double)bbox_min_orig[1]) / bbox_size_orig;
				double pz_norm = ((double)Pos[2] - (double)bbox_min_orig[2]) / bbox_size_orig;

				// Check if particle is within bounding box
				// Dense grids need extended bounds (particle radius) for proper splatting
				if (grid.type == VDBParticles::VDBParticleType::eDense) {
					// Get particle radius for boundary check with tolerance
					double radiusxyz_max = get_particle_radius(
						i,
						1,
						bbox_min_orig,
						bbox_size_orig,
						1,
						particle_fix_size,
						particle_type
					);

					if (px_norm + radiusxyz_max < bbox_x_min_norm || px_norm - radiusxyz_max > bbox_x_max_norm)
						continue;

					if (py_norm + radiusxyz_max < bbox_y_min_norm || py_norm - radiusxyz_max > bbox_y_max_norm)
						continue;

					if (pz_norm + radiusxyz_max < bbox_z_min_norm || pz_norm - radiusxyz_max > bbox_z_max_norm)
						continue;
				}
				else {
				// Simple point-in-box test for sparse grids
						continue;

					if (py_norm < bbox_y_min_norm || py_norm > bbox_y_max_norm)
						continue;

					if (pz_norm < bbox_z_min_norm || pz_norm > bbox_z_max_norm)
						continue;
				}

				// Convert normalized position to voxel coordinates
				int px = static_cast<int>((double)px_norm * (double)bbox_dim / scale_space_diagonal);
				int py = static_cast<int>((double)py_norm * (double)bbox_dim / scale_space_diagonal);
				int pz = static_cast<int>((double)pz_norm * (double)bbox_dim / scale_space_diagonal);

				// Get particle value (density, temperature, etc.)
				float v_orig = get_particle_norm_value(block_name_id, i);
				if (use_simple_density) {
					// Simple density mode: treat all particles with uniform weight
					v_orig = 1;
				}

				// Apply value range filter
				if (v_orig < filter_min || v_orig > filter_max)
					continue;

				// Apply spherical bounding filter if specified (radius > 0)
				if(bbox_sphere_r > 0.0f)
				{
					// Compute squared distance from particle to sphere center
					// Using squared distance avoids expensive sqrt operation
					double distSquared = 
					(Pos[0] - bbox_sphere_pos[0]) * (Pos[0] - bbox_sphere_pos[0]) +
					(Pos[1] - bbox_sphere_pos[1]) * (Pos[1] - bbox_sphere_pos[1]) +
					(Pos[2] - bbox_sphere_pos[2]) * (Pos[2] - bbox_sphere_pos[2]);

					// Skip particles outside the spherical region
					if (distSquared > (bbox_sphere_r * bbox_sphere_r)) {
						continue;
					}
				}
				particles_count_temp++;

				// Track min/max values for normalization
				if (v_orig < min) min = v_orig;
				if (v_orig > max) max = v_orig;

				// Rasterize particle into grid based on grid type
				float v = v_orig;

				if (grid.type == VDBParticles::VDBParticleType::eDense) {
					// Splat particle into dense grid using kernel function (SPH-like smoothing)
					fill_voxels(grid.dense_grid,
						i, v,
						bbox_dim, bbox_min_orig, bbox_size_orig,
						scale_space_diagonal,
						dense_type, dense_norm,
						particle_fix_size, particle_type, block_name_id, Pos);
				}
				else if (grid.type == VDBParticles::VDBParticleType::eNanoVDB) {
					nanovdb::Coord xyz(px, py, pz);
					auto acc = grid.nano_grid->getAccessor();
					// Accumulate value if voxel already contains data (multiple particles per voxel)
					if (acc.isValueOn(xyz)) {
						v += acc.getValue(xyz);
					}

					acc.setValue(xyz, v);
				}
#ifdef WITH_OPENVDB
				else if (grid.type == VDBParticles::VDBParticleType::eOpenVDB) {
					openvdb::Coord xyz(px, py, pz);
					auto acc = grid.vdb_grid->getAccessor();

					// Use probeValue for faster lookup (avoids double tree traversal: check + set)
					float dst_value;
					bool dst_exists = acc.probeValue(xyz, dst_value);

					if (dst_exists) {
						// Accumulate values (multiple particles contribute to same voxel)
						acc.setValue(xyz, v + dst_value);
					}
					else {
						// First particle contributing to this voxel
						acc.setValue(xyz, v);
					}
				}
#endif
				else if (grid.type == VDBParticles::VDBParticleType::eRawParticles) {
					// Store raw particle data (no rasterization)
					// Position is stored in world space coordinates
					raw_positions.values.push_back((double)px_norm * (double)object_size);
					raw_positions.values.push_back((double)py_norm * (double)object_size);
					raw_positions.values.push_back((double)pz_norm * (double)object_size);

					// Get particle attribute values (may be multi-component, e.g., velocity vector)
					int n_comp = get_particle_value_comp(block_name_id, i);
					std::vector<float> values(n_comp);
					get_particle_value(block_name_id, i, values.data());

					raw_values.num_comp = n_comp;
					raw_values.values.insert(raw_values.values.end(), values.begin(), values.end());
	
					double pr = get_particle_radius(
						i,
						object_size,
						bbox_min_orig,
						bbox_size_orig,
						1.0f,
						particle_fix_size,
						particle_type
					);

					raw_radius.values.push_back(pr);
				
					// Store animation frame number for temporal data
					raw_frame.values.push_back(frame);
				}
			}

			// Store final particle count after filtering
			particles_count = particles_count_temp;

			// Store raw particle data arrays in output structure
			if (grid.type == VDBParticles::VDBParticleType::eRawParticles) {
				grid.raw_particles.data.push_back(raw_positions);
				grid.raw_particles.data.push_back(raw_values);
				grid.raw_particles.data.push_back(raw_radius);
				grid.raw_particles.data.push_back(raw_frame);
			}

			min_value = min;
			max_value = max;

#ifdef WITH_OPENMP
			t_step = omp_get_wtime();
#endif
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
			if (grid_dst.type == VDBParticles::VDBParticleType::eDense && grid_recv.type == VDBParticles::VDBParticleType::eDense) {
				// Merge dense grids by summing voxel values (element-wise addition)
#pragma omp parallel for
				for (int z = 0; z < grid_dst.dense_grid.dims[2]; z++) {
					for (int y = 0; y < grid_dst.dense_grid.dims[1]; y++) {
						for (int x = 0; x < grid_dst.dense_grid.dims[0]; x++) {
							size_t index = grid_dst.dense_grid.get_index(x, y, z);
							grid_dst.dense_grid.data_density[index] += grid_recv.dense_grid.data_density[index];
#ifndef WITH_NO_DATA_TEMP							
							grid_dst.dense_grid.data_temp[index] += grid_recv.dense_grid.data_temp[index];
#endif							
						}
					}
				}
			}

			else if (grid_dst.type == VDBParticles::VDBParticleType::eNanoVDB && grid_recv.type == VDBParticles::VDBParticleType::eVector) {
				// Deserialize and merge sparse NanoVDB grid
				auto acc_dst = grid_dst.nano_grid->getAccessor();
				auto* grid_src_float = (nanovdb::NanoGrid<float>*)grid_recv.vector_grid.data();

				// Traverse NanoVDB sparse tree hierarchy ( Root -> Upper Internal -> Lower Internal -> Leaf )
				for (auto it2 = grid_src_float->tree().root().cbeginChild(); it2; ++it2) {
					// Iterate over upper internal nodes (level 2)
					for (auto it1 = it2->cbeginChild(); it1; ++it1) {
						// Iterate over lower internal nodes (level 1)
						for (auto it0 = it1->cbeginChild(); it0; ++it0) {
							// Iterate over active voxels in leaf nodes
							for (auto it = it0->cbeginValueOn(); it; ++it) {
								float v = *it;

								nanovdb::Coord xyz = it.getCoord();

								// Accumulate values if voxel exists in destination, otherwise set new value
								if (acc_dst.isValueOn(xyz)) {
									v += acc_dst.getValue(xyz);
								}
								acc_dst.setValue(xyz, v);
							}
						}
					}
				}
			}
#ifdef WITH_OPENVDB
			else if (grid_dst.type == VDBParticles::VDBParticleType::eOpenVDB && grid_recv.type == VDBParticles::VDBParticleType::eVector) {
				// Deserialize OpenVDB grid from vector
				auto grid_src_float = vector_to_openvdb(grid_recv.vector_grid);

				// Use OpenVDB's optimized composite operation for summing grids
				openvdb::tools::compSum(*grid_dst.vdb_grid, *grid_src_float);
			}
#endif
			else if (grid_dst.type == VDBParticles::VDBParticleType::eRawParticles && grid_recv.type == VDBParticles::VDBParticleType::eVector) {
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
			if (particle_type < 0 || particle_type >= (int)cache_manager.pos_particles_per_ptype.size()) {
				printf("Warning: Invalid particle type %d in iolib_find_bbox\n", particle_type);
				bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0f;
				bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0f;
				return;
			}

			// Get particle positions for this type
			size_t num_particles = cache_manager.pos_particles_per_ptype[particle_type].size() / 3;

			if (num_particles == 0) {
				// No particles of this type
				bbox_min[0] = bbox_min[1] = bbox_min[2] = 0.0f;
				bbox_max[0] = bbox_max[1] = bbox_max[2] = 0.0f;
				return;
			}

#ifdef WITH_GPU_CUDA
			if (cache_manager.use_gpu_cuda) {
				// Use GPU kernel with cached GPU particle positions
				const float* d_pos_particles = cache_manager.d_pos_particles_per_ptype[particle_type];
				kernel::find_bbox_gpu(d_pos_particles, num_particles, offset_position, bbox_min, bbox_max);
			}
			else
#endif
			{
				// Use CPU kernel with cached CPU particle positions
				const float* pos_particles = cache_manager.pos_particles_per_ptype[particle_type].data();
				kernel::find_bbox_cpu(pos_particles, num_particles, offset_position, bbox_min, bbox_max);
			}
		}

		/**
		 * @brief Find min/max value range for a particle data block
		 * 
		 * Scans all particles to determine value range for normalization.
		 */
		void ConvertVDBBase::iolib_find_minmax(
			int particle_type,
			int block_nr,
			float& v_min,
			float& v_max
		)
		{
			size_t no_points = get_local_num_particles();

			float min = FLT_MAX;
			float max = -FLT_MAX;

			// Parallel reduction to find min/max values
#pragma omp parallel for reduction(min : min) reduction(max : max)
			for (size_t i = 0; i < no_points; ++i) {
				if (get_particle_type(i) != particle_type)
					continue;

				float v = get_particle_norm_value(block_nr, i);

				if (v < min) min = v;
				if (v > max) max = v;
			}

			v_min = min;
			v_max = max;
		}

#ifdef WITH_NANOVDB

#if OPENVDB_VERSION == 11
		/**
		 * @brief Convert dense grid to sparse NanoVDB grid (OpenVDB v11)
		 * 
		 * Converts regular dense grid to memory-efficient sparse NanoVDB format.
		 * Applies normalization and only stores non-zero voxels.
		 */
		std::shared_ptr<nanovdb::build::FloatGrid> ConvertVDBBase::dense_to_nanovdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm)
#else
		std::shared_ptr<nanovdb::tools::build::FloatGrid> ConvertVDBBase::dense_to_nanovdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm)
#endif
		{

#if OPENVDB_VERSION == 11
			std::shared_ptr<nanovdb::build::FloatGrid> nano_grid = std::make_shared<nanovdb::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#else
			std::shared_ptr<nanovdb::tools::build::FloatGrid> nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
#endif

// Create NanoVDB transform with scale factor and translation offset

			nano_grid->setTransform(
				transform_scale,
				nanovdb::Vec3d(
					particles.offset[0] * transform_scale,
					particles.offset[1] * transform_scale,
					particles.offset[2] * transform_scale
				)
			);

			// Normalize density values using temp buffer if available
#ifndef WITH_NO_DATA_TEMP
#pragma omp parallel for
			for (int z = 0; z < particles.z(); z++) {
				for (int y = 0; y < particles.y(); y++) {
					for (int x = 0; x < particles.x(); x++) {

						// Get raw density and temp values from dense grid
						size_t index = particles.get_index(x, y, z);
						float density = particles.data_density[index];
						
						float temp = 0.0f;
						temp = particles.data_temp[index];

						// Apply normalization: divide accumulated density by accumulated weights (temp buffer)
						// This computes the weighted average for SPH-like density estimation
						if (dense_norm != common::SpaceData::DenseNorm::eNone) {
							density = density / temp;
						}

						// If the value is non-zero, set it in the grid
						if (!std::isnan(density)) {
							//accessor.setValue(openvdb::Coord(x + particles.offset[0], y + particles.offset[1], z + particles.offset[2]), density);							
							//if (dense_type == common::SpaceData::DenseType::eType2)
							//	particles.data_density[index] = std::log10(density);
							//else
							particles.data_density[index] = density;
						}
						else {
							particles.data_density[index] = 0.0f;
						}
					}
				}
			}
#endif

		/**
		 * @brief Copy normalized dense data to sparse NanoVDB grid
		 * 
		 * Iterates through dense array and copies non-zero values to sparse structure.
		 * This maintains memory efficiency by only storing occupied voxels.
		 */
			auto acc_dst = nano_grid->getAccessor();
			for (int z = 0; z < particles.z(); z++) {
				for (int y = 0; y < particles.y(); y++) {
					for (int x = 0; x < particles.x(); x++) {
						nanovdb::Coord xyz(x, y, z);
						float value = particles.data_density[(size_t)x + (size_t)y * (size_t)particles.x() + (size_t)z * (size_t)particles.x() * (size_t)particles.y()];
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
#endif

#if defined(WITH_OPENVDB)
		/**
		 * @brief Convert dense grid to sparse OpenVDB grid
		 * 
		 * Converts regular dense grid to memory-efficient sparse OpenVDB format.
		 * Applies normalization and uses OpenVDB's optimized dense-to-sparse conversion.
		 */
		openvdb::FloatGrid::Ptr ConvertVDBBase::dense_to_openvdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm)
		{
			openvdb::FloatGrid::Ptr floatgrid = openvdb::FloatGrid::create(0.0f);

			// Configure grid metadata
			floatgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
			std::string grid_name("density");
			floatgrid->setName(grid_name);

			// Create transform with scale and translation offset
			openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
			transform->postTranslate(openvdb::Vec3d(particles.offset[0] * transform_scale, particles.offset[1] * transform_scale, particles.offset[2] * transform_scale));
			floatgrid->setTransform(transform);

			// Normalize density values using temp buffer if available
#ifndef WITH_NO_DATA_TEMP
#pragma omp parallel for
			for (int z = 0; z < particles.z(); z++) {
				for (int y = 0; y < particles.y(); y++) {
					for (int x = 0; x < particles.x(); x++) {

						// Get the value from the array						
						size_t index = particles.get_index(x, y, z);
						float density = particles.data_density[index];

						float temp = 0.0f;
						temp = particles.data_temp[index];						

						// Apply normalization if enabled
						if (dense_norm != common::SpaceData::DenseNorm::eNone) {
							density = density / temp;
						}

						// Store normalized value or zero for invalid data
						if (!std::isnan(density)) {
							particles.data_density[index] = density;
						}
						else {
							particles.data_density[index] = 0.0f;
						}
					}
				}
			}
#endif

			// Use OpenVDB's optimized dense-to-sparse conversion
			// This efficiently identifies non-zero regions and builds the sparse tree
			openvdb::math::CoordBBox bbox(openvdb::Coord(0, 0, 0), openvdb::Coord(particles.x() - 1, particles.y() - 1, particles.z() - 1));
			openvdb::tools::Dense<const float, openvdb::tools::LayoutXYZ> dense(bbox, particles.data_density.data());
			openvdb::tools::copyFromDense(dense, floatgrid->tree(), 0.0f);

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
		 * @brief Get particle radius for rasterization
		 * 
		 * Computes particle smoothing radius from:
		 * - Precomputed neighbor search results (if available)
		 * - SPH smoothing length (hsml)
		 * - Constant radius setting
		 * - Particle mass and density
		 */
		double ConvertVDBBase::get_particle_radius(
			uint64_t pid,
			int bbox_dim,
			int* bbox_min_orig,
			double bbox_size_orig,
			double scale_space_diagonal,
			float particle_fix_size,
			int particle_type
		) {
			if (cache_manager.radius_particle_const > 0.0) {
				return cache_manager.radius_particle_const;
			}

			// Conversion factor from world space to voxel space
			double norm_fac = (double)bbox_dim / scale_space_diagonal;

			// Get SPH properties for this particle
			double hsml = get_particle_hsml(pid);  // Smoothing length
			double mass = get_particle_mass(pid);
			double rho = get_particle_rho(pid);    // Density

			double radius = hsml;

			// Use precomputed radius from k-NN search if available (overrides hsml)
			if (cache_manager.radius_particles_per_ptype.size() > 0 && cache_manager.radius_particles_per_ptype[particle_type].size() > 0) {
				radius = cache_manager.radius_particles_per_ptype[particle_type][pid - cache_manager.particles_ptype_offset[particle_type]];
			}

			if (particle_fix_size != 0.0f) {
				radius *= particle_fix_size;
			}

			// Convert radius to voxel units
			double len_to_pix = norm_fac / bbox_size_orig;
			double radiusxyz_max = radius * len_to_pix;

			return radiusxyz_max;
		}

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
		std::string ConvertVDBBase::format_filename(const std::string& pattern, int number) {
			// Create the formatted number as a string
			std::ostringstream formattedNumber;
			if (number < 1000) {
				// Zero-pad numbers to 3 digits for consistent file sorting
				formattedNumber << std::setw(3) << std::setfill('0') << number;
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

}//common