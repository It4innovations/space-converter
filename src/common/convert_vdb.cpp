/*
 * Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
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

#include <iostream>
#include <string>
#include <vector>
#include <float.h>

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

//#	include <nanovdb/NanoVDB.h>
//#	include <nanovdb/util/GridHandle.h>
//#	include <nanovdb/util/HostBuffer.h>

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/NanoToOpenVDB.h>
#else
#	include <nanovdb/tools/NanoToOpenVDB.h>
#endif
#endif

#include "dense_common.h"

#include <mpi.h>
#include <algorithm>

#ifdef WITH_EMBREE
#	include <embree4/rtcore.h>
#endif

#ifdef WITH_CUDAKDTREE
#	include "cudakdtree_tool.h"
#endif

#ifdef WITH_NANOFLANN
#	include "nanoflann_tool.h"
#	include <nanoflann.hpp>
#endif

//#define CONVERT_VDB_TEST
#include "utility/dense_utility.h"

namespace common {
	namespace vdb {	
		
//#ifdef CONVERT_VDB_TEST
//		void gen_particles1(std::vector<float> &particles_xyz, int &calc_radius_neigh)
//		{
//			calc_radius_neigh = 6;
//			
//			// Generate particles in a uniform grid pattern
//			const int grid_size = 32;  // Number of particles in each dimension
//			const float spacing = 1.0f; // Spacing between particles
//
//			// Calculate the total number of particles
//			const int total_particles = grid_size * grid_size * grid_size;
//
//			// Resize the vector to hold all particles
//			particles_xyz.resize(total_particles * 3);
//
//			// Generate particles in a uniform grid
//			#pragma omp parallel for
//			for (int z = 0; z < grid_size; z++) {
//				for (int y = 0; y < grid_size; y++) {
//					for (int x = 0; x < grid_size; x++) {
//						int index = (z * grid_size * grid_size + y * grid_size + x) * 3;
//						
//						// Calculate position with a small offset from origin
//						particles_xyz[index]     = x * spacing - (grid_size - 1) * spacing / 2.0f;
//						particles_xyz[index + 1] = y * spacing - (grid_size - 1) * spacing / 2.0f;
//						particles_xyz[index + 2] = z * spacing - (grid_size - 1) * spacing / 2.0f;
//					}
//				}
//			}
//		}
//
//		void reduce_particles(std::vector<float>& local_radius, std::vector<float>& radius)
//		{
//			int mpi_rank, mpi_size;
//			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//			MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
//
//			int total_particles = 0;
//			int local_count = local_radius.size();
//
//			MPI_Allreduce(&local_count, &total_particles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
//
//			//// Divide the grid along the z-axis across MPI ranks
//			//int particles_per_rank = (total_particles + mpi_size - 1) / mpi_size; // ceiling division
//			//int i_start = mpi_rank * particles_per_rank;
//			//int i_end = std::min((mpi_rank + 1) * particles_per_rank, total_particles);
//			//int local_count = i_end - i_start;
//
//			//if (local_radius_size != local_count)
//			//	printf("local_radius_size %d != local_count %d\n", local_radius_size, local_count);
//
//			// Gather the counts from all processes to know how much to receive
//			std::vector<int> counts(mpi_size);
//			MPI_Gather(&local_count, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
//
//			// Only rank 0 needs to prepare for the reduction
//			std::vector<int> displacements(mpi_size, 0);
//
//			if (mpi_rank == 0) {
//				// Calculate displacements for the gatherv operation
//				for (int i = 1; i < mpi_size; ++i) {
//					displacements[i] = displacements[i-1] + counts[i-1];
//				}
//				
//				// Resize particles_xyz to hold all particles in correct order
//				radius.resize(total_particles);
//			}
//
//			// Gather all particles to rank 0
//			MPI_Gatherv(
//				local_radius.data(), local_count, MPI_FLOAT,
//				radius.data(), counts.data(), displacements.data(), MPI_FLOAT,
//				0, MPI_COMM_WORLD
//			);
//
//			MPI_Barrier(MPI_COMM_WORLD);
//
//			if (mpi_rank == 0) {
//				// Reorder the radius based on counts and displacements
//				// This ensures that each MPI rank will have access to the correctly ordered data
//				std::vector<float> reordered_radius(total_particles);
//
//				// Copy each segment from local_radius to the correct position in reordered_radius
//				size_t offset = 0;
//				for (int i = 0; i < mpi_size; ++i) {
//					std::memcpy(
//						reordered_radius.data() + offset,
//						radius.data() + displacements[mpi_size - i - 1],
//						counts[mpi_size - i - 1] * sizeof(float)
//					);
//
//					offset += counts[mpi_size - i - 1];
//				}
//
//				// Replace the original radius vector with the reordered one
//				radius = std::move(reordered_radius);
//			}
//
//			if (mpi_rank == 0) {
//				printf("Reduced all particles to rank 0: %zu particles\n", radius.size());
//			}			
//		}
//
//		void gen_particles(std::vector<float>& particles_xyz, int& calc_radius_neigh)
//		{
//			int mpi_rank, mpi_size;
//			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//			MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
//
//			calc_radius_neigh = 6;
//
//			// Generate particles in a uniform grid pattern
//			const int grid_size = 32;  // Number of particles in each dimension
//			const float spacing = 1.0f; // Spacing between particles
//
//			// Calculate the total number of particles
//			const int total_particles = grid_size * grid_size * grid_size;
//
//			// Divide the grid along the z-axis across MPI ranks
//			int particles_per_rank = (total_particles + mpi_size - 1) / mpi_size; // ceiling division
//			int i_start = mpi_rank * particles_per_rank;
//			int i_end = std::min((mpi_rank + 1) * particles_per_rank, total_particles);
//
//			// Resize the vector to hold all particles
//			particles_xyz.resize(particles_per_rank * 3);
//
//			// Generate particles in a uniform grid
//#pragma omp parallel for
//			for (int z = 0; z < grid_size; z++) {
//				for (int y = 0; y < grid_size; y++) {
//					for (int x = 0; x < grid_size; x++) {
//						int index = z * grid_size * grid_size + y * grid_size + x;
//
//						if (index < i_start || index >= i_end)
//							continue;
//
//						int local_i = index - i_start;
//
//						// Calculate position with a small offset from origin
//						particles_xyz[local_i * 3 + 0] = x * spacing - (grid_size - 1) * spacing / 2.0f;
//						particles_xyz[local_i * 3 + 1] = y * spacing - (grid_size - 1) * spacing / 2.0f;
//						particles_xyz[local_i * 3 + 2] = z * spacing - (grid_size - 1) * spacing / 2.0f;
//					}
//				}
//			}
//		}
//#endif

		void ConvertVDBBase::read_radius_from_file(std::string& calc_radius_neigh_file)
		{
			int ptype_count = get_num_types(); // Ensure that radius_particles is resized correctly based on the number of particle types
			particles_ptype_offset.resize(ptype_count + 1, 0);

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

				radius_particles_per_ptype.push_back(radius_particles); // Store the radius particles for this type
				rho_particles_per_ptype.push_back(rho_particles); // Store the rho particles for this type

				particles_ptype_offset[ptype + 1] = particles_ptype_offset[ptype] + radius_particles.size();

				file.close();
			}
		}

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
				size_t size = radius_particles_per_ptype[ptype].size();
				if (!file.write(reinterpret_cast<const char*>(&size), sizeof(size_t))) {
					printf("Error writing size to file %s\n", calc_radius_neigh_file_ptype.c_str());
					file.close();
					return;
				}

				// Write the vector data
				if (!file.write(reinterpret_cast<const char*>(radius_particles_per_ptype[ptype].data()), size * sizeof(float))) {
					printf("Error writing data to file %s\n", calc_radius_neigh_file_ptype.c_str());
					file.close();
					return;
				}

				// Write the vector data
				if (!file.write(reinterpret_cast<const char*>(rho_particles_per_ptype[ptype].data()), size * sizeof(float))) {
					printf("Error writing data to file %s\n", calc_radius_neigh_file_ptype.c_str());
					file.close();
					return;
				}

				file.close();
			}
		}

#ifdef WITH_CUDAKDTREE
		void ConvertVDBBase::calculate_radius_by_cudakdtree(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, bool use_cudakdtree_cpu, float maxRadius, common::SpaceData::DenseType &rho_kernel)
		{
			double t_start = omp_get_wtime();

			if (calc_radius_neigh == -1) {
				read_radius_from_file(calc_radius_neigh_file);
			}
			else {
				size_t no_points = get_local_num_particles();

				int ptype_count = get_num_types();
				radius_particles_per_ptype.resize(ptype_count);
				rho_particles_per_ptype.resize(ptype_count);
				particles_ptype_offset.resize(ptype_count + 1, 0);

				int num_threads = omp_get_max_threads();

				for (int ptype = 0; ptype < ptype_count; ptype++) {

					std::vector<float> points;
					std::vector < std::vector<float> > points_thread(num_threads);

					std::vector<float> pmass;
					std::vector < std::vector<float> > pmass_thread(num_threads);

#pragma omp parallel num_threads(num_threads) 
					{
						int tid = omp_get_thread_num();

#pragma omp for
						for (size_t i = 0; i < no_points; i++) {

							if (get_particle_type(i) != ptype)
								continue;

							double Pos[3];
							get_particle_position(i, Pos);

							//points[i * 3 + 0] = Pos[0];
							//points[i * 3 + 1] = Pos[1];
							//points[i * 3 + 2] = Pos[2];

							points_thread[tid].push_back(Pos[0]);
							points_thread[tid].push_back(Pos[1]);
							points_thread[tid].push_back(Pos[2]);

							double mass = get_particle_mass(i);
							pmass_thread[tid].push_back(mass);
						}
//#pragma omp critical
//						{
//							points.insert(points.end(), points_thread.begin(), points_thread.end());
//						}
					}

					// Merge in thread order
					std::vector<float> result;
					for (int t = 0; t < num_threads; ++t) {
						points.insert(points.end(), points_thread[t].begin(), points_thread[t].end());
						pmass.insert(pmass.end(), pmass_thread[t].begin(), pmass_thread[t].end());
					}

					particles_ptype_offset[ptype + 1] = particles_ptype_offset[ptype] + points.size() / 3; // Store the offset for this particle type

					printf("cudakdtree: init: %f\n", omp_get_wtime() - t_start);
//#if CONVERT_VDB_TEST
//					//gen_particles(points, calc_radius_neigh);
//					//no_points = points.size() / 3;
//					std::vector<float> local_radius;
//					utility::cudakdtree::run_knn(points.data(), points.size(), calc_radius_neigh + 1, local_radius, false, false);
//					reduce_particles(local_radius, radius_particles_per_ptype[ptype]);
//#else

					//if (points.size() > 0) //TODO wrong for MPI
					//{
						utility::cudakdtree::run_knn(points.data(), points.size() / 3, calc_radius_neigh + 1, radius_particles_per_ptype[ptype], rho_particles_per_ptype[ptype], pmass, !use_cudakdtree_cpu, use_cycling, maxRadius, rho_kernel); // TODO: calc_radius_neigh + 1 ?
					//}
//#endif
				}

				if (calc_radius_neigh_file.length() > 0) {
					write_radius_from_file(calc_radius_neigh_file);
				}
			}

#if 0
			{
				double t_cic3d = omp_get_wtime();

				std::vector<std::vector<double>> Pos; // Matrix (3xNpart) with particle positions
				std::vector<double> HSML; // Array with particle hsml
				std::vector<double> M; // Array with particle masses
				std::vector<double> Rho; // Array with particle densities
				std::vector<double> Bin_Q; // Array with particle quantity to be mapped
				std::vector<double> Weights; // Array with weights. Defaults to density-weighted

				utility::dense::gadget::GadgetPhysical GU = utility::dense::gadget::create_gadget_physical(redshift, hubble_param);

				// test - reset values
				//GU.x_physical = 1;
				//GU.m_physical = 1;
				//GU.rho_cgs = 1;
				//GU.x_cgs = 1;

				int rvir = 15597;//100; // 2 * 500;
				double max_size = rvir * GU.x_physical;

				//int boxsize = 1000;
				int Npixels[3] = { 100, 100, 100 };
				double pixelSideLength = max_size / Npixels[0];

				//std::vector<int> npix;
				//npix.push_back(Npixels); // Number of pixels in x direction
				//npix.push_back(Npixels); // Number of pixels in y direction
				//npix.push_back(Npixels); // Number of pixels in z direction

				double len2pix = 1.0 / pixelSideLength;

				//utility::dense::cic_3d::MappingParameters param(npix, len2pix);
				utility::dense::sph_kernel::WendlandC6 kernel;
				bool show_progress = true;
				bool calc_mean = false;


				//////////////////////FILL DATA//////////////////////////////////				
				size_t no_points = get_local_num_particles();

				int ptype = 1;
				size_t no_points_ptype = radius_particles_per_ptype[ptype].size();

				//Pos.resize(no_points_ptype);
				//HSML.resize(no_points_ptype);
				//M.resize(no_points_ptype);
				//Rho.resize(no_points_ptype);
				//Bin_Q.resize(no_points_ptype);
				//Weights.resize(no_points_ptype);

				size_t id = 0;

				for (size_t i = 0; i < no_points; i++) {

					if (get_particle_type(i) != ptype)
						continue;

					double p[3];
					get_particle_position(i, p);

					p[0] *= GU.x_physical;
					p[1] *= GU.x_physical;
					p[2] *= GU.x_physical;

					//if (p[0] < -max_size || p[0] > max_size ||
					//	p[1] < -max_size || p[1] > max_size ||
					//	p[2] < -max_size || p[2] > max_size) {

					//	id++;
					//	continue; // Skip particles outside the bounding box
					//}

					Pos.push_back({ p[0], p[1], p[2] });

					HSML.push_back(radius_particles_per_ptype[ptype][id]);

					M.push_back(get_particle_mass(i) * GU.m_physical);

					Rho.push_back(rho_particles_per_ptype[ptype][id]); // Value
					Bin_Q.push_back(rho_particles_per_ptype[ptype][id] * GU.rho_cgs); // Quantity == Rho

					Weights.push_back(utility::dense::sph_kernel::part_weight_physical(pixelSideLength, GU.x_cgs));

					id++;
				}
				////////////////////////////////////////////////////////

				auto image3d = utility::dense::cic_3d::cic_mapping_3D(
					Pos, HSML, M, Rho, Bin_Q, Weights,
					Npixels, len2pix, kernel,
					show_progress,
					calc_mean
				);

				printf("cic_mapping_3D finished: %f\n", omp_get_wtime() - t_cic3d);
			}
#endif

			printf("cudakdtree: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
		}

#endif		

//#if 0
//
//		struct Particle {
//			float x, y, z;
//		};
//
//		// PointCloud Adaptor for nanoflann
//		struct PointCloud {
//			std::vector<Particle> pts;
//
//			// Must return the number of data points
//			inline size_t kdtree_get_point_count() const { return pts.size(); }
//
//			// Returns the dim'th component of the idx'th point in the class
//			inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
//				if (dim == 0) return pts[idx].x;
//				else if (dim == 1) return pts[idx].y;
//				else return pts[idx].z;
//			}
//
//			// Optional bounding-box computation
//			template <class BBOX>
//			bool kdtree_get_bbox(BBOX&) const { return false; }
//		};
//
//		using KDTreeType = nanoflann::KDTreeSingleIndexAdaptor<
//			nanoflann::L2_Simple_Adaptor<float, PointCloud>,
//			PointCloud,
//			3 // dimensionality
//		>;
//
//		void ConvertVDBBase::calculate_radius_by_nanoflann(int calc_radius_neigh, std::string& calc_radius_neigh_file)
//		{
//			double t_start = omp_get_wtime();
//
//			if (calc_radius_neigh == -1) {
//				read_radius_from_file(calc_radius_neigh_file);
//			}
//			else {
//				size_t no_points = get_local_num_particles();
//
//				int ptype_count = get_num_types();
//				radius_particles_per_ptype.resize(ptype_count);
//				rho_particles_per_ptype.resize(ptype_count);
//
//				for (int ptype = 0; ptype < ptype_count; ptype++) {
//
//					PointCloud cloud;
//					//cloud.pts.resize(no_points);
//
//#pragma omp parallel
//					{
//						std::vector<Particle> points_thread;
//#pragma omp	for
//						for (size_t i = 0; i < no_points; i++) {
//							double Pos[3];
//							get_particle_position(i, Pos);
//
//							//cloud.pts[i].x = Pos[0];
//							//cloud.pts[i].y = Pos[1];
//							//cloud.pts[i].z = Pos[2];
//
//							Particle pts;
//							pts.x = Pos[0];
//							pts.y = Pos[1];
//							pts.z = Pos[2];
//
//							points_thread.push_back(pts);
//						}
//
//#pragma omp critical
//						{
//							cloud.pts.insert(cloud.pts.end(), points_thread.begin(), points_thread.end());
//						}
//					}
//
//					printf("nanoflann: init: %f\n", omp_get_wtime() - t_start);
////#if CONVERT_VDB_TEST				
////					std::vector<float> points;
////					gen_particles(points, calc_radius_neigh);
////					no_points = points.size() / 3;
////					cloud.pts.resize(points.size() / 3);
////					memcpy(cloud.pts.data(), points.data(), points.size() * sizeof(float));
////#endif
//
//					// Build the KD-Tree
//					KDTreeType kdtree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
//					kdtree.buildIndex();
//
//					printf("nanoflann: Build the KD-Tree: %f\n", omp_get_wtime() - t_start);
//
//					radius_particles_per_ptype[ptype].resize(cloud.pts.size());
//					rho_particles_per_ptype[ptype].resize(cloud.pts.size());
//
//					// Find nearest neighbors for each particle
//#pragma omp parallel for
//					for (size_t i = 0; i < cloud.pts.size(); ++i) {
//						float query_pt[3] = { cloud.pts[i].x, cloud.pts[i].y, cloud.pts[i].z };
//
//						std::vector<size_t> ret_indexes(calc_radius_neigh);
//						std::vector<float> out_dists_sqr(calc_radius_neigh);
//
//						nanoflann::KNNResultSet<float> resultSet(calc_radius_neigh);
//						resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
//
//						kdtree.findNeighbors(resultSet, query_pt);
//
//						float max_squared_dist = 0.0f;
//						for (float d : out_dists_sqr) {
//							if (d > max_squared_dist)
//								max_squared_dist = d;
//						}
//
//						float max_distance = std::sqrt(max_squared_dist);
//						radius_particles_per_ptype[ptype][i] = max_distance;
//					}
//				}
//
//				if (calc_radius_neigh_file.length() > 0) {
//					write_radius_from_file(calc_radius_neigh_file);
//				}
//			}
//
//			printf("nanoflann: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
//		}
//
//#endif

#ifdef WITH_NANOFLANN
		void ConvertVDBBase::calculate_radius_by_nanoflann(int calc_radius_neigh, std::string& calc_radius_neigh_file, bool use_cycling, common::SpaceData::DenseType& rho_kernel)
		{
			double t_start = omp_get_wtime();

			if (calc_radius_neigh == -1) {
				read_radius_from_file(calc_radius_neigh_file);
			}
			else {
				size_t no_points = get_local_num_particles();

				int ptype_count = get_num_types();
				radius_particles_per_ptype.resize(ptype_count);
				rho_particles_per_ptype.resize(ptype_count);
				particles_ptype_offset.resize(ptype_count + 1, 0);

				int num_threads = omp_get_max_threads();

				for (int ptype = 0; ptype < ptype_count; ptype++) {

					std::vector<float> points;
					std::vector < std::vector<float> > points_thread(num_threads);

					std::vector<float> pmass;
					std::vector < std::vector<float> > pmass_thread(num_threads);

#pragma omp parallel num_threads(num_threads) 
					{
						int tid = omp_get_thread_num();

#pragma omp for
						for (size_t i = 0; i < no_points; i++) {

							if (get_particle_type(i) != ptype)
								continue;

							double Pos[3];
							get_particle_position(i, Pos);

							//points[i * 3 + 0] = Pos[0];
							//points[i * 3 + 1] = Pos[1];
							//points[i * 3 + 2] = Pos[2];

							points_thread[tid].push_back(Pos[0]);
							points_thread[tid].push_back(Pos[1]);
							points_thread[tid].push_back(Pos[2]);

							double mass = get_particle_mass(i);
							pmass_thread[tid].push_back(mass);
						}
						//#pragma omp critical
						//						{
						//							points.insert(points.end(), points_thread.begin(), points_thread.end());
						//						}
					}

					// Merge in thread order
					std::vector<float> result;
					for (int t = 0; t < num_threads; ++t) {
						points.insert(points.end(), points_thread[t].begin(), points_thread[t].end());
						pmass.insert(pmass.end(), pmass_thread[t].begin(), pmass_thread[t].end());
					}

					particles_ptype_offset[ptype + 1] = particles_ptype_offset[ptype] + points.size() / 3; // Store the offset for this particle type

					printf("nanoflann: init: %f\n", omp_get_wtime() - t_start);

					if (points.size() > 0) {
						utility::nanoflann_tool::run_knn(points.data(), points.size() / 3, calc_radius_neigh + 1, radius_particles_per_ptype[ptype], rho_particles_per_ptype[ptype], pmass, use_cycling, rho_kernel); // TODO: calc_radius_neigh + 1 ?
					}
				}

				if (calc_radius_neigh_file.length() > 0) {
					write_radius_from_file(calc_radius_neigh_file);
				}
			}

#if 0
			{
				double t_cic3d = omp_get_wtime();

				std::vector<std::vector<double>> Pos; // Matrix (3xNpart) with particle positions
				std::vector<double> HSML; // Array with particle hsml
				std::vector<double> M; // Array with particle masses
				std::vector<double> Rho; // Array with particle densities
				std::vector<double> Bin_Q; // Array with particle quantity to be mapped
				std::vector<double> Weights; // Array with weights. Defaults to density-weighted

				utility::dense::gadget::GadgetPhysical GU = utility::dense::gadget::create_gadget_physical(redshift, hubble_param);

				int rvir = 100; // 2 * 500;
				double max_size = rvir * GU.x_physical;

				//int boxsize = 1000;
				int Npixels = 100;
				double pixelSideLength = max_size / Npixels;

				std::vector<int> npix;
				npix.push_back(Npixels); // Number of pixels in x direction
				npix.push_back(Npixels); // Number of pixels in y direction
				npix.push_back(Npixels); // Number of pixels in z direction

				double len2pix = 1.0 / pixelSideLength;

				utility::dense::cic_3d::MappingParameters param(npix, len2pix);
				utility::dense::sph_kernel::WendlandC6 kernel;
				bool show_progress = true;
				bool calc_mean = false;


				//////////////////////FILL DATA//////////////////////////////////				
				size_t no_points = get_local_num_particles();

				int ptype = 1;
				size_t no_points_ptype = radius_particles_per_ptype[ptype].size();

				//Pos.resize(no_points_ptype);
				//HSML.resize(no_points_ptype);
				//M.resize(no_points_ptype);
				//Rho.resize(no_points_ptype);
				//Bin_Q.resize(no_points_ptype);
				//Weights.resize(no_points_ptype);

				size_t id = 0;

				for (size_t i = 0; i < no_points; i++) {

					if (get_particle_type(i) != ptype)
						continue;

					double p[3];
					get_particle_position(i, p);

					p[0] *= GU.x_physical;
					p[1] *= GU.x_physical;
					p[2] *= GU.x_physical;

					//if (p[0] < -max_size || p[0] > max_size ||
					//	p[1] < -max_size || p[1] > max_size ||
					//	p[2] < -max_size || p[2] > max_size) {

					//	id++;
					//	continue; // Skip particles outside the bounding box
					//}

					Pos.push_back({ p[0], p[1], p[2] });

					HSML.push_back(radius_particles_per_ptype[ptype][id]);

					M.push_back(get_particle_mass(i) * GU.m_physical);

					Rho.push_back(rho_particles_per_ptype[ptype][id]); // Value
					Bin_Q.push_back(rho_particles_per_ptype[ptype][id] * GU.rho_cgs); // Quantity == Rho

					Weights.push_back(utility::dense::sph_kernel::part_weight_physical(pixelSideLength, GU.x_cgs));

					id++;
				}
				////////////////////////////////////////////////////////

				auto image3d = utility::dense::cic_3d::cic_mapping_3D(
					Pos, HSML, M, Rho, Bin_Q, Weights,
					param, kernel,
					show_progress,
					calc_mean
				);

				printf("cic_mapping_3D finished: %f\n", omp_get_wtime() - t_cic3d);
			}
#endif

			printf("nanoflann: Find nearest neighbors: %f\n", omp_get_wtime() - t_start);
		}
#endif

#ifdef WITH_EMBREE
		void* alignedMalloc(size_t size, size_t align)
		{
			if (size == 0)
				return nullptr;

			assert((align & (align - 1)) == 0);
			void* ptr = _mm_malloc(size, align);
			if (size != 0 && ptr == nullptr)
				throw std::bad_alloc();
			return ptr;
		}

		void alignedFree(void* ptr)
		{
			if (ptr)
				_mm_free(ptr);
		}

		void ConvertVDBBase::create_embree_scene(int particle_type)
		{
			///
			size_t no_points = get_local_num_particles();
			nanovdb::Vec4f* point_vertices = (nanovdb::Vec4f*)alignedMalloc(no_points * sizeof(nanovdb::Vec4f), 16);

#pragma omp parallel for
			for (size_t i = 0; i < no_points; i++) {
				if (get_particle_type(i) != particle_type)
					continue;

				double Pos[3];
				get_particle_position(i, Pos);
				//if (offset_position[0] != 0.0f || offset_position[1] != 0.0f || offset_position[2] != 0.0f) {
				//	Pos[0] -= offset_position[0];
				//	Pos[1] -= offset_position[1];
				//	Pos[2] -= offset_position[2];
				//}

				double radius = get_particle_hsml(i) * 0.75f; //todo - factor

				point_vertices[i][0] = Pos[0];
				point_vertices[i][1] = Pos[1];
				point_vertices[i][2] = Pos[2];
				point_vertices[i][3] = radius;
			}

			rtc_device = rtcNewDevice("");
			rtc_scene = rtcNewScene((RTCDevice)rtc_device);

			rtcSetSceneBuildQuality((RTCScene)rtc_scene, RTC_BUILD_QUALITY_LOW);
			rtcSetSceneFlags((RTCScene)rtc_scene, RTC_SCENE_FLAG_DYNAMIC);

			RTCGeometry geom = rtcNewGeometry((RTCDevice)rtc_device, RTC_GEOMETRY_TYPE_USER);

			auto boundsFunc = [](const struct RTCBoundsFunctionArguments* args)
				{
					const nanovdb::Vec4f* spheres = (const nanovdb::Vec4f*)args->geometryUserPtr;
					RTCBounds* bounds_o = args->bounds_o;
					const nanovdb::Vec4f& sphere = spheres[args->primID];
					bounds_o->lower_x = sphere[0] - sphere[3];
					bounds_o->lower_y = sphere[1] - sphere[3];
					bounds_o->lower_z = sphere[2] - sphere[3];
					bounds_o->upper_x = sphere[0] + sphere[3];
					bounds_o->upper_y = sphere[1] + sphere[3];
					bounds_o->upper_z = sphere[2] + sphere[3];
				};
			auto intersectFunc = [](const RTCIntersectFunctionNArguments* args) {};
			auto occludedFunc = [](const RTCOccludedFunctionNArguments* args) {};

			rtcSetGeometryUserPrimitiveCount(geom, no_points);
			rtcSetGeometryUserData(geom, point_vertices);
			rtcSetGeometryBoundsFunction(geom, boundsFunc, nullptr);
			rtcSetGeometryIntersectFunction(geom, intersectFunc);
			rtcSetGeometryOccludedFunction(geom, occludedFunc);

			rtcCommitGeometry(geom);
			rtcAttachGeometry((RTCScene)rtc_scene, geom);
			rtcReleaseGeometry(geom);

			rtcCommitScene((RTCScene)rtc_scene);
		}
#endif


		///////////////////////	
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

			double orig_box_space_diagonal = object_size * sqrt(3.0);
			double new_box_space_diagonal = sqrt(pow(bbox_max[0] - bbox_min[0], 2) + pow(bbox_max[1] - bbox_min[1], 2) + pow(bbox_max[2] - bbox_min[2], 2));

			double scale_space_diagonal = new_box_space_diagonal / orig_box_space_diagonal;

			transform_scale = scale_space_diagonal * object_size / (double)bbox_dim;

			size_t no_points = get_local_num_particles();

			//RawParticles
			RawParticles::ParticleData raw_positions;
			RawParticles::ParticleData raw_values;
			RawParticles::ParticleData raw_radius;
			RawParticles::ParticleData raw_frame;

			if (grid.type == VDBParticles::VDBParticleType::eNanoVDB) {
				grid.nano_grid->setTransform(transform_scale);
				//auto acc = grid.nano_grid.getAccessor();
			}
#ifdef WITH_OPENVDB
			else if (grid.type == VDBParticles::VDBParticleType::eOpenVDB) {
				// Create a transform object with the scale factor
				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
				grid.vdb_grid->setTransform(transform);
				//auto acc = grid->getAccessor();
			}
#endif
			else if (grid.type == VDBParticles::VDBParticleType::eRawParticles) {
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

			float bbox_x_min_norm = bbox_min[0] / (float)object_size;
			float bbox_x_max_norm = bbox_max[0] / (float)object_size;
			float bbox_y_min_norm = bbox_min[1] / (float)object_size;
			float bbox_y_max_norm = bbox_max[1] / (float)object_size;
			float bbox_z_min_norm = bbox_min[2] / (float)object_size;
			float bbox_z_max_norm = bbox_max[2] / (float)object_size;

			if (grid.type == VDBParticles::VDBParticleType::eDense) {
				//transform_scale = object_size / (double)bbox_dim;
				grid.dense_grid.offset[0] = (size_t)(bbox_x_min_norm * (double)bbox_dim / scale_space_diagonal);
				grid.dense_grid.offset[1] = (size_t)(bbox_y_min_norm * (double)bbox_dim / scale_space_diagonal);
				grid.dense_grid.offset[2] = (size_t)(bbox_z_min_norm * (double)bbox_dim / scale_space_diagonal);

				//scale_space_diagonal = 1.0;				
			}

			float min = FLT_MAX;
			float max = -FLT_MAX;

			size_t particles_count_temp = 0;

#ifdef WITH_OPENMP			
			int nthreads = 1;
			if (grid.type == VDBParticles::VDBParticleType::eDense) {
				nthreads = omp_get_max_threads();
			}

#pragma omp parallel for reduction(min : min) reduction(max : max) reduction(+ : particles_count_temp) num_threads(nthreads) schedule(dynamic, 256) //TODO
#endif
			for (size_t i = 0; i < no_points; i++) {
				if (get_particle_type(i) != particle_type)
					continue;

				if (anim_type == common::SpaceData::AnimType::eFrameExtract) { // ONLY SELECTED FRAME
					if (frame_req != frame)
						continue;
				}

				double Pos[3];
				get_particle_position(i, Pos);

				if (offset_position[0] != 0.0f || offset_position[1] != 0.0f || offset_position[2] != 0.0f) {
					Pos[0] -= offset_position[0];
					Pos[1] -= offset_position[1];
					Pos[2] -= offset_position[2];
				}

				double px_norm = ((double)Pos[0] - (double)bbox_min_orig[0]) / bbox_size_orig;
				double py_norm = ((double)Pos[1] - (double)bbox_min_orig[1]) / bbox_size_orig;
				double pz_norm = ((double)Pos[2] - (double)bbox_min_orig[2]) / bbox_size_orig;

				// is in the bbox?
				if (grid.type == VDBParticles::VDBParticleType::eDense) {
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
					if (px_norm < bbox_x_min_norm || px_norm > bbox_x_max_norm)
						continue;

					if (py_norm < bbox_y_min_norm || py_norm > bbox_y_max_norm)
						continue;

					if (pz_norm < bbox_z_min_norm || pz_norm > bbox_z_max_norm)
						continue;
				}
				////////////////////

				int px = static_cast<int>((double)px_norm * (double)bbox_dim / scale_space_diagonal);
				int py = static_cast<int>((double)py_norm * (double)bbox_dim / scale_space_diagonal);
				int pz = static_cast<int>((double)pz_norm * (double)bbox_dim / scale_space_diagonal);

				float v_orig = get_particle_norm_value(block_name_id, i);
				if (use_simple_density) {
					v_orig = 1;
				}

				if (v_orig < filter_min || v_orig > filter_max)
					continue;

				/////////////////////////
				if(bbox_sphere_r > 0.0f)
				{
				    // Compute squared distance from the point to the sphere center
					double distSquared = 
					(Pos[0] - bbox_sphere_pos[0]) * (Pos[0] - bbox_sphere_pos[0]) +
					(Pos[1] - bbox_sphere_pos[1]) * (Pos[1] - bbox_sphere_pos[1]) +
					(Pos[2] - bbox_sphere_pos[2]) * (Pos[2] - bbox_sphere_pos[2]);

					// Compare with the squared radius
					if (distSquared > (bbox_sphere_r * bbox_sphere_r)) {
						continue;
					}
				}
				/////////////////////////
				particles_count_temp++;

				if (v_orig < min) min = v_orig;
				if (v_orig > max) max = v_orig;

				/////////////////////////

				float v = v_orig;

				if (grid.type == VDBParticles::VDBParticleType::eDense) {
					fill_voxels_v5(grid.dense_grid,
						i, v,
						bbox_dim, bbox_min_orig, bbox_size_orig,
						scale_space_diagonal,
						dense_type, dense_norm,
						particle_fix_size, particle_type, block_name_id, Pos);
				}
				else if (grid.type == VDBParticles::VDBParticleType::eNanoVDB) {
					nanovdb::Coord xyz(px, py, pz);
					auto acc = grid.nano_grid->getAccessor();
					if (acc.isValueOn(xyz)) {
						v += acc.getValue(xyz); //ADD
					}

					acc.setValue(xyz, v);
				}
#ifdef WITH_OPENVDB
				else if (grid.type == VDBParticles::VDBParticleType::eOpenVDB) {
					openvdb::Coord xyz(px, py, pz);
					auto acc = grid.vdb_grid->getAccessor();

					//if (acc.isValueOn(xyz)) {
					//	v += acc.getValue(xyz); //ADD
					//}

					//acc.setValue(xyz, v);
					// Optimization 2: Use probeValue for faster lookup (avoids double tree traversal)
					float dst_value;
					bool dst_exists = acc.probeValue(xyz, dst_value);

					if (dst_exists) {
						acc.setValue(xyz, v + dst_value);
					}
					else {
						acc.setValue(xyz, v);
					}
				}
#endif
				else if (grid.type == VDBParticles::VDBParticleType::eRawParticles) {
#if 0					
					raw_positions.values.push_back(Pos[0]);
					raw_positions.values.push_back(Pos[1]);
					raw_positions.values.push_back(Pos[2]);					
#else
					raw_positions.values.push_back((double)px_norm * (double)object_size);
					raw_positions.values.push_back((double)py_norm * (double)object_size);
					raw_positions.values.push_back((double)pz_norm * (double)object_size);

					//printf("%lld: Pos: %f, %f, %f\n", i, Pos[0], Pos[1], Pos[2]);
					//printf("%lld: px_norm: %f, py_norm: %f, pz_norm: %f, bbox_dim: %d, scale_space_diagonal: %f\n", i, px_norm, py_norm, pz_norm, bbox_dim, scale_space_diagonal);
#endif
					int n_comp = get_particle_value_comp(block_name_id, i);
					std::vector<float> values(n_comp);
					get_particle_value(block_name_id, i, values.data());

					raw_values.num_comp = n_comp;
					raw_values.values.insert(raw_values.values.end(), values.begin(), values.end());

//#if 0					
//					float radiusxyz_max = get_particle_radius(i);
//#else
//					double hsml = 0.0;
//					double rho = 0.0;
//					double mass = 0.0;
//					if (particle_type == 0) {  //only sph
//						hsml = get_particle_hsml(i);
//						rho = get_particle_rho(i);
//						mass = get_particle_mass(i);
//					}
//
//					if (particle_fix_size != 0.0f) {
//						hsml = particle_fix_size * pow((mass / rho), 1.0 / 3.0); // todo
//					}
//
//#if 0					
//					printf("particles_ptype_offset.size(): %zu\n", particles_ptype_offset.size());
//					for (size_t j = 0; j < particles_ptype_offset.size(); j++) {
//						printf("  particles_ptype_offset[%zu]: %zu\n", j, particles_ptype_offset[j]);
//					}
//					printf("radius_particles_per_ptype.size(): %zu\n", radius_particles_per_ptype.size());
//					for (size_t j = 0; j < radius_particles_per_ptype[particle_type].size(); j++) {
//						printf("  radius_particles_per_ptype[%zu][%zu]: %f\n", particle_type, j, radius_particles_per_ptype[particle_type][j]);
//					}
//#endif
//
//					double radius = 2.0 * hsml;
//					if (radius_particles_per_ptype.size() > 0 && radius_particles_per_ptype[particle_type].size() > 0) {
//						radius = radius_particles_per_ptype[particle_type][i - particles_ptype_offset[particle_type]];
//					}
//
//					double radiusx_norm = ((double)radius) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
//					double radiusy_norm = ((double)radius) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
//					double radiusz_norm = ((double)radius) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);
//
//					double radiusx = (double)radiusx_norm * (double)object_size;
//					double radiusy = (double)radiusy_norm * (double)object_size;
//					double radiusz = (double)radiusz_norm * (double)object_size;
//
//					double radiusxyz_max = std::max(radiusx, std::max(radiusy, radiusz));
//#endif					
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
				
					raw_frame.values.push_back(frame);
				}
			}

			particles_count = particles_count_temp;

			if (grid.type == VDBParticles::VDBParticleType::eRawParticles) {
				grid.raw_particles.data.push_back(raw_positions);
				grid.raw_particles.data.push_back(raw_values);
				grid.raw_particles.data.push_back(raw_radius);
				grid.raw_particles.data.push_back(raw_frame);
			}

			min_value = min;
			max_value = max;

#ifdef WITH_OPENMP
			//printf("convert_iolib_to_nvdb: %s, time: %f, %f\n", grid_name.c_str(), omp_get_wtime() - t_step, omp_get_wtime() - t);
			t_step = omp_get_wtime();
#endif

			//fflush(stdout);
		}


		void ConvertVDBBase::merge_grid(
			VDBParticles& grid_dst,
			VDBParticles& grid_recv
		)
		{
			if (grid_dst.type == VDBParticles::VDBParticleType::eDense && grid_recv.type == VDBParticles::VDBParticleType::eDense) {

#pragma omp parallel for
				for (int z = 0; z < grid_dst.dense_grid.dims[2]; z++) {
					for (int y = 0; z < grid_dst.dense_grid.dims[1]; y++) {
						for (int x = 0; z < grid_dst.dense_grid.dims[0]; x++) {
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
				auto acc_dst = grid_dst.nano_grid->getAccessor();
				auto* grid_src_float = (nanovdb::NanoGrid<float>*)grid_recv.vector_grid.data();

				// loop over child nodes of the root node
				for (auto it2 = grid_src_float->tree().root().cbeginChild(); it2; ++it2) {
					// loop over child nodes of the upper internal node
					for (auto it1 = it2->cbeginChild(); it1; ++it1) {
						// loop over child nodes of the lower internal node
						for (auto it0 = it1->cbeginChild(); it0; ++it0) {
							// loop over values
							for (auto it = it0->cbeginValueOn(); it; ++it) {
								float v = *it;

								nanovdb::Coord xyz = it.getCoord();

								if (acc_dst.isValueOn(xyz)) {
									v += acc_dst.getValue(xyz); //ADD
								}
								acc_dst.setValue(xyz, v);
							}
						}
					}
				}
			}
#ifdef WITH_OPENVDB
			else if (grid_dst.type == VDBParticles::VDBParticleType::eOpenVDB && grid_recv.type == VDBParticles::VDBParticleType::eVector) {
#if 0
				auto acc_dst = grid_dst.vdb_grid->getAccessor();
				auto grid_src_float = vector_to_openvdb(grid_recv.vector_grid);

				auto it0 = grid_src_float;
				for (auto it = it0->cbeginValueOn(); it; ++it) {
					float v = *it;

					openvdb::Coord xyz = it.getCoord();

					if (acc_dst.isValueOn(xyz)) {
						v += acc_dst.getValue(xyz); //ADD
					}
					acc_dst.setValue(xyz, v);
				}
#endif
#if 1
				auto grid_src_float = vector_to_openvdb(grid_recv.vector_grid);

				// Optimization 3: Use OpenVDB's composite operation for maximum performance
				openvdb::tools::compSum(*grid_dst.vdb_grid, *grid_src_float);
#endif

#if 0
				auto grid_src_float = vector_to_openvdb(grid_recv.vector_grid);

				// Optimization 4: Use tree combining for efficient parallel processing
				using TreeType = openvdb::FloatTree;
				auto combine_op = [](const float& dst, const float& src, bool /*dst_active*/, bool /*src_active*/) -> float {
					return dst + src;
					};

				grid_dst.vdb_grid->tree().combine2(grid_src_float->tree(), grid_dst.vdb_grid->tree(), combine_op);

#endif
			}
#endif
			else if (grid_dst.type == VDBParticles::VDBParticleType::eRawParticles && grid_recv.type == VDBParticles::VDBParticleType::eVector) {
				RawParticles temp;
				temp.deserialize(grid_recv.vector_grid);
				grid_dst.raw_particles.merge(temp);
			}
		}

		void ConvertVDBBase::iolib_find_bbox(
			int particle_type,
			float* bbox_min,
			float* bbox_max,
			float* offset_position
		)
		{
			size_t no_points = get_local_num_particles();// All.TotNumPart;

			float min_x = FLT_MAX, min_y = FLT_MAX, min_z = FLT_MAX;
			float max_x = -FLT_MAX, max_y = -FLT_MAX, max_z = -FLT_MAX;

#pragma omp parallel for reduction(min : min_x, min_y, min_z) reduction(max : max_x, max_y, max_z)
			for (size_t i = 0; i < no_points; ++i) {

				if (particle_type != -1 && get_particle_type(i) != particle_type)
					continue;

				double Pos[3];
				get_particle_position(i, Pos);

				if (offset_position[0] != 0.0f || offset_position[1] != 0.0f || offset_position[2] != 0.0f) {
					Pos[0] -= offset_position[0];
					Pos[1] -= offset_position[1];
					Pos[2] -= offset_position[2];
				}

				if (Pos[0] < min_x) 
					min_x = Pos[0];
				if (Pos[1] < min_y) 
					min_y = Pos[1];
				if (Pos[2] < min_z) 
					min_z = Pos[2];

				if (Pos[0] > max_x) 
					max_x = Pos[0];
				if (Pos[1] > max_y) 
					max_y = Pos[1];
				if (Pos[2] > max_z) 
					max_z = Pos[2];
			}

			bbox_min[0] = min_x;
			bbox_min[1] = min_y;
			bbox_min[2] = min_z;

			bbox_max[0] = max_x;
			bbox_max[1] = max_y;
			bbox_max[2] = max_z;
		}

		void ConvertVDBBase::iolib_find_minmax(
			int particle_type,
			int block_nr,
			float& v_min,
			float& v_max
		)
		{
			size_t no_points = get_local_num_particles();// All.TotNumPart;

			float min = FLT_MAX;
			float max = -FLT_MAX;

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

			//// Create a transform object with the scale factor
			//openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
			//transform->postTranslate(openvdb::Vec3d(particles.offset[0] * transform_scale, particles.offset[1] * transform_scale, particles.offset[2] * transform_scale));
			//floatgrid->setTransform(transform);

			// Create NanoVDB transform with scale factor and translation

			nano_grid->setTransform(
				transform_scale,
				nanovdb::Vec3d(
					particles.offset[0] * transform_scale,
					particles.offset[1] * transform_scale,
					particles.offset[2] * transform_scale
				)
			);

			// Populate
			//common::vdb::exptable xexp(-20.0);
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

						//if (fabs(temp) > FDATA_EPSILON) {
						//	density = density / temp;
						//}else{
						//	density = 0.0f; // Avoid division by zero
						//}

						if (dense_norm != common::SpaceData::DenseNorm::eNone) {
							density = density / temp;
						}

						// if (particles.type == 5) {
						// 	density = particles.data_density[index];
						// }

						//if (particles.type == 2) {
						//	density = -xexp.expm1(density);
						//}

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
			//nanovdb::GridBuilder<float> builder(0.0f); // background value = 0.0f
			//builder.setGridClass(nanovdb::GridClass::FogVolume);

			auto acc_dst = nano_grid->getAccessor();
			// Fill the builder with your data
			for (int z = 0; z < particles.z(); z++) {
				for (int y = 0; y < particles.y(); y++) {
					for (int x = 0; x < particles.x(); x++) {
						nanovdb::Coord xyz(x, y, z);
						float value = particles.data_density[(size_t)x + (size_t)y * (size_t)particles.x() + (size_t)z * (size_t)particles.x() * (size_t)particles.y()];
						// Only store non-zero values to maintain sparsity
						if (value != 0.0f) {
							acc_dst.setValue(xyz, value);
						}
					}
				}
			}

			//openvdb::math::CoordBBox bbox(openvdb::Coord(0, 0, 0), openvdb::Coord(particles.x() - 1, particles.y() - 1, particles.z() - 1));
			//openvdb::tools::Dense<const float, openvdb::tools::LayoutXYZ> dense(bbox, particles.data_density.data());
			//openvdb::tools::copyFromDense(dense, floatgrid->tree(), 0.0f);

			return nano_grid;
		}
#endif

#if defined(WITH_OPENVDB)
		openvdb::FloatGrid::Ptr ConvertVDBBase::dense_to_openvdb(DenseParticles& particles, double transform_scale, common::SpaceData::DenseType dense_type, common::SpaceData::DenseNorm dense_norm)
		{
			/////convert
			openvdb::FloatGrid::Ptr floatgrid = openvdb::FloatGrid::create(0.0f);

			// Identify the grid as a level set.
			floatgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
			// Name the grid "LevelSetSphere".
			std::string grid_name("density");
			floatgrid->setName(grid_name);


			// Create a transform object with the scale factor
			openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(transform_scale);
			transform->postTranslate(openvdb::Vec3d(particles.offset[0] * transform_scale, particles.offset[1] * transform_scale, particles.offset[2] * transform_scale));
			floatgrid->setTransform(transform);


			// Populate
			//common::vdb::exptable xexp(-20.0);

			// float minvalue = FLT_MAX; 
			// float maxvalue = -FLT_MAX;

#ifndef WITH_NO_DATA_TEMP
#pragma omp parallel for //reduction(min : minvalue) reduction(max : maxvalue)
			for (int z = 0; z < particles.z(); z++) {
				for (int y = 0; y < particles.y(); y++) {
					for (int x = 0; x < particles.x(); x++) {

						// Get the value from the array						
						size_t index = particles.get_index(x, y, z);
						float density = particles.data_density[index];

						float temp = 0.0f;
						temp = particles.data_temp[index];						

						//if (fabs(temp) > FDATA_EPSILON) {
						//	density = density / temp;
						//}else{
						//	density = 0.0f; // Avoid division by zero
						//}

						if (dense_norm != common::SpaceData::DenseNorm::eNone) {
							density = density / temp;
						}

						// if (particles.type == 5) {
						// 	density = particles.data_density[index];
						// }

						//if (particles.type == 2) {
						//	density = -xexp.expm1(density);
						//}

						// If the value is non-zero, set it in the grid
						if (!std::isnan(density)) {
							//accessor.setValue(openvdb::Coord(x + particles.offset[0], y + particles.offset[1], z + particles.offset[2]), density);							

							// if (dense_type == common::SpaceData::DenseType::eType2)
							// 	particles.data_density[index] = std::log10(density);
							// else
							particles.data_density[index] = density;
						}
						else {
							particles.data_density[index] = 0.0f;
						}
// #if 1
// 						minvalue = std::min(minvalue, particles.data_density[index]);
// 						maxvalue = std::max(maxvalue, particles.data_density[index]);
// #endif
					}
				}
			}
#endif

// #if 1
// #pragma omp parallel for
// 			for (int z = 0; z < particles.z(); z++) {
// 				for (int y = 0; y < particles.y(); y++) {
// 					for (int x = 0; x < particles.x(); x++) {

// 						// Get the value from the array						
// 						size_t index = particles.get_index(x, y, z);
// 						particles.data_density[index] = (particles.data_density[index] - minvalue + 1.0f);// / (maxvalue - minvalue);
// 					}
// 				}
// 			}
// #endif

// #if 0
			openvdb::math::CoordBBox bbox(openvdb::Coord(0, 0, 0), openvdb::Coord(particles.x() - 1, particles.y() - 1, particles.z() - 1));
			openvdb::tools::Dense<const float, openvdb::tools::LayoutXYZ> dense(bbox, particles.data_density.data());
			openvdb::tools::copyFromDense(dense, floatgrid->tree(), 0.0f);
// #else
// 			auto acc_dst = floatgrid->getAccessor();
// 			// Fill the builder with your data
// 			for (int z = 0; z < particles.z(); z++) {
// 				for (int y = 0; y < particles.y(); y++) {
// 					for (int x = 0; x < particles.x(); x++) {
// 						openvdb::Coord xyz(x, y, z);
// 						float value = particles.data_density[(size_t)x + (size_t)y * (size_t)particles.x() + (size_t)z * (size_t)particles.x() * (size_t)particles.y()];
// 						// Only store non-zero values to maintain sparsity
// 						if (value != 0.0f) {
// 							acc_dst.setValue(xyz, value);
// 						}
// 					}
// 				}
// 			}
// #endif

			return floatgrid;
		}

		void ConvertVDBBase::openvdb_to_vector(openvdb::FloatGrid::Ptr grid, std::vector<uint8_t>& file_content)
		{
			// Convert the stringstream to a vector<char>			
			std::ostringstream stream(std::ios_base::binary);
			openvdb::io::Stream(stream).write({ grid });
			stream.flush();
			const std::string& str = stream.str();
			file_content.assign(str.begin(), str.end());
		}

		void ConvertVDBBase::openvdb_to_vector2(openvdb::FloatGrid::Ptr grid1, openvdb::FloatGrid::Ptr grid2, std::vector<uint8_t>& file_content)
		{
			// Convert the stringstream to a vector<char>			
			std::ostringstream stream(std::ios_base::binary);
			openvdb::io::Stream(stream).write({ grid1, grid2 });
			stream.flush();
			const std::string& str = stream.str();
			file_content.assign(str.begin(), str.end());
		}

		// Converts a vector<char> to an openvdb::FloatGrid::Ptr
		openvdb::FloatGrid::Ptr ConvertVDBBase::vector_to_openvdb(std::vector<uint8_t>& file_content)
		{
			// Convert the vector<char> to a stringstream
			std::istringstream stream(std::string(file_content.begin(), file_content.end()), std::ios_base::binary);

			// Create a VDB input stream from the stringstream
			openvdb::io::Stream vdbStream(stream);

			// Read the grid from the stream
			openvdb::GridPtrVecPtr grids = vdbStream.getGrids();

			// Find the first FloatGrid in the grids vector and return it
			for (auto& grid : *grids) {
				if (grid->isType<openvdb::FloatGrid>()) {
					return openvdb::gridPtrCast<openvdb::FloatGrid>(grid);
				}
			}

			// If no FloatGrid is found, return a null pointer
			return nullptr;
		}

#endif

		double ConvertVDBBase::get_particle_radius(
			uint64_t pid,
			int bbox_dim,
			int* bbox_min_orig,
			double bbox_size_orig,
			double scale_space_diagonal,
			float particle_fix_size,
			int particle_type
		) {

			double norm_fac = (double)bbox_dim / scale_space_diagonal;

			double hsml = get_particle_hsml(pid);
			double mass = get_particle_mass(pid);
			double rho = get_particle_rho(pid);

			//if (particle_fix_size != 0.0f) {
			//	if (mass != 0.0 && rho != 0)
			//		hsml = particle_fix_size * pow((mass / rho), 1.0 / 3.0); // TODO
			//	else
			//		hsml = particle_fix_size; //TODO: check
			//}

			double radius = hsml; // TODO? 2.0 * ?

			if (radius_particles_per_ptype.size() > 0 && radius_particles_per_ptype[particle_type].size() > 0) {
				radius = radius_particles_per_ptype[particle_type][pid - particles_ptype_offset[particle_type]];
			}

			if (particle_fix_size != 0.0f) {
				radius *= particle_fix_size; // TODO: check
			}

			//double radiusx_norm = ((double)radius) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			//double radiusy_norm = ((double)radius) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			//double radiusz_norm = ((double)radius) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double len_to_pix = norm_fac / bbox_size_orig;

			//double radiusx = (double)radiusx_norm * norm_fac;
			//double radiusy = (double)radiusy_norm * norm_fac;
			//double radiusz = (double)radiusz_norm * norm_fac;

			double radiusxyz_max = radius * len_to_pix;// std::max(radiusx, std::max(radiusy, radiusz));

			return radiusxyz_max;
		}

		double ConvertVDBBase::get_particle_rho(uint64_t id) {
			double mass = get_particle_mass(id);
			if (mass != 0.0) {
				int particle_type = get_particle_type(id);
				if (rho_particles_per_ptype.size() > 0 && rho_particles_per_ptype[particle_type].size() > 0) {
					return rho_particles_per_ptype[particle_type][id - particles_ptype_offset[particle_type]];
				}
			}

			return get_particle_rho_internal(id);
		}

		void ConvertVDBBase::get_types_and_blocks(std::vector<int>& types_and_blocks) {
			get_types_and_blocks_internal(types_and_blocks);

			int num_types = get_num_types();
			int rho_blocknr = get_particle_rho_blocknr();

			for (int type = 0; type < num_types; type++) {
				if (rho_particles_per_ptype.size() > 0 && rho_particles_per_ptype[type].size() > 0) {
					types_and_blocks[num_types * rho_blocknr + type]++;
				}
			}
		}

		float ConvertVDBBase::get_particle_norm_value(int blocknr, uint64_t id) {
			if (blocknr == get_particle_rho_blocknr()) {
				return get_particle_rho(id);
			}

			return get_particle_norm_value_internal(blocknr, id);
		}

		int ConvertVDBBase::get_particle_value(int blocknr, uint64_t id, float* out_value) {
			if (blocknr == get_particle_rho_blocknr()) {
				*out_value = get_particle_rho(id);
				return 1; // 1 component
			}

			return get_particle_value_internal(blocknr, id, out_value);
		}
		int ConvertVDBBase::get_particle_value_comp(int blocknr, uint64_t id) {
			return get_particle_value_comp_internal(blocknr, id);
		}

	}//vdb

}//common