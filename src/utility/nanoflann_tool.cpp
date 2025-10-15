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

#include "nanoflann_tool.h"

#include <mpi.h>
#include <nanoflann.hpp>
#include <iostream>

//#undef WITH_OPENMP

#ifdef WITH_OPENMP
#	include <omp.h>
#endif

#define NANOFLANN_MPI_CALL(fctCall)                                          \
  { int rc = MPI_##fctCall;                                             \
    if (rc != MPI_SUCCESS) {                                             \
       std::cerr << (std::string(__FUNCTION__) + " " + #fctCall) << std::endl; }}

#include <climits>  // For INT_MAX
// Use a macro for max bytes per transfer (less than INT_MAX to be safe)
#define MAX_MPI_BYTES (INT_MAX - 1024)	   

#include "utility/dense_utility.h"

namespace utility {
	namespace nanoflann_tool {

		// Create point cloud adaptor for nanoflann
		struct PointCloud {
			std::vector<float> pts;
			size_t N = 0;

			// Must return the number of data points
			inline size_t kdtree_get_point_count() const { return pts.size() / 3; }

			// Returns the dim'th component of the idx'th point in the class
			inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
				return pts[idx * 3 + dim];
			}

			// Optional bounding-box computation
			template <class BBOX>
			bool kdtree_get_bbox(BBOX&) const { return false; }

			void resize(size_t n) {
				pts.resize(n * 3); // Each point has 3 dimensions (x, y, z)
			};

			void copy_from(float* data, size_t count) {
				if (count > pts.size() / 3) {
					std::cerr << "Error: Not enough space in PointCloud to copy data." << std::endl;
					return;
				}

				N = count;
				memcpy(pts.data(), data, N * sizeof(float) * 3);
			}
		};

        void mpi_cycling(int recvPeer, char* recvPtr, size_t totalBytesRecv, int sendPeer, char* sendPtr, size_t totalBytesSend) 
        {
            // char* recvPtr = reinterpret_cast<char*>(d_tree_recv);
            // char* sendPtr = reinterpret_cast<char*>(d_tree);
            // size_t totalBytesRecv = recvCount * sizeof(*d_tree);
            // size_t totalBytesSend = sendCount * sizeof(*d_tree);

            MPI_Request requests[2];
            size_t offsetRecv = 0;
            size_t offsetSend = 0;

            while (offsetRecv < totalBytesRecv || offsetSend < totalBytesSend) {
                int chunkRecv = static_cast<int>(std::min(totalBytesRecv - offsetRecv, static_cast<size_t>(MAX_MPI_BYTES)));
                int chunkSend = static_cast<int>(std::min(totalBytesSend - offsetSend, static_cast<size_t>(MAX_MPI_BYTES)));

                if (offsetRecv < totalBytesRecv) {
                    NANOFLANN_MPI_CALL(Irecv(recvPtr + offsetRecv, chunkRecv, MPI_BYTE, recvPeer, 0,
                                        MPI_COMM_WORLD, &requests[0]));
                } else {
                    requests[0] = MPI_REQUEST_NULL;
                }

                if (offsetSend < totalBytesSend) {
                    NANOFLANN_MPI_CALL(Isend(sendPtr + offsetSend, chunkSend, MPI_BYTE, sendPeer, 0,
                                        MPI_COMM_WORLD, &requests[1]));
                } else {
                    requests[1] = MPI_REQUEST_NULL;
                }

                NANOFLANN_MPI_CALL(Waitall(2, requests, MPI_STATUSES_IGNORE));

                offsetRecv += chunkRecv;
                offsetSend += chunkSend;
            }

        }		

		void run_knn_cpu(float* points, size_t numPointsThatIHave, int k, std::vector<float>& radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_cycling, common::SpaceData::DenseType& rho_kernel)
		{
			int mpi_rank, mpi_size;
			MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
			MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

			// -----------------------------------------------------------------------------
			// find out max num points anybody has, so we can allocate
			// -----------------------------------------------------------------------------
			int N = numPointsThatIHave;
			int maxNumPointsAnybodyHas = 0;

			NANOFLANN_MPI_CALL(Allreduce(&N, &maxNumPointsAnybodyHas, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));

			// Create point clouds for tree data and received data
			PointCloud tree_cloud, tree_recv_cloud;
			tree_cloud.resize(maxNumPointsAnybodyHas + 1);
			tree_recv_cloud.resize(maxNumPointsAnybodyHas + 1);

			// Copy input points to tree cloud with OpenMP
//#ifdef WITH_OPENMP
//#pragma omp parallel for schedule(static)
//#endif
//			for (int i = 0; i < N; i++) {
//				tree_cloud.pts[i][0] = points[i * 3 + 0];
//				tree_cloud.pts[i][1] = points[i * 3 + 1];
//				tree_cloud.pts[i][2] = points[i * 3 + 2];
//			}
//			tree_cloud.pts.resize(N);
			tree_cloud.copy_from(points, N);

			// Define KD-tree type
			using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
				nanoflann::L2_Simple_Adaptor<float, PointCloud>,
				PointCloud,
				3 /* dimensionality */
			>;

			// Add timing
			double start_time, end_time;
			MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

			// Build initial KD-tree with OpenMP support
			std::unique_ptr<KDTree> index = std::make_unique<KDTree>(
				3 /* dimensionality */, tree_cloud,
				nanoflann::KDTreeSingleIndexAdaptorParams(
					10 /* max leaf */,
					nanoflann::KDTreeSingleIndexAdaptorFlags::None,
#ifdef WITH_OPENMP
					omp_get_max_threads() /* n_thread_build */
#else
					1 /* n_thread_build */
#endif
				)
			);

			MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

			if (mpi_rank == 0) {
				printf("Total execution time (buildTree nanoflann): %.6f seconds\n", end_time - start_time);
			}

			// Create query points with OpenMP
			size_t numQueries = N;
			PointCloud queries_cloud;
			queries_cloud.resize(numQueries);

//#ifdef WITH_OPENMP
//#pragma omp parallel for schedule(static)
//#endif
//			for (size_t i = 0; i < numQueries; i++) {
//				queries_cloud.pts[i][0] = points[i * 3 + 0];
//				queries_cloud.pts[i][1] = points[i * 3 + 1];
//				queries_cloud.pts[i][2] = points[i * 3 + 2];
//			}
			queries_cloud.copy_from(points, numQueries);

			// Storage for candidates
			std::vector<std::vector<std::pair<size_t, float>>> all_candidates(numQueries);
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (size_t i = 0; i < numQueries; i++) {
				all_candidates[i].reserve(k * mpi_size); // Reserve space for all rounds
			}

			// -----------------------------------------------------------------------------
			// now, do the queries and cycling:
			// -----------------------------------------------------------------------------
			MPI_Barrier(MPI_COMM_WORLD);
			start_time = MPI_Wtime();

			int round_size = mpi_size;
			if (!use_cycling) {
				round_size = 1;
			}

			for (int round = 0; round < round_size; round++) {
				if (mpi_rank == 0) {
					printf("Starting round (nanoflann cycling) %d\n", round);
				}

				if (round == 0) {
					// Use the initial tree
				}
				else {
					MPI_Request requests[2];
					int sendCount = N;
					int recvCount = 0;
					int sendPeer = (mpi_rank + 1) % mpi_size;
					int recvPeer = (mpi_rank + mpi_size - 1) % mpi_size;

					// Exchange counts
					NANOFLANN_MPI_CALL(Irecv(&recvCount, sizeof(int), MPI_BYTE, recvPeer, 0,
						MPI_COMM_WORLD, &requests[0]));
					NANOFLANN_MPI_CALL(Isend(&sendCount, sizeof(int), MPI_BYTE, sendPeer, 0,
						MPI_COMM_WORLD, &requests[1]));
					NANOFLANN_MPI_CALL(Waitall(2, requests, MPI_STATUSES_IGNORE));

					// Prepare receive buffer
					//tree_recv_cloud.pts.resize(recvCount);
					tree_recv_cloud.N = recvCount;

#if 0					
					// Exchange point data
					NANOFLANN_MPI_CALL(Irecv(tree_recv_cloud.pts.data(), recvCount * sizeof(std::array<float, 3>),
						MPI_BYTE, recvPeer, 0, MPI_COMM_WORLD, &requests[0]));
					NANOFLANN_MPI_CALL(Isend(tree_cloud.pts.data(), sendCount * sizeof(std::array<float, 3>),
						MPI_BYTE, sendPeer, 0, MPI_COMM_WORLD, &requests[1]));
					NANOFLANN_MPI_CALL(Waitall(2, requests, MPI_STATUSES_IGNORE));
#else
                    char* recvPtr = reinterpret_cast<char*>(tree_recv_cloud.pts.data());
                    char* sendPtr = reinterpret_cast<char*>(tree_cloud.pts.data());
                    size_t totalBytesRecv = recvCount * sizeof(std::array<float, 3>);
                    size_t totalBytesSend = sendCount * sizeof(std::array<float, 3>);
                    mpi_cycling(recvPeer, recvPtr, totalBytesRecv, sendPeer, sendPtr, totalBytesSend);
#endif 					

					N = recvCount;
					std::swap(tree_cloud, tree_recv_cloud);

					// Rebuild KD-tree with new data
					index = std::make_unique<KDTree>(
						3 /* dimensionality */, tree_cloud,
						nanoflann::KDTreeSingleIndexAdaptorParams(
							10 /* max leaf */,
							nanoflann::KDTreeSingleIndexAdaptorFlags::None,
#ifdef WITH_OPENMP
							omp_get_max_threads() /* n_thread_build */
#else
							1 /* n_thread_build */
#endif
						)
					);
				}

				// Run queries for this round with OpenMP
				//std::vector<std::vector<size_t>> indices_per_thread(omp_get_max_threads());
				//std::vector<std::vector<float>> distances_per_thread(omp_get_max_threads());

#ifdef WITH_OPENMP
#pragma omp parallel
#endif
				{
#ifdef WITH_OPENMP
					int thread_id = omp_get_thread_num();
#else
					int thread_id = 0;
#endif
					std::vector<uint32_t> indices(k);
					std::vector<float> distances(k);

#ifdef WITH_OPENMP
#pragma omp for //schedule(dynamic, 100) //TODO
#endif
					for (size_t q = 0; q < numQueries; q++) {
						const float query_pt[3] = {
							queries_cloud.pts[q * 3 + 0],
							queries_cloud.pts[q * 3 + 1],
							queries_cloud.pts[q * 3 + 2]
						};

						size_t num_results = index->knnSearch(
							&query_pt[0], k, indices.data(), distances.data() 
						);

						// Store results with global indices adjusted for current round
						for (size_t i = 0; i < num_results; i++) {
							size_t global_idx = indices[i] + round * maxNumPointsAnybodyHas;
#ifdef WITH_OPENMP
#pragma omp critical
#endif
							{
								all_candidates[q].emplace_back(global_idx, distances[i]);
							}
						}
					}
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
			end_time = MPI_Wtime();

			if (mpi_rank == 0) {
				printf("Total execution time (queries and cycling are done): %.6f seconds\n", end_time - start_time);
			}

			// Extract final results - get k closest from all rounds with OpenMP
			radius_particles.resize(numPointsThatIHave);
			rho_particles.resize(numPointsThatIHave);

			//utility::dense::sph_kernel::WendlandC6 kernel_wendland;

// TODO uncomment
#ifdef WITH_OPENMP
#pragma omp parallel for schedule(static)
#endif
			for (size_t q = 0; q < numPointsThatIHave; q++) {
				// Sort all candidates for this query by distance
				std::sort(all_candidates[q].begin(), all_candidates[q].end(),
					[](const std::pair<size_t, float>& a, const std::pair<size_t, float>& b) {
						return a.second < b.second;
					});

				// Take the k-th closest distance (0-indexed, so k-1)
				if (all_candidates[q].size() >= k) {
					radius_particles[q] = std::sqrt(all_candidates[q][k - 1].second); // Convert squared distance to distance

					double h_inv = 1.0 / radius_particles[q];
					for (size_t i = 0; i < k; i++) {
						// SPH density estimate by SPHtoGrid						
						rho_particles[q] += mass_particles[q] * utility::dense::sph_kernel::W(rho_kernel, std::sqrt(all_candidates[q][i].second) * h_inv, h_inv);
					}

					//Corrects the density estimate for the kernel bias
					//See Dehnen & Aly 2012, eq. 18 + 19
					rho_particles[q] = utility::dense::sph_kernel::bias_correction(rho_kernel, rho_particles[q], mass_particles[q], h_inv, k);
				}
				else {
					radius_particles[q] = std::numeric_limits<float>::infinity();
					rho_particles[q] = std::numeric_limits<float>::infinity();
				}
			}

			MPI_Barrier(MPI_COMM_WORLD);
		}

		void run_knn(float* points, size_t N, int k, std::vector<float>& radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_cycling, common::SpaceData::DenseType& rho_kernel)
		{
			run_knn_cpu(points, N, k, radius_particles, rho_particles, mass_particles, use_cycling, rho_kernel);
		}

	}// nanoflann
} //utility