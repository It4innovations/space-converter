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

#include "cudakdtree_tool.h"

#include "cukd/cukd-math.h"
#include "cukd/traverse-stack-free.h"
#include "cukd/knn.h"
#include <mpi.h>
#include <stdexcept>
#include <cuda_runtime.h>

#define CUKD_MPI_CALL(fctCall)                                          \
  { int rc = MPI_##fctCall;                                             \
    if (rc != MPI_SUCCESS)                                              \
      throw std::runtime_error(std::string(__PRETTY_FUNCTION__)+#fctCall); }

using cukd::divRoundUp;

#include <climits>  // For INT_MAX
// Use a macro for max bytes per transfer (less than INT_MAX to be safe)
#define MAX_MPI_BYTES (INT_MAX - 1024)

#include "utility/dense_utility.h"

namespace utility {
	namespace cudakdtree {

        __global__ void runQuery(float3* tree, int N,
            uint64_t* candidateLists, int k, float maxRadius,
            float3* queries, int numQueries,
            int round)
        {
            size_t tid = threadIdx.x + (size_t)blockIdx.x * blockDim.x;
            if (tid >= numQueries) return;

            float3 qp = queries[tid];
            cukd::FlexHeapCandidateList cl(candidateLists + (size_t)k * tid, k,
                round == 0 ? maxRadius : -1.f);
            cukd::stackFree::knn(cl, qp, tree, N);
        }

        __global__ void extractFinalResult(
            float* d_radius_particles,
            float* d_rho_particles,
            float* d_mass_particles,
            int numPoints,
            int k,
            uint64_t* candidateLists,
            common::SpaceData::DenseType rho_kernel
            )
        {
            size_t tid = threadIdx.x + (size_t)blockIdx.x * blockDim.x;
            if (tid >= numPoints) return;

            cukd::FlexHeapCandidateList cl(candidateLists + (size_t)k * tid, k, -1.f);
            float result = cl.returnValue();
            if (!isinf(result)) {
                result = sqrtf(result);
                d_radius_particles[tid] = result;

                //utility::dense::sph_kernel::WendlandC6 kernel_wendland;

                double h_inv = 1.0 / d_radius_particles[tid];

                for (int i = 0; i < k; i++) {
                    float result_i = cukd::uint_as_float(cl.entry[i] >> 32); //cl.get_dist2(i);
                    if (!isinf(result_i) && !isinf(result)) {
                        result_i = sqrtf(result_i);

                        // SPH density estimate by SPHtoGrid
                        d_rho_particles[tid] += d_mass_particles[tid] * utility::dense::sph_kernel::W(rho_kernel, result_i * h_inv, h_inv);
                    }
                    else {
                        d_rho_particles[tid] = result_i;
                        break;
                    }
                }

                //SPH density estimate by SPHtoGrid
                //Corrects the density estimate for the kernel bias
                //See Dehnen & Aly 2012, eq. 18 + 19
                if (!isinf(d_rho_particles[tid])) {
                    d_rho_particles[tid] = utility::dense::sph_kernel::bias_correction(rho_kernel, d_rho_particles[tid], d_mass_particles[tid], h_inv, k);
                }

            }
            else {
                d_radius_particles[tid] = result;
                d_rho_particles[tid] = result;
            }
        }

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
                    CUKD_MPI_CALL(Irecv(recvPtr + offsetRecv, chunkRecv, MPI_BYTE, recvPeer, 0,
                                        MPI_COMM_WORLD, &requests[0]));
                } else {
                    requests[0] = MPI_REQUEST_NULL;
                }

                if (offsetSend < totalBytesSend) {
                    CUKD_MPI_CALL(Isend(sendPtr + offsetSend, chunkSend, MPI_BYTE, sendPeer, 0,
                                        MPI_COMM_WORLD, &requests[1]));
                } else {
                    requests[1] = MPI_REQUEST_NULL;
                }

                CUKD_MPI_CALL(Waitall(2, requests, MPI_STATUSES_IGNORE));

                offsetRecv += chunkRecv;
                offsetSend += chunkSend;
            }

        }

        void run_knn_gpu(float* points, size_t numPointsThatIHave, int k, std::vector<float>& radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_cycling, float maxRadius, common::SpaceData::DenseType& rho_kernel)
        {
            int mpi_rank, mpi_size;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

            //float maxRadius = std::numeric_limits<float>::infinity();

            // -----------------------------------------------------------------------------
            // find out max num points anybody has, so we can allocate
            // -----------------------------------------------------------------------------
            int N = numPointsThatIHave;
            int maxNumPointsAnybodyHas = 0;

            CUKD_MPI_CALL(Allreduce(&N, &maxNumPointsAnybodyHas, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));

            float3* d_tree = 0;
            float3* d_tree_recv = 0;
            //int N = myPoints.size();
            // alloc N+1 so we can store one more if anytoher rank gets oen more point
            CUKD_CUDA_CALL(MallocManaged((void**)&d_tree, (maxNumPointsAnybodyHas + 1) * sizeof(float3)));
            CUKD_CUDA_CALL(MallocManaged((void**)&d_tree_recv, (maxNumPointsAnybodyHas + 1) * sizeof(float3)));
            CUKD_CUDA_CALL(Memcpy(d_tree, points, (size_t)N * sizeof(float3),
                cudaMemcpyDefault));


            double start_time, end_time;
            // Start timing before your main computation
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();

            cukd::buildTree(d_tree, N);

            // End timing after computation
            MPI_Barrier(MPI_COMM_WORLD);
            end_time = MPI_Wtime();

            // Print results on rank 0
            if (mpi_rank == 0) {
                printf("Total execution time (buildTree): %.6f seconds\n", end_time - start_time);
            }

            int numQueries = N;// myPoints.size();

#if 0 //nosplit
            float3* d_queries;            
            uint64_t* d_cand;
            CUKD_CUDA_CALL(MallocManaged((void**)&d_queries, (size_t)N * sizeof(float3)));
            CUKD_CUDA_CALL(Memcpy(d_queries, points, N * sizeof(float3), cudaMemcpyDefault));
            CUKD_CUDA_CALL(MallocManaged((void**)&d_cand, (size_t)N * k * sizeof(uint64_t)));

#else
            //TODO!!!
            // Choose S (number of splits/batches)
            int S = 1;//k; // or make this a parameter

            if (S > N) {
                S = 1;
            }

            size_t batch_size = numQueries / S;

            // Allocate output buffers for all queries
            radius_particles.resize(numPointsThatIHave);
            rho_particles.resize(numPointsThatIHave);

            //float* d_radius_particles = 0;
            //float* d_rho_particles = 0;
            //float* d_mass_particles = 0;
            //CUKD_CUDA_CALL(MallocManaged((void**)&d_radius_particles, numPointsThatIHave * sizeof(float)));
            //CUKD_CUDA_CALL(MallocManaged((void**)&d_rho_particles, numPointsThatIHave * sizeof(float)));
            //CUKD_CUDA_CALL(MallocManaged((void**)&d_mass_particles, numPointsThatIHave * sizeof(float)));
            //CUKD_CUDA_CALL(Memcpy(d_mass_particles, mass_particles.data(), numPointsThatIHave * sizeof(float), cudaMemcpyHostToDevice));
#endif            

            // -----------------------------------------------------------------------------
            // now, do the queries and cycling:
            // -----------------------------------------------------------------------------
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();

			int round_size = mpi_size;
			if (!use_cycling) {
				round_size = 1;
			}

			// Allocate
            uint64_t* d_cand = nullptr;
            size_t d_cand_size = 0;

            float3* d_batch_queries = nullptr;
			size_t d_batch_queries_size = 0;

            float* d_batch_mass = nullptr;
			size_t d_batch_mass_size = 0;

            float* d_batch_radius = nullptr;
			size_t d_batch_radius_size = 0;

            float* d_batch_rho = nullptr;
			size_t d_batch_rho_size = 0;

            for (int round = 0; round < round_size; round++) {
				if (mpi_rank == 0) {
					printf("Starting round (cudakdtree-gpu cycling) %d\n", round);
				}  

                if (round == 0) {
                    // nothing to do , we already have our own tree
                }
                else {
                    MPI_Request requests[2];
                    int sendCount = N;
                    int recvCount = 0;
                    int sendPeer = (mpi_rank + 1) % mpi_size;
                    int recvPeer = (mpi_rank + mpi_size - 1) % mpi_size;
                    CUKD_MPI_CALL(Irecv(&recvCount, 1 * sizeof(int), MPI_BYTE, recvPeer, 0,
                        MPI_COMM_WORLD, &requests[0]));
                    CUKD_MPI_CALL(Isend(&sendCount, 1 * sizeof(int), MPI_BYTE, sendPeer, 0,
                        MPI_COMM_WORLD, &requests[1]));
                    CUKD_MPI_CALL(Waitall(2, requests, MPI_STATUSES_IGNORE));

                    char* recvPtr = reinterpret_cast<char*>(d_tree_recv);
                    char* sendPtr = reinterpret_cast<char*>(d_tree);
                    size_t totalBytesRecv = recvCount * sizeof(*d_tree);
                    size_t totalBytesSend = sendCount * sizeof(*d_tree);
                    mpi_cycling(recvPeer, recvPtr, totalBytesRecv, sendPeer, sendPtr, totalBytesSend);

                    N = recvCount;
                    std::swap(d_tree, d_tree_recv);
                }

#if 0 //nosplit                
                // -----------------------------------------------------------------------------
                runQuery << <divRoundUp(numQueries, 1024), 1024 >> >
                    (/* tree */d_tree, N,
                        /* query params */d_cand, k, maxRadius,
                        /* query points */d_queries, numQueries,
                        round);
                CUKD_CUDA_CALL(DeviceSynchronize());
#else
                // Process queries in S batches
                for (int s = 0; s < S; ++s) {
                    size_t start = s * batch_size;
                    size_t this_batch = (s == S - 1) ? numQueries - start : batch_size;
                    if (this_batch == 0) continue;

                    // Allocate only for this batch
                    if (d_cand_size < this_batch * (size_t)k * sizeof(uint64_t)) {

                        if (d_cand != nullptr) {
							CUKD_CUDA_CALL(Free(d_cand));
                        }

                        d_cand_size = this_batch * (size_t)k * sizeof(uint64_t);
                        CUKD_CUDA_CALL(MallocManaged((void**)&d_cand, d_cand_size));
                    }
					CUKD_CUDA_CALL(Memset(d_cand, 0, d_cand_size));

                    // Prepare batch queries
					if (d_batch_queries_size < this_batch * sizeof(float3)) {
						if (d_batch_queries != nullptr) {
							CUKD_CUDA_CALL(Free(d_batch_queries));
						}
						d_batch_queries_size = this_batch * sizeof(float3);
                        CUKD_CUDA_CALL(MallocManaged((void**)&d_batch_queries, d_batch_queries_size));
					}
                    CUKD_CUDA_CALL(Memcpy(d_batch_queries, ((float3*)points) + start, this_batch * sizeof(float3), cudaMemcpyDefault));

                    // Prepare batch mass
					if (d_batch_mass_size < this_batch * sizeof(float)) {
						if (d_batch_mass != nullptr) {
							CUKD_CUDA_CALL(Free(d_batch_mass));
						}
						d_batch_mass_size = this_batch * sizeof(float);
                        CUKD_CUDA_CALL(MallocManaged((void**)&d_batch_mass, d_batch_mass_size));
					}                   
                    CUKD_CUDA_CALL(Memcpy(d_batch_mass, mass_particles.data() + start, this_batch * sizeof(float), cudaMemcpyHostToDevice));

                    // Output for this batch
					if (d_batch_radius_size < this_batch * sizeof(float)) {
						if (d_batch_radius != nullptr) {
							CUKD_CUDA_CALL(Free(d_batch_radius));
						}
						d_batch_radius_size = this_batch * sizeof(float);
                        CUKD_CUDA_CALL(MallocManaged((void**)&d_batch_radius, d_batch_radius_size));
					}
					CUKD_CUDA_CALL(Memset(d_batch_radius, 0, d_batch_radius_size));

					if (d_batch_rho_size < this_batch * sizeof(float)) {
						if (d_batch_rho != nullptr) {
							CUKD_CUDA_CALL(Free(d_batch_rho));
						}
						d_batch_rho_size = this_batch * sizeof(float);
                        CUKD_CUDA_CALL(MallocManaged((void**)&d_batch_rho, d_batch_rho_size));
					}
					CUKD_CUDA_CALL(Memset(d_batch_rho, 0, d_batch_rho_size));

                    // Run query for this batch
                    runQuery << <divRoundUp(this_batch, 1024ULL), 1024ULL >> >(
                        d_tree, N,
                        d_cand, k, maxRadius,
                        d_batch_queries, this_batch,
                        round
                    );
                    CUKD_CUDA_CALL(DeviceSynchronize());

                    // Extract results for this batch
                    extractFinalResult<<<divRoundUp(this_batch, 1024ULL), 1024ULL >>>(
						d_batch_radius, d_batch_rho, d_batch_mass, this_batch, k, d_cand, rho_kernel
                    );
                    CUKD_CUDA_CALL(DeviceSynchronize());

                    // Copy results to the full output arrays
                    CUKD_CUDA_CALL(Memcpy(radius_particles.data() + start, d_batch_radius, this_batch * sizeof(float), cudaMemcpyDeviceToHost));
                    CUKD_CUDA_CALL(Memcpy(rho_particles.data() + start, d_batch_rho, this_batch * sizeof(float), cudaMemcpyDeviceToHost));
                }
#endif
            }

            // Free batch memory
            CUKD_CUDA_CALL(Free(d_cand));
            CUKD_CUDA_CALL(Free(d_batch_queries));
            CUKD_CUDA_CALL(Free(d_batch_mass));
            CUKD_CUDA_CALL(Free(d_batch_radius));
            CUKD_CUDA_CALL(Free(d_batch_rho));

            // End timing after computation
            MPI_Barrier(MPI_COMM_WORLD);
            end_time = MPI_Wtime();

            // Print results on rank 0
            if (mpi_rank == 0) {
                printf("Total execution time (queries and cycling are done): %.6f seconds\n", end_time - start_time);
            }

#if 0  //nosplit            
            std::cout << "done all queries..." << std::endl;
            float* d_radius_particles = 0;
            float* d_rho_particles = 0;
            float* d_mass_particles = 0;
            CUKD_CUDA_CALL(MallocManaged((void**)&d_radius_particles, numPointsThatIHave * sizeof(float)));
            CUKD_CUDA_CALL(MallocManaged((void**)&d_rho_particles, numPointsThatIHave * sizeof(float)));
            
            CUKD_CUDA_CALL(MallocManaged((void**)&d_mass_particles, numPointsThatIHave * sizeof(float)));
            CUKD_CUDA_CALL(Memcpy(d_mass_particles, mass_particles.data(), numPointsThatIHave * sizeof(float), cudaMemcpyHostToDevice));

            extractFinalResult << <divRoundUp(numQueries, 1024), 1024 >> >
                (d_radius_particles, d_rho_particles, d_mass_particles, numQueries, k, d_cand);

            CUKD_CUDA_CALL(DeviceSynchronize());

            radius_particles.resize(numPointsThatIHave);
            rho_particles.resize(numPointsThatIHave);
            CUKD_CUDA_CALL(Memcpy(radius_particles.data(), d_radius_particles, numPointsThatIHave * sizeof(float), cudaMemcpyDeviceToHost));
            CUKD_CUDA_CALL(Memcpy(rho_particles.data(), d_rho_particles, numPointsThatIHave * sizeof(float), cudaMemcpyDeviceToHost));

            MPI_Barrier(MPI_COMM_WORLD);
#endif            
        }

        void runQuery_host(
            float3 *tree, size_t N,
            uint64_t *candidateLists, int k, float maxRadius,
            float3 *queries, size_t numQueries,
            int round)
        {
#pragma omp parallel for
            for (size_t tid = 0; tid < numQueries; tid++) {
                float3 qp = queries[tid];
                cukd::FlexHeapCandidateList cl(&candidateLists[(size_t)k * tid], k,
                    round == 0 ? maxRadius : -1.f);
                cukd::stackFree::knn(cl, qp, tree, N);
            }
        }

        void extractFinalResult_host(
            float *radius_particles,
            float *rho_particles,
            float *mass_particles,
            size_t numPoints,
            int k,
            uint64_t *candidateLists,
            common::SpaceData::DenseType rho_kernel
            )
        {
            //utility::dense::sph_kernel::WendlandC6 kernel_wendland;

#pragma omp parallel for
            for (size_t tid = 0; tid < numPoints; tid++) {
                cukd::FlexHeapCandidateList cl(&candidateLists[(size_t)k * tid], k, -1.f);
                float result = cl.returnValue();
                if (!isinf(result)) {
                    result = sqrtf(result);

                    radius_particles[tid] = result;

                    double h_inv = 1.0 / radius_particles[tid];

                    for (int i = 0; i < k; i++) {
                        float result_i = cukd::uint_as_float(cl.entry[i] >> 32); //cl.get_dist2(i);
                        if (!isinf(result_i) && !isinf(result)) {
                            result_i = sqrtf(result_i);

                            // SPH density estimate by SPHtoGrid						
                            rho_particles[tid] += mass_particles[tid] * utility::dense::sph_kernel::W(rho_kernel, result_i * h_inv, h_inv);
                        }
                        else {
                            rho_particles[tid] = result_i;
                            break;
                        }
                    }

                    //SPH density estimate by SPHtoGrid
                    //Corrects the density estimate for the kernel bias
                    //See Dehnen & Aly 2012, eq. 18 + 19
                    if (!isinf(rho_particles[tid])) {
                        rho_particles[tid] = utility::dense::sph_kernel::bias_correction(rho_kernel, rho_particles[tid], mass_particles[tid], h_inv, k);
                    }
                }
                else {
                    radius_particles[tid] = result;
					rho_particles[tid] = result;
                }
            }
        }

        void run_knn_cpu(float *points, size_t numPointsThatIHave, int k, std::vector<float> &radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_cycling, float maxRadius, common::SpaceData::DenseType& rho_kernel)
        {
            int mpi_rank, mpi_size;
            MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
            MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

            //float maxRadius = std::numeric_limits<float>::infinity();

            // -----------------------------------------------------------------------------
            // find out max num points anybody has, so we can allocate
            // -----------------------------------------------------------------------------
            int N = numPointsThatIHave;
            int maxNumPointsAnybodyHas = 0;

            CUKD_MPI_CALL(Allreduce(&N, &maxNumPointsAnybodyHas, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD));

            std::vector<float3> tree((maxNumPointsAnybodyHas + 1));
            std::vector<float3> tree_recv((maxNumPointsAnybodyHas + 1));
            memcpy(tree.data(), points, N * sizeof(float3));

            // Add timing to your mpiHugeQuery.cu
            double start_time, end_time;
            // Start timing before your main computation
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();

            cukd::buildTree_host(tree.data(), N);

            // End timing after computation
            MPI_Barrier(MPI_COMM_WORLD);
            end_time = MPI_Wtime();

            // Print results on rank 0
            if (mpi_rank == 0) {
                printf("Total execution time (buildTree_host): %.6f seconds\n", end_time - start_time);
            }

            size_t numQueries = N;

#if 0  //nosplit            
            std::vector<float3>  queries(N);
            memcpy(queries.data(), points, (size_t)N * sizeof(float3));
            std::vector<uint64_t>  cand((size_t)N * k);
#else
            //TODO!!!
            int S = 1;//k; // for example, or make it a parameter

            if (S > N) {
                S = 1;
            }

            size_t batch_size = numQueries / S;

            // Allocate output buffers for all queries
            radius_particles.resize(numPointsThatIHave);
            rho_particles.resize(numPointsThatIHave);            
#endif            

            // -----------------------------------------------------------------------------
            // now, do the queries and cycling:
            // -----------------------------------------------------------------------------
            MPI_Barrier(MPI_COMM_WORLD);
            start_time = MPI_Wtime();

			int round_size = mpi_size;
			if (!use_cycling) {
				round_size = 1;
			}

            std::vector<uint64_t> cand;
            std::vector<float> batch_radius;
            std::vector<float> batch_rho;

            for (int round = 0; round < round_size; round++) {
				if (mpi_rank == 0) {
					printf("Starting round (cudakdtree-cpu cycling) %d\n", round);
				}                

                if (round == 0) {
                    // nothing to do , we already have our own tree
                }
                else {
                    MPI_Request requests[2];
                    int sendCount = N;
                    int recvCount = 0;
                    int sendPeer = (mpi_rank + 1) % mpi_size;
                    int recvPeer = (mpi_rank + mpi_size - 1) % mpi_size;
                    CUKD_MPI_CALL(Irecv(&recvCount, 1 * sizeof(int), MPI_BYTE, recvPeer, 0,
                        MPI_COMM_WORLD, &requests[0]));
                    CUKD_MPI_CALL(Isend(&sendCount, 1 * sizeof(int), MPI_BYTE, sendPeer, 0,
                        MPI_COMM_WORLD, &requests[1]));
                    CUKD_MPI_CALL(Waitall(2, requests, MPI_STATUSES_IGNORE));

                    char* recvPtr = reinterpret_cast<char*>(tree_recv.data());
                    char* sendPtr = reinterpret_cast<char*>(tree.data());
                    size_t totalBytesRecv = recvCount * sizeof(float3);
                    size_t totalBytesSend = sendCount * sizeof(float3);
                    mpi_cycling(recvPeer, recvPtr, totalBytesRecv, sendPeer, sendPtr, totalBytesSend);               

                    N = recvCount;
                    std::swap(tree, tree_recv);
                }
#if 0  //nosplit                
                // -----------------------------------------------------------------------------
                runQuery_host(tree, N,
                    cand, k, maxRadius,
                    queries, numQueries,
                    round);
#else
                // Process queries in S batches
                for (int s = 0; s < S; ++s) {
                    size_t start = s * batch_size;
                    //size_t end = std::min(start + batch_size, numQueries);
                    //size_t this_batch = end - start;
                    size_t this_batch = (s == S - 1) ? numQueries - start : batch_size;

                    if (this_batch == 0) continue;

                    // Allocate only for this batch
                    if (this_batch * (size_t)k > cand.size()) {
                        cand.resize(this_batch * (size_t)k);
                    }
					memset(cand.data(), 0, this_batch * (size_t)k * sizeof(uint64_t));

                    if (this_batch > batch_radius.size()) {
						batch_radius.resize(this_batch);
                    }
					memset(batch_radius.data(), 0, this_batch * sizeof(float));

                    if (this_batch > batch_rho.size()) {
                        batch_rho.resize(this_batch);
                    }
					memset(batch_rho.data(), 0, this_batch * sizeof(float));

                    // Prepare batch queries
                    float3 *batch_queries = (float3*)points + start;

                    // Run query for this batch
                    runQuery_host(tree.data(), N, cand.data(), k, maxRadius, batch_queries, this_batch, round);

                    // Extract results for this batch
                    float *batch_mass = mass_particles.data() + start;

                    extractFinalResult_host(batch_radius.data(), batch_rho.data(), batch_mass, this_batch, k, cand.data(), rho_kernel);

                    // Copy results to the full output arrays
                    std::copy(batch_radius.begin(), batch_radius.end(), radius_particles.begin() + start);
                    std::copy(batch_rho.begin(), batch_rho.end(), rho_particles.begin() + start);
                }
#endif                    
            }

            // End timing after computation
            MPI_Barrier(MPI_COMM_WORLD);
            end_time = MPI_Wtime();

            // Print results on rank 0
            if (mpi_rank == 0) {
                printf("Total execution time (queries and cycling are done): %.6f seconds\n", end_time - start_time);
            }

#if 0  //nosplit
            radius_particles.resize(numPointsThatIHave);
            rho_particles.resize(numPointsThatIHave);
            extractFinalResult_host(radius_particles, rho_particles, mass_particles, numPointsThatIHave, k, cand);

            MPI_Barrier(MPI_COMM_WORLD);
#endif           
        }

        void run_knn(float* points, size_t N, int k, std::vector<float>& radius_particles, std::vector<float>& rho_particles, std::vector<float>& mass_particles, bool use_gpu, bool use_cycling, float max_radius, common::SpaceData::DenseType& rho_kernel)
        {
            if (use_gpu)
                run_knn_gpu(points, N, k, radius_particles, rho_particles, mass_particles, use_cycling, max_radius, rho_kernel);
            else
                run_knn_cpu(points, N, k, radius_particles, rho_particles, mass_particles, use_cycling, max_radius, rho_kernel);
        }

	}// cudakdtree
} //utility