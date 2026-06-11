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

#ifdef _WIN32
#	define NOMINMAX
#endif

#include "data_communication.h"
#include "data_common.h"
#include "convert_common.h"

#ifdef WITH_GPU_CUDA
#include "utility/gpu_logging.h"
#include <cuda_runtime.h>
#endif

#ifdef WITH_OPENMP
# include <omp.h>
#endif

#include "utility/logging.h"
#include <iostream>

namespace space_converter {

	// ============================================================================
	// MPI Communication Functions
	// ============================================================================

	/**
	 * @brief Send data via MPI with automatic chunking for large transfers.
	 * 
	 * This function handles sending arbitrarily large data by splitting it into
	 * chunks if necessary. This avoids MPI message size limitations and improves
	 * reliability for large data transfers.
	 * 
	 * @param data Pointer to the data buffer to send
	 * @param size Total size of data in bytes
	 * @param data_type MPI datatype of the data
	 * @param from Rank of the sending process
	 * @param to Rank of the receiving process
	 */
	void mpi_send(void* data, size_t size, MPI_Datatype data_type, int from, int to)
	{
		// Maximum chunk size: 2 GB (2 billion bytes)
		const size_t unit_giga = 1000L * 1000L * 1000L * 2L;

		// Send data in chunks if it exceeds the maximum size
		size_t size_sended = 0;
		while (size - size_sended > unit_giga) {
			MPI_Send((char*)data + size_sended, unit_giga, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD);
			size_sended += unit_giga;
		}

		// Send remaining data
		MPI_Send((char*)data + size_sended, size - size_sended, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD);
	}

	/**
	 * @brief Receive data via MPI with automatic chunking for large transfers.
	 * 
	 * This function receives data sent by mpi_send(), handling large messages
	 * by receiving them in chunks. It must be paired with a corresponding mpi_send().
	 * 
	 * @param data Pointer to the buffer to store received data
	 * @param size Total size of data to receive in bytes
	 * @param data_type MPI datatype of the data
	 * @param from Rank of the sending process
	 * @param to Rank of the receiving process
	 */
	void mpi_recv(void* data, size_t size, MPI_Datatype data_type, int from, int to)
	{
		// Maximum chunk size: 2 GB (2 billion bytes)
		const size_t unit_giga = 1000L * 1000L * 1000L * 2L;

		// Receive data in chunks if it exceeds the maximum size
		size_t size_sended = 0;
		while (size - size_sended > unit_giga) {
			MPI_Recv((char*)data + size_sended, unit_giga, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			size_sended += unit_giga;
		}

		// Receive remaining data
		MPI_Recv((char*)data + size_sended, size - size_sended, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	/**
	 * @brief Perform MPI reduction (sum) across all processes with chunking.
	 * 
	 * This function aggregates float data from all processes using the SUM operation.
	 * It handles large arrays by processing them in chunks. The result is stored
	 * on the root process (rank 0).
	 * 
	 * @param ldata Pointer to local float data array
	 * @param gdata Pointer to store global summed result (valid on root)
	 * @param size Number of float elements in the arrays
	 */
	void mpi_reduce(float* ldata, float* gdata, size_t size)
	{
		// Maximum chunk size: 2 billion floats
		const size_t unit_giga = 1000L * 1000L * 1000L * 2L;

		// Reduce data in chunks if it exceeds the maximum size
		size_t size_sended = 0;
		while (size - size_sended > unit_giga) {
			MPI_Reduce(ldata + size_sended, gdata + size_sended, unit_giga, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			size_sended += unit_giga;
		}

		// Reduce remaining data
		MPI_Reduce(ldata + size_sended, gdata + size_sended, size - size_sended, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	// ============================================================================
	// MPI Initialization and Finalization
	// ============================================================================

	/// Global start time for performance measurement
	double t_start = 0;

	/**
	 * @brief Initialize MPI environment and configure process information.
	 * 
	 * This function initializes the MPI runtime, optionally starts performance
	 * profiling with MERIC, and queries process rank and size information.
	 * The root process (rank 0) prints initialization information.
	 * 
	 * @param argc Number of command-line arguments
	 * @param argv Array of command-line argument strings
	 * @param from_cl Reference to FromCL structure to populate with MPI info
	 */
	void init_mpi(int argc, char** argv, space_converter::FromCL& from_cl)
	{
#ifdef WITH_GPU_CUDA
		// Initialize CUDA before MPI to avoid CUDA_ERROR_INVALID_CONTEXT errors
		// from UCX/MPI trying to query CUDA devices during MPI_Init
		int deviceCount = 0;
		cudaError_t err = cudaGetDeviceCount(&deviceCount);
		if (err == cudaSuccess && deviceCount > 0) {
			// Set device 0 as default and establish CUDA context
			cudaSetDevice(0);
			cudaFree(0);  // Force context creation
			//printf("CUDA initialized before MPI: %d device(s) found\n", deviceCount);
		} else {
			printf("Warning: CUDA initialization before MPI failed or no devices: %s\n", 
			       cudaGetErrorString(err));
		}
#endif
		// Initialize MPI
		MPI_Init(&argc, &argv);
		
#ifdef WITH_MERIC
		// Initialize performance monitoring
		MERIC_Init();		
#endif		

#ifdef WITH_OPENMP
		// Record start time for performance measurement
		t_start = omp_get_wtime();
#endif

		// Get the total number of MPI processes
		MPI_Comm_size(MPI_COMM_WORLD, &from_cl.world_size);

		// Get the rank of this process
		MPI_Comm_rank(MPI_COMM_WORLD, &from_cl.world_rank);

		// LOG Init 
		LOG_Init(from_cl.world_rank);
		LOG_MeasureStart("Main");

#ifdef WITH_GPU_CUDA
		// Initialize GPU logging with MPI rank
		gpu_logging::set_gpu_log_rank(from_cl.world_rank);
#endif

		// Get the processor name
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int name_len;
		MPI_Get_processor_name(processor_name, &name_len);

		// Root process prints initialization message
		if (from_cl.world_rank < 1) {
			printf(
				"Start from processor %s, rank %d"
				" out of %d processors\n",
				processor_name,
				from_cl.world_rank,
				from_cl.world_size);
		}

#ifdef _WIN32
		// Debug
		//wait_on_attach_process(from_cl);
#endif
	}

	/**
	 * @brief Finalize MPI environment and print statistics.
	 * 
	 * This function finalizes the MPI runtime, stops performance profiling,
	 * and prints execution time statistics from the root process.
	 * 
	 * @param from_cl Reference to FromCL structure with MPI configuration
	 */
	void close_mpi(space_converter::FromCL& from_cl)
	{
#ifdef WITH_OPENMP
		// Root process prints execution statistics
		if (from_cl.world_rank == 0) {
			printf(
				"End from processor, rank %d"
				" out of %d processors, total time: %f\n",
				from_cl.world_rank,
				from_cl.world_size,
				omp_get_wtime() - t_start);

			fflush(0);
		}
#endif

		LOG_MeasureStop("Main");
#ifdef WITH_MERIC
		// Stop performance monitoring		
		MERIC_Close();
#endif		

		// Finalize MPI
		MPI_Finalize();
	}

	// ============================================================================
	// TCP Communication Functions
	// ============================================================================

	/**
	 * @brief Initialize TCP communication on the root process.
	 * 
	 * Only the root process (rank 0) establishes a TCP connection if remote
	 * communication is enabled. This centralizes network communication.
	 * 
	 * @param tcp_connection Reference to TcpConnection object to initialize
	 * @param from_cl Reference to FromCL structure with server configuration
	 */
	void init_communication(TcpConnection& tcp_connection, space_converter::FromCL& from_cl)
	{
		if (from_cl.remote && from_cl.world_rank == 0) {
			tcp_connection.init_sockets_data(from_cl.server.c_str(), from_cl.port);
		}
	}

	/**
	 * @brief Close TCP communication on the root process.
	 * 
	 * Cleanly shuts down the TCP connection established by init_communication().
	 * 
	 * @param tcp_connection Reference to TcpConnection object to close
	 * @param from_cl Reference to FromCL structure with server configuration
	 */
	void close_communication(TcpConnection& tcp_connection, space_converter::FromCL& from_cl)
	{
		if (from_cl.remote && from_cl.world_rank == 0) {
			tcp_connection.client_close();
			tcp_connection.server_close();
		}
	}

	// ============================================================================
	// Message Handling Functions
	// ============================================================================

	/**
	 * @brief Wait for and receive a message type from the TCP connection.
	 * 
	 * The root process receives a message type from the remote client. If remote
	 * mode is disabled or a connection error occurs, it sets the message to exit.
	 * The message type is then broadcast to all MPI processes.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData to store the received message type
	 */
	void wait_on_message(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				// Receive message type from remote client
				tcp_connection.recv_data_data((char*)&space_data.message_type, sizeof(space_data.message_type));

				// Set exit message on connection error
				if (tcp_connection.is_error()) {
					space_data.message_type = common::SpaceData::MessageType::eExit;
				}
			}

			// printf("message_type: %d\n", space_data.message_type); fflush(0);
		}

		// Broadcast message type to all processes
		MPI_Bcast((char*)&space_data.message_type, sizeof(space_data.message_type), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	/**
	 * @brief Send dataset information to the remote client.
	 * 
	 * Transmits animation parameters and particle data type information to the
	 * remote client. This allows the client to understand the available data
	 * before requesting specific subsets.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData containing animation parameters
	 * @param particle_data_types String describing available particle data types
	 */
	void send_info(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, std::string particle_data_types)
	{
		if (from_cl.world_rank == 0) {
			// Send animation information
			int anim_type = (int)space_data.anim_type;
			tcp_connection.send_data_data((char*)&anim_type, sizeof(int));
			tcp_connection.send_data_data((char*)&space_data.anim_start, sizeof(int));
			tcp_connection.send_data_data((char*)&space_data.anim_end, sizeof(int));

			// Send particle data type names
			int s = particle_data_types.length();
			tcp_connection.send_data_data((char*)&s, sizeof(int));
			tcp_connection.send_data_data((char*)particle_data_types.c_str(), sizeof(char) * s);

			// Wait for acknowledgment
			int ack;
			tcp_connection.recv_data_data((char*)&ack, sizeof(ack));
			printf("sended: particle and data types\n");
		}
	}
	/**
	 * @brief Receive bounding box request parameters from remote client.
	 * 
	 * The root process receives which particle type and data block the client
	 * wants to query for bounding box information, then broadcasts these
	 * parameters to all MPI processes.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData to store requested parameters
	 */
	void recv_requested_bbox(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				// Receive particle type and block ID
				tcp_connection.recv_data_data((char*)&space_data.particle_type, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.block_name_id, sizeof(int));
			}
		}

		// Broadcast to all processes
		MPI_Bcast((char*)&space_data.particle_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.block_name_id, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	/**
	 * @brief Send bounding box information to the remote client.
	 * 
	 * Calculates and transmits the normalized bounding box coordinates based on
	 * the original data extent. The bounding box is transformed to match the
	 * grid resolution and coordinate system used for visualization.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData with overall dataset bounds
	 * @param space_data_bbox Reference to SpaceData with specific bounding box
	 * @param convert_vdb_base Pointer to VDB converter (currently unused)
	 */
	void send_bbox(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::SpaceData& space_data_bbox, common::vdb::ConvertVDBBase* convert_vdb_base)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				// Initialize bounding box coordinates
				float bbox_min[3] = { 0.0f, 0.0f, 0.0f };
				float bbox_max[3] = { 0.0f, 0.0f, 0.0f };

				float bbox_dim = (float)space_data.object_size;

				// Calculate original data dimensions
				float dims_orig[3];
				dims_orig[0] = space_data.bbox_max_orig[0] - space_data.bbox_min_orig[0];
				dims_orig[1] = space_data.bbox_max_orig[1] - space_data.bbox_min_orig[1];
				dims_orig[2] = space_data.bbox_max_orig[2] - space_data.bbox_min_orig[2];

				float dims_orig_max = STDMAX(dims_orig[0], STDMAX(dims_orig[1], dims_orig[2]));

				// Transform bounding box to normalized grid coordinates
				bbox_min[0] = (space_data_bbox.bbox_min_orig[0] - space_data.bbox_min_orig[0]) * bbox_dim / dims_orig[0];
				bbox_min[1] = (space_data_bbox.bbox_min_orig[1] - space_data.bbox_min_orig[1]) * bbox_dim / dims_orig[1];
				bbox_min[2] = (space_data_bbox.bbox_min_orig[2] - space_data.bbox_min_orig[2]) * bbox_dim / dims_orig[2];

				bbox_max[0] = (space_data_bbox.bbox_max_orig[0] - space_data.bbox_min_orig[0]) * bbox_dim / dims_orig[0];
				bbox_max[1] = (space_data_bbox.bbox_max_orig[1] - space_data.bbox_min_orig[1]) * bbox_dim / dims_orig[1];
				bbox_max[2] = (space_data_bbox.bbox_max_orig[2] - space_data.bbox_min_orig[2]) * bbox_dim / dims_orig[2];

				// Send bounding box to client
				tcp_connection.send_data_data((char*)&bbox_min[0], sizeof(float) * 3);
				tcp_connection.send_data_data((char*)&bbox_max[0], sizeof(float) * 3);

				// Wait for acknowledgment
				int ack;
				tcp_connection.recv_data_data((char*)&ack, sizeof(ack));
				printf("sended: BBOX\n");
			}
		}
	}

	/**
	 * @brief Receive data request parameters from the remote client.
	 * 
	 * The root process receives all the parameters specifying what data the client
	 * wants (bounding box, resolution, particle type, etc.) and broadcasts them
	 * to all MPI processes for coordinated processing.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData to store all request parameters
	 */
	void recv_requested_data(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				// Receive spatial parameters
				tcp_connection.recv_data_data((char*)&space_data.bbox_min[0], sizeof(float) * 3);
				tcp_connection.recv_data_data((char*)&space_data.bbox_max[0], sizeof(float) * 3);
				tcp_connection.recv_data_data((char*)&space_data.bbox_dim, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.grid_transform, sizeof(float));
				
				// Receive data selection parameters
				tcp_connection.recv_data_data((char*)&space_data.particle_type, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.block_name_id, sizeof(int));
				
				// Receive processing parameters
				tcp_connection.recv_data_data((char*)&space_data.extracted_type, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.dense_type, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.dense_norm, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.object_size, sizeof(float));
				tcp_connection.recv_data_data((char*)&space_data.particle_radius_multiplier, sizeof(float));
				
				// Receive filtering parameters
				tcp_connection.recv_data_data((char*)&space_data.filter_min, sizeof(float));
				tcp_connection.recv_data_data((char*)&space_data.filter_max, sizeof(float));
				
				// Receive animation parameters
				tcp_connection.recv_data_data((char*)&space_data.frame, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.anim_type, sizeof(int));
				tcp_connection.recv_data_data((char*)&space_data.anim_task_counter, sizeof(int));
			}
		}

		// Broadcast all parameters to all MPI processes
		MPI_Bcast((char*)&space_data.bbox_min[0], sizeof(float) * 3, MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.bbox_max[0], sizeof(float) * 3, MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.bbox_dim, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.grid_transform, sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.particle_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.block_name_id, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.extracted_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.dense_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.dense_norm, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.object_size, sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.filter_min, sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.filter_max, sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.frame, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.anim_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.anim_task_counter, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	// ============================================================================
	// VDB Data Transmission Functions
	// ============================================================================

	/**
	 * @brief Send VDB volumetric data to the remote client.
	 * 
	 * Transmits the processed volumetric data in OpenVDB, NanoVDB, or raw particle
	 * format to the remote client. Includes metadata about value ranges and animation
	 * frame information.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData containing metadata
	 * @param grid_handle Reference to VDBParticles containing the grid data
	 */
	void send_vdb(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_handle)
	{
		if (from_cl.world_rank == 0 && grid_handle.type == common::vdb::VDBParticleType::eVector) {
			double t = omp_get_wtime();
			
			// Get grid data size
			std::size_t size = grid_handle.vector_grid.size();

			// Determine output file type based on extraction mode
			int file_type = common::FTI_OPENVDB;
			if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
				file_type = common::FTI_RAW_PART; // Raw particle data
			}
			else if (from_cl.use_nanovdb) {
				file_type = common::FTI_NANOVDB; // GPU-friendly NanoVDB format
			}
			else if (from_cl.use_cub) {
				file_type = common::FTI_CUB; // GPU-friendly CUB format
			}		

			// Send file type identifier
			tcp_connection.send_data_data((char*)&file_type, sizeof(file_type));

			// Send VDB grid data
			tcp_connection.send_data_data((char*)&size, sizeof(size));
			tcp_connection.send_data_data((char*)grid_handle.vector_grid.data(), size);

			// Send value range metadata
			tcp_connection.send_data_data((char*)&space_data.min_value, sizeof(space_data.min_value));
			tcp_connection.send_data_data((char*)&space_data.max_value, sizeof(space_data.max_value));
			tcp_connection.send_data_data((char*)&space_data.min_value_reduced, sizeof(space_data.min_value_reduced));
			tcp_connection.send_data_data((char*)&space_data.max_value_reduced, sizeof(space_data.max_value_reduced));
						
			// Send animation frame count
			int frames = 1;
			if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
				frames = from_cl.world_size; // One frame per MPI process
			}
			tcp_connection.send_data_data((char*)&frames, sizeof(frames));

			int ack;
			tcp_connection.recv_data_data((char*)&ack, sizeof(ack));
			printf("rank #%d: sendOpenVDB time: %f\n", from_cl.world_rank, omp_get_wtime() - t);
			printf("sended: vdb\n");
		}
	}

	/**
	 * @brief Send VDB file path reference to the remote client.
	 * 
	 * Instead of sending the entire VDB data, transmits only the file path.
	 * Useful when the client has direct access to the file system where VDB
	 * files are stored. Also includes metadata about value ranges.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData containing file path and metadata
	 * @param grid_handle Reference to VDBParticles (used for type checking)
	 */
	void send_path(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_handle)
	{
		if (from_cl.world_rank == 0 && grid_handle.type == common::vdb::VDBParticleType::eVector) {
			double t = omp_get_wtime();

			// Send file path type identifier
			int file_type = common::FTI_PATH;
			tcp_connection.send_data_data((char*)&file_type, sizeof(file_type));

			// Send VDB file path
			std::string full_filepath = space_data.full_filepath;
			const char* cf = full_filepath.c_str();

			size_t size = strlen(cf);
			tcp_connection.send_data_data((char*)&size, sizeof(size));
			tcp_connection.send_data_data((char*)cf, size);

			// Send value range metadata
			tcp_connection.send_data_data((char*)&space_data.min_value, sizeof(space_data.min_value));
			tcp_connection.send_data_data((char*)&space_data.max_value, sizeof(space_data.max_value));
			tcp_connection.send_data_data((char*)&space_data.min_value_reduced, sizeof(space_data.min_value_reduced));
			tcp_connection.send_data_data((char*)&space_data.max_value_reduced, sizeof(space_data.max_value_reduced));

			// Send animation frame count
			int frames = 1;
			if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
				frames = from_cl.world_size;
			}
			tcp_connection.send_data_data((char*)&frames, sizeof(frames));

			int ack;
			tcp_connection.recv_data_data((char*)&ack, sizeof(ack));
			printf("rank #%d: sendOpenVDB time: %f\n", from_cl.world_rank, omp_get_wtime() - t);
			printf("sended: vdb\n");
		}
	}

	// ============================================================================
	// Debugging and Utility Functions
	// ============================================================================

	/**
	 * @brief Wait for debugger attachment (debugging utility).
	 * 
	 * This function creates a waiting loop that allows a debugger to be attached
	 * to the process. The root process waits in a loop that can be broken by
	 * setting attach=0 in the debugger. Useful for debugging MPI applications.
	 * 
	 * @param from_cl Reference to FromCL structure with configuration
	 * 
	 * @note Set the 'attach' variable to 0 in the debugger to continue execution.
	 */
	void wait_on_attach_process(space_converter::FromCL& from_cl)
	{
		if (from_cl.world_rank == 0) {
			// Wait loop for debugger attachment
			int attach = 1;			
			while (attach) {
#ifdef _WIN32
				Sleep(1000);
#endif
				printf("wait_on_attach_process\n"); fflush(0);
			}			
		}

		// Synchronize all processes
		int attach_process = 1;
		MPI_Bcast((char*)&attach_process, sizeof(attach_process), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	/**
	 * @brief Read a large test file for benchmarking or testing purposes.
	 * 
	 * This utility function reads a large file specified by the environment variable
	 * CONVERTER_LARGE_TEST_FILE. All MPI processes synchronize before and after
	 * reading to ensure coordinated I/O testing. Useful for testing file system
	 * performance and MPI synchronization.
	 * 
	 * @note Set the CONVERTER_LARGE_TEST_FILE environment variable to enable.
	 */
	void read_large_test_file() {
		const char* filename = getenv("CONVERTER_LARGE_TEST_FILE");
		if (filename) {
			// Synchronize all processes before file I/O
			MPI_Barrier(MPI_COMM_WORLD);

			// Open file at end to determine size
			std::ifstream file(filename, std::ios::binary | std::ios::ate);
			if (!file) {
				std::cerr << "Failed to open file: " << filename << std::endl;
				exit(-1);
			}

			// Get file size and rewind
			std::streamsize size = file.tellg();
			file.seekg(0, std::ios::beg);

			// Allocate buffer and read entire file
			char* buffer = new char[size];
			if (!file.read(buffer, size)) {
				std::cerr << "Failed to read file: " << filename << std::endl;
				exit(-1);
			}

			std::cout << "Read " << size << " bytes from " << filename << std::endl;

			// Clean up buffer
			delete [] buffer;
			
			// Synchronize all processes after file I/O
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}

} // namespace space_converter