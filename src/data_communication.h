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

#pragma once

#include <mpi.h>
#include <vector>

#include "args_processing.h"
#include "convert_vdb.h"
#include "data_common.h"
#include "renderengine_tcp.h"

/**
 * @namespace space_converter
 * @brief Provides MPI communication and data exchange functionalities.
 * 
 * This namespace contains functions for managing parallel communication between
 * MPI processes, TCP client-server communication, and data transmission for
 * volumetric simulation data conversion and visualization.
 */
namespace space_converter {

	// ========================================================================
	// MPI Communication Functions
	// ========================================================================

	/**
	 * @brief Send data via MPI with automatic chunking for large transfers.
	 * 
	 * Handles sending arbitrarily large data by splitting into chunks automatically.
	 * This avoids MPI message size limitations for large data transfers.
	 * 
	 * @param data Pointer to the data buffer to send
	 * @param size Total size of data in bytes
	 * @param data_type MPI datatype of the data
	 * @param from Rank of the sending process
	 * @param to Rank of the receiving process
	 */
	void mpi_send(void* data, size_t size, MPI_Datatype data_type, int from, int to);

	/**
	 * @brief Receive data via MPI with automatic chunking for large transfers.
	 * 
	 * Receives data sent by mpi_send(), handling large messages by receiving
	 * them in chunks. Must be paired with a corresponding mpi_send() call.
	 * 
	 * @param data Pointer to the buffer to store received data
	 * @param size Total size of data to receive in bytes
	 * @param data_type MPI datatype of the data
	 * @param from Rank of the sending process
	 * @param to Rank of the receiving process
	 */
	void mpi_recv(void* data, size_t size, MPI_Datatype data_type, int from, int to);

	/**
	 * @brief Perform MPI reduction (sum) across all processes with chunking.
	 * 
	 * Aggregates float data from all processes using the SUM operation.
	 * Handles large arrays by processing them in chunks. Result stored on root.
	 * 
	 * @param ldata Pointer to local float data array
	 * @param gdata Pointer to store global summed result (valid on root process)
	 * @param size Number of float elements in the arrays
	 */
	void mpi_reduce(float* ldata, float* gdata, size_t size);
    // @param argc: Number of command-line arguments.
    // @param argv: Array of command-line argument strings.
    // @param from_cl: Reference to the FromCL struct to populate with MPI-related configuration.
    void init_mpi(int argc, char** argv, space_converter::FromCL& from_cl);

    // Finalize the MPI environment and clean up resources.
    // @param from_cl: Reference to the FromCL struct for MPI-related configuration.
    void close_mpi(space_converter::FromCL& from_cl);

	// ========================================================================
	// TCP Communication Functions
	// ========================================================================

	/**
	 * @brief Initialize TCP communication on root process.
	 * 
	 * Establishes TCP connection on root process if remote communication is enabled.
	 * Centralizes network communication through the root process.
	 * 
	 * @param tcp_connection Reference to TcpConnection object to initialize
	 * @param from_cl Reference to FromCL structure with server configuration
	 */
	void init_communication(TcpConnection& tcp_connection, space_converter::FromCL& from_cl);

	/**
	 * @brief Close TCP communication on root process.
	 * 
	 * Cleanly shuts down the TCP connection established by init_communication().
	 * 
	 * @param tcp_connection Reference to TcpConnection object to close
	 * @param from_cl Reference to FromCL structure with server configuration
	 */
	void close_communication(TcpConnection& tcp_connection, space_converter::FromCL& from_cl);
    // @param tcp_connection: TCP connection object to listen on.
    // @param from_cl: Reference to the FromCL struct for communication configuration.
    // @param space_data: Reference to the common::SpaceData object to populate based on received message.
    void wait_on_message(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data);

    // Send informational data via the TCP connection.
    // @param tcp_connection: TCP connection object to send data through.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    // @param space_data: Reference to the common::SpaceData object containing data to send.
    // @param particle_data_types: Description of particle data types to include in the message.
    void send_info(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, std::string particle_data_types);

    // Receive requested data via the TCP connection.
    // @param tcp_connection: TCP connection object to receive data from.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    // @param space_data: Reference to the common::SpaceData object to populate with received data.
    void recv_requested_data(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data);

	// ========================================================================
	// VDB Data Transmission Functions
	// ========================================================================

	/**
	 * @brief Send VDB volumetric data to remote client.
	 * 
	 * Transmits processed volumetric data in OpenVDB, NanoVDB, or raw particle
	 * format, including metadata about value ranges and animation information.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData containing metadata
	 * @param grid_handle Reference to VDBParticles containing grid data
	 */
	void send_vdb(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_handle);

	/**
	 * @brief Send VDB file path reference to remote client.
	 * 
	 * Transmits only the file path instead of entire VDB data. Useful when
	 * client has direct file system access. Includes value range metadata.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData with file path and metadata
	 * @param grid_handle Reference to VDBParticles for type checking
	 */
	void send_path(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_handle);

	// ========================================================================
	// Debugging and Utility Functions
	// ========================================================================

	/**
	 * @brief Wait for debugger attachment (debugging utility).
	 * 
	 * Creates a waiting loop allowing debugger attachment. Root process waits
	 * until 'attach' variable is set to 0 in debugger. Useful for debugging MPI apps.
	 * 
	 * @param from_cl Reference to FromCL structure with configuration
	 * 
	 * @note Set the 'attach' variable to 0 in debugger to continue execution.
	 */
	void wait_on_attach_process(space_converter::FromCL& from_cl);

	/**
	 * @brief Receive bounding box request parameters from remote client.
	 * 
	 * Root process receives particle type and data block to query for bounding
	 * box information, then broadcasts to all MPI processes.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData to store requested parameters
	 */
	void recv_requested_bbox(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data);
	
	/**
	 * @brief Send bounding box information to remote client.
	 * 
	 * Calculates and transmits normalized bounding box coordinates based on
	 * original data extent, transformed for grid resolution and visualization.
	 * 
	 * @param tcp_connection Reference to TcpConnection object
	 * @param from_cl Reference to FromCL structure with configuration
	 * @param space_data Reference to SpaceData with overall dataset bounds
	 * @param space_data_bbox Reference to SpaceData with specific bounding box
	 * @param convert_vdb_base Pointer to VDB converter
	 */
	void send_bbox(TcpConnection& tcp_connection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::SpaceData& space_data_bbox, common::vdb::ConvertVDBBase* convert_vdb_base);

	/**
	 * @brief Read large test file for benchmarking or testing.
	 * 
	 * Reads a large file specified by CONVERTER_LARGE_TEST_FILE environment variable.
	 * All MPI processes synchronize before and after reading. Useful for testing
	 * file system performance and MPI synchronization.
	 * 
	 * @note Set CONVERTER_LARGE_TEST_FILE environment variable to enable.
	 */
	void read_large_test_file();

} // namespace space_converter


