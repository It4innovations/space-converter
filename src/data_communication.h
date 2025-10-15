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

#pragma once

#include <mpi.h>
#include <vector>

#include "args_processing.h"
#include "convert_vdb.h"
#include "data_common.h"
#include "connection_tcp.h"

 // Namespace space_converter provides MPI communication and data exchange functionalities.
namespace space_converter {

    // Send data via MPI from one process to another.
    // @param data: Pointer to the data to send.
    // @param size: Size of the data in bytes.
    // @param data_type: MPI data type of the data.
    // @param from: Rank of the sending process.
    // @param to: Rank of the receiving process.
    void mpi_send(void* data, size_t size, MPI_Datatype data_type, int from, int to);

    // Receive data via MPI from one process to another.
    // @param data: Pointer to the buffer to store received data.
    // @param size: Size of the data in bytes.
    // @param data_type: MPI data type of the data.
    // @param from: Rank of the sending process.
    // @param to: Rank of the receiving process.
    void mpi_recv(void* data, size_t size, MPI_Datatype data_type, int from, int to);

    // Perform an MPI reduction operation to aggregate data from all processes.
    // @param ldata: Pointer to the local data to reduce.
    // @param gdata: Pointer to store the global reduced data.
    // @param size: Number of elements in the data arrays.
    void mpi_reduce(float* ldata, float* gdata, size_t size);

    // Initialize the MPI environment and configure the process with parsed command-line options.
    // @param argc: Number of command-line arguments.
    // @param argv: Array of command-line argument strings.
    // @param from_cl: Reference to the FromCL struct to populate with MPI-related configuration.
    void init_mpi(int argc, char** argv, space_converter::FromCL& from_cl);

    // Finalize the MPI environment and clean up resources.
    // @param from_cl: Reference to the FromCL struct for MPI-related configuration.
    void close_mpi(space_converter::FromCL& from_cl);

    // Initialize communication via TCP.
    // @param tcpConnection: TCP connection object to initialize.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    void init_communication(TcpConnection& tcpConnection, space_converter::FromCL& from_cl);

    // Close communication via TCP.
    // @param tcpConnection: TCP connection object to close.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    void close_communication(TcpConnection& tcpConnection, space_converter::FromCL& from_cl);

    // Wait for a message to be received from the TCP connection.
    // @param tcpConnection: TCP connection object to listen on.
    // @param from_cl: Reference to the FromCL struct for communication configuration.
    // @param spaceData: Reference to the common::SpaceData object to populate based on received message.
    void wait_on_message(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData);

    // Send informational data via the TCP connection.
    // @param tcpConnection: TCP connection object to send data through.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    // @param spaceData: Reference to the common::SpaceData object containing data to send.
    // @param particle_data_types: Description of particle data types to include in the message.
    void send_info(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData, std::string particle_data_types);

    // Receive requested data via the TCP connection.
    // @param tcpConnection: TCP connection object to receive data from.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    // @param spaceData: Reference to the common::SpaceData object to populate with received data.
    void recv_requested_data(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData);

    // Send VDB data via the TCP connection.
    // @param tcpConnection: TCP connection object to send data through.
    // @param from_cl: Reference to the FromCL struct with communication configuration.
    // @param spaceData: Reference to the common::SpaceData object containing metadata for the VDB data.
    // @param grid_handle: Reference to the VDBParticles grid handle to send.
    void send_vdb(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_handle);

    void send_path(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_handle);

    // Wait for an attach.
    // @param from_cl: Reference to the FromCL struct for communication configuration.
    void wait_on_attach_process(space_converter::FromCL& from_cl);

    void recv_requested_bbox(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData);
    void send_bbox(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& spaceData, common::SpaceData& spaceDataBBox, common::vdb::ConvertVDBBase* convert_vdb_base);

    void read_large_test_file();

} // namespace space_converter


