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

#include "data_communication.h"
#include "data_processing.h"
#include "args_processing.h"

/**
 * @brief Main entry point for the space converter application.
 * 
 * Orchestrates the entire conversion pipeline:
 * 1. Initializes MPI and parses command-line arguments
 * 2. Sets up VDB converter for the specified simulation format
 * 3. In remote mode: Enters server loop waiting for client requests
 * 4. In local mode: Processes data directly and saves output
 * 5. Handles data extraction, grid conversion, reduction, and output
 * 
 * Supports both standalone batch processing and remote client-server operation
 * for interactive visualization workflows.
 * 
 * @param argc Number of command-line arguments
 * @param argv Array of command-line argument strings
 * @return 0 on successful completion
 */
int main(int argc, char** argv)
{
	// ========================================================================
	// Initialization Phase
	// ========================================================================

	space_converter::FromCL from_cl;
	space_converter::init_mpi(argc, argv, from_cl);

#ifdef _WIN32
	// Disable stdout buffering on Windows for immediate output
	setbuf(stdout, NULL);
#endif

	// Read large test file if CONVERTER_LARGE_TEST_FILE environment variable is set
	space_converter::read_large_test_file();

#ifdef WITH_OPENVDB
	openvdb::initialize();
#endif

	common::SpaceData space_data;
	space_converter::parse_args(from_cl, space_data, argc, argv);

	// Initialize the VDB converter for the specified simulation format
	common::vdb::ConvertVDBBase* convert_vdb_base = space_converter::init_converter(argc, argv, from_cl, space_data);

	// Query and print available particle types and data blocks
	std::vector<int> types_and_blocks_global;
	space_converter::print_info(convert_vdb_base, from_cl, types_and_blocks_global);

	// ========================================================================
	// Main Processing Loop (Remote Server Mode)
	// ========================================================================

	while (true) {
		// Initialize TCP communication if in remote mode
		TcpConnection tcp_connection;
		space_converter::init_communication(tcp_connection, from_cl);

		// Calculate the bounding box for the entire dataset
		space_converter::find_bbox(convert_vdb_base, from_cl, space_data);

		// ====================================================================
		// Message Handling Loop (Process Client Requests)
		// ====================================================================

		while (true) {
			// Wait for a message from the TCP connection or proceed in local mode
			space_converter::wait_on_message(tcp_connection, from_cl, space_data);

			// Handle exit message from client
			if (from_cl.remote && space_data.message_type == common::SpaceData::MessageType::eExit) {
				break;
			}

			// Handle info request - send dataset metadata to client
			if (from_cl.remote && space_data.message_type == common::SpaceData::MessageType::eInfo) {
				space_converter::send_info(
					tcp_connection,
					from_cl,
					space_data,
					convert_vdb_base->get_particle_data_type_names(types_and_blocks_global)
				);
			}

			// Handle bounding box query for specific particle type
			if (from_cl.remote && space_data.message_type == common::SpaceData::MessageType::eBBOX) {
				common::SpaceData space_data_bbox(space_data);
				space_converter::recv_requested_bbox(
					tcp_connection,
					from_cl,
					space_data_bbox
				);
				space_converter::find_bbox(convert_vdb_base, from_cl, space_data_bbox, space_data_bbox.particle_type);
				space_converter::send_bbox(
					tcp_connection,
					from_cl,
					space_data,
					space_data_bbox,
					convert_vdb_base
				);
			}

			// ================================================================
			// Data Conversion Pipeline
			// ================================================================

			if (space_data.message_type == common::SpaceData::MessageType::eData || !from_cl.remote) {
				// Receive conversion parameters from client (particle type, resolution, etc.)
				space_converter::recv_requested_data(tcp_connection, from_cl, space_data);

				// Create empty VDB grid with specified dimensions and transform
				common::vdb::VDBParticles grid_main;
				space_converter::create_grid(grid_main, from_cl, space_data);
				
				// Convert raw particle data to volumetric grid using density estimation
				space_converter::convert_to_grid(convert_vdb_base, from_cl, space_data, grid_main);

				// Find local min/max values for this MPI rank
				space_converter::find_minmax_value(from_cl, space_data);

				// Save per-rank intermediate results if requested
				if (from_cl.use_save_mpirank) {
					space_converter::save_vdb(convert_vdb_base, from_cl, space_data, grid_main, grid_main.type, false);
				}

				// Combine grids from all MPI ranks using reduction operations
				common::vdb::VDBParticles grid_main_sum;
				space_converter::reduction(convert_vdb_base, from_cl, space_data, grid_main, grid_main_sum);

				// Apply final transformations, filtering, and optimizations
				common::vdb::VDBParticles grid_main_final;
				space_converter::finalize_grid(convert_vdb_base, from_cl, space_data, grid_main_sum, grid_main_final);

				// Find global min/max values after reduction
				space_converter::find_minmax_reduced_value(from_cl, space_data);

				// Output results - either send to client or save to file
				if ((space_data.anim_type == common::SpaceData::AnimType::eNone || 
				     space_data.anim_type == common::SpaceData::AnimType::eAllMerge || 
				     space_data.anim_type == common::SpaceData::AnimType::eFrameExtract) && from_cl.remote) {
					// Send VDB data directly to remote client
					space_converter::send_vdb(tcp_connection, from_cl, space_data, grid_main_final);
				}
				else {
					// Save VDB file to disk
					space_converter::save_vdb(convert_vdb_base, from_cl, space_data, grid_main_final, grid_main_final.type);
					
					// For animation sequences, send file path instead of full data
					if ((space_data.anim_type != common::SpaceData::AnimType::eNone && 
					     space_data.anim_type != common::SpaceData::AnimType::eAllMerge) && from_cl.remote) {
						space_converter::send_path(tcp_connection, from_cl, space_data, grid_main_final);
					}
				}
			}

			// Exit inner loop if not in remote mode (single batch processing)
			if (!from_cl.remote) {
				break;
			}
		}

		// Close TCP connection and exit outer loop based on mode
		if (from_cl.remote) {
			space_converter::close_communication(tcp_connection, from_cl);
		}
		else {
			break;
		}
	}

	// ========================================================================
	// Cleanup Phase
	// ========================================================================

	space_converter::deinit_converter(convert_vdb_base);

#ifdef WITH_OPENVDB
	openvdb::uninitialize();
#endif

	space_converter::close_mpi(from_cl);

	return 0;
}
