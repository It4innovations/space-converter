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

#include "data_communication.h"
#include "data_processing.h"
#include "args_processing.h"

int main(int argc, char** argv)
{
	// Initialize configuration from command-line arguments.
	space_converter::FromCL from_cl;
	space_converter::init_mpi(argc, argv, from_cl);  // Initialize MPI for parallel processing.

#ifdef _WIN32
	//setvbuf(stdout, NULL, _IONBF, 0);  // No buffering
	//setvbuf(stdout, NULL, _IOLBF, 0);  // Line buffering
	setbuf(stdout, NULL);
#endif

#if 1
	space_converter::read_large_test_file();
#endif

#if 0 //debug
	space_converter::wait_on_attach_process(from_cl);
#endif

#ifdef WITH_OPENVDB
	openvdb::initialize();  // Initialize OpenVDB if available.
#endif

	// Create a SpaceData object for managing current processing state.
	common::SpaceData space_data;

	// Parse command-line arguments to populate configuration.
	space_converter::parse_args(from_cl, space_data, argc, argv);

#if 0
	space_converter::test_converter(argc, argv, from_cl);
#endif

	// Initialize the VDB converter.
	common::vdb::ConvertVDBBase* convert_vdb_base = space_converter::init_converter(argc, argv, from_cl, space_data);

	// Print global particle types and blocks information.
	std::vector<int> types_and_blocks_global;
	space_converter::print_info(convert_vdb_base, from_cl, types_and_blocks_global);

	while (true) {
		// Create a SpaceData object for managing current processing state.
		//common::SpaceData space_data(
		//	from_cl.export_type, 
		//	from_cl.export_dataset, 
		//	from_cl.export_extracted_type,
		//	from_cl.export_dense_type,
		//	from_cl.export_dense_norm,
		//	from_cl.grid_dim, 
		//	from_cl.use_bbox, 
		//	from_cl.bbox_pos);

		//TODO
		//space_data.dense_type = 1;

		// Initialize TCP communication if in remote mode.
		TcpConnection tcp_connection;
		space_converter::init_communication(tcp_connection, from_cl);

#if 1 //fixing BBOX
		// Calculate the bounding box for the data.
		space_converter::find_bbox(convert_vdb_base, from_cl, space_data);
#endif

		while (true) {
			// Wait for a message from the TCP connection or control logic.
			space_converter::wait_on_message(tcp_connection, from_cl, space_data);

			// Handle exit message.
			if (from_cl.remote && space_data.message_type == common::SpaceData::MessageType::eExit) {
				break;
			}

			// Handle informational messages.
			if (from_cl.remote && space_data.message_type == common::SpaceData::MessageType::eInfo) {
				space_converter::send_info(
					tcp_connection,
					from_cl,
					space_data,
					convert_vdb_base->get_particle_data_type_names(types_and_blocks_global)
				);
			}

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

			// Handle data extraction and conversion.
			if (space_data.message_type == common::SpaceData::MessageType::eData || !from_cl.remote) {
				// Receive requested data (if applicable).
				space_converter::recv_requested_data(tcp_connection, from_cl, space_data);

#if 0 //fixing BBOX
				// Calculate the bounding box for the data.
				space_converter::find_bbox(convert_vdb_base, from_cl, space_data);
#endif
				// Create the VDB grid and convert particle data.
				common::vdb::VDBParticles grid_main;
				space_converter::create_grid(grid_main, from_cl, space_data);
				space_converter::convert_to_grid(convert_vdb_base, from_cl, space_data, grid_main);

				// Find the minimum and maximum values in the data.
				space_converter::find_minmax_value(from_cl, space_data);

				// Save file per rank
				if (from_cl.use_save_mpirank) {
					space_converter::save_vdb(convert_vdb_base, from_cl, space_data, grid_main, grid_main.type, false);
				}

				// Perform reduction to combine data across processes.
				common::vdb::VDBParticles grid_main_sum;
				space_converter::reduction(convert_vdb_base, from_cl, space_data, grid_main, grid_main_sum);

				// Finalize the grid with transformations and optimizations.
				common::vdb::VDBParticles grid_main_final;
				space_converter::finalize_grid(convert_vdb_base, from_cl, space_data, grid_main_sum, grid_main_final);

				// Find reduced min/max
				space_converter::find_minmax_reduced_value(from_cl, space_data);

				// Send or save the finalized VDB data.
				if ((space_data.anim_type == common::SpaceData::AnimType::eNone || space_data.anim_type == common::SpaceData::AnimType::eAllMerge || space_data.anim_type == common::SpaceData::AnimType::eFrameExtract) && from_cl.remote) {
					space_converter::send_vdb(tcp_connection, from_cl, space_data, grid_main_final);
				}
				else {
					space_converter::save_vdb(convert_vdb_base, from_cl, space_data, grid_main_final, grid_main_final.type);
					if ((space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) && from_cl.remote) {
						space_converter::send_path(tcp_connection, from_cl, space_data, grid_main_final);
					}

				}
			}

			// Exit inner loop if not in remote mode.
			if (!from_cl.remote) {
				break;
			}
		}

		// Close communication in remote mode; exit outer loop if not remote.
		if (from_cl.remote) {
			space_converter::close_communication(tcp_connection, from_cl);
		}
		else {
			break;
		}
	}

	// Clean up resources.
	space_converter::deinit_converter(convert_vdb_base);

#ifdef WITH_OPENVDB
	openvdb::uninitialize();  // Uninitialize OpenVDB if initialized.
#endif

	space_converter::close_mpi(from_cl);  // Finalize MPI.

	return 0;  // Exit the program.
}
