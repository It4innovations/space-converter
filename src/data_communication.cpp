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
#include "data_common.h"

#ifdef WITH_OPENMP
# include <omp.h>
#endif

#include <iostream>

 // file_type_items
#define	FTI_NONE 0
#define	FTI_OPENVDB 1
#define	FTI_NANOVDB 2
#define	FTI_PATH 3
#define	FTI_RAW_PART 4

namespace space_converter {

	//MPI_Send(&ns, 1, MPI_INT, world_rank - step, 0, MPI_COMM_WORLD);
	void mpi_send(void* data, size_t size, MPI_Datatype data_type, int from, int to)
	{
		const size_t unit_giga = 1000L * 1000L * 1000L * 2L;

		size_t size_sended = 0;
		while (size - size_sended > unit_giga) {
			MPI_Send((char*)data + size_sended, unit_giga, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD);
			size_sended += unit_giga;
		}

		MPI_Send((char*)data + size_sended, size - size_sended, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD);
	}

	//MPI_Recv(&ns, 1, MPI_INT, world_rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	void mpi_recv(void* data, size_t size, MPI_Datatype data_type, int from, int to)
	{
		const size_t unit_giga = 1000L * 1000L * 1000L * 2L;

		size_t size_sended = 0;
		while (size - size_sended > unit_giga) {
			MPI_Recv((char*)data + size_sended, unit_giga, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			size_sended += unit_giga;
		}

		MPI_Recv((char*)data + size_sended, size - size_sended, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	void mpi_reduce(float* ldata, float* gdata, size_t size)
	{
		const size_t unit_giga = 1000L * 1000L * 1000L * 2L;

		size_t size_sended = 0;
		while (size - size_sended > unit_giga) {
			//MPI_Send((char*)data + size_sended, unit_giga, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD);
			MPI_Reduce(ldata + size_sended, gdata + size_sended, unit_giga, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			size_sended += unit_giga;
		}

		//MPI_Send((char*)data + size_sended, size - size_sended, (MPI_Datatype)data_type, from, to, MPI_COMM_WORLD);
		MPI_Reduce(ldata + size_sended, gdata + size_sended, size - size_sended, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	double t_start = 0;
	void init_mpi(int argc, char** argv, space_converter::FromCL& from_cl)
	{
		MPI_Init(&argc, &argv);

#ifdef WITH_OPENMP
		t_start = omp_get_wtime();
#endif

		// Get the number of processes
		MPI_Comm_size(MPI_COMM_WORLD, &from_cl.world_size);

		// Get the rank of the process
		MPI_Comm_rank(MPI_COMM_WORLD, &from_cl.world_rank);

		// Get the name of the processor
		char processor_name[MPI_MAX_PROCESSOR_NAME];
		int name_len;
		MPI_Get_processor_name(processor_name, &name_len);

		if (from_cl.world_rank < 1) {
			// Print a hello mpi message
			printf(
				"Start from processor %s, rank %d"
				" out of %d processors\n",
				processor_name,
				from_cl.world_rank,
				from_cl.world_size);
		}
	}

	void close_mpi(space_converter::FromCL& from_cl)
	{
#ifdef WITH_OPENMP
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

		MPI_Finalize();
	}

	void init_communication(TcpConnection& tcpConnection, space_converter::FromCL& from_cl)
	{
		if (from_cl.world_rank == 0 && from_cl.remote) {
			tcpConnection.init_sockets_data(from_cl.server.c_str(), from_cl.port);
		}
	}

	void close_communication(TcpConnection& tcpConnection, space_converter::FromCL& from_cl)
	{
		if (from_cl.world_rank == 0 && from_cl.remote) {
			tcpConnection.client_close();
			tcpConnection.server_close();
		}
	}

	void wait_on_message(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				tcpConnection.recv_data_data((char*)&space_data.message_type, sizeof(space_data.message_type));

				if (tcpConnection.is_error()) {
					space_data.message_type = common::SpaceData::MessageType::eExit;
				}
			}

			printf("message_type: %d\n", space_data.message_type); fflush(0);
		}

		MPI_Bcast((char*)&space_data.message_type, sizeof(space_data.message_type), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	void send_info(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data, std::string particle_data_types)
	{
		if (from_cl.world_rank == 0) {
			//std::string particle_data_types = convert_vdb_base->get_particle_data_type_names(types_and_blocks_global);

			int anim_type = (int)space_data.anim_type;
			tcpConnection.send_data_data((char*)&anim_type, sizeof(int));
			tcpConnection.send_data_data((char*)&space_data.anim_start, sizeof(int));
			tcpConnection.send_data_data((char*)&space_data.anim_end, sizeof(int));

			int s = particle_data_types.length();
			tcpConnection.send_data_data((char*)&s, sizeof(int));
			tcpConnection.send_data_data((char*)particle_data_types.c_str(), sizeof(char) * s);

			int ack;
			tcpConnection.recv_data_data((char*)&ack, sizeof(ack));
			printf("sended: particle and data types\n");
		}
	}
	void recv_requested_bbox(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				tcpConnection.recv_data_data((char*)&space_data.particle_type, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.block_name_id, sizeof(int));
			}
		}

		MPI_Bcast((char*)&space_data.particle_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.block_name_id, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	void send_bbox(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::SpaceData& spaceDataBBox, common::vdb::ConvertVDBBase* convert_vdb_base)
	{
		//TODO
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				float bbox_min[3] = { 0.0f, 0.0f, 0.0f };
				float bbox_max[3] = { 0.0f, 0.0f, 0.0f };

				float bbox_dim = (float)space_data.object_size;

				float dims_orig[3];
				//float dims_bbox[3];

				dims_orig[0] = space_data.bbox_max_orig[0] - space_data.bbox_min_orig[0];
				dims_orig[1] = space_data.bbox_max_orig[1] - space_data.bbox_min_orig[1];
				dims_orig[2] = space_data.bbox_max_orig[2] - space_data.bbox_min_orig[2];

				float dims_orig_max = std::max(dims_orig[0], std::max(dims_orig[1], dims_orig[2]));

				//dims_bbox[0] = spaceDataBBox.bbox_max_orig[0] - spaceDataBBox.bbox_min_orig[0];
				//dims_bbox[1] = spaceDataBBox.bbox_max_orig[1] - spaceDataBBox.bbox_min_orig[1];
				//dims_bbox[2] = spaceDataBBox.bbox_max_orig[2] - spaceDataBBox.bbox_min_orig[2];

				bbox_min[0] = (spaceDataBBox.bbox_min_orig[0] - space_data.bbox_min_orig[0]) * bbox_dim / dims_orig[0];
				bbox_min[1] = (spaceDataBBox.bbox_min_orig[1] - space_data.bbox_min_orig[1]) * bbox_dim / dims_orig[1];
				bbox_min[2] = (spaceDataBBox.bbox_min_orig[2] - space_data.bbox_min_orig[2]) * bbox_dim / dims_orig[2];

				bbox_max[0] = (spaceDataBBox.bbox_max_orig[0] - space_data.bbox_min_orig[0]) * bbox_dim / dims_orig[0];
				bbox_max[1] = (spaceDataBBox.bbox_max_orig[1] - space_data.bbox_min_orig[1]) * bbox_dim / dims_orig[1];
				bbox_max[2] = (spaceDataBBox.bbox_max_orig[2] - space_data.bbox_min_orig[2]) * bbox_dim / dims_orig[2];

				// send bbox
				tcpConnection.send_data_data((char*)&bbox_min[0], sizeof(float) * 3);
				tcpConnection.send_data_data((char*)&bbox_max[0], sizeof(float) * 3);

				int ack;
				tcpConnection.recv_data_data((char*)&ack, sizeof(ack));
				printf("sended: BBOX\n");
			}
		}
	}

	void recv_requested_data(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		if (from_cl.world_rank == 0) {
			if (from_cl.remote) {
				tcpConnection.recv_data_data((char*)&space_data.bbox_min[0], sizeof(float) * 3);
				tcpConnection.recv_data_data((char*)&space_data.bbox_max[0], sizeof(float) * 3);
				tcpConnection.recv_data_data((char*)&space_data.bbox_dim, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.grid_transform, sizeof(float));
				tcpConnection.recv_data_data((char*)&space_data.particle_type, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.block_name_id, sizeof(int));
				//tcpConnection.recv_data_data((char*)&space_data.value_convert, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.extracted_type, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.dense_type, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.dense_norm, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.object_size, sizeof(float));
				tcpConnection.recv_data_data((char*)&space_data.particle_fix_size, sizeof(float));
				tcpConnection.recv_data_data((char*)&space_data.filter_min, sizeof(float));
				tcpConnection.recv_data_data((char*)&space_data.filter_max, sizeof(float));
				tcpConnection.recv_data_data((char*)&space_data.frame, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.anim_type, sizeof(int));
				tcpConnection.recv_data_data((char*)&space_data.anim_task_counter, sizeof(int));
			}
		}

		MPI_Bcast((char*)&space_data.bbox_min[0], sizeof(float) * 3, MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.bbox_max[0], sizeof(float) * 3, MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.bbox_dim, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.grid_transform, sizeof(float), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.particle_type, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		MPI_Bcast((char*)&space_data.block_name_id, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
		//MPI_Bcast((char*)&space_data.value_convert, sizeof(int), MPI_BYTE, 0, MPI_COMM_WORLD);
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

	void send_vdb(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_handle)
	{
		if (from_cl.world_rank == 0 && grid_handle.type == common::vdb::VDBParticles::VDBParticleType::eVector) {
			double t = omp_get_wtime();
			//send vdb        
			std::size_t size = grid_handle.vector_grid.size();

			// file type
			int file_type = FTI_OPENVDB;
			// Particle
			if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
				file_type = FTI_RAW_PART;				
			}
			else if (from_cl.use_nanovdb) {
				file_type = FTI_NANOVDB;
			}			

			tcpConnection.send_data_data((char*)&file_type, sizeof(file_type));

			//vdb file
			tcpConnection.send_data_data((char*)&size, sizeof(size));
			tcpConnection.send_data_data((char*)grid_handle.vector_grid.data(), size);

			//vdb info
			tcpConnection.send_data_data((char*)&space_data.min_value, sizeof(space_data.min_value));
			tcpConnection.send_data_data((char*)&space_data.max_value, sizeof(space_data.max_value));
			tcpConnection.send_data_data((char*)&space_data.min_value_reduced, sizeof(space_data.min_value_reduced));
			tcpConnection.send_data_data((char*)&space_data.max_value_reduced, sizeof(space_data.max_value_reduced));
						
			//anim
			int frames = 1;
			if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
				frames = from_cl.world_size;
			}
			tcpConnection.send_data_data((char*)&frames, sizeof(frames));

			int ack;
			tcpConnection.recv_data_data((char*)&ack, sizeof(ack));
			printf("rank: %d: sendOpenVDB time: %f\n", from_cl.world_rank, omp_get_wtime() - t);
			printf("sended: vdb\n");
		}
	}

	void send_path(TcpConnection& tcpConnection, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_handle)
	{
		if (from_cl.world_rank == 0 && grid_handle.type == common::vdb::VDBParticles::VDBParticleType::eVector) {
			double t = omp_get_wtime();

			// file type
			int file_type = FTI_PATH;
			tcpConnection.send_data_data((char*)&file_type, sizeof(file_type));

			//vdb file path
			std::string full_filepath = space_data.full_filepath;
			const char* cf = full_filepath.c_str();

			size_t size = strlen(cf);
			tcpConnection.send_data_data((char*)&size, sizeof(size));
			tcpConnection.send_data_data((char*)cf, size);

			//vdb info
			tcpConnection.send_data_data((char*)&space_data.min_value, sizeof(space_data.min_value));
			tcpConnection.send_data_data((char*)&space_data.max_value, sizeof(space_data.max_value));
			tcpConnection.send_data_data((char*)&space_data.min_value_reduced, sizeof(space_data.min_value_reduced));
			tcpConnection.send_data_data((char*)&space_data.max_value_reduced, sizeof(space_data.max_value_reduced));

			//anim
			int frames = 1;
			if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
				frames = from_cl.world_size;
			}
			tcpConnection.send_data_data((char*)&frames, sizeof(frames));

			int ack;
			tcpConnection.recv_data_data((char*)&ack, sizeof(ack));
			printf("rank: %d: sendOpenVDB time: %f\n", from_cl.world_rank, omp_get_wtime() - t);
			printf("sended: vdb\n");
		}
	}

	void wait_on_attach_process(space_converter::FromCL& from_cl)
	{
		if (from_cl.world_rank == 0) {
			int attach = 1;			
			while (attach) {
#ifdef _WIN32
				Sleep(1);
#endif
				printf("wait_on_attach_process\n"); fflush(0);
			}			
		}

		int attach_process = 1;
		MPI_Bcast((char*)&attach_process, sizeof(attach_process), MPI_BYTE, 0, MPI_COMM_WORLD);
	}

	void read_large_test_file() {
		const char* filename = getenv("CONVERTER_LARGE_TEST_FILE");
		if (filename) {
			MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes are synchronized before reading the file

			std::ifstream file(filename, std::ios::binary | std::ios::ate); // open at end to get size
			if (!file) {
				std::cerr << "Failed to open file: " << filename << std::endl;
				exit(-1);
			}

			std::streamsize size = file.tellg(); // get size from end position
			file.seekg(0, std::ios::beg);        // rewind to beginning

			char* buffer = new char[size];
			if (!file.read(buffer, size)) {
				std::cerr << "Failed to read file: " << filename << std::endl;
				exit(-1);
			}

			std::cout << "Read " << size << " bytes from " << filename << std::endl;

			delete [] buffer;
			MPI_Barrier(MPI_COMM_WORLD); // Ensure all processes are synchronized before reading the file
		}
	}

} //space_converter