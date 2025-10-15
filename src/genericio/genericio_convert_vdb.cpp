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

#include "genericio_convert_vdb.h"
#include "genericio_extract_iolib.h"

namespace genericio {
	void ConvertVDBGenericIO::print_CPU_steps() {
		genericio::io::print_CPU_steps();
	}

	float ConvertVDBGenericIO::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return genericio::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBGenericIO::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return genericio::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBGenericIO::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return genericio::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBGenericIO::get_particle_type(uint64_t id) {
		return genericio::io::get_particle_type(id);
	}
	void ConvertVDBGenericIO::get_particle_position(uint64_t id, double* pos) const {
		genericio::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBGenericIO::get_local_num_particles() const {
		return genericio::io::get_local_num_particles();
	}
	size_t ConvertVDBGenericIO::get_global_num_particles() const {
		return genericio::io::get_global_num_particles();
	}

	//double ConvertVDBGenericIO::get_particle_radius(uint64_t id) {
	//	return genericio::io::get_particle_radius(id);
	//}

	double ConvertVDBGenericIO::get_particle_hsml(uint64_t id) {
		return genericio::io::get_particle_hsml(id);
	}

	double ConvertVDBGenericIO::get_particle_mass(uint64_t id) {
		return genericio::io::get_particle_mass(id);
	}

	double ConvertVDBGenericIO::get_particle_rho_internal(uint64_t id) {
		return genericio::io::get_particle_rho(id);
	}

	int ConvertVDBGenericIO::get_particle_rho_blocknr() {
		return genericio::io::get_particle_rho_blocknr();
	}

	void ConvertVDBGenericIO::init_lib(int argc, char** argv, int world_rank, int world_size) {	
		std::string genericio_file;
		// std::string pos_names = "x,y,z";
		// std::string vel_names = "vx,vy,vz";
		std::vector<std::string> pos_names_vec;
		std::vector<std::string> vel_names_vec;

		std::string vel_name_mass;
		std::string vel_name_rho;
		std::string vel_name_hsml;		

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--genericio-file") {
				genericio_file = argv[++i];
			}
			else if (arg == "--pos-names") {
				pos_names_vec.push_back(argv[++i]);
				pos_names_vec.push_back(argv[++i]);
				pos_names_vec.push_back(argv[++i]);
			}
			else if (arg == "--vel-names") {
				vel_names_vec.push_back(argv[++i]);
				vel_names_vec.push_back(argv[++i]);
				vel_names_vec.push_back(argv[++i]);
			}
			else if (arg == "--mass-name") {
				vel_name_mass = argv[++i];
			}
			else if (arg == "--rho-name") {
				vel_name_rho = argv[++i];
			}
			else if (arg == "--hsml-name") {
				vel_name_hsml = argv[++i];
			}			
		}

		genericio::io::init_lib(genericio_file, world_rank, world_size, pos_names_vec, vel_names_vec, vel_name_mass, vel_name_rho, vel_name_hsml);

		print_CPU_steps();
	}

	void ConvertVDBGenericIO::finish_lib()
	{
		genericio::io::finish_lib();
	}

	void ConvertVDBGenericIO::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		genericio::io::get_types_and_blocks(types_and_blocks);
	}
	void ConvertVDBGenericIO::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		genericio::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		//char buf[500];
		for (int type = 0; type < genericio::GenericIOParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / genericio::GenericIOParticleType::PTMax; blocknr++) {

				if (types_and_blocks[genericio::GenericIOParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBGenericIO::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		//char buf[500];
		for (int type = 0; type < genericio::GenericIOParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / genericio::GenericIOParticleType::PTMax; blocknr++) {

				if (types_and_blocks[genericio::GenericIOParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBGenericIO::get_type_name(int type) {		

		genericio::GenericIOParticleType pt = (genericio::GenericIOParticleType)type;
		switch (pt) {
		case genericio::GenericIOParticleType::HACC:
			return "HACC";
		//case genericio::GenericIOParticleType::Dark:
		//	return "Dark";
		//case genericio::GenericIOParticleType::Stars:
		//	return "Stars";
		}
		return "Unknown";
	}
	std::string ConvertVDBGenericIO::get_dataset_name(int blocknr) {
		return genericio::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBGenericIO::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < genericio::GenericIOParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / genericio::GenericIOParticleType::PTMax; bnr++) {
				if (types_and_blocks[genericio::GenericIOParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBGenericIO::get_num_types() {
		return genericio::GenericIOParticleType::PTMax;
	}
	int ConvertVDBGenericIO::get_num_blocks() {
		std::vector<int> types_and_blocks;
		genericio::io::get_types_and_blocks(types_and_blocks);
		return types_and_blocks.size() / genericio::GenericIOParticleType::PTMax;
	}
}