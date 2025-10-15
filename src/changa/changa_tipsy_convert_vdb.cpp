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

#include "changa_tipsy_convert_vdb.h"
#include "changa_tipsy_extract_iolib.h"

namespace changa {
	void ConvertVDBChangaTipsy::print_CPU_steps() {
		changa::tipsy::io::print_CPU_steps();
	}

	float ConvertVDBChangaTipsy::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return changa::tipsy::io::get_particle_norm_value(blocknr, id);
	}

	int ConvertVDBChangaTipsy::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return changa::tipsy::io::get_particle_value(blocknr, id, value);
	}

	int ConvertVDBChangaTipsy::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return changa::tipsy::io::get_particle_value_comp(blocknr, id);
	}

	int ConvertVDBChangaTipsy::get_particle_type(uint64_t id) {
		return changa::tipsy::io::get_particle_type(id);
	}

	void ConvertVDBChangaTipsy::get_particle_position(uint64_t id, double* pos) const {
		changa::tipsy::io::get_particle_position(id, pos);
	}

	size_t ConvertVDBChangaTipsy::get_local_num_particles() const {
		return changa::tipsy::io::get_local_num_particles();
	}

	size_t ConvertVDBChangaTipsy::get_global_num_particles() const {
		return changa::tipsy::io::get_global_num_particles();
	}

	//double ConvertVDBChangaTipsy::get_particle_radius(uint64_t id) {
	//	return changa::tipsy::io::get_particle_radius(id);
	//}
	double ConvertVDBChangaTipsy::get_particle_hsml(uint64_t id) {
		return changa::tipsy::io::get_particle_hsml(id);
	}
	double ConvertVDBChangaTipsy::get_particle_mass(uint64_t id) {
		return changa::tipsy::io::get_particle_mass(id);
	}
	double ConvertVDBChangaTipsy::get_particle_rho_internal(uint64_t id) {
		return changa::tipsy::io::get_particle_rho(id);
	}
	int ConvertVDBChangaTipsy::get_particle_rho_blocknr() {
		return changa::tipsy::io::get_particle_rho_blocknr();
	}

	void ConvertVDBChangaTipsy::init_lib(int argc, char** argv, int world_rank, int world_size) {
		std::string base_file;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--tipsy-file") {
				base_file = argv[++i]; //std::stof(av[++i]);
			}
		}

		changa::tipsy::io::init_lib(base_file, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBChangaTipsy::finish_lib() {
		changa::tipsy::io::finish_lib();
	}

	void ConvertVDBChangaTipsy::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		changa::tipsy::io::get_types_and_blocks(types_and_blocks);
	}

	void ConvertVDBChangaTipsy::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		changa::tipsy::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		//char buf[500];
		for (int type = 0; type < changa::tipsy::ParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / changa::tipsy::ParticleType::PTMax; blocknr++) {

				if (types_and_blocks[changa::tipsy::ParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBChangaTipsy::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		//char buf[500];
		for (int type = 0; type < changa::tipsy::ParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / changa::tipsy::ParticleType::PTMax; blocknr++) {

				if (types_and_blocks[changa::tipsy::ParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBChangaTipsy::get_type_name(int type) {

		changa::tipsy::ParticleType pt = (changa::tipsy::ParticleType)type;
		switch (pt) {
		case changa::tipsy::ParticleType::Gas:
			return "Gas";
		case changa::tipsy::ParticleType::Dark:
			return "Dark";
		case changa::tipsy::ParticleType::Stars:
			return "Stars";
		}
		return "Unknown";
	}

	std::string ConvertVDBChangaTipsy::get_dataset_name(int blocknr) {
		return changa::tipsy::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBChangaTipsy::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < changa::tipsy::ParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / changa::tipsy::ParticleType::PTMax; bnr++) {
				if (types_and_blocks[changa::tipsy::ParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBChangaTipsy::get_num_types() {
		return changa::tipsy::ParticleType::PTMax;
	}
	int ConvertVDBChangaTipsy::get_num_blocks() {

		std::vector<int> types_and_blocks;
		changa::tipsy::io::get_types_and_blocks(types_and_blocks);
		return types_and_blocks.size() / changa::tipsy::ParticleType::PTMax;
	}
}