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

#include "changa_nchilada_convert_vdb.h"
#include "changa_nchilada_extract_iolib.h"

namespace changa {
	void ConvertVDBChangaNchilada::print_CPU_steps() {
		changa::nchilada::io::print_CPU_steps();
	}

	float ConvertVDBChangaNchilada::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return changa::nchilada::io::get_particle_norm_value(blocknr, id);
	}

	int ConvertVDBChangaNchilada::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return changa::nchilada::io::get_particle_value(blocknr, id, value);
	}
	
	int ConvertVDBChangaNchilada::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return changa::nchilada::io::get_particle_value_comp(blocknr, id);
	}

	int ConvertVDBChangaNchilada::get_particle_type(uint64_t id) {
		return changa::nchilada::io::get_particle_type(id);
	}

	void ConvertVDBChangaNchilada::get_particle_position(uint64_t id, double* pos) const {
		changa::nchilada::io::get_particle_position(id, pos);
	}

	size_t ConvertVDBChangaNchilada::get_local_num_particles() const {
		return changa::nchilada::io::get_local_num_particles();
	}

	size_t ConvertVDBChangaNchilada::get_global_num_particles() const {
		return changa::nchilada::io::get_global_num_particles();
	}

	//double ConvertVDBChangaNchilada::get_particle_radius(uint64_t id) {
	//	return changa::nchilada::io::get_particle_radius(id);
	//}
	double ConvertVDBChangaNchilada::get_particle_hsml(uint64_t id) {
		return changa::nchilada::io::get_particle_hsml(id);
	}
	double ConvertVDBChangaNchilada::get_particle_mass(uint64_t id) {
		return changa::nchilada::io::get_particle_mass(id);
	}
	double ConvertVDBChangaNchilada::get_particle_rho_internal(uint64_t id) {
		return changa::nchilada::io::get_particle_rho(id);
	}
	int ConvertVDBChangaNchilada::get_particle_rho_blocknr() {
		return changa::nchilada::io::get_particle_rho_blocknr();
	}

	void ConvertVDBChangaNchilada::init_lib(int argc, char** argv, int world_rank, int world_size) {
		std::string base_dir;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--nc-dir") {
				base_dir = argv[++i]; //std::stof(av[++i]);
			}
		}

		changa::nchilada::io::init_lib(base_dir, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBChangaNchilada::finish_lib() {
		changa::nchilada::io::finish_lib();
	}

	void ConvertVDBChangaNchilada::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		changa::nchilada::io::get_types_and_blocks(types_and_blocks);
	}

	void ConvertVDBChangaNchilada::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		changa::nchilada::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		//char buf[500];
		for (int type = 0; type < changa::nchilada::ParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < changa::nchilada::BlockType::BTMax; blocknr++) {

				if (types_and_blocks[changa::nchilada::ParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBChangaNchilada::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		//char buf[500];
		for (int type = 0; type < changa::nchilada::ParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < changa::nchilada::BlockType::BTMax; blocknr++) {

				if (types_and_blocks[changa::nchilada::ParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBChangaNchilada::get_type_name(int type) {

		changa::nchilada::ParticleType pt = (changa::nchilada::ParticleType)type;
		switch (pt) {
		case changa::nchilada::ParticleType::Gas:
			return "Gas";
		case changa::nchilada::ParticleType::Dark:
			return "Dark";
		case changa::nchilada::ParticleType::Stars:
			return "Stars";
		}
		return "Unknown";
	}

	std::string ConvertVDBChangaNchilada::get_dataset_name(int blocknr) {

		changa::nchilada::BlockType bt = (changa::nchilada::BlockType)blocknr;
		switch (bt) {
		case changa::nchilada::BlockType::Pos:
			return "Pos";
		case changa::nchilada::BlockType::Mass:
			return "Mass";
		case changa::nchilada::BlockType::Vel:
			return "Vel";
		case changa::nchilada::BlockType::Soft:
			return "Soft";
		case changa::nchilada::BlockType::Phi:
			return "Phi";
		case changa::nchilada::BlockType::HSmooth:
			return "HSmooth";
		case changa::nchilada::BlockType::Rho:
			return "Rho";
		case changa::nchilada::BlockType::Temp:
			return "Temp";
		case changa::nchilada::BlockType::MetalsOx:
			return "MetalsOx";
		case changa::nchilada::BlockType::MetalsFe:
			return "MetalsFe";
		case changa::nchilada::BlockType::TForm:
			return "TForm";

		}
		return "Unknown";
	}

	std::string ConvertVDBChangaNchilada::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < changa::nchilada::ParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < changa::nchilada::BlockType::BTMax; bnr++) {
				if (types_and_blocks[changa::nchilada::ParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBChangaNchilada::get_num_types() {
		return changa::nchilada::ParticleType::PTMax;
	}
	int ConvertVDBChangaNchilada::get_num_blocks() {
		return changa::nchilada::BlockType::BTMax;
	}
}