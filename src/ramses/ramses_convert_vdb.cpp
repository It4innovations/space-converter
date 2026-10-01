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

#include "ramses_convert_vdb.h"
#include "ramses_extract_iolib.h"

#include <iostream>

namespace ramses {
	void ConvertVDBRamses::print_CPU_steps() {
		ramses::io::print_CPU_steps();
	}

	float ConvertVDBRamses::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return ramses::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBRamses::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return ramses::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBRamses::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return ramses::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBRamses::get_particle_type(uint64_t id) {
		return ramses::io::get_particle_type(id);
	}
	void ConvertVDBRamses::get_particle_position(uint64_t id, double* pos) const {
		ramses::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBRamses::get_local_num_particles() const {
		return ramses::io::get_local_num_particles();
	}
	size_t ConvertVDBRamses::get_global_num_particles() const {
		return ramses::io::get_global_num_particles();
	}

	double ConvertVDBRamses::get_particle_hsml(uint64_t id) {
		return ramses::io::get_particle_hsml(id);
	}

	double ConvertVDBRamses::get_particle_mass(uint64_t id) {
		return ramses::io::get_particle_mass(id);
	}

	double ConvertVDBRamses::get_particle_rho_internal(uint64_t id) {
		return ramses::io::get_particle_rho(id);
	}

	int ConvertVDBRamses::get_particle_rho_blocknr() {
		return ramses::io::get_particle_rho_blocknr();
	}

	void ConvertVDBRamses::init_lib(int argc, char** argv, int world_rank, int world_size) {
		std::string output_dir;
		bool read_gas = true;
		bool read_particles = true;
		int level_max = 0;
		std::vector<ramses::io::FieldFilter> filters;

		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;
		int anim_step = -1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--ramses-output") {
				output_dir = argv[++i];
			}
			else if (arg == "--no-gas") {
				read_gas = false;
			}
			else if (arg == "--no-particles") {
				read_particles = false;
			}
			else if (arg == "--ramses-levelmax") {
				level_max = std::stoi(argv[++i]);
			}
			else if (arg == "--ramses-filter") {
				ramses::io::FieldFilter f;
				f.name = argv[++i];
				f.min = std::stod(argv[++i]);
				f.max = std::stod(argv[++i]);
				filters.push_back(f);
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
				anim_step = std::stoi(argv[++i]);
			}
		}

		// Anim: one output per rank (the pattern holds the output number)
		if (use_anim) {
			// Clamp to anim_end so trailing ranks do not address outputs past the animation range
			int anim_frame = anim_start + anim_step * world_rank;
			if (anim_frame > anim_end)
				anim_frame = anim_end;
			output_dir = format_filename(output_dir, anim_frame);
			std::cout << "Reading RAMSES output: " << output_dir << std::endl;
			world_rank = 0;
			world_size = 1;
		}

		ramses::io::init_lib(output_dir, world_rank, world_size, read_gas, read_particles, level_max, filters);

		print_CPU_steps();
	}

	void ConvertVDBRamses::finish_lib()
	{
		ramses::io::finish_lib();
	}

	void ConvertVDBRamses::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		ramses::io::get_types_and_blocks(types_and_blocks);
	}

	void ConvertVDBRamses::print_types_and_blocks_local() {
		ramses::io::print_types_and_blocks_local();
	}

	void ConvertVDBRamses::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");
		ramses::io::print_types_and_blocks(types_and_blocks);
	}

	std::string ConvertVDBRamses::get_type_name(int type) {
		switch ((ramses::RamsesParticleType)type) {
		case ramses::RamsesParticleType::Gas: return "Gas";
		case ramses::RamsesParticleType::DM: return "DM";
		case ramses::RamsesParticleType::Star: return "Star";
		case ramses::RamsesParticleType::Cloud: return "Cloud";
		case ramses::RamsesParticleType::Debris: return "Debris";
		case ramses::RamsesParticleType::Other: return "Other";
		case ramses::RamsesParticleType::Sink: return "Sink";
		default: break;
		}
		return "Unknown";
	}

	std::string ConvertVDBRamses::get_dataset_name(int blocknr) {
		return ramses::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBRamses::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < ramses::RamsesParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < (int)types_and_blocks.size() / ramses::RamsesParticleType::PTMax; bnr++) {
				if (types_and_blocks[ramses::RamsesParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBRamses::get_num_types() {
		return ramses::RamsesParticleType::PTMax;
	}

	int ConvertVDBRamses::get_num_blocks() {
		std::vector<int> types_and_blocks;
		ramses::io::get_types_and_blocks(types_and_blocks);
		return (int)types_and_blocks.size() / ramses::RamsesParticleType::PTMax;
	}
}
