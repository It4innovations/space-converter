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

#include "fil_convert_vdb.h"
#include "fil_extract_iolib.h"

#include <iostream>
#include <sstream>
#include <stdexcept>

namespace fil {
	void ConvertVDBFil::print_CPU_steps() {
		fil::io::print_CPU_steps();
	}

	float ConvertVDBFil::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return fil::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBFil::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return fil::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBFil::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return fil::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBFil::get_particle_type(uint64_t id) {
		return fil::io::get_particle_type(id);
	}
	void ConvertVDBFil::get_particle_position(uint64_t id, double* pos) const {
		fil::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBFil::get_local_num_particles() const {
		return fil::io::get_local_num_particles();
	}
	size_t ConvertVDBFil::get_global_num_particles() const {
		return fil::io::get_global_num_particles();
	}

	double ConvertVDBFil::get_particle_hsml(uint64_t id) {
		return fil::io::get_particle_hsml(id);
	}

	double ConvertVDBFil::get_particle_mass(uint64_t id) {
		return fil::io::get_particle_mass(id);
	}

	double ConvertVDBFil::get_particle_rho_internal(uint64_t id) {
		return fil::io::get_particle_rho(id);
	}

	int ConvertVDBFil::get_particle_rho_blocknr() {
		return fil::io::get_particle_rho_blocknr();
	}

	void ConvertVDBFil::init_lib(int argc, char** argv, int world_rank, int world_size) {
		fil::io::Options options;

		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;
		int anim_step = -1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--fil-file") {
				options.files.push_back(argv[++i]);
			}
			else if (arg == "--fil-dir") {
				options.dirs.push_back(argv[++i]);
			}
			else if (arg == "--fil-vars") {
				// comma- or space-separated list in one argument
				std::string list = argv[++i];
				for (char& c : list)
					if (c == ',') c = ' ';
				std::stringstream ss(list);
				std::string v;
				while (ss >> v)
					options.vars.push_back(v);
			}
			else if (arg == "--fil-iteration") {
				options.iteration = std::stoll(argv[++i]);
			}
			else if (arg == "--fil-levels") {
				options.level_min = std::stoi(argv[++i]);
				options.level_max = std::stoi(argv[++i]);
			}
			else if (arg == "--fil-ghosts") {
				options.keep_ghosts = true;
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
				anim_step = std::stoi(argv[++i]);
			}
		}

		// Anim: one iteration per rank (the frame number is the Cactus iteration;
		// "{}" in the file and directory names is replaced by it as well)
		if (use_anim) {
			// Clamp to anim_end so trailing ranks do not address iterations past the animation range
			int anim_frame = anim_start + anim_step * world_rank;
			if (anim_frame > anim_end)
				anim_frame = anim_end;
			options.iteration = anim_frame;
			for (std::string& f : options.files)
				f = format_filename(f, anim_frame);
			for (std::string& d : options.dirs)
				d = format_filename(d, anim_frame);
			std::cout << "Reading FIL iteration: " << anim_frame << std::endl;
			world_rank = 0;
			world_size = 1;
		}

		fil::io::init_lib(options, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBFil::finish_lib()
	{
		fil::io::finish_lib();
	}

	void ConvertVDBFil::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		fil::io::get_types_and_blocks(types_and_blocks);
	}

	void ConvertVDBFil::print_types_and_blocks_local() {
		fil::io::print_types_and_blocks_local();
	}

	void ConvertVDBFil::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");
		fil::io::print_types_and_blocks(types_and_blocks);
	}

	std::string ConvertVDBFil::get_type_name(int type) {
		switch ((fil::FilParticleType)type) {
		case fil::FilParticleType::Cell: return "Cell";
		default: break;
		}
		return "Unknown";
	}

	std::string ConvertVDBFil::get_dataset_name(int blocknr) {
		return fil::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBFil::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < fil::FilParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < (int)types_and_blocks.size() / fil::FilParticleType::PTMax; bnr++) {
				if (types_and_blocks[fil::FilParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBFil::get_num_types() {
		return fil::FilParticleType::PTMax;
	}

	int ConvertVDBFil::get_num_blocks() {
		std::vector<int> types_and_blocks;
		fil::io::get_types_and_blocks(types_and_blocks);
		return (int)types_and_blocks.size() / fil::FilParticleType::PTMax;
	}
}
