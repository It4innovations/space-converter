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

#include "gadget_convert_vdb.h"
#include "gadget_extract_iolib.h"

#include <iomanip>

namespace gadget {
	//std::string format_filename(const std::string& pattern, int number) {
	//	// Find the formatting marker ({}) in the pattern
	//	size_t start = pattern.find('{');
	//	size_t end = pattern.find('}', start);

	//	if (start == std::string::npos || end == std::string::npos) {
	//		std::cerr << "Pattern must contain a format specifier enclosed in curly braces" << std::endl;
	//	}

	//	char temp[1024];
	//	if (number < 999) {
	//		sprintf(temp, "%03d", number);
	//	}
	//	else {
	//		sprintf(temp, "%d", number);
	//	}

	//	// Construct the final filename by replacing the format specifier in the pattern
	//	std::string result = pattern;
	//	result.replace(start, end - start + 1, temp);

	//	return result;
	//}

	std::string format_filename(const std::string& pattern, int number) {
		// Create the formatted number as a string
		std::ostringstream formattedNumber;
		if (number < 1000) {
			formattedNumber << std::setw(3) << std::setfill('0') << number;
		}
		else {
			formattedNumber << number;
		}

		std::string result = pattern;
		std::string placeholder = "{}";
		size_t pos = result.find(placeholder);

		while (pos != std::string::npos) {
			// Replace the first occurrence of the placeholder with the formatted number
			result.replace(pos, placeholder.length(), formattedNumber.str());

			// Find the next occurrence of the placeholder
			pos = result.find(placeholder, pos + formattedNumber.str().length());
		}

		return result;
	}

	void ConvertVDBGadget::print_CPU_steps() {
		gadget::io::print_CPU_steps();
	}

	float ConvertVDBGadget::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return gadget::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBGadget::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return gadget::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBGadget::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return gadget::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBGadget::get_particle_type(uint64_t id) {
		return gadget::io::get_particle_type(id);
	}
	void ConvertVDBGadget::get_particle_position(uint64_t id, double* pos) const {
		gadget::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBGadget::get_local_num_particles() const {
		return gadget::io::get_local_num_particles();
	}
	size_t ConvertVDBGadget::get_global_num_particles() const {
		return gadget::io::get_global_num_particles();
	}

	//double ConvertVDBGadget::get_particle_radius(uint64_t id) {
	//	return gadget::io::get_particle_radius(id);
	//}

	double ConvertVDBGadget::get_particle_hsml(uint64_t id) {
		return gadget::io::get_particle_hsml(id);
	}

	double ConvertVDBGadget::get_particle_mass(uint64_t id) {
		return gadget::io::get_particle_mass(id);
	}

	double ConvertVDBGadget::get_particle_rho_internal(uint64_t id) {
		return gadget::io::get_particle_rho(id);
	}

	int ConvertVDBGadget::get_particle_rho_blocknr() {
		return gadget::io::get_particle_rho_blocknr();
	}

	void ConvertVDBGadget::init_lib(int argc, char** argv, int world_rank, int world_size) {	
		std::string param_file;
		std::string snap_file;
		int MaxMemSize = 1000;
		double BufferSize = 100;
		double PartAllocFactor = 1.2;
		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;
		int TotBHs = 1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--param-file") {
				param_file = argv[++i]; //std::stof(av[++i]);
			}
			else if (arg == "--gadget-file") {
				snap_file = argv[++i]; //std::stoi(argv[++i]);
			}
			else if (arg == "--max-mem-size") {
				MaxMemSize = std::stoi(argv[++i]);
			}
			else if (arg == "--buffer-size") {
				BufferSize = std::stof(argv[++i]);
			}
			else if (arg == "--part-alloc-factor") {
				PartAllocFactor = std::stof(argv[++i]);
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
			}
			else if (arg == "--bh-count") {
				TotBHs = std::stoi(argv[++i]);
			}
		}

		// Anim
		if (use_anim) {
			snap_file = format_filename(snap_file, anim_start + world_rank);
			std::cout << "Reading step file: " << snap_file << std::endl;
		}

		if (use_anim) {
			world_rank = 0;
			world_size = 1;
		}

		gadget::io::gadget_init_lib(world_rank, world_size);

		//read_parameter_file(param_file, NULL, NULL, NULL, 0);
		//void set_parameter(int ICFormat, int SnapFormat, int NumFilesWrittenInParallel, int MaxMemSize, double BufferSize, double PartAllocFactor)
		gadget::io::gadget_set_parameter(2, 2, world_size, MaxMemSize, BufferSize, PartAllocFactor, TotBHs);
		gadget::io::gadget_mymalloc_init();
		gadget::io::gadget_set_units();

		// Allocate memory for the char array
		std::vector<char> csnap_file(snap_file.length() + 1);

		// Copy the contents of the std::string to the char array
		std::strcpy(csnap_file.data(), snap_file.c_str());

		gadget::io::gadget_read_ic(csnap_file.data());

		print_CPU_steps();
	}

	void ConvertVDBGadget::finish_lib()
	{
		gadget::io::gadget_finish_lib();
	}

	void ConvertVDBGadget::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		gadget::io::get_types_and_blocks(types_and_blocks);
	}
	void ConvertVDBGadget::print_types_and_blocks_local() {
		gadget::io::print_types_and_blocks_local();
	}
	void ConvertVDBGadget::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		gadget::io::print_types_and_blocks(types_and_blocks);
	}

	std::string ConvertVDBGadget::get_type_name(int type) {
		return gadget::io::gadget_get_type_name(type);
	}
	std::string ConvertVDBGadget::get_dataset_name(int blocknr) {
		return gadget::io::gadget_get_dataset_name(blocknr);
	}

	std::string ConvertVDBGadget::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < 6; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / 6; bnr++) {
				if (types_and_blocks[6 * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}
	int ConvertVDBGadget::get_num_types() {
		return gadget::io::get_num_types();
	}
	int ConvertVDBGadget::get_num_blocks() {
		return gadget::io::get_num_blocks();
	}
}