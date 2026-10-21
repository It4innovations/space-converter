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

#include "pluto_vtk_convert_vdb.h"
#include "pluto_vtk_extract_iolib.h"

namespace plutovtk {
	void ConvertVDBPlutoVTK::print_CPU_steps() {
		plutovtk::io::print_CPU_steps();
	}

	float ConvertVDBPlutoVTK::get_particle_norm_value_internal(int blocknr, uint64_t id) {
		return plutovtk::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBPlutoVTK::get_particle_value_internal(int blocknr, uint64_t id, float* value) {
		return plutovtk::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBPlutoVTK::get_particle_value_comp_internal(int blocknr, uint64_t id) {
		return plutovtk::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBPlutoVTK::get_particle_type(uint64_t id) {
		return plutovtk::io::get_particle_type(id);
	}
	void ConvertVDBPlutoVTK::get_particle_position(uint64_t id, double* pos) const {
		plutovtk::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBPlutoVTK::get_local_num_particles() const {
		return plutovtk::io::get_local_num_particles();
	}
	size_t ConvertVDBPlutoVTK::get_global_num_particles() const {
		return plutovtk::io::get_global_num_particles();
	}

	double ConvertVDBPlutoVTK::get_particle_hsml(uint64_t id) {
		return plutovtk::io::get_particle_hsml(id);
	}

	double ConvertVDBPlutoVTK::get_particle_mass(uint64_t id) {
		return plutovtk::io::get_particle_mass(id);
	}

	double ConvertVDBPlutoVTK::get_particle_rho_internal(uint64_t id) {
		return plutovtk::io::get_particle_rho(id);
	}

	int ConvertVDBPlutoVTK::get_particle_rho_blocknr() {
		return plutovtk::io::get_particle_rho_blocknr();
	}

	void ConvertVDBPlutoVTK::init_lib(int argc, char** argv, int world_rank, int world_size) {	
		std::string vtk_file;
		std::vector<std::string> scalar_names_vec;

		bool use_anim = false;
		int anim_start = -1;
		int anim_end = -1;
		int anim_step = -1;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--vtk-file") {
				vtk_file = argv[++i];
			}
			else if (arg == "--scalar-names") {
				// Read scalar field names (can be multiple)
				while (i + 1 < argc && argv[i + 1][0] != '-') {
					scalar_names_vec.push_back(argv[++i]);
				}
			}
			else if (arg == "--anim") {
				use_anim = true;
				anim_start = std::stoi(argv[++i]);
				anim_end = std::stoi(argv[++i]);
				anim_step = std::stoi(argv[++i]);
			}
		}

		// Anim
		if (use_anim) {
			vtk_file = format_filename(vtk_file, anim_start + anim_step * world_rank);
			std::cout << "Reading timestep file: " << vtk_file << std::endl;
		}

		if (use_anim) {
			world_rank = 0;
			world_size = 1;
		}

		plutovtk::io::init_lib(vtk_file, world_rank, world_size, scalar_names_vec);

		print_CPU_steps();
	}

	void ConvertVDBPlutoVTK::finish_lib()
	{
		plutovtk::io::finish_lib();
	}

	void ConvertVDBPlutoVTK::get_types_and_blocks_internal(std::vector<int>& types_and_blocks) {
		plutovtk::io::get_types_and_blocks(types_and_blocks);
	}
	void ConvertVDBPlutoVTK::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		plutovtk::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		for (int type = 0; type < plutovtk::PlutoVTKParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / plutovtk::PlutoVTKParticleType::PTMax; blocknr++) {

				if (types_and_blocks[plutovtk::PlutoVTKParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBPlutoVTK::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		for (int type = 0; type < plutovtk::PlutoVTKParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / plutovtk::PlutoVTKParticleType::PTMax; blocknr++) {

				if (types_and_blocks[plutovtk::PlutoVTKParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBPlutoVTK::get_type_name(int type) {		
		plutovtk::PlutoVTKParticleType pt = (plutovtk::PlutoVTKParticleType)type;
		switch (pt) {
		case plutovtk::PlutoVTKParticleType::PLUTO:
			return "PLUTO";
		}
		return "Unknown";
	}
	
	std::string ConvertVDBPlutoVTK::get_dataset_name(int blocknr) {
		return plutovtk::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBPlutoVTK::get_particle_data_type_names(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < plutovtk::PlutoVTKParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / plutovtk::PlutoVTKParticleType::PTMax; bnr++) {
				if (types_and_blocks[plutovtk::PlutoVTKParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}

	int ConvertVDBPlutoVTK::get_num_types() {
		return plutovtk::PlutoVTKParticleType::PTMax;
	}
	
	int ConvertVDBPlutoVTK::get_num_blocks() {
		std::vector<int> types_and_blocks;
		plutovtk::io::get_types_and_blocks(types_and_blocks);
		return types_and_blocks.size() / plutovtk::PlutoVTKParticleType::PTMax;
	}
}
