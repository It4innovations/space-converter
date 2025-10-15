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

#include "hdf5_convert_vdb.h"
#include "hdf5_extract_iolib.h"

namespace hdf5 {
	void ConvertVDBHDF5::print_CPU_steps() {
		hdf5::io::print_CPU_steps();
	}

	float ConvertVDBHDF5::get_particle_norm_value(int blocknr, uint64_t id) {
		return hdf5::io::get_particle_norm_value(blocknr, id);
	}
	int ConvertVDBHDF5::get_particle_value(int blocknr, uint64_t id, float* value) {
		return hdf5::io::get_particle_value(blocknr, id, value);
	}
	int ConvertVDBHDF5::get_particle_value_comp(int blocknr, uint64_t id) {
		return hdf5::io::get_particle_value_comp(blocknr, id);
	}
	int ConvertVDBHDF5::get_particle_type(uint64_t id) {
		return hdf5::io::get_particle_type(id);
	}
	void ConvertVDBHDF5::get_particle_position(uint64_t id, double* pos) const {
		hdf5::io::get_particle_position(id, pos);
	}
	size_t ConvertVDBHDF5::get_local_num_particles() const {
		return hdf5::io::get_local_num_particles();
	}
	size_t ConvertVDBHDF5::get_global_num_particles() const {
		return hdf5::io::get_global_num_particles();
	}

	double ConvertVDBHDF5::get_particle_radius(uint64_t id) {
		return hdf5::io::get_particle_radius(id);
	}

	double ConvertVDBHDF5::get_particle_hsml(uint64_t id) {
		return hdf5::io::get_particle_hsml(id);
	}

	double ConvertVDBHDF5::get_particle_mass(uint64_t id) {
		return hdf5::io::get_particle_mass(id);
	}

	double ConvertVDBHDF5::get_particle_rho(uint64_t id) {
		return hdf5::io::get_particle_rho(id);
	}

	void ConvertVDBHDF5::init_lib(int argc, char** argv, int world_rank, int world_size) {	
		std::string hdf5_file;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];
			if (arg == "--hdf5-file") {
				hdf5_file = argv[++i]; //std::stoi(argv[++i]);
			}
		}

		hdf5::io::init_lib(hdf5_file, world_rank, world_size);

		print_CPU_steps();
	}

	void ConvertVDBHDF5::finish_lib()
	{
		hdf5::io::finish_lib();
	}

	void ConvertVDBHDF5::get_types_and_blocks(std::vector<int>& types_and_blocks) {
		hdf5::io::get_types_and_blocks(types_and_blocks);
	}
	void ConvertVDBHDF5::print_types_and_blocks_local() {
		std::vector<int> types_and_blocks;
		hdf5::io::get_types_and_blocks(types_and_blocks);

		printf("\n");

		//char buf[500];
		for (int type = 0; type < hdf5::HDF5ParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / hdf5::HDF5ParticleType::PTMax; blocknr++) {

				if (types_and_blocks[hdf5::HDF5ParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	void ConvertVDBHDF5::print_types_and_blocks(std::vector<int>& types_and_blocks) {
		printf("\nAll snapshots contain:\n");

		//char buf[500];
		for (int type = 0; type < hdf5::HDF5ParticleType::PTMax; type++) {
			printf("Type: %s (%d)\n", get_type_name(type).c_str(), type);
			for (int blocknr = 0; blocknr < types_and_blocks.size() / hdf5::HDF5ParticleType::PTMax; blocknr++) {

				if (types_and_blocks[hdf5::HDF5ParticleType::PTMax * blocknr + type] > 0) {
					std::string buf = get_dataset_name(blocknr);
					printf("\t%s (%d)\n", buf.c_str(), blocknr);
				}
			}
		}
	}

	std::string ConvertVDBHDF5::get_type_name(int type) {		

		hdf5::HDF5ParticleType pt = (hdf5::HDF5ParticleType)type;
		switch (pt) {
		case hdf5::HDF5ParticleType::Gas:
			return "Gas";
		case hdf5::HDF5ParticleType::Dark:
			return "Dark";
		case hdf5::HDF5ParticleType::Stars:
			return "Stars";
		}
		return "Unknown";
	}
	std::string ConvertVDBHDF5::get_dataset_name(int blocknr) {
		return hdf5::io::get_dataset_name(blocknr);
	}

	std::string ConvertVDBHDF5::get_particle_data_types(std::vector<int>& types_and_blocks) {
		std::string particle_data_types = "";

		for (int t = 0; t < hdf5::HDF5ParticleType::PTMax; t++) {
			for (int bnr = 0; bnr < types_and_blocks.size() / hdf5::HDF5ParticleType::PTMax; bnr++) {
				if (types_and_blocks[hdf5::HDF5ParticleType::PTMax * bnr + t] == 0)
					continue;

				particle_data_types = particle_data_types + get_type_name(t) + ";" + std::to_string(t) + ";" + get_dataset_name(bnr) + ";" + std::to_string(bnr) + "\n";
			}
		}

		return particle_data_types;
	}
}