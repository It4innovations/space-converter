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

#include "haccbin_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>

#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#include <omp.h>

#include "convert_common.h"

#define RETURN_NORM_EMPTY return 0; //return std::numeric_limits<float>::quiet_NaN();//0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)common::calculate_dmagnitude3(v[0], v[1], v[2]);
#define RETURN_NORM_DVECTORN(v,n) return (float)common::calculate_dmagnituden(v,n);
#define RETURN_NORM_FVECTORN(v,n) return (float)common::calculate_fmagnituden(v,n);

#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;
#define RETURN_COMP_DVECTORN(v,n) return n;
#define RETURN_COMP_FVECTORN(v,n) return n;

#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)
#define RETURN_ORIG_DVECTORN(v,n) for(int iv=0;iv<n;iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_DVECTORN(v,n)
#define RETURN_ORIG_FVECTORN(v,n) for(int iv=0;iv<n;iv++) out_value[iv] = (float)v[iv]; RETURN_COMP_FVECTORN(v,n)

#ifdef _WIN32
#	define STDMAX max
#	define STDMIN min
#else
#	define STDMAX std::max
#	define STDMIN std::min
#endif

namespace haccbin {
	namespace io {

		// AUX
		//struct HACCBinAuxValue {
		//	/// The array of gas particles
		//	std::vector<float> gas;
		//	/// The array of dark matter particles
		//	std::vector<float> darks;
		//	/// The array of star particles
		//	std::vector<float> stars;

		//	int num_components = 0;
		//	std::string name;
		//};

		//struct HACCBinHeader {
		//	/// The time of the output
		//	double time;
		//	/// The number of particles of all types in this file
		//	unsigned int nbodies;
		//	/// The number of dimensions, must be equal to MAXDIM
		//	int ndim;
		//	/// The number of SPH (gas) particles in this file
		//	unsigned int nsph;
		//	/// The number of dark matter particles in this file
		//	unsigned int ndark;
		//	/// The number of star particles in this file
		//	unsigned int nstar;
		//	//int pad; //unused on x86
		//};

		//struct HACCBinReader {
		//	std::string basefile;
		//	//HACCBinHeader header;
		//	size_t line_count;
		//	std::vector<std::string> headers;


		//	HACCBinReader(std::string basefile_) :
		//		basefile(basefile_), line_count(0)
		//	{
		//	}

		//	void read_header() {
		//		// Open the HACCBin file
		//		std::ifstream file(basefile);

		//		if (!file.is_open()) {
		//			std::cerr << "Error: Could not open the file " << basefile << std::endl;
		//			return;
		//		}

		//		line_count = 0;
		//		headers.clear();

		//		std::string line;
		//		// Read the first line (header)
		//		if (std::getline(file, line)) {
		//			std::cout << "Header row: " << line << std::endl;

		//			// Parse the header row

		//			std::stringstream ss(line);
		//			std::string column;

		//			while (std::getline(ss, column, ',')) {
		//				headers.push_back(column);
		//			}

		//			// Print out the header columns
		//			//std::cout << "Header columns:" << std::endl;
		//			//for (const auto& header : headers) {
		//			//	std::cout << header << std::endl;
		//			//}
		//		}
		//		else {
		//			std::cerr << "Error: The file is empty or could not read the header row." << std::endl;
		//		}

		//		// Count the remaining lines
		//		while (std::getline(file, line)) {
		//			line_count++;
		//		}

		//		file.close();
		//	}
		//};

#pragma pack(push, 1)  // Set structure packing to 1-byte alignment
		struct HACCBinParticle {
			// position and velocity vectors
			float x; //4
			float vx;//8
			float y; //12
			float vy;
			float z;
			float vz;//24
			// other fields
			float mass;//28
			float uu;//32
			float hh;
			float mu;
			float rho;
			float phi;//48
			// particle ID and mask
			int64_t id;//56
			uint16_t mask;//58
		};
#pragma pack(pop)  // Restore default packing

		struct HACCBinFile {
			std::string basefile;
			size_t start_particles;
			size_t end_particles;

			/// The array of gas particles
			std::vector<HACCBinParticle> particles;
			///// The array of dark matter particles
			//std::vector<HACCBinParticle> darks;
			///// The array of star particles
			//std::vector<HACCBinParticle> stars;

			//std::vector<HACCBinAuxValue> aux_values;

			/// The header for the full file
			//HACCBinHeader fullHeader;
			size_t nbodies;
			/// The header for the part of the file we hold
			//HACCBinHeader h;

			HACCBinFile(std::string basefile_, size_t start_particles_, size_t end_particles_, size_t nbodies_) :
				basefile(basefile_), start_particles(start_particles_), end_particles(end_particles_), nbodies(nbodies_)
			{
			}

			//int get_dataset(const std::string& name) {
			//	// Check predefined names
			//	if (name == "Pos") return HACCBinBlockType::Pos;
			//	if (name == "Mass") return HACCBinBlockType::Mass;
			//	if (name == "Vel") return HACCBinBlockType::Vel;
			//	if (name == "Soft") return HACCBinBlockType::Soft;
			//	if (name == "Phi") return HACCBinBlockType::Phi;
			//	if (name == "HSmooth") return HACCBinBlockType::HSmooth;
			//	if (name == "Rho") return HACCBinBlockType::Rho;
			//	if (name == "Temp") return HACCBinBlockType::Temp;
			//	if (name == "Metals") return HACCBinBlockType::Metals;
			//	if (name == "TForm") return HACCBinBlockType::TForm;

			//	//// Check auxiliary data names
			//	//if (tf_data) {
			//	//	for (size_t i = 0; i < tf_data->aux_values.size(); ++i) {
			//	//		if (tf_data->aux_values[i].name == name) {
			//	//			return static_cast<HACCBinBlockType>(static_cast<int>(HACCBinBlockType::BTMax) + i);
			//	//		}
			//	//	}
			//	//}

			//	// If no match is found, print an error
			//	std::cerr << "Error: Unknown dataset name: " << name << std::endl;

			//	return -1;
			//}

			static size_t get_count(std::string filename) {
				std::ifstream file(filename, std::ios::binary | std::ios::ate); // Open in binary mode and seek to end

				if (!file) {
					std::cerr << "Error: Unable to open file " << filename << std::endl;
					return 0;
				}

				std::streamsize fileSize = file.tellg(); // Get total file size
				file.close();

				return fileSize / sizeof(HACCBinParticle); // Compute number of particles
			}

			void read_data() {
				// Open the HACCBin file
				std::ifstream file(basefile, std::ios::binary);

				if (!file) {
					std::cerr << "Error: Unable to open file " << basefile << std::endl;
					return;
				}

				//// Compute the total number of particles in file
				//file.seekg(0, std::ios::end);
				//size_t totalParticles = file.tellg() / sizeof(HACCBinParticle);

				//// Validate the requested range
				//if (start_particles >= totalParticles || end_particles > totalParticles || start_particles > end_particles) {
				//	std::cerr << "Error: Invalid particle range (" << start_particles << " to " << end_particles << ")" << std::endl;
				//	file.close();
				//	return;
				//}

				// Seek to the starting particle position
				file.seekg(start_particles * sizeof(HACCBinParticle), std::ios::beg);

				size_t count = end_particles - start_particles;
				particles.resize(count); // Preallocate memory for efficiency

				// Read only the required particles
				file.read((char*)particles.data(), count * sizeof(HACCBinParticle));

				//std::vector<char> temp(count * sizeof(HACCBinParticle));
				//file.read(reinterpret_cast<char*>(temp.data()), count * sizeof(HACCBinParticle));
				
				//for (int i = 0; i < temp.size() - 58; i+=58) {
				//	float *t = (float*)&temp[i];
				//	//if (t[0] > -100000.0f && t[0] < 100000.0f && fabs(t[0]) > 0.00001f &&
				//	//	t[2] > -100000.0f && t[2] < 100000.0f && fabs(t[2]) > 0.00001f &&
				//	//	t[4] > -100000.0f && t[4] < 100000.0f && fabs(t[4]) > 0.00001f
				//	//	) {
				//	//	printf("id: %d: %f, %f, %f: check: %f\n", i, t[0], t[2], t[4], (double)count * sizeof(HACCBinParticle) / (double)i );
				//	//}
				//	printf("id: %d: %f, %f, %f: check: %f\n", i, t[0], t[2], t[4], (double)count * sizeof(HACCBinParticle) / (double)i);
				//	printf("id: %d: %f, %f, %f\n", i, t[1], t[3], t[5]);
				//}

				//file.read(reinterpret_cast<char*>(&particle), sizeof(HACCBinParticle));

				//for (int i = 0; i < 10; i++) {
				//	HACCBinParticle particle;
				//	file.read(reinterpret_cast<char*>(&particle), sizeof(HACCBinParticle));

				//	particles.push_back(particle);
				//}

				//if (file.eof()) {
				//	std::cout << "End of file reached successfully.\n";
				//}
				//else if (file.fail()) {
				//	std::cerr << "Error: Read operation failed (possible format issue).\n";
				//}
				//else if (file.bad()) {
				//	std::cerr << "Error: Serious I/O error occurred.\n";
				//}

				file.close();
			}
		};

		//static double fh_time; // gross, but quick way to get time
		double steps_time[2];

		//Tipsy::TipsyFile* tf_data = nullptr;
		HACCBinFile* tf_data = nullptr;
		size_t g_total_particles = 0;
		size_t g_start_particles = 0;
		size_t g_num_particles = 0;

		int g_smoothlength_blocknr = -1;

		void init_lib(std::string basefile, int world_rank, int world_size) {
			steps_time[0] = omp_get_wtime();

			//tf_data = new Tipsy::TipsyFile(basefile);
			//HACCBinReader treader(basefile);
			//treader.read_header();
			//Tipsy::header th = treader.getHeader();

			g_total_particles = HACCBinFile::get_count(basefile);// .line_count;//STDMAX(th.nsph, STDMAX(th.ndark, th.nstar));

			// Each process calculates its range of particles to read
			int particles_per_process = g_total_particles / world_size;
			g_start_particles = world_rank * particles_per_process;
			g_num_particles = (world_rank == world_size - 1) ? g_total_particles - g_start_particles : particles_per_process;

			tf_data = new HACCBinFile(basefile, g_start_particles, g_start_particles + g_num_particles, g_total_particles);
			tf_data->read_data();

			//check smoothlength
			//for (int i = 0; i < tf_data->aux_values.size(); i++) {
			//	if (tf_data->aux_values[i].name == "smoothlength") {
			//		g_smoothlength_blocknr = i + HACCBinBlockType::BTMax;
			//		break;
			//	}
			//}

			steps_time[1] = omp_get_wtime();
		}

		void finish_lib() {
			if (tf_data != nullptr)
				delete tf_data;

			tf_data = nullptr;
		}

		//uint64_t get_particle_type_offset(HACCBinParticleType pt) {
		//	switch (pt) {
		//	case HACCBinParticleType::Gas:
		//		return 0;
		//	case HACCBinParticleType::Dark:
		//		return tf_data->gas.size();
		//	case HACCBinParticleType::Stars:
		//		return tf_data->gas.size() + tf_data->darks.size();
		//	}

		//	return 0;
		//}

		void print_CPU_steps() {
			printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		float get_particle_norm_value(int blocknr, uint64_t id) {
			HACCBinBlockType bt = (HACCBinBlockType)blocknr;

			switch (bt) {
			case HACCBinBlockType::Pos:
				float pos[3];
				pos[0] = tf_data->particles[id].x;
				pos[1] = tf_data->particles[id].y;
				pos[2] = tf_data->particles[id].z;
				RETURN_NORM_VECTOR3(pos);
			case HACCBinBlockType::Vel:
				float vel[3];
				vel[0] = tf_data->particles[id].vx;
				vel[1] = tf_data->particles[id].vy;
				vel[2] = tf_data->particles[id].vz;
				RETURN_NORM_VECTOR3(vel);
			case HACCBinBlockType::Mass:
				RETURN_NORM_VALUE(tf_data->particles[id].mass);
			case HACCBinBlockType::UU:
				RETURN_NORM_VALUE(tf_data->particles[id].uu);
			case HACCBinBlockType::HH:
				RETURN_NORM_VALUE(tf_data->particles[id].hh);
			case HACCBinBlockType::MU:
				RETURN_NORM_VALUE(tf_data->particles[id].mu);
			case HACCBinBlockType::Rho:
				RETURN_NORM_VALUE(tf_data->particles[id].rho);
			case HACCBinBlockType::Phi:
				RETURN_NORM_VALUE(tf_data->particles[id].phi);
			case HACCBinBlockType::Id:
				RETURN_NORM_VALUE(tf_data->particles[id].id);
			case HACCBinBlockType::Mask:
				RETURN_NORM_VALUE(tf_data->particles[id].mask);
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value) {

			HACCBinBlockType bt = (HACCBinBlockType)blocknr;

			switch (bt) {
			case HACCBinBlockType::Pos:
				float pos[3];
				pos[0] = tf_data->particles[id].x;
				pos[1] = tf_data->particles[id].y;
				pos[2] = tf_data->particles[id].z;
				RETURN_ORIG_VECTOR3(pos);
			case HACCBinBlockType::Vel:
				float vel[3];
				vel[0] = tf_data->particles[id].vx;
				vel[1] = tf_data->particles[id].vy;
				vel[2] = tf_data->particles[id].vz;
				RETURN_ORIG_VECTOR3(vel);
			case HACCBinBlockType::Mass:
				RETURN_ORIG_VALUE(tf_data->particles[id].mass);
			case HACCBinBlockType::UU:
				RETURN_ORIG_VALUE(tf_data->particles[id].uu);
			case HACCBinBlockType::HH:
				RETURN_ORIG_VALUE(tf_data->particles[id].hh);
			case HACCBinBlockType::MU:
				RETURN_ORIG_VALUE(tf_data->particles[id].mu);
			case HACCBinBlockType::Rho:
				RETURN_ORIG_VALUE(tf_data->particles[id].rho);
			case HACCBinBlockType::Phi:
				RETURN_ORIG_VALUE(tf_data->particles[id].phi);
			case HACCBinBlockType::Id:
				RETURN_ORIG_VALUE(tf_data->particles[id].id);
			case HACCBinBlockType::Mask:
				RETURN_ORIG_VALUE(tf_data->particles[id].mask);
			}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id) {
			HACCBinBlockType bt = (HACCBinBlockType)blocknr;

			switch (bt) {
			case HACCBinBlockType::Pos:
				float pos[3];
				RETURN_COMP_VECTOR3(pos);
			case HACCBinBlockType::Vel:
				float vel[3];
				RETURN_COMP_VECTOR3(vel);
			case HACCBinBlockType::Mass:
				RETURN_COMP_VALUE(tf_data->particles[id].mass);
			case HACCBinBlockType::UU:
				RETURN_COMP_VALUE(tf_data->particles[id].uu);
			case HACCBinBlockType::HH:
				RETURN_COMP_VALUE(tf_data->particles[id].hh);
			case HACCBinBlockType::MU:
				RETURN_COMP_VALUE(tf_data->particles[id].mu);
			case HACCBinBlockType::Rho:
				RETURN_COMP_VALUE(tf_data->particles[id].rho);
			case HACCBinBlockType::Phi:
				RETURN_COMP_VALUE(tf_data->particles[id].phi);
			case HACCBinBlockType::Id:
				RETURN_COMP_VALUE(tf_data->particles[id].id);
			case HACCBinBlockType::Mask:
				RETURN_COMP_VALUE(tf_data->particles[id].mask);
			}

			RETURN_COMP_EMPTY;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			// 3 * 10
			types_and_blocks.resize(HACCBinParticleType::PTMax * (HACCBinBlockType::BTMax));
			memset(types_and_blocks.data(), 0, types_and_blocks.size() * sizeof(int));
			//memcpy(types_and_blocks.data(), all_particles_blocks, sizeof(all_particles_blocks));

			int* tab = types_and_blocks.data();

			for (int ipt = 0; ipt < HACCBinParticleType::PTMax; ipt++) {
				// tipsy default
				for (int ibt = 0; ibt < HACCBinBlockType::BTMax; ibt++) {
					tab[ibt * HACCBinParticleType::PTMax + ipt]++;
				}
			}
		}

		std::string get_dataset_name(int blocknr) {
			HACCBinBlockType bt = (HACCBinBlockType)blocknr;

			switch (bt) {
			case HACCBinBlockType::Pos: return "Positions";
			case HACCBinBlockType::Vel: return "Velocities";
			case HACCBinBlockType::Mass: return "Mass";
			case HACCBinBlockType::UU: return "Internal energy";
			case HACCBinBlockType::HH: return "SPH smoothing length";
			case HACCBinBlockType::MU: return "Molecular weight";
			case HACCBinBlockType::Rho: return "Density";
			case HACCBinBlockType::Phi: return "Gravitational potential";
			case HACCBinBlockType::Id: return "Id";
			case HACCBinBlockType::Mask: return "Mask";
			}

			return "Unknown";
		}

		int get_particle_type(uint64_t id) {

			static constexpr uint16_t BARYON_BIT = 1u << 1;  // 2nd bit
			static constexpr uint16_t STAR_BIT = 1u << 5;  // 6th bit
			static constexpr uint16_t WIND_BIT = 1u << 6;  // 7th bit
			static constexpr uint16_t STAR_FORM_BIT = 1u << 7;  // 8th bit
			static constexpr uint16_t AGN_FLAG_BIT = 1u << 8;  // 9th bit

			uint16_t mask = tf_data->particles[id].mask;

			//6th bit : Denotes if a baryon particle is also a star particle
			if ((mask & STAR_BIT) != 0)
				return HACCBinParticleType::BaryonStart;

			//7th bit : Denotes if a baryon particle is also a wind particle
			if ((mask & WIND_BIT) != 0)
				return HACCBinParticleType::BaryonWind;

			//8th bit : Denotes if a baryon particle is also a star - forming gas particle
			if ((mask & STAR_FORM_BIT) != 0)
				return HACCBinParticleType::BaryonGas;

			//9th bit : Denotes if a dark matter particle has been flagged as an active galactic nuclei(AGN)
			if ((mask & AGN_FLAG_BIT) != 0)
				return HACCBinParticleType::DarkMatterAGN;

			//2nd bit : Denotes whether a particle is dark matter(0) or a baryon(1)
			if ((mask & BARYON_BIT) != 0)
				return HACCBinParticleType::Baryon;
			else
				return HACCBinParticleType::DarkMatter;
		}

		void get_particle_position(uint64_t id, double* pos) {
			pos[0] = tf_data->particles[id].x;
			pos[1] = tf_data->particles[id].y;
			pos[2] = tf_data->particles[id].z;
		}

		size_t get_local_num_particles() {
			return tf_data->particles.size();
		}
		size_t get_global_num_particles() {
			return tf_data->nbodies;
		}

		//double get_particle_radius(uint64_t id) {
		//	return get_particle_hsml(id) * 2.0;
		//}

		double get_particle_hsml(uint64_t id) {
			return tf_data->particles[id].hh;
		}

		double get_particle_mass(uint64_t id) {
			return get_particle_norm_value(HACCBinBlockType::Mass, id);
		}

		int get_particle_rho_blocknr() {
			return HACCBinBlockType::Rho;
		}

		double get_particle_rho(uint64_t id) {
			return get_particle_norm_value(HACCBinBlockType::Rho, id);
		}

		std::string get_particle_unit(int blocknr) {
			return "";
		}
	}
}