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

#include "hdf5_extract_iolib.h"

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

#include <H5Cpp.h>

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

namespace hdf5 {
	namespace io {

		// AUX
		struct HDF5AuxValue {
			/// The array of gas particles
			std::vector<float> gas;
			/// The array of dark matter particles
			std::vector<float> darks;
			/// The array of star particles
			std::vector<float> stars;

			int num_components = 0;
			std::string name;
		};

		struct HDF5Header {
			/// The time of the output
			double time;
			/// The number of particles of all types in this file
			unsigned int nbodies;
			/// The number of dimensions, must be equal to MAXDIM
			int ndim;
			/// The number of SPH (gas) particles in this file
			unsigned int nsph;
			/// The number of dark matter particles in this file
			unsigned int ndark;
			/// The number of star particles in this file
			unsigned int nstar;
			//int pad; //unused on x86
		};

		struct HDF5Reader {
			std::string basefile;
			//HDF5Header header;
			size_t line_count;
			std::vector<std::string> headers;


			HDF5Reader(std::string basefile_) :
				basefile(basefile_), line_count(0)
			{
			}

			void read_header() {
				// Open the HDF5 file
				std::ifstream file(basefile);

				if (!file.is_open()) {
					std::cerr << "Error: Could not open the file " << basefile << std::endl;
					return;
				}

				line_count = 0;
				headers.clear();

				std::string line;
				// Read the first line (header)
				if (std::getline(file, line)) {
					std::cout << "Header row: " << line << std::endl;

					// Parse the header row

					std::stringstream ss(line);
					std::string column;

					while (std::getline(ss, column, ',')) {
						headers.push_back(column);
					}

					// Print out the header columns
					//std::cout << "Header columns:" << std::endl;
					//for (const auto& header : headers) {
					//	std::cout << header << std::endl;
					//}
				}
				else {
					std::cerr << "Error: The file is empty or could not read the header row." << std::endl;
				}

				// Count the remaining lines
				while (std::getline(file, line)) {
					line_count++;
				}

				file.close();
			}
		};

		struct HDF5Particle {
			/** The position of this particle. */
			float pos[3];

			/** The mass of the particle. */
			float mass;
			/** The velocity vector of this particle. */
			float vel[3];

			// gass
			/** The local density of gas at this particle's location. */
			float rho;
			/** The temperature of the gas in this particle. */
			float temp;
			/** The gravitational softening length of this particle. */
			float hsmooth;
			/** The metal content of this gas particle. */
			float metals;
			/** The gravitational potential at this particle. */
			float phi;

			// dark
			/** The gravitational softening length of this particle. */
			float eps;

			// stars
			/** The time of formation of this star particle. */
			float tform;
		};

		struct HDF5File {
			std::string basefile;
			size_t start_particles;
			size_t end_particles;

			/// The array of gas particles
			std::vector<HDF5Particle> gas;
			/// The array of dark matter particles
			std::vector<HDF5Particle> darks;
			/// The array of star particles
			std::vector<HDF5Particle> stars;

			std::vector<HDF5AuxValue> aux_values;

			/// The header for the full file
			//HDF5Header fullHeader;
			size_t nbodies;
			/// The header for the part of the file we hold
			HDF5Header h;

			HDF5File(std::string basefile_, size_t start_particles_, size_t end_particles_, size_t nbodies_) :
				basefile(basefile_), start_particles(start_particles_), end_particles(end_particles_), nbodies(nbodies_)
			{
			}

			//int get_dataset(const std::string& name) {
			//	// Check predefined names
			//	if (name == "Pos") return HDF5BlockType::Pos;
			//	if (name == "Mass") return HDF5BlockType::Mass;
			//	if (name == "Vel") return HDF5BlockType::Vel;
			//	if (name == "Soft") return HDF5BlockType::Soft;
			//	if (name == "Phi") return HDF5BlockType::Phi;
			//	if (name == "HSmooth") return HDF5BlockType::HSmooth;
			//	if (name == "Rho") return HDF5BlockType::Rho;
			//	if (name == "Temp") return HDF5BlockType::Temp;
			//	if (name == "Metals") return HDF5BlockType::Metals;
			//	if (name == "TForm") return HDF5BlockType::TForm;

			//	//// Check auxiliary data names
			//	//if (tf_data) {
			//	//	for (size_t i = 0; i < tf_data->aux_values.size(); ++i) {
			//	//		if (tf_data->aux_values[i].name == name) {
			//	//			return static_cast<HDF5BlockType>(static_cast<int>(HDF5BlockType::BTMax) + i);
			//	//		}
			//	//	}
			//	//}

			//	// If no match is found, print an error
			//	std::cerr << "Error: Unknown dataset name: " << name << std::endl;

			//	return -1;
			//}

			void read_data() {
				const char line_del = ',';
				const char value_del = ';';

				// Open the HDF5 file
				std::ifstream file(basefile);

				if (!file.is_open()) {
					std::cerr << "Error: Could not open the file " << basefile << std::endl;
					return;
				}

				// read headers
				std::vector<std::string> headers;
				std::string line;
				// Read the first line (header)
				if (std::getline(file, line)) {
					//std::cout << "Header row: " << line << std::endl;

					// Parse the header row
					std::stringstream ss(line);
					std::string column;

					while (std::getline(ss, column, line_del)) {
						headers.push_back(column);
					}

					// Print out the header columns
					//std::cout << "Header columns:" << std::endl;
					//for (const auto& header : headers) {
					//	std::cout << header << std::endl;
					//}
				}
				else {
					std::cerr << "Error: The file is empty or could not read the header row." << std::endl;
				}

				for (size_t l = 0; l < end_particles; l++) {
					if (std::getline(file, line)) {
						if (l < start_particles) //skip lines
							break;

						// Parse the row
						std::stringstream ss(line);
						std::string column;

						std::vector<std::string> line_data;
						while (std::getline(ss, column, line_del)) {
							line_data.push_back(column);
						}

						if (line_data.size() != headers.size()) {
							std::cerr << "Error: The file has wrong number of columuns for the row: " << l << std::endl;
							continue;
						}

						HDF5Particle hdf5_particle;
						memset(&hdf5_particle, 0, sizeof(hdf5_particle));

						HDF5ParticleType ptype;
						std::vector<std::string> aux_names_tmp;
						std::vector<float> aux_values_tmp;
						std::vector<int> aux_comp_tmp;

						for (int h = 0; h < headers.size(); h++) {
							if (headers[h] == "Type") {
								ptype = (HDF5ParticleType)atoi(line_data[h].c_str());
							}
							else if (headers[h] == "Pos") {
								std::stringstream ss_values(line_data[h]);
								std::string value;
								std::vector<std::string> values;
								while (std::getline(ss_values, value, value_del)) {
									values.push_back(value);
								}

								if (values.size() != 3) {
									std::cerr << "Error: The file has wrong number of values for the row (expected 3): " << l << std::endl;
									continue;
								}

								hdf5_particle.pos[0] = atof(values[0].c_str());
								hdf5_particle.pos[1] = atof(values[1].c_str());
								hdf5_particle.pos[2] = atof(values[2].c_str());
							}
							else if (headers[h] == "Mass") {
								hdf5_particle.mass = atof(line_data[h].c_str());
							}
							else if (headers[h] == "Vel") {
								std::stringstream ss_values(line_data[h]);
								std::string value;
								std::vector<std::string> values;
								while (std::getline(ss_values, value, value_del)) {
									values.push_back(value);
								}

								if (values.size() != 3) {
									std::cerr << "Error: The file has wrong number of values for the row (expected 3): " << l << std::endl;
									continue;
								}

								hdf5_particle.vel[0] = atof(values[0].c_str());
								hdf5_particle.vel[1] = atof(values[1].c_str());
								hdf5_particle.vel[2] = atof(values[2].c_str());
							}
							else if (headers[h] == "Soft") {
								hdf5_particle.eps = atof(line_data[h].c_str());
							}
							else if (headers[h] == "Phi") {
								hdf5_particle.phi = atof(line_data[h].c_str());
							}
							else if (headers[h] == "HSmooth") {
								hdf5_particle.hsmooth = atof(line_data[h].c_str());
							}
							else if (headers[h] == "Rho") {
								hdf5_particle.rho = atof(line_data[h].c_str());
							}
							else if (headers[h] == "Temp") {
								hdf5_particle.temp = atof(line_data[h].c_str());
							}
							else if (headers[h] == "Metals") {
								hdf5_particle.metals = atof(line_data[h].c_str());
							}
							else if (headers[h] == "TForm") {
								hdf5_particle.tform = atof(line_data[h].c_str());
							}
							else {
								//std::cerr << "Error: The file has unknown column: " << headers[h] << std::endl;
								aux_names_tmp.push_back(headers[h]);
								//aux_values_tmp.push_back(atof(line_data[h].c_str()));

								std::stringstream ss_values(line_data[h]);
								std::string value;
								std::vector<std::string> values;
								while (std::getline(ss_values, value, value_del)) {
									values.push_back(value);
								}

								aux_comp_tmp.push_back(values.size());

								for (int c = 0; c < values.size(); c++) {
									aux_values_tmp.push_back(atof(values[c].c_str()));
								}
							}
						}

						if (aux_values.size() == 0 && aux_names_tmp.size() != 0) {
							aux_values.resize(aux_names_tmp.size());
						}
						else if (aux_values.size() != aux_names_tmp.size() && aux_names_tmp.size() > 0) {
							std::cerr << "Error: The file has wrong number of aux_values!" << std::endl;
						}

						for (int a = 0; a < aux_values.size(); a++) {
							if (aux_values[a].num_components == 0) {
								aux_values[a].name = aux_names_tmp[a];
								aux_values[a].num_components = aux_comp_tmp[a];

								if (aux_values[a].num_components == 0) {
									std::cerr << "Error: The file has wrong data of aux_values (num_components == 0)!" << std::endl;
								}
							}
							else if (aux_values[a].num_components != aux_comp_tmp[a] ||
								aux_values[a].name != aux_names_tmp[a]) {
								std::cerr << "Error: The file has wrong data of aux_values!" << std::endl;
							}						
						}

						switch (ptype) {
						case HDF5ParticleType::Gas: {
							for (int a = 0; a < aux_values.size(); a++) {
								aux_values[a].gas.push_back(aux_values_tmp[a]);
							}

							gas.push_back(hdf5_particle);
							break;
						}
						case HDF5ParticleType::Dark: {
							for (int a = 0; a < aux_values.size(); a++) {
								aux_values[a].darks.push_back(aux_values_tmp[a]);
							}

							darks.push_back(hdf5_particle);
							break;
						}
						case HDF5ParticleType::Stars: {
							for (int a = 0; a < aux_values.size(); a++) {
								aux_values[a].stars.push_back(aux_values_tmp[a]);
							}

							stars.push_back(hdf5_particle);
							break;
						}
						}
					}
					else {
						std::cerr << "Error: The file is empty or could not read the row: " << l << std::endl;
					}
				}

				file.close();

				//set header
				h.time = 0;
				h.nbodies = nbodies;
				h.ndim = 0;
				h.nsph = gas.size();
				h.ndark = darks.size();
				h.nstar = stars.size();
			}
		};

		//static double fh_time; // gross, but quick way to get time
		double steps_time[2];

		//Tipsy::TipsyFile* tf_data = nullptr;
		HDF5File* tf_data = nullptr;
		size_t g_total_particles = 0;
		size_t g_start_particles = 0;
		size_t g_num_particles = 0;

		int g_smoothlength_blocknr = -1;

		void init_lib(std::string basefile, int world_rank, int world_size) {
			steps_time[0] = omp_get_wtime();

			H5::H5File file(basefile, H5F_ACC_RDONLY);

			//tf_data = new Tipsy::TipsyFile(basefile);
			HDF5Reader treader(basefile);
			treader.read_header();
			//Tipsy::header th = treader.getHeader();

			g_total_particles = treader.line_count;//STDMAX(th.nsph, STDMAX(th.ndark, th.nstar));

			// Each process calculates its range of particles to read
			int particles_per_process = g_total_particles / world_size;
			g_start_particles = world_rank * particles_per_process;
			g_num_particles = (world_rank == world_size - 1) ? g_total_particles - g_start_particles : particles_per_process;

			tf_data = new HDF5File(basefile, g_start_particles, g_start_particles + g_num_particles, treader.line_count);
			tf_data->read_data();

			//check smoothlength
			for (int i = 0; i < tf_data->aux_values.size(); i++) {
				if (tf_data->aux_values[i].name == "smoothlength") {
					g_smoothlength_blocknr = i + HDF5BlockType::BTMax;
					break;
				}
			}

			steps_time[1] = omp_get_wtime();
		}

		void finish_lib() {
			if (tf_data != nullptr)
				delete tf_data;

			tf_data = nullptr;
		}

		uint64_t get_particle_type_offset(HDF5ParticleType pt) {
			switch (pt) {
			case HDF5ParticleType::Gas:
				return 0;
			case HDF5ParticleType::Dark:
				return tf_data->gas.size();
			case HDF5ParticleType::Stars:
				return tf_data->gas.size() + tf_data->darks.size();
			}

			return 0;
		}

		void print_CPU_steps() {
			printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		float get_particle_norm_value(int blocknr, uint64_t id) {
			HDF5ParticleType pt = (HDF5ParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			HDF5BlockType bt = (HDF5BlockType)blocknr;

			switch (pt) {
			case HDF5ParticleType::Gas: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_NORM_VECTOR3(tf_data->gas[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_NORM_VALUE(tf_data->gas[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_NORM_VECTOR3(tf_data->gas[id - offset].vel);
				case HDF5BlockType::Phi:
					RETURN_NORM_VALUE(tf_data->gas[id - offset].phi);
				case HDF5BlockType::HSmooth:
					RETURN_NORM_VALUE(tf_data->gas[id - offset].hsmooth);
				case HDF5BlockType::Rho:
					RETURN_NORM_VALUE(tf_data->gas[id - offset].rho);
				case HDF5BlockType::Temp:
					RETURN_NORM_VALUE(tf_data->gas[id - offset].temp);
				case HDF5BlockType::Metals:
					RETURN_NORM_VALUE(tf_data->gas[id - offset].metals);
				}
				break;
			}
			case HDF5ParticleType::Dark: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_NORM_VECTOR3(tf_data->darks[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_NORM_VALUE(tf_data->darks[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_NORM_VECTOR3(tf_data->darks[id - offset].vel);
				case HDF5BlockType::Soft:
					RETURN_NORM_VALUE(tf_data->darks[id - offset].eps);
				case HDF5BlockType::Phi:
					RETURN_NORM_VALUE(tf_data->darks[id - offset].phi);
				}
				break;
			}
			case HDF5ParticleType::Stars: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_NORM_VECTOR3(tf_data->stars[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_NORM_VALUE(tf_data->stars[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_NORM_VECTOR3(tf_data->stars[id - offset].vel);
				case HDF5BlockType::Soft:
					RETURN_NORM_VALUE(tf_data->stars[id - offset].eps);
				case HDF5BlockType::Phi:
					RETURN_NORM_VALUE(tf_data->stars[id - offset].phi);
				case HDF5BlockType::Metals:
					RETURN_NORM_VALUE(tf_data->stars[id - offset].metals);
				case HDF5BlockType::TForm:
					RETURN_NORM_VALUE(tf_data->stars[id - offset].tform);
				}
				break;
			}
			}

			if (blocknr >= HDF5BlockType::BTMax && tf_data->aux_values.size() > blocknr - HDF5BlockType::BTMax) {

				HDF5AuxValue& av = tf_data->aux_values[blocknr - HDF5BlockType::BTMax];

				if (pt == HDF5ParticleType::Gas) {
					if (av.num_components == 1) {
						RETURN_NORM_VALUE(av.gas[id - offset]);
					}

					if (av.num_components == 3) {
						float v[3];
						v[0] = av.gas[(id - offset) * 3 + 0];
						v[1] = av.gas[(id - offset) * 3 + 1];
						v[2] = av.gas[(id - offset) * 3 + 2];
						RETURN_NORM_VECTOR3(v);
					}
				}

				if (pt == HDF5ParticleType::Dark) {
					if (av.num_components == 1) {
						RETURN_NORM_VALUE(av.darks[id - offset]);
					}

					if (av.num_components == 3) {
						float v[3];
						v[0] = av.darks[(id - offset) * 3 + 0];
						v[1] = av.darks[(id - offset) * 3 + 1];
						v[2] = av.darks[(id - offset) * 3 + 2];
						RETURN_NORM_VECTOR3(v);
					}
				}

				if (pt == HDF5ParticleType::Stars) {
					if (av.num_components == 1) {
						RETURN_NORM_VALUE(av.stars[id - offset]);
					}

					if (av.num_components == 3) {
						float v[3];
						v[0] = av.stars[(id - offset) * 3 + 0];
						v[1] = av.stars[(id - offset) * 3 + 1];
						v[2] = av.stars[(id - offset) * 3 + 2];
						RETURN_NORM_VECTOR3(v);
					}
				}
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value) {
			HDF5ParticleType pt = (HDF5ParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			HDF5BlockType bt = (HDF5BlockType)blocknr;

			switch (pt) {
			case HDF5ParticleType::Gas: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_ORIG_VECTOR3(tf_data->gas[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_ORIG_VALUE(tf_data->gas[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_ORIG_VECTOR3(tf_data->gas[id - offset].vel);
				case HDF5BlockType::Phi:
					RETURN_ORIG_VALUE(tf_data->gas[id - offset].phi);
				case HDF5BlockType::HSmooth:
					RETURN_ORIG_VALUE(tf_data->gas[id - offset].hsmooth);
				case HDF5BlockType::Rho:
					RETURN_ORIG_VALUE(tf_data->gas[id - offset].rho);
				case HDF5BlockType::Temp:
					RETURN_ORIG_VALUE(tf_data->gas[id - offset].temp);
				case HDF5BlockType::Metals:
					RETURN_ORIG_VALUE(tf_data->gas[id - offset].metals);
				}
				break;
			}
			case HDF5ParticleType::Dark: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_ORIG_VECTOR3(tf_data->darks[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_ORIG_VALUE(tf_data->darks[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_ORIG_VECTOR3(tf_data->darks[id - offset].vel);
				case HDF5BlockType::Soft:
					RETURN_ORIG_VALUE(tf_data->darks[id - offset].eps);
				case HDF5BlockType::Phi:
					RETURN_ORIG_VALUE(tf_data->darks[id - offset].phi);
				}
				break;
			}
			case HDF5ParticleType::Stars: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_ORIG_VECTOR3(tf_data->stars[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_ORIG_VALUE(tf_data->stars[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_ORIG_VECTOR3(tf_data->stars[id - offset].vel);
				case HDF5BlockType::Soft:
					RETURN_ORIG_VALUE(tf_data->stars[id - offset].eps);
				case HDF5BlockType::Phi:
					RETURN_ORIG_VALUE(tf_data->stars[id - offset].phi);
				case HDF5BlockType::Metals:
					RETURN_ORIG_VALUE(tf_data->stars[id - offset].metals);
				case HDF5BlockType::TForm:
					RETURN_ORIG_VALUE(tf_data->stars[id - offset].tform);
				}
				break;
			}
			}

			if (blocknr >= HDF5BlockType::BTMax && tf_data->aux_values.size() > blocknr - HDF5BlockType::BTMax) {

				HDF5AuxValue& av = tf_data->aux_values[blocknr - HDF5BlockType::BTMax];

				if (pt == HDF5ParticleType::Gas) {
					if (av.num_components == 1) {
						RETURN_ORIG_VALUE(av.gas[id - offset]);
					}

					if (av.num_components == 3) {
						float v[3];
						v[0] = av.gas[(id - offset) * 3 + 0];
						v[1] = av.gas[(id - offset) * 3 + 1];
						v[2] = av.gas[(id - offset) * 3 + 2];
						RETURN_ORIG_VECTOR3(v);
					}
				}

				if (pt == HDF5ParticleType::Dark) {
					if (av.num_components == 1) {
						RETURN_ORIG_VALUE(av.darks[id - offset]);
					}

					if (av.num_components == 3) {
						float v[3];
						v[0] = av.darks[(id - offset) * 3 + 0];
						v[1] = av.darks[(id - offset) * 3 + 1];
						v[2] = av.darks[(id - offset) * 3 + 2];
						RETURN_ORIG_VECTOR3(v);
					}
				}

				if (pt == HDF5ParticleType::Stars) {
					if (av.num_components == 1) {
						RETURN_ORIG_VALUE(av.stars[id - offset]);
					}

					if (av.num_components == 3) {
						float v[3];
						v[0] = av.stars[(id - offset) * 3 + 0];
						v[1] = av.stars[(id - offset) * 3 + 1];
						v[2] = av.stars[(id - offset) * 3 + 2];
						RETURN_ORIG_VECTOR3(v);
					}
				}
			}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id) {
			HDF5ParticleType pt = (HDF5ParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			HDF5BlockType bt = (HDF5BlockType)blocknr;

			switch (pt) {
			case HDF5ParticleType::Gas: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_COMP_VECTOR3(tf_data->gas[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_COMP_VALUE(tf_data->gas[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_COMP_VECTOR3(tf_data->gas[id - offset].vel);
				case HDF5BlockType::Phi:
					RETURN_COMP_VALUE(tf_data->gas[id - offset].phi);
				case HDF5BlockType::HSmooth:
					RETURN_COMP_VALUE(tf_data->gas[id - offset].hsmooth);
				case HDF5BlockType::Rho:
					RETURN_COMP_VALUE(tf_data->gas[id - offset].rho);
				case HDF5BlockType::Temp:
					RETURN_COMP_VALUE(tf_data->gas[id - offset].temp);
				case HDF5BlockType::Metals:
					RETURN_COMP_VALUE(tf_data->gas[id - offset].metals);
				}
				break;
			}
			case HDF5ParticleType::Dark: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_COMP_VECTOR3(tf_data->darks[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_COMP_VALUE(tf_data->darks[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_COMP_VECTOR3(tf_data->darks[id - offset].vel);
				case HDF5BlockType::Soft:
					RETURN_COMP_VALUE(tf_data->darks[id - offset].eps);
				case HDF5BlockType::Phi:
					RETURN_COMP_VALUE(tf_data->darks[id - offset].phi);
				}
				break;
			}
			case HDF5ParticleType::Stars: {

				switch (bt) {
				case HDF5BlockType::Pos:
					RETURN_COMP_VECTOR3(tf_data->stars[id - offset].pos);
				case HDF5BlockType::Mass:
					RETURN_COMP_VALUE(tf_data->stars[id - offset].mass);
				case HDF5BlockType::Vel:
					RETURN_COMP_VECTOR3(tf_data->stars[id - offset].vel);
				case HDF5BlockType::Soft:
					RETURN_COMP_VALUE(tf_data->stars[id - offset].eps);
				case HDF5BlockType::Phi:
					RETURN_COMP_VALUE(tf_data->stars[id - offset].phi);
				case HDF5BlockType::Metals:
					RETURN_COMP_VALUE(tf_data->stars[id - offset].metals);
				case HDF5BlockType::TForm:
					RETURN_COMP_VALUE(tf_data->stars[id - offset].tform);
				}
				break;
			}
			}

			if (blocknr >= HDF5BlockType::BTMax && tf_data->aux_values.size() > blocknr - HDF5BlockType::BTMax) {

				HDF5AuxValue& av = tf_data->aux_values[blocknr - HDF5BlockType::BTMax];

				if (pt == HDF5ParticleType::Gas) {
					if (av.num_components == 1) {
						RETURN_COMP_VALUE(av.gas[id - offset]);
					}

					if (av.num_components == 3) {
						RETURN_COMP_VECTOR3(v);
					}
				}

				if (pt == HDF5ParticleType::Dark) {
					if (av.num_components == 1) {
						RETURN_COMP_VALUE(av.darks[id - offset]);
					}

					if (av.num_components == 3) {
						RETURN_COMP_VECTOR3(v);
					}
				}

				if (pt == HDF5ParticleType::Stars) {
					if (av.num_components == 1) {
						RETURN_COMP_VALUE(av.stars[id - offset]);
					}

					if (av.num_components == 3) {
						RETURN_COMP_VECTOR3(v);
					}
				}
			}

			RETURN_COMP_EMPTY;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			// 3 * 10
			types_and_blocks.resize(HDF5ParticleType::PTMax * (HDF5BlockType::BTMax + tf_data->aux_values.size()));
			memset(types_and_blocks.data(), 0, types_and_blocks.size() * sizeof(int));
			//memcpy(types_and_blocks.data(), all_particles_blocks, sizeof(all_particles_blocks));

			int* tab = types_and_blocks.data();

			for (int ipt = 0; ipt < HDF5ParticleType::PTMax; ipt++) {
				// tipsy default
				for (int ibt = 0; ibt < HDF5BlockType::BTMax; ibt++) {

					HDF5ParticleType pt = (HDF5ParticleType)ipt;
					HDF5BlockType bt = (HDF5BlockType)ibt;

					switch (pt) {
					case HDF5ParticleType::Gas: {

						switch (bt) {
						case HDF5BlockType::Pos:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Mass:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Vel:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Phi:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::HSmooth:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Rho:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Temp:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Metals:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						}
						break;
					}
					case HDF5ParticleType::Dark: {

						switch (bt) {
						case HDF5BlockType::Pos:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Mass:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Vel:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Soft:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Phi:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						}
						break;
					}
					case HDF5ParticleType::Stars: {

						switch (bt) {
						case HDF5BlockType::Pos:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Mass:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Vel:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Soft:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Phi:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::Metals:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						case HDF5BlockType::TForm:
							tab[ibt * HDF5ParticleType::PTMax + ipt]++;
							break;
						}
						break;
					}
					}
				}

				// tipsy aux
				for (int ibt = HDF5BlockType::BTMax; ibt < HDF5BlockType::BTMax + tf_data->aux_values.size(); ibt++) {
					HDF5ParticleType pt = (HDF5ParticleType)ipt;

					switch (pt) {
					case HDF5ParticleType::Gas: {
						tab[ibt * HDF5ParticleType::PTMax + ipt]++;
						break;
					}
					case HDF5ParticleType::Dark: {
						tab[ibt * HDF5ParticleType::PTMax + ipt]++;
						break;
					}
					case HDF5ParticleType::Stars: {
						tab[ibt * HDF5ParticleType::PTMax + ipt]++;
						break;
					}
					}
				}
			}
		}

		std::string get_dataset_name(int blocknr) {
			HDF5BlockType bt = (HDF5BlockType)blocknr;
			switch (bt) {
			case HDF5BlockType::Pos:
				return "Pos";
			case HDF5BlockType::Mass:
				return "Mass";
			case HDF5BlockType::Vel:
				return "Vel";
			case HDF5BlockType::Soft:
				return "Soft";
			case HDF5BlockType::Phi:
				return "Phi";
			case HDF5BlockType::HSmooth:
				return "HSmooth";
			case HDF5BlockType::Rho:
				return "Rho";
			case HDF5BlockType::Temp:
				return "Temp";
			case HDF5BlockType::Metals:
				return "Metals";
			case HDF5BlockType::TForm:
				return "TForm";
			}

			if (blocknr >= HDF5BlockType::BTMax && tf_data->aux_values.size() > blocknr - HDF5BlockType::BTMax) {
				return tf_data->aux_values[blocknr - HDF5BlockType::BTMax].name;
			}

			return "Unknown";
		}

		int get_particle_type(uint64_t id) {
			if (id < tf_data->gas.size()) {
				return HDF5ParticleType::Gas;
			}
			else if (id < tf_data->gas.size() + tf_data->darks.size()) {
				return HDF5ParticleType::Dark;
			}
			else {
				return HDF5ParticleType::Stars;
			}
		}

		void get_particle_position(uint64_t id, double* pos) {
			HDF5ParticleType pt = (HDF5ParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			switch (pt) {
			case HDF5ParticleType::Gas: {
				pos[0] = tf_data->gas[id - offset].pos[0];
				pos[1] = tf_data->gas[id - offset].pos[1];
				pos[2] = tf_data->gas[id - offset].pos[2];
				break;
			}
			case HDF5ParticleType::Dark: {
				pos[0] = tf_data->darks[id - offset].pos[0];
				pos[1] = tf_data->darks[id - offset].pos[1];
				pos[2] = tf_data->darks[id - offset].pos[2];
				break;
			}
			case HDF5ParticleType::Stars: {
				pos[0] = tf_data->stars[id - offset].pos[0];
				pos[1] = tf_data->stars[id - offset].pos[1];
				pos[2] = tf_data->stars[id - offset].pos[2];
				break;
			}
			}
		}

		size_t get_local_num_particles() {
			return tf_data->h.nbodies; //tf_data->gas.size() + tf_data->darks.size() + tf_data->stars.size();
		}
		size_t get_global_num_particles() {
			return tf_data->nbodies; //tf_data->gas.size() + tf_data->darks.size() + tf_data->stars.size();
		}

		double get_particle_radius(uint64_t id) {
			return get_particle_hsml(id) * 2.0;
		}

		double get_particle_hsml(uint64_t id) {
			if (g_smoothlength_blocknr != -1)
				return get_particle_norm_value(g_smoothlength_blocknr, id);

			return 0.0;
		}

		double get_particle_mass(uint64_t id) {
			return get_particle_norm_value(HDF5BlockType::Mass, id);
		}

		double get_particle_rho(uint64_t id) {
			return get_particle_norm_value(HDF5BlockType::Rho, id);
		}

		std::string get_particle_unit(int blocknr) {
			return "";
		}
	}
}