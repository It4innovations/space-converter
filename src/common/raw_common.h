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

#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <algorithm>

namespace space_converter {
	namespace common {
		namespace vdb {

			/**
			 * @brief Structure to hold raw particle data for serialization and transmission
			 *
			 * RawParticles stores particle data in a compact format suitable for serialization,
			 * network transmission, and file I/O operations.
			 */
			struct RawParticles {
				/**
				 * @brief Single particle data attribute (position, velocity, etc.)
				 */
				struct ParticleData {
					std::string name;              ///< Name of the data attribute (e.g., "position", "velocity")
					int num_comp = 0;              ///< Number of components per particle (3 for position, 1 for scalar)
					std::vector<float> values;     ///< Flat array of values (size = num_particles * num_comp)
				};

				std::vector<ParticleData> data;  ///< Collection of all particle data attributes

				/**
				 * @brief Serialize particle data to a binary file
				 * @param filename Path to output file
				 */
				void serialize(const std::string& filename) const {
					std::ofstream out(filename, std::ios::binary);
					if (!out) {
						throw std::runtime_error("Failed to open file for writing");
					}

					// Write the size of the data vector
					size_t dataSize = data.size();
					out.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));

					// Write each ParticleData
					for (const auto& particle : data) {
						// Write the size of the name string
						size_t nameSize = particle.name.size();
						out.write(reinterpret_cast<const char*>(&nameSize), sizeof(nameSize));

						// Write the characters of the name string
						out.write(particle.name.data(), nameSize);

						// Write num_comp
						out.write(reinterpret_cast<const char*>(&particle.num_comp), sizeof(particle.num_comp));

						// Write the size of the values vector
						size_t valuesSize = particle.values.size();
						out.write(reinterpret_cast<const char*>(&valuesSize), sizeof(valuesSize));

						// Write the values
						out.write(reinterpret_cast<const char*>(particle.values.data()), valuesSize * sizeof(float));
					}

					out.close();
				}

				/**
				 * @brief Deserialize particle data from a binary file
				 * @param filename Path to input file
				 */
				void deserialize(const std::string& filename) {
					std::ifstream in(filename, std::ios::binary);
					if (!in) {
						throw std::runtime_error("Failed to open file for reading");
					}

					// Read the size of the data vector
					size_t dataSize;
					in.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));

					data.resize(dataSize);

					// Read each ParticleData
					for (auto& particle : data) {
						// Read the size of the name string
						size_t nameSize;
						in.read(reinterpret_cast<char*>(&nameSize), sizeof(nameSize));

						// Read the characters of the name string
						particle.name.resize(nameSize);
						in.read(&particle.name[0], nameSize);

						// Read num_comp
						in.read(reinterpret_cast<char*>(&particle.num_comp), sizeof(particle.num_comp));

						// Read the size of the values vector
						size_t valuesSize;
						in.read(reinterpret_cast<char*>(&valuesSize), sizeof(valuesSize));

						// Read the values
						particle.values.resize(valuesSize);
						in.read(reinterpret_cast<char*>(particle.values.data()), valuesSize * sizeof(float));
					}

					in.close();
				}

				/**
				 * @brief Serialize particle data to a binary buffer (for MPI/network transfer)
				 * @param bin_data Output buffer to store serialized data
				 */
				void serialize(std::vector<uint8_t>& bin_data) const {
					std::ostringstream oss(std::ios::binary);

					// Write the size of the data vector
					size_t dataSize = data.size();
					oss.write(reinterpret_cast<const char*>(&dataSize), sizeof(dataSize));

					// Write each ParticleData
					for (const auto& particle : data) {
						// Write the size of the name string
						size_t nameSize = particle.name.size();
						oss.write(reinterpret_cast<const char*>(&nameSize), sizeof(nameSize));

						// Write the characters of the name string
						oss.write(particle.name.data(), nameSize);

						// Write num_comp
						oss.write(reinterpret_cast<const char*>(&particle.num_comp), sizeof(particle.num_comp));

						// Write the size of the values vector
						size_t valuesSize = particle.values.size();
						oss.write(reinterpret_cast<const char*>(&valuesSize), sizeof(valuesSize));

						// Write the values
						oss.write(reinterpret_cast<const char*>(particle.values.data()), valuesSize * sizeof(float));
					}

					const std::string& str = oss.str();
					bin_data.assign(str.begin(), str.end());
				}
				/**
				 * @brief Deserialize particle data from a binary buffer
				 * @param bin_data Input buffer containing serialized data
				 */
				void deserialize(const std::vector<uint8_t>& bin_data) {
					std::istringstream iss(std::string(bin_data.begin(), bin_data.end()), std::ios::binary);

					// Read the size of the data vector
					size_t dataSize;
					iss.read(reinterpret_cast<char*>(&dataSize), sizeof(dataSize));

					data.resize(dataSize);

					// Read each ParticleData
					for (auto& particle : data) {
						// Read the size of the name string
						size_t nameSize;
						iss.read(reinterpret_cast<char*>(&nameSize), sizeof(nameSize));

						// Read the characters of the name string
						particle.name.resize(nameSize);
						iss.read(&particle.name[0], nameSize);

						// Read num_comp
						iss.read(reinterpret_cast<char*>(&particle.num_comp), sizeof(particle.num_comp));

						// Read the size of the values vector
						size_t valuesSize;
						iss.read(reinterpret_cast<char*>(&valuesSize), sizeof(valuesSize));

						// Read the values
						particle.values.resize(valuesSize);
						iss.read(reinterpret_cast<char*>(particle.values.data()), valuesSize * sizeof(float));
					}
				}

				/**
				 * @brief Merge another RawParticles object into this one
				 * @param other Source particle data to merge
				 *
				 * Appends values for matching particle data names, or adds new data attributes.
				 */
				void merge(const RawParticles& other) {
					for (const auto& otherParticle : other.data) {
						auto it = std::find_if(data.begin(), data.end(), [&](const ParticleData& particle) {
							return particle.name == otherParticle.name;
							});

						if (it != data.end()) {
							// If a particle with the same name exists, append values
							it->values.insert(it->values.end(), otherParticle.values.begin(), otherParticle.values.end());
						}
						else {
							data.push_back(otherParticle);
						}
					}
				}
			};

		}// vdb
	} //common
} // namespace space_converter