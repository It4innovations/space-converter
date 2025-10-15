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

#include <stdlib.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>

#include <string>
#include <sstream>

#if defined(WITH_OPENVDB)
#	include <openvdb/openvdb.h>
#	include <openvdb/io/Stream.h>
#endif

 // stb_image
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

std::string get_prefix(const std::string& input) {
	std::stringstream ss(input);
	std::string item;
	std::vector<std::string> tokens;

	while (std::getline(ss, item, '-')) {
		tokens.push_back(item);
		if (tokens.size() == 2) break; // We only need the first two
	}

	if (tokens.size() >= 2) {
		return tokens[0] + "-" + tokens[1];
	}
	else {
		return input; // Handle unexpected input
	}
}

static void calc_histo(openvdb::FloatGrid::Ptr grid, std::string &histogram_filename)
{
	float min = 0.0f, max = 0.0f;
	grid->tree().evalMinMax(min, max);
	std::cout << "Computed value range: " << min << "," << max << '\n';

	openvdb::CoordBBox bbox;
	grid->tree().evalActiveVoxelBoundingBox(bbox);
	openvdb::Coord dim =  bbox.dim();
	std::cout << "Computed dim: " << dim.x() << "," << dim.y() << "," << dim.z() << '\n';

	auto acc = grid->getAccessor();

	// compute histogram
	//#define N (1<<10)
	constexpr int N = 1024;

	std::vector<uint64_t> counts(N);
	std::memset(counts.data(), 0, sizeof(uint64_t) * N);

#ifdef WITH_OPENMP
#pragma omp parallel
	{
		std::vector<uint64_t> local_counts(N);
		std::memset(local_counts.data(), 0, sizeof(uint64_t) * N);

#pragma omp for
		for (int z = 0; z < dim.z(); ++z) {
			for (int y = 0; y < dim.y(); ++y) {
				for (int x = 0; x < dim.x(); ++x) {
					openvdb::Coord xyz(x, y, z);

					if (!acc.isValueOn(xyz))
						continue; // Skip inactive voxels

					float value = acc.getValue(xyz);
					value -= min;
					value /= (max - min);
					value *= (N - 1);
					int index = int(value);

					if (index < 0) {
						std::cerr << "Warning: Value out of bounds: " << value << " at (" << x << ", " << y << ", " << z << ")\n";
						index = 0;
					}
					else if (index >= N) {
						std::cerr << "Warning: Value out of bounds: " << value << " at (" << x << ", " << y << ", " << z << ")\n";
						index = N - 1;
					}

					// Update the thread-local counts
					local_counts[index]++;
				}
			}
		}

		// Combine local counts into the global counts array
#pragma omp critical
		{
			for (int i = 0; i < N; ++i) {
				counts[i] += local_counts[i];
			}
		}
	}
#else
	for (int z = 0; z < dim.z(); ++z) {
		for (int y = 0; y < dim.y(); ++y) {
			for (int x = 0; x < dim.x(); ++x) {
				openvdb::Coord xyz(x, y, z);
				float value = acc.getValue(xyz);
				value -= min;
				value /= max - min;
				value *= N - 1;
				counts[int(value)]++;
			}
		}
	}
#endif  

	uint64_t max_count = 0;
	for (int i = 0; i < N; ++i) {
		max_count = std::max(max_count, counts[i]);
	}
	std::vector<float> normalizedCounts(N);
	for (int i = 0; i < N; ++i) {
		if (1) {
			normalizedCounts[i] = logf(counts[i]) / logf(max_count);
		}
		else {
			normalizedCounts[i] = counts[i] / float(max_count);
		}
	}
	std::cout << max_count << '\n';

	constexpr int W = 1024, H = 768;
	std::vector<unsigned int> img(W * H);
	for (int y = 0; y < H; ++y) {
		for (int x = 0; x < W; ++x) {
			int index = x + W * y;
			img[index] = 0xffffffff;
			float binf = x / float(W - 1);
			binf *= normalizedCounts.size() - 1;
			int bin = binf;

			int yy = H - y - 1;
			if (yy <= normalizedCounts[bin] * H) {
				img[index] = 0xff000000;
			}
		}
	}

#if 1
	stbi_write_png(histogram_filename.c_str(), W, H, 4, img.data(), W * 4);
#endif
}

int main(int argc, char** argv) {
	
	if (argc < 3) {
		printf("usage: vdb2histo file.vdb -g grid-name -o out.png\n");
		exit(0);
	}

	std::string out_file_png = "";
	std::string grid_name = "";
	std::string vdb_file = "";

	for (int i = 1; i < argc; i++) {
		const std::string arg = argv[i];
		if (arg == "-o") {
			out_file_png = argv[++i];
		}
		else if (arg == "-g") {
			grid_name = argv[++i];
		}
		else {
			vdb_file = arg;
		}
	}

#if defined(WITH_OPENVDB)
	if (out_file_png.empty()) {
		std::cerr << "Error: Output filename not specified. Use -o option.\n";
		return 1;
	}

	if (grid_name.empty()) {
		std::cerr << "Error: Grid name not specified. Use -g option.\n";
		return 1;
	}

	if (vdb_file.empty()) {
		std::cerr << "Error: No input file specified.\n";
		return 1;
	}

	openvdb::initialize();  // Initialize OpenVDB if available.

	// Process each input file
	try {
		// Open the VDB file
		openvdb::io::File input_file(vdb_file);
		input_file.open();
			
		// Get all grids from the file
		openvdb::GridPtrVecPtr file_grids = input_file.getGrids();
		for (auto& base_grid : *file_grids) {
			if (base_grid->getName() == grid_name) {
				openvdb::FloatGrid::Ptr grid = openvdb::gridPtrCast<openvdb::FloatGrid>(base_grid);
				calc_histo(grid, out_file_png);

				break;
			}
		}
			
		input_file.close();
		std::cout << "Processed file: " << vdb_file << " (" << file_grids->size() << " grids)\n";
	}
	catch (const std::exception& e) {
		std::cerr << "Error processing file " << vdb_file << ": " << e.what() << std::endl;
	}

#else
	std::cerr << "Error: OpenVDB support is not enabled.\n";
	return 1;
#endif

	return 0;
}