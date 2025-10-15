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

#include <float.h>

#if defined(WITH_OPENVDB)
#	include <openvdb/openvdb.h>
#	include <openvdb/io/Stream.h>
#	include <openvdb/tools/Interpolation.h>
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

void parse_comma_separated_ints(const std::string& input, std::vector<int>& result) {
	result.clear();
	std::stringstream ss(input);
	std::string token;

	while (std::getline(ss, token, ',')) {
		// Trim whitespace if needed
		token.erase(0, token.find_first_not_of(" \t"));
		token.erase(token.find_last_not_of(" \t") + 1);

		if (!token.empty()) {
			result.push_back(std::stoi(token));
		}
	}
}

static void calc_slice(openvdb::FloatGrid::Ptr grid, std::string& slice_filename, std::string slice_axis, int slice_nr, std::string ramp_values, bool use_log)
{
	std::vector<int> ramp_rgb_colors;
	parse_comma_separated_ints(ramp_values, ramp_rgb_colors);

	//float min = 0.0f, max = 0.0f;
	//grid->tree().evalMinMax(min, max);
	//std::cout << "Computed value range: " << min << "," << max << '\n';

	openvdb::CoordBBox bbox;
	grid->tree().evalActiveVoxelBoundingBox(bbox);

	int width, height;
	if (slice_axis == "x") {
		width = bbox.dim().y();
		height = bbox.dim().z();
	}
	else if (slice_axis == "y") {
		width = bbox.dim().x();
		height = bbox.dim().z();
	}
	else if (slice_axis == "z") {
		width = bbox.dim().x();
		height = bbox.dim().y();
	}
	else {
		std::cerr << "Invalid slice axis." << std::endl;
		return;
	}

	// RGB image (3 channels)
	std::vector<float> image_float(width * height, 0.0);
	std::vector<unsigned char> image_rgb(width * height * 3, 0);

	openvdb::tools::GridSampler<openvdb::FloatTree, openvdb::tools::BoxSampler> sampler(grid->constTree(), grid->transform());

	// Calculate the number of color entries (RGB triplets)
	int ramp_size = ramp_rgb_colors.size() / 3;

#ifdef WITH_OPENMP
#	pragma omp parallel for
#endif
	for (int j = 0; j < height; ++j) {
		for (int i = 0; i < width; ++i) {
			openvdb::Vec3d coord;
			if (slice_axis == "x") coord = grid->indexToWorld(openvdb::Coord(slice_nr, bbox.min().y() + i, bbox.min().z() + j));
			else if (slice_axis == "y") coord = grid->indexToWorld(openvdb::Coord(bbox.min().x() + i, slice_nr, bbox.min().z() + j));
			else coord = grid->indexToWorld(openvdb::Coord(bbox.min().x() + i, bbox.min().y() + j, slice_nr));

			float value = sampler.wsSample(coord);  // Fixed normalization
			int pixel_index = j * width + i;

			if (use_log) {
				if (value < 1.0f)
					value = 0.0f;
				else
					value = log10(value);
			}

			image_float[pixel_index] = value;     // Red
		}
	}

	float min = FLT_MAX;
	float max = -FLT_MAX;

#ifdef WITH_OPENMP
#	pragma omp parallel for reduction(min:min) reduction(max:max)
#endif
	for (int j = 0; j < height; ++j) {
		for (int i = 0; i < width; ++i) {
			int pixel_index = j * width + i;

			min = std::min(min, image_float[pixel_index]);
			max = std::max(max, image_float[pixel_index]);
		}
	}

#ifdef WITH_OPENMP
#	pragma omp parallel for
#endif
	for (int j = 0; j < height; ++j) {
		for (int i = 0; i < width; ++i) {
			// Map value to color ramp
			int pixel_index = j * width + i;

			float value = (image_float[pixel_index] - min) / (max - min);  // Fixed normalization
			value = std::clamp(value, 0.0f, 1.0f); // clamp between 0-1

			if (fabs(value) < 0.000001f) // skip, black color
				continue;

			if (ramp_size > 0) {
				// Map normalized value (0-1) to ramp index
				float ramp_pos = value * (ramp_size - 1);
				int ramp_index = static_cast<int>(ramp_pos);
				float frac = ramp_pos - ramp_index;

				// Clamp ramp_index to valid range
				ramp_index = std::clamp(ramp_index, 0, ramp_size - 1);
				int next_index = std::min(ramp_index + 1, ramp_size - 1);

				// Get colors from ramp (RGB triplets)
				int r1 = ramp_rgb_colors[ramp_index * 3 + 0];
				int g1 = ramp_rgb_colors[ramp_index * 3 + 1];
				int b1 = ramp_rgb_colors[ramp_index * 3 + 2];

				int r2 = ramp_rgb_colors[next_index * 3 + 0];
				int g2 = ramp_rgb_colors[next_index * 3 + 1];
				int b2 = ramp_rgb_colors[next_index * 3 + 2];

				// Linear interpolation between colors
				unsigned char r = static_cast<unsigned char>(r1 + frac * (r2 - r1));
				unsigned char g = static_cast<unsigned char>(g1 + frac * (g2 - g1));
				unsigned char b = static_cast<unsigned char>(b1 + frac * (b2 - b1));

				image_rgb[pixel_index * 3 + 0] = r;     // Red
				image_rgb[pixel_index * 3 + 1] = g;     // Green
				image_rgb[pixel_index * 3 + 2] = b;     // Blue
			}
			else {
				// Fallback: grayscale if no ramp provided
				unsigned char gray = static_cast<unsigned char>(value * 255);
				image_rgb[pixel_index * 3 + 0] = gray;
				image_rgb[pixel_index * 3 + 1] = gray;
				image_rgb[pixel_index * 3 + 2] = gray;
			}
		}
	}

	// Save as RGB PNG (3 channels)
	stbi_write_png(slice_filename.c_str(), width, height, 3, image_rgb.data(), width * 3);

	std::cout << "Saved slice to " << slice_filename << std::endl;
}

int main(int argc, char** argv) {
	
	if (argc < 3) {
		printf("usage: vdb2histo file.vdb -g grid-name -o out.png\n");
		exit(0);
	}

	std::string out_file_png = "";
	std::string grid_name = "";
	std::string vdb_file = "";
	std::string slice_axis = "x";
	std::string ramp_values = "";
	int slice_nr = 0;
	bool use_log = false;

	for (int i = 1; i < argc; i++) {
		const std::string arg = argv[i];
		if (arg == "-o") {
			out_file_png = argv[++i];
		}
		else if (arg == "-g") {
			grid_name = argv[++i];
		}
		else if (arg == "--axis") {
			slice_axis = argv[++i];
		}
		else if (arg == "--slice") {
			slice_nr = std::atoi(argv[++i]);
		}
		else if (arg == "--ramp") {
			ramp_values = argv[++i];
		}
		else if (arg == "--use_log") {
			use_log = true;
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
				calc_slice(grid, out_file_png, slice_axis, slice_nr, ramp_values, use_log);

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