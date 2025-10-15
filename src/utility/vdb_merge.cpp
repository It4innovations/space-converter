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

int main(int argc, char** argv) {
	
	if (argc < 3) {
		printf("usage: vdb_merge file1.vdb file2.vdb ... -o file_merged.vdb\n");
		exit(0);
	}

	std::string out_file = "";
	std::vector<std::string> files;

	for (int i = 1; i < argc; i++) {
		const std::string arg = argv[i];
		if (arg == "-o") {
			out_file = argv[++i];
		}
		else {
			files.push_back(arg);
		}
	}

#if defined(WITH_OPENVDB)
	if (out_file.empty()) {
		std::cerr << "Error: Output file not specified. Use -o option.\n";
		return 1;
	}

	if (files.empty()) {
		std::cerr << "Error: No input files specified.\n";
		return 1;
	}

	openvdb::initialize();  // Initialize OpenVDB if available.

	// Create a new grid container for merged grids
	openvdb::GridPtrVec all_grids;

	// Process each input file
	for (const auto& file : files) {
		try {
			// Open the VDB file
			openvdb::io::File input_file(file);
			input_file.open();
			
			// Get all grids from the file
			openvdb::GridPtrVecPtr file_grids = input_file.getGrids();
			for (auto& grid : *file_grids) {
				if (grid) {
					// Use std::filesystem to extract the filename without extension
					std::filesystem::path filePath(input_file.filename());
					std::string attribute_name = get_prefix(filePath.stem().string());

					grid->setName(attribute_name);
					all_grids.push_back(grid);
				}
			}
			
			input_file.close();
			std::cout << "Processed file: " << file << " (" << file_grids->size() << " grids)\n";
		}
		catch (const std::exception& e) {
			std::cerr << "Error processing file " << file << ": " << e.what() << std::endl;
		}
	}

	// Write all grids to the output file
	try {
		openvdb::io::File output_file(out_file);
		output_file.write(all_grids);
		output_file.close();
		std::cout << "Successfully merged " << all_grids.size() << " grids into " << out_file << std::endl;
	}
	catch (const std::exception& e) {
		std::cerr << "Error writing to output file " << out_file << ": " << e.what() << std::endl;
		return 1;
	}
#else
	std::cerr << "Error: OpenVDB support is not enabled.\n";
	return 1;
#endif

	return 0;
}