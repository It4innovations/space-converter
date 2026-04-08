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

#include <cstdint>
#include <cstring>
#include <array>
#include <iostream>
#include <fstream>
#include <string>

#ifdef WITH_NANOVDB

#ifdef WITH_TBB
#define NANOVDB_USE_TBB
#define NANOVDB_USE_INTRINSICS
#endif

#if OPENVDB_VERSION == 11
#	include <nanovdb/util/GridBuilder.h>
#	include <nanovdb/util/IO.h>
#else
#	include <nanovdb/tools/GridBuilder.h>
#	include <nanovdb/io/IO.h>
#endif

#	include <nanovdb/NanoVDB.h>
#endif

namespace {

enum class CombineMode {
	Mesh,
	Zip,
};

struct Args {
	std::string input_nvdb;
	std::array<std::string, 3> channel_inputs;
	std::string output_bin;
	bool use_channel_files = false;
	CombineMode mode = CombineMode::Mesh;
	uint32_t grid_index = 0;
};

void print_usage() {
	printf("usage:\n");
	printf("  nano2float3 file.nvdb file.bin [--mode mesh|zip] [--grid-index N]\n");
	printf("  nano2float3 --channel-files x.nvdb y.nvdb z.nvdb file.bin [--mode mesh|zip] [--grid-index N]\n");
	printf("  mesh: dense export over indexBBox (includes inactive voxels/background values)\n");
	printf("  zip : export only active voxels (single input) or union of active voxels (channel-files)\n");
}

bool parse_u32(const std::string& s, uint32_t& out) {
	char* end_ptr = nullptr;
	const unsigned long value = std::strtoul(s.c_str(), &end_ptr, 10);
	if (end_ptr == nullptr || *end_ptr != '\0') {
		return false;
	}
	if (value > static_cast<unsigned long>(UINT32_MAX)) {
		return false;
	}
	out = static_cast<uint32_t>(value);
	return true;
}

bool parse_args(int argc, char** argv, Args& out) {
	if (argc < 2) {
		print_usage();
		return false;
	}

	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--mode") {
			if (i + 1 >= argc) {
				fprintf(stderr, "Error: --mode expects mesh or zip\n");
				return false;
			}
			const std::string mode = argv[++i];
			if (mode == "mesh") {
				out.mode = CombineMode::Mesh;
			} else if (mode == "zip") {
				out.mode = CombineMode::Zip;
			} else {
				fprintf(stderr, "Error: unknown mode '%s'\n", mode.c_str());
				return false;
			}
		} else if (arg == "--grid-index") {
			if (i + 1 >= argc) {
				fprintf(stderr, "Error: --grid-index expects an unsigned integer\n");
				return false;
			}
			if (!parse_u32(argv[++i], out.grid_index)) {
				fprintf(stderr, "Error: invalid --grid-index value\n");
				return false;
			}
		} else if (arg == "--channel-files") {
			if (!out.input_nvdb.empty()) {
				fprintf(stderr, "Error: --channel-files cannot be used together with a positional input file\n");
				return false;
			}
			if (out.use_channel_files) {
				fprintf(stderr, "Error: --channel-files can be specified only once\n");
				return false;
			}
			if (i + 3 >= argc) {
				fprintf(stderr, "Error: --channel-files expects exactly 3 input files\n");
				return false;
			}
			out.use_channel_files = true;
			out.channel_inputs[0] = argv[++i];
			out.channel_inputs[1] = argv[++i];
			out.channel_inputs[2] = argv[++i];
		} else if (arg == "-h" || arg == "--help") {
			print_usage();
			exit(0);
		} else {
			if (!arg.empty() && arg[0] == '-') {
				fprintf(stderr, "Error: unknown argument '%s'\n", arg.c_str());
				print_usage();
				return false;
			}

			if (out.use_channel_files) {
				if (out.output_bin.empty()) {
					out.output_bin = arg;
				} else {
					fprintf(stderr, "Error: unexpected extra positional argument '%s'\n", arg.c_str());
					print_usage();
					return false;
				}
			} else {
				if (out.input_nvdb.empty()) {
					out.input_nvdb = arg;
				} else if (out.output_bin.empty()) {
					out.output_bin = arg;
				} else {
					fprintf(stderr, "Error: unexpected extra positional argument '%s'\n", arg.c_str());
					print_usage();
					return false;
				}
			}
		}
	}

	if (out.use_channel_files) {
		if (out.output_bin.empty()) {
			fprintf(stderr, "Error: missing output .bin file\n");
			print_usage();
			return false;
		}
	} else {
		if (out.input_nvdb.empty() || out.output_bin.empty()) {
			fprintf(stderr, "Error: missing input .nvdb and/or output .bin file\n");
			print_usage();
			return false;
		}
	}

	return true;
}

#ifdef WITH_NANOVDB
bool bbox_equal(const nanovdb::CoordBBox& a, const nanovdb::CoordBBox& b) {
	const nanovdb::Coord amin = a.min();
	const nanovdb::Coord amax = a.max();
	const nanovdb::Coord bmin = b.min();
	const nanovdb::Coord bmax = b.max();
	return amin[0] == bmin[0] && amin[1] == bmin[1] && amin[2] == bmin[2] &&
		amax[0] == bmax[0] && amax[1] == bmax[1] && amax[2] == bmax[2];
}
#endif

uint32_t bswap32(uint32_t x) {
	return ((x & 0x000000FFu) << 24u) |
			   ((x & 0x0000FF00u) << 8u) |
			   ((x & 0x00FF0000u) >> 8u) |
			   ((x & 0xFF000000u) >> 24u);
}

bool host_is_little_endian() {
	const uint32_t value = 1u;
	return *reinterpret_cast<const uint8_t*>(&value) == 1u;
}

void write_float3_le(std::ofstream& out, float x, float y, float z) {
	float packed[3] = { x, y, z };
	if (!host_is_little_endian()) {
		for (int i = 0; i < 3; ++i) {
			uint32_t raw = 0u;
			std::memcpy(&raw, &packed[i], sizeof(uint32_t));
			raw = bswap32(raw);
			std::memcpy(&packed[i], &raw, sizeof(uint32_t));
		}
	}
	out.write(reinterpret_cast<const char*>(packed), sizeof(packed));
}

} // namespace

int main(int argc, char** argv) {
#ifdef WITH_NANOVDB
	Args args;
	if (!parse_args(argc, argv, args)) {
		return 1;
	}

	try {
		std::ofstream out(args.output_bin, std::ios::out | std::ios::binary | std::ios::trunc);
		if (!out.is_open()) {
			fprintf(stderr, "Error: cannot open output file '%s'\n", args.output_bin.c_str());
			return 1;
		}

		nanovdb::Coord min;
		nanovdb::Coord max;
		uint64_t entries = 0u;

		if (!args.use_channel_files) {
			auto handle = nanovdb::io::readGrid<nanovdb::HostBuffer>(args.input_nvdb);

			if (handle.gridCount() == 0) {
				fprintf(stderr, "Error: no grids found in '%s'\n", args.input_nvdb.c_str());
				return 1;
			}
			if (args.grid_index >= handle.gridCount()) {
				fprintf(stderr, "Error: --grid-index %u out of range [0, %u)\n", args.grid_index, handle.gridCount());
				return 1;
			}

			const auto* meta = handle.gridMetaData(args.grid_index);
			if (!meta) {
				fprintf(stderr, "Error: failed to access grid metadata\n");
				return 1;
			}

			const auto* grid = handle.grid<nanovdb::Vec3f>(args.grid_index);
			if (!grid) {
				fprintf(stderr,
					"Error: selected grid #%u ('%s') is not a Vec3f grid\n",
					args.grid_index,
					meta->shortGridName());
				return 1;
			}

			if (meta->isEmpty()) {
				fprintf(stderr, "Error: selected grid is empty\n");
				return 1;
			}

			const auto& bbox = meta->indexBBox();
			min = bbox.min();
			max = bbox.max();

			const auto& tree = grid->tree();
			for (int z = min[2]; z <= max[2]; ++z) {
				for (int y = min[1]; y <= max[1]; ++y) {
					for (int x = min[0]; x <= max[0]; ++x) {
						const nanovdb::Coord ijk(x, y, z);
						if (args.mode == CombineMode::Zip && !tree.isActive(ijk)) {
							continue;
						}
						const auto v = tree.getValue(ijk);
						write_float3_le(out, v[0], v[1], v[2]);
						++entries;
					}
				}
			}

			printf("Grid: #%u (%s)\n", args.grid_index, meta->shortGridName());
		} else {
			auto handle_x = nanovdb::io::readGrid<nanovdb::HostBuffer>(args.channel_inputs[0]);
			auto handle_y = nanovdb::io::readGrid<nanovdb::HostBuffer>(args.channel_inputs[1]);
			auto handle_z = nanovdb::io::readGrid<nanovdb::HostBuffer>(args.channel_inputs[2]);

			if (handle_x.gridCount() == 0 || handle_y.gridCount() == 0 || handle_z.gridCount() == 0) {
				fprintf(stderr, "Error: at least one --channel-files input has no grids\n");
				return 1;
			}
			if (args.grid_index >= handle_x.gridCount() || args.grid_index >= handle_y.gridCount() || args.grid_index >= handle_z.gridCount()) {
				fprintf(stderr, "Error: --grid-index %u is out of range for at least one channel grid\n", args.grid_index);
				return 1;
			}

			const auto* meta_x = handle_x.gridMetaData(args.grid_index);
			const auto* meta_y = handle_y.gridMetaData(args.grid_index);
			const auto* meta_z = handle_z.gridMetaData(args.grid_index);
			if (!meta_x || !meta_y || !meta_z) {
				fprintf(stderr, "Error: failed to access grid metadata for --channel-files\n");
				return 1;
			}

			if (meta_x->isEmpty() || meta_y->isEmpty() || meta_z->isEmpty()) {
				fprintf(stderr, "Error: at least one selected channel grid is empty\n");
				return 1;
			}

			const auto* grid_x = handle_x.grid<float>(args.grid_index);
			const auto* grid_y = handle_y.grid<float>(args.grid_index);
			const auto* grid_z = handle_z.grid<float>(args.grid_index);
			if (!grid_x || !grid_y || !grid_z) {
				fprintf(stderr, "Error: selected channel grid #%u must be float in all three input files\n", args.grid_index);
				return 1;
			}

			const auto& bbox_x = meta_x->indexBBox();
			const auto& bbox_y = meta_y->indexBBox();
			const auto& bbox_z = meta_z->indexBBox();
			if (!bbox_equal(bbox_x, bbox_y) || !bbox_equal(bbox_x, bbox_z)) {
				fprintf(stderr, "Error: indexBBox mismatch across --channel-files inputs\n");
				return 1;
			}

			min = bbox_x.min();
			max = bbox_x.max();

			const auto& tree_x = grid_x->tree();
			const auto& tree_y = grid_y->tree();
			const auto& tree_z = grid_z->tree();
			for (int z = min[2]; z <= max[2]; ++z) {
				for (int y = min[1]; y <= max[1]; ++y) {
					for (int x = min[0]; x <= max[0]; ++x) {
						const nanovdb::Coord ijk(x, y, z);
						if (args.mode == CombineMode::Zip) {
							const bool active = tree_x.isActive(ijk) || tree_y.isActive(ijk) || tree_z.isActive(ijk);
							if (!active) {
								continue;
							}
						}

						write_float3_le(out,
							tree_x.getValue(ijk),
							tree_y.getValue(ijk),
							tree_z.getValue(ijk));
						++entries;
					}
				}
			}

			printf("Grid: #%u (x=%s, y=%s, z=%s)\n",
				args.grid_index,
				meta_x->shortGridName(),
				meta_y->shortGridName(),
				meta_z->shortGridName());
			printf("Channel files: x=%s y=%s z=%s\n",
				args.channel_inputs[0].c_str(),
				args.channel_inputs[1].c_str(),
				args.channel_inputs[2].c_str());
		}

		if (!out.good()) {
			fprintf(stderr, "Error: failed while writing '%s'\n", args.output_bin.c_str());
			return 1;
		}

		const uint64_t byte_size = entries * 3u * sizeof(float);
		printf("Mode: %s\n", args.mode == CombineMode::Mesh ? "mesh" : "zip");
		printf("Input mode: %s\n", args.use_channel_files ? "channel-files" : "single-vec3");
		printf("Index bbox: min=(%d,%d,%d) max=(%d,%d,%d)\n", min[0], min[1], min[2], max[0], max[1], max[2]);
		printf("Wrote %llu float3 entries to: %s\n", static_cast<unsigned long long>(entries), args.output_bin.c_str());
		printf("Byte size: %llu\n", static_cast<unsigned long long>(byte_size));
		return 0;
	} catch (const std::exception& exc) {
		fprintf(stderr, "Error: %s\n", exc.what());
		return 1;
	}
#else
	(void)argc;
	(void)argv;
	fprintf(stderr, "Error: nano2float3 was built without WITH_NANOVDB\n");
	return 1;
#endif
}