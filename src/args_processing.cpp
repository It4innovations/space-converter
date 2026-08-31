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

#include "args_processing.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

namespace space_converter {

	/**
	 * @brief Display usage information and command-line options, then exit.
	 *
	 * @param exit_code  Process exit code (0 for -h/--help, 1 for argument errors).
	 * @param print      When false the text is suppressed (used so that only MPI
	 *                   rank 0 prints while all ranks still exit consistently).
	 */
	static void usage(int exit_code = 0, bool print = true)
	{
		if (print) {
		std::cout << "./space_converter --data-type [GADGET, GADGET_SIMPLE, CHANGA_TIPSY, CHANGA_NCHILADA, HACC_GENERICIO, HACC_BIN, IPIC3D_HDF5, PLUTO_VTK] <options> <args>" << std::endl;

		// === General Options ===
		std::cout << "\noptions (defaults in brackets):" << std::endl;
		std::cout << "\t--grid-dim N                       : Grid resolution per axis [100]" << std::endl;
		std::cout << "\t-o, --output-path DIR              : Output directory for generated files" << std::endl;
		std::cout << "\t-f, --output-file FILE             : Output file path (base name without extension; overrides the automatic <type>_<dataset> naming, frame/rank suffixes and the format extension are still appended; a relative FILE is placed inside --output-path)" << std::endl;
		std::cout << "\t--port N                           : TCP port the server listens on [5000]" << std::endl;
		std::cout << "\t--info                             : Print dataset info and exit" << std::endl;
		std::cout << "\t--nanovdb                          : Output NanoVDB instead of OpenVDB" << std::endl;
		std::cout << "\t--cub                              : Output CUB format (GPU builds)" << std::endl;
		std::cout << "\t--dense-file N                     : Also dump the dense grid: 0=off, 1=RAW, 2=VTI [0]" << std::endl;
		std::cout << "\t--anim START END STEP              : Process an animation range of snapshots" << std::endl;
		std::cout << "\t--anim-merge                       : Merge all animation frames into one output" << std::endl;
		std::cout << "\t--raw-particles N                  : Export raw particles: 0=none, 1=RAW, 2=VDB points, 3=VTP [0]" << std::endl;
		std::cout << "\t--export, --export-data TYPE BLOCK : Batch mode: extract particle TYPE / data BLOCK and exit" << std::endl;
		std::cout << "\t--dense-type N                     : SPH splat kernel: 0=off(sparse), 1=Cubic, 2=Quintic, 3..6=WendlandC2..C8 [0]" << std::endl;
		std::cout << "\t--dense-norm N                     : Density normalization: 0=none, 1=count, 2=SPH, 3=voxel volume [0]" << std::endl;
		std::cout << "\t--simple-density                   : Deposit weight 1 per particle (implies dense mode)" << std::endl;
		std::cout << "\t--bbox x1 y1 z1 x2 y2 z2           : Axis-aligned zoom box in object space [full box]" << std::endl;
		std::cout << "\t--bbox-sphere x y z r              : Keep only particles inside this sphere (world space)" << std::endl;
		std::cout << "\t--offset-position X Y Z            : Subtract this offset from all particle positions" << std::endl;
		std::cout << "\t--radius-const R                   : Fixed particle radius in voxel units (0 = use per-particle radius) [0]" << std::endl;
		std::cout << "\t--no-norm-value                    : Store vector blocks as 3-component grids instead of magnitudes" << std::endl;

		// === Neighbor Search Options ===
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
		std::cout << "\t--calc-radius-neigh N              : Compute per-particle radius from N nearest neighbors" << std::endl;
		std::cout << "\t--calc-radius-neigh-rho-kernel N   : SPH kernel for the k-NN density estimate [6=WendlandC6]" << std::endl;
#endif
		std::cout << "\t--calc-radius-neigh-file FILE      : Read/write precomputed radii (per-type .N.bin files)" << std::endl;

#ifdef WITH_CUDAKDTREE
		std::cout << "\t--cudakdtree                       : Use cudaKDTree (GPU) for the k-NN search" << std::endl;
		std::cout << "\t--cudakdtree-cpu                   : Use cudaKDTree host build for the k-NN search" << std::endl;
		std::cout << "\t--dense-loop-over-voxels           : Dense grid: loop over voxels, query particles via KD-tree (needs --calc-radius-neigh)" << std::endl;
#endif
#ifdef WITH_NANOFLANN
		std::cout << "\t--nanoflann                        : Use nanoflann (CPU) for the k-NN search" << std::endl;
#endif

		// === Advanced Options ===
		std::cout << "\t--skip-cache-manager               : Stream particles from the reader instead of caching them (CPU only, lower memory)" << std::endl;
		std::cout << "\t--sort-by-radius                   : Sort particles by radius before splatting" << std::endl;
		std::cout << "\t--sort-by-non-overlap              : Sort particles spatially (Morton) to reduce atomic contention" << std::endl;
		std::cout << "\t--save-mpi-rank                    : Additionally save each rank's partial grid" << std::endl;
#ifdef WITH_GPU_CUDA
		std::cout << "\t--gpu                              : Enable GPU acceleration (CUDA or HIP, depending on the build)" << std::endl;
#endif
		std::cout << "\t-h, --help                         : Display this usage information" << std::endl;

		std::cout << "\nNote: without --export/--info the converter starts as a TCP server for the bspace" << std::endl;
		std::cout << "Blender addon. Some parameters exist only in that protocol and cannot be set here:" << std::endl;
		std::cout << "particle radius multiplier, value filters, object size (fixed 1000 in batch mode)." << std::endl;

		// === Format-Specific Arguments ===
		std::cout << "\nGADGET args:" << std::endl;
		std::cout << "\t--param-file FILE           : GADGET parameter file" << std::endl;
		std::cout << "\t--gadget-file FILE          : GADGET snapshot file path" << std::endl;
		std::cout << "\t--max-mem-size N            : Maximum memory size (MB)" << std::endl;
		std::cout << "\t--buffer-size N             : Buffer size for particle loading" << std::endl;
		std::cout << "\t--part-alloc-factor N       : Particle allocation factor" << std::endl;
		std::cout << "\t--bh-count N                : Number of black holes in simulation" << std::endl;

		std::cout << "\nGADGET_SIMPLE args:" << std::endl;
		std::cout << "\t--gadget-file FILE          : GADGET snapshot file path (simplified reader)" << std::endl;

		std::cout << "\nCHANGA_TIPSY args:" << std::endl;
		std::cout << "\t--tipsy-file FILE           : TIPSY format file path" << std::endl;
		std::cout << "\t--filter-in FILE            : Aux-file prefix for additional per-particle fields" << std::endl;

		std::cout << "\nCHANGA_NCHILADA args:" << std::endl;
		std::cout << "\t--nc-dir DIR                : NChilada data directory" << std::endl;

		std::cout << "\nHACC_GENERICIO args:" << std::endl;
		std::cout << "\t--genericio-file FILE       : GenericIO format file path" << std::endl;
		std::cout << "\t--pos-names x y z           : Position field names in the GenericIO file" << std::endl;
		std::cout << "\t--vel-names vx vy vz        : Velocity field names in the GenericIO file" << std::endl;
		std::cout << "\t--mass-name NAME            : Mass field name (otherwise mass = 0)" << std::endl;
		std::cout << "\t--rho-name NAME             : Density field name (otherwise rho = 0)" << std::endl;
		std::cout << "\t--hsml-name NAME            : Smoothing-length field name (otherwise hsml = 0)" << std::endl;

		std::cout << "\nHACC_BIN args:" << std::endl;
		std::cout << "\t--haccbin-file FILE         : HACC binary format file path" << std::endl;

		std::cout << "\nIPIC3D_HDF5 args:" << std::endl;
		std::cout << "\t--hdf5-file FILE            : iPIC3D HDF5 file path (all species, fields and moments at the latest cycle are loaded automatically)" << std::endl;
		std::cout << "\t--num-files N               : Number of files to load" << std::endl;
		std::cout << "\t--settings-file FILE        : Companion settings.hdf (grid geometry); defaults to settings.hdf next to --hdf5-file" << std::endl;

		std::cout << "\nPLUTO_VTK args:" << std::endl;
		std::cout << "\t--vtk-file FILE             : PLUTO VTK rectilinear grid file path" << std::endl;
		std::cout << "\t--scalar-names NAME...      : Cell arrays to load (default: all)" << std::endl;

		// === Examples ===
		std::cout << "\nexamples:" << std::endl;
		std::cout << "  # Batch extraction of block 1 for particle type 0 into OpenVDB:" << std::endl;
		std::cout << "  space_converter --data-type CHANGA_TIPSY --tipsy-file snap.00995 \\" << std::endl;
		std::cout << "      --output-path /tmp/out --grid-dim 500 --export-data 0 1" << std::endl;
		std::cout << "  # Dense SPH grid (WendlandC6) of a zoom region:" << std::endl;
		std::cout << "  space_converter --data-type GADGET_SIMPLE --gadget-file snap_091 --output-path /tmp/out \\" << std::endl;
		std::cout << "      --grid-dim 500 --export-data 0 1 --dense-type 6 --bbox 450 500 470 510 560 530" << std::endl;
		std::cout << "  # Remote server for the bspace Blender addon:" << std::endl;
		std::cout << "  space_converter --data-type GADGET_SIMPLE --gadget-file snap_091 --output-path /tmp/out --port 5000" << std::endl;
		}

		exit(exit_code);
	}

	namespace {

		/// Reader-specific flags. They are parsed a second time inside each reader's
		/// init_lib(); parse_args only needs to know them (and their arity) so that
		/// unknown arguments can be rejected. value_count == -1 means "consume values
		/// until the next token starting with '-'" (e.g. --scalar-names).
		struct ReaderFlag {
			const char* name;
			int value_count;
		};

		const ReaderFlag kReaderFlags[] = {
			// GADGET / GADGET_SIMPLE
			{ "--param-file", 1 },
			{ "--gadget-file", 1 },
			{ "--max-mem-size", 1 },
			{ "--buffer-size", 1 },
			{ "--part-alloc-factor", 1 },
			{ "--bh-count", 1 },
			// CHANGA
			{ "--tipsy-file", 1 },
			{ "--filter-in", 1 },
			{ "--nc-dir", 1 },
			// HACC
			{ "--genericio-file", 1 },
			{ "--pos-names", 3 },
			{ "--vel-names", 3 },
			{ "--mass-name", 1 },
			{ "--rho-name", 1 },
			{ "--hsml-name", 1 },
			{ "--haccbin-file", 1 },
			// iPIC3D
			{ "--hdf5-file", 1 },
			{ "--num-files", 1 },
			{ "--settings-file", 1 },
			// PLUTO
			{ "--vtk-file", 1 },
			{ "--scalar-names", -1 },
		};

		int g_parse_rank = 0;  ///< MPI rank, so parse errors are printed once

		/// Print an argument error (rank 0 only) and exit on all ranks.
		[[noreturn]] void arg_error(const std::string& message)
		{
			if (g_parse_rank == 0) {
				fprintf(stderr, "space_converter: %s\n", message.c_str());
				fprintf(stderr, "space_converter: use -h for usage\n");
			}
			exit(1);
		}

		/// Fetch the next value for `flag`, erroring out if the command line ends.
		const char* next_arg(int& i, int argc, char** argv, const char* flag)
		{
			if (i + 1 >= argc) {
				arg_error(std::string("missing value for ") + flag);
			}
			return argv[++i];
		}

		int parse_int(int& i, int argc, char** argv, const char* flag)
		{
			const char* v = next_arg(i, argc, argv, flag);
			try {
				return std::stoi(v);
			}
			catch (const std::exception&) {
				arg_error(std::string("invalid integer '") + v + "' for " + flag);
			}
		}

		float parse_float(int& i, int argc, char** argv, const char* flag)
		{
			const char* v = next_arg(i, argc, argv, flag);
			try {
				return std::stof(v);
			}
			catch (const std::exception&) {
				arg_error(std::string("invalid number '") + v + "' for " + flag);
			}
		}

		double parse_double(int& i, int argc, char** argv, const char* flag)
		{
			const char* v = next_arg(i, argc, argv, flag);
			try {
				return std::stod(v);
			}
			catch (const std::exception&) {
				arg_error(std::string("invalid number '") + v + "' for " + flag);
			}
		}

		int parse_int_range(int& i, int argc, char** argv, const char* flag, int lo, int hi)
		{
			int v = parse_int(i, argc, argv, flag);
			if (v < lo || v > hi) {
				arg_error(std::string(flag) + " must be in [" + std::to_string(lo) + ", " + std::to_string(hi) + "], got " + std::to_string(v));
			}
			return v;
		}

		/// Warn (rank 0 only) about a flag that is accepted but has no effect in this build.
		void warn_unsupported(const char* flag, const char* reason)
		{
			if (g_parse_rank == 0) {
				fprintf(stderr, "space_converter: warning: %s ignored (%s)\n", flag, reason);
			}
		}

	} // namespace

	void parse_args(FromCL& from_cl, common::SpaceData &space_data, int argc, char** argv)
	{
		g_parse_rank = from_cl.world_rank;

		// Require at least the program name + --data-type + its value
		if (argc < 3) {
			usage(1, from_cl.world_rank == 0);
		}

		// Flags that decide the extraction mode; the mode itself is derived once at
		// the end so that the result does not depend on the flag order.
		bool flag_raw_particles = false;
		bool flag_simple_density = false;

		for (int i = 1; i < argc; i++) {
			const std::string arg = argv[i];

			// === Core Configuration ===
			if (arg == "--data-type") {
				from_cl.data_type = next_arg(i, argc, argv, "--data-type");
			}
			else if (arg == "--grid-dim") {
				space_data.bbox_dim = parse_int(i, argc, argv, "--grid-dim");
				if (space_data.bbox_dim < 1) {
					arg_error("--grid-dim must be >= 1");
				}
				if (space_data.bbox_dim > 4096 && from_cl.world_rank == 0) {
					fprintf(stderr, "space_converter: warning: --grid-dim %d is very large (dense grids need dim^3 floats)\n", space_data.bbox_dim);
				}
			}
			else if (arg == "-o" || arg == "--output-path") {
				from_cl.output_path = next_arg(i, argc, argv, "--output-path");
			}
			else if (arg == "-f" || arg == "--output-file") {
				from_cl.output_file = next_arg(i, argc, argv, "--output-file");
			}
			else if (arg == "--port") {
				from_cl.port = parse_int_range(i, argc, argv, "--port", 1, 65535);
			}
			else if (arg == "--info") {
				// Display dataset information only (no conversion)
				from_cl.info = true;
				from_cl.remote = false;
			}
			// === Output Format Options ===
			else if (arg == "--nanovdb") {
				from_cl.use_nanovdb = true;
			}
			else if (arg == "--cub") {
				from_cl.use_cub = true;
			}
			else if (arg == "--dense-file") {
				space_data.extracted_dense_type = (common::SpaceData::ExtractedDenseType)parse_int_range(i, argc, argv, "--dense-file", 0, 2);
			}

			// === Animation Processing ===
			else if (arg == "--anim") {
				if (space_data.anim_type == common::SpaceData::AnimType::eNone)
					space_data.anim_type = common::SpaceData::AnimType::eAllPath;

				space_data.anim_start = parse_int(i, argc, argv, "--anim");
				space_data.anim_end = parse_int(i, argc, argv, "--anim");
				space_data.anim_step = parse_int(i, argc, argv, "--anim");
				if (space_data.anim_step < 1) {
					arg_error("--anim STEP must be >= 1");
				}
				if (space_data.anim_start > space_data.anim_end) {
					arg_error("--anim START must be <= END");
				}
			}
			else if (arg == "--anim-merge") {
				space_data.anim_type = common::SpaceData::AnimType::eAllMerge;
			}

			// === Particle Export Options ===
			else if (arg == "--raw-particles") {
				flag_raw_particles = true;
				space_data.extracted_particle_type = (common::SpaceData::ExtractedParticleType)parse_int_range(i, argc, argv, "--raw-particles", 0, 3);
			}
			else if (arg == "--save-mpi-rank") {
				from_cl.use_save_mpirank = true;
			}
			else if (arg == "--export-data" || arg == "--export") {
				space_data.particle_type = parse_int(i, argc, argv, "--export-data");
				space_data.block_name_id = parse_int(i, argc, argv, "--export-data");
				if (space_data.particle_type < 0 || space_data.block_name_id < 0) {
					arg_error("--export-data TYPE and BLOCK must be >= 0");
				}
				from_cl.remote = false;
			}
			// === Density Computation Options ===
			else if (arg == "--dense-type") {
				space_data.dense_type = (common::SpaceData::DenseType)parse_int_range(i, argc, argv, "--dense-type", 0, 6);
			}
			else if (arg == "--dense-norm") {
				space_data.dense_norm = (common::SpaceData::DenseNorm)parse_int_range(i, argc, argv, "--dense-norm", 0, 3);
			}
			else if (arg == "--simple-density") {
				flag_simple_density = true;
				space_data.use_simple_density = true;
			}

			// === Spatial Filtering Options ===
			else if (arg == "--bbox-sphere") {
				space_data.use_bbox_sphere = true;
				space_data.bbox_sphere_pos[0] = parse_float(i, argc, argv, "--bbox-sphere");
				space_data.bbox_sphere_pos[1] = parse_float(i, argc, argv, "--bbox-sphere");
				space_data.bbox_sphere_pos[2] = parse_float(i, argc, argv, "--bbox-sphere");
				space_data.bbox_sphere_r = parse_float(i, argc, argv, "--bbox-sphere");
				if (space_data.bbox_sphere_r <= 0.0f) {
					arg_error("--bbox-sphere radius must be > 0");
				}
			}
			else if (arg == "--bbox") {
				space_data.bbox_min[0] = parse_float(i, argc, argv, "--bbox");
				space_data.bbox_min[1] = parse_float(i, argc, argv, "--bbox");
				space_data.bbox_min[2] = parse_float(i, argc, argv, "--bbox");

				space_data.bbox_max[0] = parse_float(i, argc, argv, "--bbox");
				space_data.bbox_max[1] = parse_float(i, argc, argv, "--bbox");
				space_data.bbox_max[2] = parse_float(i, argc, argv, "--bbox");

				for (int a = 0; a < 3; a++) {
					if (space_data.bbox_min[a] >= space_data.bbox_max[a]) {
						arg_error("--bbox: min must be < max on every axis");
					}
				}
			}
			else if (arg == "--no-norm-value") {
				space_data.use_norm_value = false;
			}
			else if (arg == "--offset-position") {
				space_data.offset_position[0] = parse_float(i, argc, argv, "--offset-position");
				space_data.offset_position[1] = parse_float(i, argc, argv, "--offset-position");
				space_data.offset_position[2] = parse_float(i, argc, argv, "--offset-position");
			}
			else if (arg == "--radius-const") {
				space_data.particle_radius_const = parse_double(i, argc, argv, "--radius-const");
				if (space_data.particle_radius_const < 0.0) {
					arg_error("--radius-const must be >= 0");
				}
			}
			// === Neighbor Search Configuration ===
			else if (arg == "--calc-radius-neigh") {
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
				space_data.calc_radius_neigh = parse_int(i, argc, argv, "--calc-radius-neigh");
				if (space_data.calc_radius_neigh < 1) {
					arg_error("--calc-radius-neigh must be >= 1");
				}
#else
				(void)parse_int(i, argc, argv, "--calc-radius-neigh");
				warn_unsupported("--calc-radius-neigh", "built without cudaKDTree/nanoflann");
#endif
			}
			else if (arg == "--calc-radius-neigh-rho-kernel") {
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
				space_data.calc_radius_neigh_rho_kernel = (common::SpaceData::DenseType)parse_int_range(i, argc, argv, "--calc-radius-neigh-rho-kernel", 1, 6);
#else
				(void)parse_int(i, argc, argv, "--calc-radius-neigh-rho-kernel");
				warn_unsupported("--calc-radius-neigh-rho-kernel", "built without cudaKDTree/nanoflann");
#endif
			}
			else if (arg == "--calc-radius-neigh-file") {
				space_data.calc_radius_neigh_file = next_arg(i, argc, argv, "--calc-radius-neigh-file");
			}

			// === Neighbor Search Implementation ===
			else if (arg == "--cudakdtree") {
#ifdef WITH_CUDAKDTREE
				from_cl.use_cudakdtree = true;
#else
				warn_unsupported("--cudakdtree", "built without cudaKDTree");
#endif
			}
			else if (arg == "--cudakdtree-cpu") {
#ifdef WITH_CUDAKDTREE
				from_cl.use_cudakdtree_cpu = true;
#else
				warn_unsupported("--cudakdtree-cpu", "built without cudaKDTree");
#endif
			}
			else if (arg == "--dense-loop-over-voxels") {
#ifdef WITH_CUDAKDTREE
				from_cl.use_dense_loop_over_voxels = true;
#else
				warn_unsupported("--dense-loop-over-voxels", "built without cudaKDTree");
#endif
			}
			else if (arg == "--nanoflann") {
#ifdef WITH_NANOFLANN
				from_cl.use_nanoflann = true;
#else
				warn_unsupported("--nanoflann", "built without nanoflann");
#endif
			}
			else if (arg == "--gpu") {
#ifdef WITH_GPU_CUDA
				from_cl.use_gpu = true;
#else
				warn_unsupported("--gpu", "built without GPU support");
#endif
			}
			else if (arg == "--skip-cache-manager") {
				from_cl.skip_cache_manager = true;
			}
			else if (arg == "--sort-by-radius") {
				from_cl.use_sort_by_radius = true;
			}
			else if (arg == "--sort-by-non-overlap") {
				from_cl.use_sort_by_non_overlap = true;
			}
			else if (arg == "-h" || arg == "--help") {
				usage(0, from_cl.world_rank == 0);
			}
			else {
				// Reader-specific flags: recognized here (with their arity) so that the
				// unknown-argument check below stays valid; parsed again in the reader.
				bool is_reader_flag = false;
				for (const ReaderFlag& rf : kReaderFlags) {
					if (arg == rf.name) {
						is_reader_flag = true;
						if (rf.value_count >= 0) {
							for (int v = 0; v < rf.value_count; v++) {
								next_arg(i, argc, argv, rf.name);
							}
						}
						else {
							// Variable arity: consume values until the next flag
							while (i + 1 < argc && argv[i + 1][0] != '-') {
								++i;
							}
						}
						break;
					}
				}

				if (!is_reader_flag) {
					arg_error(std::string("unknown argument '") + arg + "'");
				}
			}
		}

		// === Derive the extraction mode (independent of flag order) ===
		bool wants_dense = flag_simple_density || space_data.dense_type != common::SpaceData::DenseType::eNone;
		if (flag_raw_particles && wants_dense) {
			arg_error("--raw-particles cannot be combined with --dense-type/--simple-density");
		}
		if (flag_raw_particles) {
			space_data.extracted_type = common::SpaceData::ExtractedType::eParticle;
		}
		else if (wants_dense) {
			space_data.extracted_type = common::SpaceData::ExtractedType::eDense;
			// --simple-density alone still needs a splat kernel
			if (space_data.dense_type == common::SpaceData::DenseType::eNone) {
				space_data.dense_type = common::SpaceData::DenseType::eCubic;
			}
		}
		else {
			space_data.extracted_type = common::SpaceData::ExtractedType::eSparse;
		}
	}

	void FromCL::print_info() const
	{
		std::cout << "\n=== FromCL Configuration ==="<< std::endl;
		std::cout << "data_type: " << data_type << std::endl;
		std::cout << "output_path: " << output_path << std::endl;
		std::cout << "output_file: " << (output_file.empty() ? "(auto)" : output_file) << std::endl;
		std::cout << "port: " << port << std::endl;
		std::cout << "info: " << (info ? "true" : "false") << std::endl;
		std::cout << "remote: " << (remote ? "true" : "false") << std::endl;
		std::cout << "use_nanovdb: " << (use_nanovdb ? "true" : "false") << std::endl;
		std::cout << "use_cub: " << (use_cub ? "true" : "false") << std::endl;
		std::cout << "use_save_mpirank: " << (use_save_mpirank ? "true" : "false") << std::endl;
#ifdef WITH_CUDAKDTREE
		std::cout << "use_cudakdtree: " << (use_cudakdtree ? "true" : "false") << std::endl;
		std::cout << "use_cudakdtree_cpu: " << (use_cudakdtree_cpu ? "true" : "false") << std::endl;
		std::cout << "use_dense_loop_over_voxels: " << (use_dense_loop_over_voxels ? "true" : "false") << std::endl;
#endif
#ifdef WITH_NANOFLANN
		std::cout << "use_nanoflann: " << (use_nanoflann ? "true" : "false") << std::endl;
#endif
#ifdef WITH_GPU_CUDA
		std::cout << "use_gpu: " << (use_gpu ? "true" : "false") << std::endl;
#endif
		std::cout << "skip_cache_manager: " << (skip_cache_manager ? "true" : "false") << std::endl;
		std::cout << "use_sort_by_radius: " << (use_sort_by_radius ? "true" : "false") << std::endl;
		std::cout << "use_sort_by_non_overlap: " << (use_sort_by_non_overlap ? "true" : "false") << std::endl;
		std::cout << "world_rank: " << world_rank << std::endl;
		std::cout << "world_size: " << world_size << std::endl;
		std::cout << "============================\n" << std::endl;
	}

} // namespace space_converter
