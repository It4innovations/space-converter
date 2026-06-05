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

// ============================================================================
// Header includes
// ============================================================================

#include "data_processing.h"
#include "data_communication.h"
#include "sparse_common.h"

// Data format converters
#include "gadget/gadget_convert_vdb.h"
#include "gadget/gadget_simple_convert_vdb.h"
#include "changa/changa_nchilada_convert_vdb.h"
#include "changa/changa_tipsy_convert_vdb.h"
#include "common/convert_vdb.h"

#ifdef WITH_HACC
#	include "hacc/haccbin_convert_vdb.h"
#	include "hacc/hacc_genericio_convert_vdb.h"
#endif

#ifdef WITH_IPIC3D
#	include "ipic3d/ipic3d_hdf5_convert_vdb.h"
#endif

#ifdef WITH_PLUTO
#	include "pluto/pluto_vtk_convert_vdb.h"
#endif

#ifdef WITH_OPENMP
#	include <omp.h>
#endif

#ifdef WITH_OPENVDB

#include <nanovdb/tools/NanoToOpenVDB.h>

#include <openvdb/openvdb.h>
#include <openvdb/points/PointConversion.h>
#include <openvdb/points/PointDataGrid.h>
#include <openvdb/tools/Count.h>

#include <fstream>
#include <filesystem>

#endif

#ifdef WITH_VTK
#	include <vtkCellArray.h>
#	include <vtkFloatArray.h>
#	include <vtkImageData.h>
#	include <vtkPointData.h>
#	include <vtkPoints.h>
#	include <vtkPolyData.h>
#	include <vtkSmartPointer.h>
#	include <vtkXMLImageDataWriter.h>
#	include <vtkXMLPPolyDataWriter.h>
#	include <vtkXMLPolyDataWriter.h>
#endif

#include "utility/logging.h"

#ifdef WITH_GPU_CUDA
#	include "utility/gpu_utility.h"
#endif

#ifdef _WIN32
#	undef max              // Disable Windows macro to avoid conflicts with std::max
#	undef min              // Disable Windows macro to avoid conflicts with std::min
#	include <algorithm>    // For std::max and std::min
#endif

// ============================================================================
// Macros
// ============================================================================

// /**
//  * @brief Synchronization barrier for all MPI processes
//  * @details Ensures all processes reach this point before any can proceed
//  */
// #define CALL_MPI_BARRIER \
// 	{ \
// 		MPI_Barrier(MPI_COMM_WORLD); \
// 	}

namespace space_converter {
	/**
	 * @brief Test function for converter development and debugging
	 * @param argc Command line argument count
	 * @param argv Command line arguments
	 * @param from_cl Command line parameters structure
	 * @note This function is currently disabled and contains various test code
	 *       for VDB conversion, GenericIO reading, and format testing
	 */
	void test_converter(int argc, char** argv, space_converter::FromCL& from_cl)
	{
		//TODO: Implement test cases for converter functionality
	}

	/**
	 * @brief Initialize the appropriate data converter based on input data type
	 * @param argc Command line argument count
	 * @param argv Command line arguments
	 * @param from_cl Command line parameters structure
	 * @param space_data Spatial data configuration
	 * @return Pointer to the initialized converter (caller must delete)
	 * @throws std::runtime_error if data type is unknown
	 * 
	 * @details Supports multiple data formats:
	 *   - GADGET: Gadget simulation format
	 *   - GADGET_SIMPLE: Simplified Gadget format
	 *   - CHANGA_TIPSY: ChaNGa Tipsy format
	 *   - CHANGA_NCHILADA: ChaNGa NChilada format
	 *   - HACC_GENERICIO: HACC GenericIO format (if compiled with support)
	 *   - HDF5: HDF5 format (if compiled with support)
	 *   - HACC_BIN: HACC binary format
	 *   - IPIC3D_HDF5: iPIC3D HDF5 format (if compiled with support)
	 *   - PLUTO_VTK: PLUTO VTK rectilinear grid format (if compiled with support)
	 */
	common::vdb::ConvertVDBBase* init_converter(int argc, char** argv, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		common::vdb::ConvertVDBBase* convert_vdb_base = nullptr;

		if (from_cl.world_rank == 0)
			printf("Data: %s\n", from_cl.data_type.c_str());

		// Initialize the appropriate converter based on data type
		if (from_cl.data_type == "GADGET") {
			convert_vdb_base = new gadget::ConvertVDBGadget();
		}
		else if (from_cl.data_type == "GADGET_SIMPLE") {
			convert_vdb_base = new gadget_simple::ConvertVDBGadgetSimple();
		}
		else if (from_cl.data_type == "CHANGA_TIPSY") {
			convert_vdb_base = new changa::ConvertVDBChangaTipsy();
		}
		else if (from_cl.data_type == "CHANGA_NCHILADA") {
			convert_vdb_base = new changa::ConvertVDBChangaNchilada();
		}
#ifdef WITH_HACC
		else if (from_cl.data_type == "HACC_GENERICIO") {
			convert_vdb_base = new genericio::ConvertVDBHACCGenericIO();
		}
		else if (from_cl.data_type == "HACC_BIN") {
			convert_vdb_base = new haccbin::ConvertVDBHACCBin();
		}
#endif
#ifdef WITH_IPIC3D
		else if (from_cl.data_type == "IPIC3D_HDF5") {
			convert_vdb_base = new ipic3d::ConvertVDBIPIC3DHDF5();
		}
#endif		
#ifdef WITH_PLUTO
		else if (from_cl.data_type == "PLUTO_VTK") {
			convert_vdb_base = new plutovtk::ConvertVDBPlutoVTK();
		}
#endif
		else {
			throw std::runtime_error("Unknown data type [GADGET, CHANGA_TIPSY, CHANGA_NCHILADA, HACC_GENERICIO, HACC_BIN, IPIC3D_HDF5, PLUTO_VTK]");
		}

		// Other params
#ifdef WITH_GPU_CUDA
		if (from_cl.use_gpu_cuda) {
			printf("Using GPU CUDA for computations\n");
		}		
		convert_vdb_base->cache_manager.use_gpu_cuda = from_cl.use_gpu_cuda;
#endif
		// MPI
		convert_vdb_base->cache_manager.world_rank = from_cl.world_rank;
		convert_vdb_base->cache_manager.world_size = from_cl.world_size;

		// CALL_MPI_BARRIER;
		// Initialize the converter library with MPI rank and size
		LOG_MeasureStart("init_lib");

		const char* converter_split_init_count = getenv("CONVERTER_SPLIT_INIT_COUNT");
		if (converter_split_init_count) {
			printf("WARNING (use CONVERTER_SPLIT_INIT_COUNT): This is a test version of the converter!\n");

			// Split initialization across process groups to reduce memory pressure
			// This is useful for very large scale runs where simultaneous initialization
			// can cause memory or I/O bottlenecks
			
			//double t_init = omp_get_wtime();

			int group_size = atoi(converter_split_init_count);
			int num_groups = (from_cl.world_size + group_size - 1) / group_size;

			for (int g = 0; g < num_groups; ++g) {
				int start = g * group_size;
				int end = std::min(start + group_size, from_cl.world_size);

				// Initialize library for processes in current group
				if (from_cl.world_rank >= start && from_cl.world_rank < end) {
					convert_vdb_base->init_lib(argc, argv, from_cl.world_rank, from_cl.world_size);
				}

				// Synchronize before next group proceeds
				// CALL_MPI_BARRIER;
				MPI_Barrier(MPI_COMM_WORLD);
			}				
			// CALL_MPI_BARRIER;
			// double t_init_end = omp_get_wtime();
			// if (from_cl.world_rank == 0) {
			// 	printf("rank #%d: CONVERTER_SPLIT_INIT time: %f\n", from_cl.world_rank, t_init_end - t_init);
			// }

			//printf("rank #%d: init_lib done\n", from_cl.world_rank);			
		}
		else {
			// CALL_MPI_BARRIER;
			// double t_init = omp_get_wtime();

			convert_vdb_base->init_lib(argc, argv, from_cl.world_rank, from_cl.world_size);

			// CALL_MPI_BARRIER;
			// double t_init_end = omp_get_wtime();
			// if (from_cl.world_rank == 0) {
			// 	printf("rank #%d: FINAL INIT time: %f\n", from_cl.world_rank, t_init_end - t_init);
			// }
		}

		// CALL_MPI_BARRIER;
		// Initialize the converter library with MPI rank and size
		LOG_MeasureStop("init_lib");

		if (from_cl.info) {
			return convert_vdb_base;
		}

		LOG_MeasureStart("find_particle_positions");
		convert_vdb_base->find_particle_positions();

#ifdef WITH_GPU_CUDA
		if (convert_vdb_base->cache_manager.use_gpu_cuda) {
			convert_vdb_base->cache_manager.copy_particles_to_gpu();
		}
#endif		

		LOG_MeasureStop("find_particle_positions");

		// Calculate particle radius using KD-tree algorithms
		// Radius is needed for proper density estimation and particle size
		LOG_MeasureStart("calculate_radius");

		// Determine if we should use cycling mode (no animation or using merged animation)
		bool use_cycling = space_data.anim_type == common::SpaceData::AnimType::eNone;		

		std::string calc_radius_neigh_file = space_data.calc_radius_neigh_file;

		// Calculate particle radii using CUDA KD-tree (GPU or CPU fallback)
#ifdef WITH_CUDAKDTREE
		if ((from_cl.use_cudakdtree || from_cl.use_cudakdtree_cpu) && (space_data.calc_radius_neigh > 0 || calc_radius_neigh_file.length() > 0)) {
			// Maximum radius for KNN search (can be set via environment variable)
			float max_radius = std::numeric_limits<float>::infinity();
			const char* converter_knn_max_radius = getenv("CONVERTER_KNN_MAX_RADIUS");
			if (converter_knn_max_radius) {
				max_radius = static_cast<float>(atof(converter_knn_max_radius));
				printf("Using max_radius: %f\n", max_radius);
			}
			common::SpaceData::DenseType rho_kernel = space_data.calc_radius_neigh_rho_kernel;
			convert_vdb_base->calculate_radius_by_cudakdtree(space_data.calc_radius_neigh, calc_radius_neigh_file, use_cycling, from_cl.use_cudakdtree_cpu, max_radius, rho_kernel);
		}
		else
#endif

		// Calculate particle radii using Nanoflann KD-tree (CPU)
#ifdef WITH_NANOFLANN		
		if (from_cl.use_nanoflann && (space_data.calc_radius_neigh > 0 || calc_radius_neigh_file.length() > 0)) {
			common::SpaceData::DenseType rho_kernel = space_data.calc_radius_neigh_rho_kernel;
			convert_vdb_base->calculate_radius_by_nanoflann(space_data.calc_radius_neigh, calc_radius_neigh_file, use_cycling, rho_kernel);
		}
		else
#endif
		// Load pre-computed radii from file if no KD-tree calculation requested
		if (calc_radius_neigh_file.length() > 0) {
			convert_vdb_base->read_radius_from_file(calc_radius_neigh_file);
		}

		// Override with constant radius if specified via environment variable
		//const char* converter_radius_particle_const = getenv("CONVERTER_RADIUS_PARTICLE_CONST");
		//if (converter_radius_particle_const) {
		//	convert_vdb_base->cache_manager.particle_radius_const = atof(converter_radius_particle_const);
		//}
		convert_vdb_base->cache_manager.particle_radius_const = space_data.particle_radius_const;

		// Enable sorting by radius if specified via command-line argument
		if (from_cl.use_sort_by_radius) {
			// Sorting methods
			//void sort_particles_by_radius_cpu();                        ///< Sort particle IDs by radius on CPU (ascending)
			//void sort_particles_by_radius_gpu();                        ///< Sort particle IDs by radius on GPU (ascending)
			//void sort_particles_by_radius_gpu_inplace();                ///< Sort particle IDs by radius on GPU using device pointers (no CPU<->GPU copy)

			if (convert_vdb_base->cache_manager.use_gpu_cuda) {
				convert_vdb_base->cache_manager.sort_particles_by_radius_gpu_inplace();
			}
			else {
				convert_vdb_base->cache_manager.sort_particles_by_radius_cpu();
			}			
		}

		// Enable spatial sorting to avoid atomic operations if specified via command-line argument
		if (from_cl.use_sort_by_non_overlap) {
			// Sort particles using Morton codes (Z-order curve) for spatial coherence
			// This minimizes overlapping voxel regions, significantly reducing atomic operation overhead
			// Can reduce GPU processing time from 20s to ~0.005s by eliminating atomic contention
			
			if (convert_vdb_base->cache_manager.use_gpu_cuda) {
				convert_vdb_base->cache_manager.sort_particles_by_nonoverlap_gpu_inplace();
			}
			else {
				convert_vdb_base->cache_manager.sort_particles_by_nonoverlap_cpu();
			}
		}

		LOG_MeasureStop("calculate_radius");

		return convert_vdb_base;
	}

	/**
	 * @brief Clean up and destroy the converter instance
	 * @param convert_vdb_base Pointer to converter to be destroyed
	 */
	void deinit_converter(common::vdb::ConvertVDBBase* convert_vdb_base)
	{
		convert_vdb_base->finish_lib();
		delete convert_vdb_base;
	}

#ifdef WITH_OPENVDB
	/**
	 * @brief Calculate minimum and maximum values in an OpenVDB grid
	 * @param grid OpenVDB float grid to analyze
	 * @param min_val Output parameter for minimum value
	 * @param max_val Output parameter for maximum value
	 */
	void eval_min_max(openvdb::FloatGrid::Ptr grid, float &min_val, float &max_val)
	{
		openvdb::math::MinMax<float> extrema = openvdb::tools::minMax(grid->tree());
		min_val = extrema.min();
		max_val = extrema.max();
	}
#endif

	/**
	 * @brief Print information about particle types and blocks across all MPI ranks
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param types_and_blocks_global Output vector for aggregated particle type/block counts
	 * @note If info flag is set, program will exit after printing
	 */
	void print_info(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, std::vector<int>& types_and_blocks_global)
	{
		std::vector<int> types_and_blocks_local;
		convert_vdb_base->get_types_and_blocks(types_and_blocks_local);

		types_and_blocks_global.resize(types_and_blocks_local.size());

		// Aggregate counts from all MPI ranks
		MPI_Allreduce(types_and_blocks_local.data(), types_and_blocks_global.data(), types_and_blocks_local.size(), MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if (from_cl.world_rank == 0) {
			convert_vdb_base->print_types_and_blocks(types_and_blocks_global);
		}

		if (from_cl.info) {
			convert_vdb_base->finish_lib();

			MPI_Finalize();
			exit(0);
		}

	}

	/**
	 * @brief Find bounding box of particles across all MPI ranks
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration (updated with bbox info)
	 * @param particle_type Type of particles to process
	 * 
	 * @details Computes local bounding box on each rank, then uses MPI reduction
	 *          to find global min/max. The bbox is expanded by 1 unit in each
	 *          direction and made symmetric around the center point.
	 */
	void find_bbox(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, int particle_type)
	{
		//CALL_MPI_BARRIER;
		LOG_MeasureStart("find_bbox");

		// double t_message_type2 = omp_get_wtime();
		// double t_find_bbox = omp_get_wtime();

		float bbox_min_orig_local[3], bbox_max_orig_local[3], bbox_min_orig[3], bbox_max_orig[3];

		// Find local bounding box for this rank's particles
		convert_vdb_base->iolib_find_bbox(
			particle_type,
			bbox_min_orig_local,
			bbox_max_orig_local,
			space_data.offset_position
		);

		//CALL_MPI_BARRIER;

		// if (from_cl.world_rank == 0)
		// 	printf("rank #%d: find bbox local: %f\n", from_cl.world_rank, omp_get_wtime() - t_message_type2);

		// t_find_bbox = omp_get_wtime();

		// Find global min/max across all MPI ranks
		MPI_Allreduce(bbox_min_orig_local, bbox_min_orig, 3, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(bbox_max_orig_local, bbox_max_orig, 3, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

		// Expand bounding box by 1 unit in each direction
		space_data.bbox_min_orig[0] = static_cast<int>(bbox_min_orig[0] - 1.0f);
		space_data.bbox_min_orig[1] = static_cast<int>(bbox_min_orig[1] - 1.0f);
		space_data.bbox_min_orig[2] = static_cast<int>(bbox_min_orig[2] - 1.0f);

		space_data.bbox_max_orig[0] = static_cast<int>(bbox_max_orig[0] + 1.0f);
		space_data.bbox_max_orig[1] = static_cast<int>(bbox_max_orig[1] + 1.0f);
		space_data.bbox_max_orig[2] = static_cast<int>(bbox_max_orig[2] + 1.0f);

		// Calculate the largest dimension to create a symmetric bounding box
		space_data.bbox_size_orig = std::max((double)space_data.bbox_max_orig[0] - (double)space_data.bbox_min_orig[0], (double)space_data.bbox_max_orig[1] - (double)space_data.bbox_min_orig[1]);
		space_data.bbox_size_orig = std::max(space_data.bbox_size_orig, (double)space_data.bbox_max_orig[2] - (double)space_data.bbox_min_orig[2]);

		// Make bounding box symmetric around center point
		space_data.bbox_min_orig[0] = (space_data.bbox_min_orig[0] + space_data.bbox_max_orig[0]) / 2.0 - space_data.bbox_size_orig / 2.0;
		space_data.bbox_min_orig[1] = (space_data.bbox_min_orig[1] + space_data.bbox_max_orig[1]) / 2.0 - space_data.bbox_size_orig / 2.0;
		space_data.bbox_min_orig[2] = (space_data.bbox_min_orig[2] + space_data.bbox_max_orig[2]) / 2.0 - space_data.bbox_size_orig / 2.0;

		space_data.bbox_max_orig[0] = (space_data.bbox_min_orig[0] + space_data.bbox_max_orig[0]) / 2.0 + space_data.bbox_size_orig / 2.0;
		space_data.bbox_max_orig[1] = (space_data.bbox_min_orig[1] + space_data.bbox_max_orig[1]) / 2.0 + space_data.bbox_size_orig / 2.0;
		space_data.bbox_max_orig[2] = (space_data.bbox_min_orig[2] + space_data.bbox_max_orig[2]) / 2.0 + space_data.bbox_size_orig / 2.0;

		//CALL_MPI_BARRIER;

		if (from_cl.world_rank == 0)
			printf("rank #%d: find bbox mpi - box_size: %f\n", from_cl.world_rank, space_data.bbox_size_orig);

		LOG_MeasureStop("find_bbox");
	}

	/**
	 * @brief Create an empty VDB grid structure based on extraction type
	 * @param grid_main Grid structure to be initialized
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * 
	 * @details Creates different grid types based on extraction mode:
	 *   - Dense: Regular 3D grid for volume data
	 *   - Particle: Raw particle data storage
	 *   - NanoVDB: Compact sparse grid format
	 *   - OpenVDB: Full-featured sparse grid format
	 */
	void create_grid(common::vdb::ConvertVDBBase* convert_vdb_base, common::vdb::VDBParticles& grid_main, space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		LOG_MeasureStart("create_grid");

		//Dense
		if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
			grid_main.type = common::vdb::VDBParticleType::eDense;

#ifdef WITH_GPU_CUDA
			if (from_cl.use_gpu_cuda) 
			{
				grid_main.dense_grid = std::make_shared<common::vdb::dense::VoxelGPUDenseManager>();
			}
			else 
#endif			
			{
				grid_main.dense_grid = std::make_shared<common::vdb::dense::VoxelCPUDenseManager>();
			}

			grid_main.dense_grid->create(space_data.bbox_dim, space_data.bbox_dim, space_data.bbox_dim);
		}
		//Raw particle data
		else if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
			grid_main.type = common::vdb::VDBParticleType::eRawParticles;
		}
		// NanoVDB (compact sparse format)
		// else if (from_cl.use_nanovdb) {
		// 	grid_main.type = common::vdb::VDBParticleType::eNanoVDB;		
// 			grid_main.nano_grid = std::make_shared<nanovdb::tools::build::FloatGrid>(0.0f, "density", nanovdb::GridClass::FogVolume);
		// }
		// Sparse
		else {
// #ifdef WITH_OPENVDB
			// grid_main.type = common::vdb::VDBParticleType::eOpenVDB;
// 			grid_main.vdb_grid = openvdb::FloatGrid::create(0.0f);
// 			grid_main.vdb_grid->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);
// 			grid_main.vdb_grid->setName("density");
// #endif			
			if (from_cl.use_nanovdb) {
				grid_main.type = common::vdb::VDBParticleType::eNanoVDB;				
			} else {
				grid_main.type = common::vdb::VDBParticleType::eOpenVDB;
			}

#ifdef WITH_GPU_CUDA			
			if (from_cl.use_gpu_cuda) {
				grid_main.sparse_grid = std::make_shared<common::vdb::sparse::VoxelGPUManagerSortReduce>();
				grid_main.sparse_grid->init(convert_vdb_base->get_local_num_particles()); // TODO: convert_vdb_base->get_local_num_particles()
			} 
			else 
#endif			
			{
				//grid_main.sparse_grid = std::make_shared<common::vdb::sparse::VoxelOpenMPManager>();
				grid_main.sparse_grid = std::make_shared<common::vdb::sparse::VoxelNanoVDBManager>();				
				grid_main.sparse_grid->init(convert_vdb_base->get_local_num_particles()); // TODO: convert_vdb_base->get_local_num_particles()
			}
		}

		LOG_MeasureStop("create_grid");
	}

	/**
	 * @brief Convert particle data to VDB grid format
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main Grid to populate with converted data
	 * 
	 * @details Reads particle data from the converter's I/O library and
	 *          rasterizes it into the VDB grid based on the specified
	 *          density estimation method and filtering parameters.
	 */
	void convert_to_grid(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main)
	{
		LOG_MeasureStart("convert_to_grid");

		//CALL_MPI_BARRIER;

		// double t_convert = omp_get_wtime();

		// Convert particles to grid using specified parameters
		convert_vdb_base->convert_iolib_to_grid(
			space_data.particle_type,
			space_data.particle_radius_multiplier,
			"density",
			space_data.grid_transform,
			space_data.bbox_min,
			space_data.bbox_max,
			space_data.bbox_dim,
			space_data.bbox_min_orig,
			space_data.bbox_size_orig,
			space_data.extracted_type,
			space_data.extracted_particle_type,
			space_data.dense_type,
			space_data.dense_norm,
			space_data.block_name_id,
			space_data.object_size,
			space_data.min_value_local,
			space_data.max_value_local,
			space_data.min_value,
			space_data.max_value,
			space_data.particles_count_local,
			grid_main,
			space_data.transform_scale,
			space_data.filter_min,
			space_data.filter_max,
			space_data.min_rho,
			space_data.max_rho,
			space_data.anim_type,
			space_data.frame,
			space_data.anim_start + from_cl.world_rank,
			space_data.bbox_sphere_pos,
			space_data.bbox_sphere_r,
			space_data.use_simple_density,
			space_data.use_norm_value,
			space_data.offset_position
		);

		//CALL_MPI_BARRIER;

		if (from_cl.world_rank == 0) {
			// printf("rank #%d: grid convert time: %f\n", from_cl.world_rank, omp_get_wtime() - t_convert);

			// Calculate bounding sphere center in original coordinate space
			float bbox_x_mid_orig = ((((space_data.bbox_min[0] / (float)space_data.object_size) + 
				(space_data.bbox_max[0] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[0] - (double)space_data.bbox_min_orig[0])) + (double)space_data.bbox_min_orig[0];

			float bbox_y_mid_orig = ((((space_data.bbox_min[1] / (float)space_data.object_size) +
				(space_data.bbox_max[1] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[1] - (double)space_data.bbox_min_orig[1])) + (double)space_data.bbox_min_orig[1];

			float bbox_z_mid_orig = ((((space_data.bbox_min[2] / (float)space_data.object_size) +
				(space_data.bbox_max[2] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[2] - (double)space_data.bbox_min_orig[2])) + (double)space_data.bbox_min_orig[2];

			float bbox_x_r_orig = ((((-space_data.bbox_min[0] / (float)space_data.object_size) +
				(space_data.bbox_max[0] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[0] - (double)space_data.bbox_min_orig[0]));// +(double)space_data.bbox_min_orig[0];

			float bbox_y_r_orig = ((((-space_data.bbox_min[1] / (float)space_data.object_size) +
				(space_data.bbox_max[1] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[1] - (double)space_data.bbox_min_orig[1]));// + (double)space_data.bbox_min_orig[1];

			float bbox_z_r_orig = ((((-space_data.bbox_min[2] / (float)space_data.object_size) +
				(space_data.bbox_max[2] / (float)space_data.object_size)) / 2.0) * ((double)space_data.bbox_max_orig[2] - (double)space_data.bbox_min_orig[2]));// + (double)space_data.bbox_min_orig[2];

			printf("rank #%d: bbox-sphere coord: %f, %f, %f, bbox-sphere radius: %f\n", 
				from_cl.world_rank, bbox_x_mid_orig, bbox_y_mid_orig, bbox_z_mid_orig,
				std::max(bbox_x_r_orig, std::max(bbox_y_r_orig, bbox_z_r_orig)));

			printf("rank #%d: bbox coord: %f, %f, %f, %f, %f, %f\n",
				from_cl.world_rank, space_data.bbox_min[0], space_data.bbox_min[1], space_data.bbox_min[2],
				space_data.bbox_max[0], space_data.bbox_max[1], space_data.bbox_max[2]);
		}

		LOG_MeasureStop("convert_to_grid");
	}

	/**
	 * @brief Find global minimum and maximum density values across all MPI ranks
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration (updated with global min/max)
	 * 
	 * @details Uses MPI_Allreduce to aggregate local min/max values from all
	 *          ranks and also sums up the total particle count.
	 */
	void find_minmax_value(space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{
		LOG_MeasureStart("find_minmax_value");
		// double t_find_minmax = omp_get_wtime();

		// Find global min/max across all MPI ranks
		MPI_Allreduce(&space_data.min_value_local, &space_data.min_value, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
		MPI_Allreduce(&space_data.max_value_local, &space_data.max_value, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

		// Sum up total particle count from all ranks
		MPI_Allreduce(&space_data.particles_count_local, &space_data.particles_count, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

		if (from_cl.world_rank == 0)
			printf("rank #%d: find minmax mpi - particles_count: %lld\n", from_cl.world_rank, space_data.particles_count);

		LOG_MeasureStop("find_minmax_value");
	}

	/**
	 * @brief Find min/max values after grid reduction (for animation frames)
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * 
	 * @details Only applicable when using animation mode (not None or AllMerge)
	 */
	void find_minmax_reduced_value(space_converter::FromCL& from_cl, common::SpaceData& space_data)
	{	
		if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
			LOG_MeasureStart("find_minmax_reduced_value");
			// double t_find_minmax = omp_get_wtime();

			float min_value_reduced = space_data.min_value_reduced;
			float max_value_reduced = space_data.max_value_reduced;

			// Find global min/max for reduced grids
			MPI_Allreduce(&min_value_reduced, &space_data.min_value_reduced, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
			MPI_Allreduce(&max_value_reduced, &space_data.max_value_reduced, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

			// if (from_cl.world_rank == 0)
			// 	printf("rank #%d: find reduced minmax mpi: %f\n", from_cl.world_rank, omp_get_wtime() - t_find_minmax);
			LOG_MeasureStop("find_minmax_reduced_value");
		}		
	}

	/**
	 * @brief Reduce/merge VDB grids from all MPI ranks into a single grid
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main Local grid from this rank
	 * @param grid_main_sum Output merged grid (valid on rank 0 or all ranks for animation)
	 * 
	 * @details For dense grids, uses MPI_Reduce to sum values. For sparse grids,
	 *          uses a logarithmic tree reduction pattern to merge grids efficiently.
	 *          In animation mode, each rank keeps its own grid.
	 */
	void reduction(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, common::vdb::VDBParticles& grid_main_sum)
	{
		LOG_MeasureStart("reduction");
		//CALL_MPI_BARRIER;

		// double t_grid = omp_get_wtime();

		// Dense grid reduction
		if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
			grid_main_sum.type = common::vdb::VDBParticleType::eDense;

#ifdef WITH_GPU_CUDA
			if (from_cl.use_gpu_cuda)
			{
				grid_main_sum.dense_grid = std::make_shared<common::vdb::dense::VoxelGPUDenseManager>();
			}
			else
#endif			
			{
				grid_main_sum.dense_grid = std::make_shared<common::vdb::dense::VoxelCPUDenseManager>();
			}

			if (from_cl.world_rank == 0 || space_data.anim_type != common::SpaceData::AnimType::eNone) {
				// Allocate data_temp only if normalization is needed
				bool allocate_data_temp = (space_data.dense_norm != common::SpaceData::DenseNorm::eNone);
				grid_main_sum.dense_grid->create(space_data.bbox_dim, space_data.bbox_dim, space_data.bbox_dim, allocate_data_temp);
				memcpy(grid_main_sum.dense_grid->offset, grid_main.dense_grid->offset, sizeof(grid_main.dense_grid->offset));
			}

#ifdef WITH_GPU_CUDA
			common::vdb::dense::VoxelGPUDenseManager* grid_main_gpu_sum = dynamic_cast<common::vdb::dense::VoxelGPUDenseManager*>(grid_main_sum.dense_grid.get());
			common::vdb::dense::VoxelGPUDenseManager* grid_main_gpu = dynamic_cast<common::vdb::dense::VoxelGPUDenseManager*>(grid_main.dense_grid.get());
#endif

			if (space_data.anim_type != common::SpaceData::AnimType::eNone) {
#ifdef WITH_GPU_CUDA
				if (grid_main_gpu && grid_main_gpu_sum) {
					// Animation mode: each rank keeps its own grid (no reduction)
					CUDA_CHECK_ERROR(cudaMemcpy(grid_main_gpu_sum->d_data_density, grid_main_gpu->d_data_density, grid_main_gpu->memsize(), cudaMemcpyDeviceToDevice));
					// Copy temp buffer if both grids have it allocated
					if (grid_main_gpu->d_data_temp != nullptr && grid_main_gpu_sum->d_data_temp != nullptr) {
						CUDA_CHECK_ERROR(cudaMemcpy(grid_main_gpu_sum->d_data_temp, grid_main_gpu->d_data_temp, grid_main_gpu->memsize(), cudaMemcpyDeviceToDevice));
					}

				}
				else 
#endif
				{
					// Animation mode: each rank keeps its own grid (no reduction)
					memcpy(grid_main_sum.dense_grid->data_density.data(), grid_main.dense_grid->data_density.data(), grid_main.dense_grid->memsize());
					// Copy temp buffer if both grids have it allocated
					if (!grid_main.dense_grid->data_temp.empty() && !grid_main_sum.dense_grid->data_temp.empty()) {
						memcpy(grid_main_sum.dense_grid->data_temp.data(), grid_main.dense_grid->data_temp.data(), grid_main.dense_grid->memsize());
					}
				}
			}
			else {
#if defined(WITH_GPU_CUDA) && defined(WITH_CUDA_AWARE_MPI)
				if (grid_main_gpu && grid_main_gpu_sum) {
					// Standard mode: sum all grids to rank 0
					mpi_reduce(grid_main_gpu->d_data_density, grid_main_gpu_sum->d_data_density, grid_main.dense_grid->size());
					
					// Debug logging: print first 10 values from both arrays
					DEBUG_PRINT_GPU_ARRAY(grid_main_gpu->d_data_density, grid_main.dense_grid->size(), "grid_main_gpu->d_data_density after mpi_reduce");
					DEBUG_PRINT_GPU_ARRAY(grid_main_gpu_sum->d_data_density, grid_main.dense_grid->size(), "grid_main_gpu_sum->d_data_density after mpi_reduce");
					
					// Reduce temp buffer if both grids have it allocated
					if (grid_main_gpu->d_data_temp != nullptr && grid_main_gpu_sum->d_data_temp != nullptr) {
						mpi_reduce(grid_main_gpu->d_data_temp, grid_main_gpu_sum->d_data_temp, grid_main.dense_grid->size());
					}
				}
				else
#endif
				{
#if defined(WITH_GPU_CUDA)
					if (grid_main_gpu && grid_main_gpu_sum) {
						grid_main_gpu->from_device();
						grid_main_gpu_sum->from_device();
					}
#endif
					// Standard mode: sum all grids to rank 0
					mpi_reduce(grid_main.dense_grid->data_density.data(), grid_main_sum.dense_grid->data_density.data(), grid_main.dense_grid->size());
					// Reduce temp buffer if both grids have it allocated
					if (!grid_main.dense_grid->data_temp.empty() && !grid_main_sum.dense_grid->data_temp.empty()) {
						mpi_reduce(grid_main.dense_grid->data_temp.data(), grid_main_sum.dense_grid->data_temp.data(), grid_main.dense_grid->size());
					}		

#if defined(WITH_GPU_CUDA)
					if (grid_main_gpu_sum) {
						grid_main_gpu_sum->to_device();
					}
#endif
				}

			}

			grid_main.dense_grid->clear();
		}
		// Sparse grid reduction (VDB/NanoVDB/Particles)
		else {
			grid_main_sum.type = grid_main.type;

			if (grid_main_sum.type == common::vdb::VDBParticleType::eOpenVDB ||
				grid_main_sum.type == common::vdb::VDBParticleType::eNanoVDB) {
				grid_main_sum.sparse_grid = grid_main.sparse_grid;
// #ifdef WITH_OPENVDB
// 				grid_main_sum.vdb_grid = grid_main.vdb_grid;
// #endif
 			}
// 			else if (grid_main_sum.type == common::vdb::VDBParticleType::eNanoVDB)
// 				grid_main_sum.nano_grid = grid_main.nano_grid;
			else if (grid_main_sum.type == common::vdb::VDBParticleType::eRawParticles) {
				grid_main_sum.raw_particles = grid_main.raw_particles;
			}

			// Perform logarithmic tree reduction for sparse grids
			if (space_data.anim_type == common::SpaceData::AnimType::eNone || space_data.anim_type == common::SpaceData::AnimType::eAllMerge || space_data.anim_type == common::SpaceData::AnimType::eFrameExtract) {
				// Logarithmic reduction: pair-wise merging in a tree structure
				for (int step = 1; step < from_cl.world_size; step *= 2) {
					if (from_cl.world_rank % (2 * step) == 0) {
						// Receiving rank: merge data from partner
						if (from_cl.world_rank + step < from_cl.world_size) {

							size_t ns = 0;
							mpi_recv(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank + step, 0);
#if defined(WITH_GPU_CUDA) && defined(WITH_CUDA_AWARE_MPI)
							if (auto* gpu_mgr = dynamic_cast<common::vdb::sparse::VoxelGPUManagerSortReduce*>(grid_main_sum.sparse_grid.get())) {
								// RDMA path: receive directly into device buffer, merge on GPU
								uint8_t* d_buf = nullptr;
								cudaMalloc(&d_buf, ns);
								mpi_recv(d_buf, ns, MPI_BYTE, from_cl.world_rank + step, 0);
								common::vdb::sparse::VoxelGPUManagerSortReduce recv_mgr;
								recv_mgr.deserializeGPU(d_buf);
								cudaFree(d_buf);
								gpu_mgr->merge(&recv_mgr);
							}
							else
#endif
							{
								common::vdb::VDBParticles nanogrid_recv;
								nanogrid_recv.type = common::vdb::VDBParticleType::eVector;
								nanogrid_recv.vector_grid.resize(ns);
								mpi_recv(nanogrid_recv.vector_grid.data(), ns, MPI_BYTE, from_cl.world_rank + step, 0);
								convert_vdb_base->merge_grid(grid_main_sum, nanogrid_recv);
							}
						}
					}
					else if (from_cl.world_rank % (2 * step) == step) {
						// Sending rank: serialize and send grid to partner
						if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
							std::vector<uint8_t> grid_handle_main;
							grid_main_sum.raw_particles.serialize(grid_handle_main);
							size_t ns = grid_handle_main.size();
							mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank - step, 0);

							// TODO: GPU
							mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank - step, 0);
						}
						else {
							// Sparse particle data
							size_t ns = grid_main_sum.sparse_grid->mem_size();
							mpi_send(&ns, sizeof(ns), MPI_BYTE, from_cl.world_rank - step, 0);
#if defined(WITH_GPU_CUDA) && defined(WITH_CUDA_AWARE_MPI)
							if (auto* gpu_mgr = dynamic_cast<common::vdb::sparse::VoxelGPUManagerSortReduce*>(grid_main_sum.sparse_grid.get())) {
								// RDMA path: serialize directly into device buffer, send via GPU-direct
								uint8_t* d_buf = nullptr;
								cudaMalloc(&d_buf, ns);
								gpu_mgr->serializeGPU(d_buf);
								mpi_send(d_buf, ns, MPI_BYTE, from_cl.world_rank - step, 0);
								cudaFree(d_buf);
							}
							else
#endif
							{
								std::vector<uint8_t> grid_handle_main(ns);
								grid_main_sum.sparse_grid->serialize(grid_handle_main.data());
								mpi_send(grid_handle_main.data(), ns, MPI_BYTE, from_cl.world_rank - step, 0);
							}
						}

						break; // This rank is done participating
					}
				}
			}
		}

		//CALL_MPI_BARRIER;

		// if (from_cl.world_rank == 0) {
		// 	printf("rank #%d: merged time: %f\n", from_cl.world_rank, omp_get_wtime() - t_grid);
		// }

		LOG_MeasureStop("reduction");
	}
	
	void finalize_grid(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main_sum, common::vdb::VDBParticles& grid_main_final)
	{
		LOG_MeasureStart("finalize_grid");

		if (from_cl.world_rank == 0 || space_data.anim_type != common::SpaceData::AnimType::eNone) {
			//Dense
			if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {

				// double t = omp_get_wtime();

				grid_main_final.type = common::vdb::VDBParticleType::eVector;

				// Convert dense grid to sparse NanoVDB format
				if (from_cl.use_nanovdb) {
					auto nanovdb_handleI = convert_vdb_base->dense_to_nanovdb(grid_main_sum.dense_grid.get(), space_data.transform_scale, space_data.dense_type, space_data.dense_norm);

					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*nanovdb_handleI);

					// Serialize NanoVDB to binary vector
					grid_main_final.vector_grid.resize(grid_handle_final.bufferSize());
					memcpy(grid_main_final.vector_grid.data(), grid_handle_final.data(), grid_handle_final.bufferSize());

					nanovdb::NanoGrid<float>* nanogrid = (nanovdb::NanoGrid<float>*) grid_handle_final.data();

					// Extract min/max from NanoVDB tree
					nanogrid->tree().extrema(space_data.min_value_reduced, space_data.max_value_reduced);
				}
				// Convert dense grid to sparse OpenVDB format
				else {

#ifdef WITH_OPENVDB
					auto openvdb_handleI = convert_vdb_base->dense_to_openvdb(grid_main_sum.dense_grid.get(), space_data.transform_scale, space_data.dense_type, space_data.dense_norm);
					
					convert_vdb_base->openvdb_to_vector(openvdb_handleI, grid_main_final.vector_grid);
					printf("rank #%d: Active voxels in input grid: %lld\n", from_cl.world_rank, openvdb_handleI->activeVoxelCount());

					eval_min_max(openvdb_handleI, space_data.min_value_reduced, space_data.max_value_reduced);
#endif
				}

				printf("rank #%d: minI: %e, maxI: %e, reduced: minI: %e, maxI: %e\n", from_cl.world_rank, space_data.min_value, space_data.max_value, space_data.min_value_reduced, space_data.max_value_reduced);
				//printf("rank #%d: final grid time (Dense): %f\n", from_cl.world_rank, omp_get_wtime() - t);

				// Optionally save raw volume data
				if (space_data.extracted_dense_type == common::SpaceData::ExtractedDenseType::eRAW)
					save_raw_volume(convert_vdb_base, from_cl, space_data, grid_main_sum);
				else if (space_data.extracted_dense_type == common::SpaceData::ExtractedDenseType::eVTI)
					save_vti_volume(convert_vdb_base, from_cl, space_data, grid_main_sum);
			}
			// Sparse grid finalization
			else {
				// double t = omp_get_wtime();

				grid_main_final.type = common::vdb::VDBParticleType::eVector;

				// Grid already in vector format (from reduction)
				if (grid_main_sum.type == common::vdb::VDBParticleType::eVector) {
					grid_main_final.vector_grid = std::move(grid_main_sum.vector_grid);
				}
				// Raw particle data
				else if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
					grid_main_sum.raw_particles.serialize(grid_main_final.vector_grid);

					space_data.min_value_reduced = space_data.min_value;
					space_data.max_value_reduced = space_data.max_value;

					// Optionally save as VDB point cloud
					if (space_data.extracted_particle_type == common::SpaceData::ExtractedParticleType::eVDB)
						save_raw_particles_to_vdb(convert_vdb_base, from_cl, space_data, grid_main_sum);
					else if (space_data.extracted_particle_type == common::SpaceData::ExtractedParticleType::eVTP)
						save_raw_particles_to_vtp(convert_vdb_base, from_cl, space_data, grid_main_sum);

					printf("rank #%d: Particles count: %lld\n", from_cl.world_rank, (size_t)(grid_main_sum.raw_particles.data[0].values.size() / 3));
				}
				// Convert NanoVDB to binary format
				else if (from_cl.use_nanovdb) {
//					convert_vdb_base->sparse_to_nanovdb(grid_main_sum.sparse_grid.get());
//
//					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*grid_main_sum.nano_grid);
//
//					// Serialize to binary
//					grid_main_final.vector_grid.resize(grid_handle_final.size());
//					memcpy(grid_main_final.vector_grid.data(), grid_handle_final.data(), grid_handle_final.size());
//
//					nanovdb::NanoGrid<float>* nanogrid = (nanovdb::NanoGrid<float>*) grid_handle_final.data();
//
//					nanogrid->tree().extrema(space_data.min_value_reduced, space_data.max_value_reduced);

					auto nanovdb_handleI = convert_vdb_base->sparse_to_nanovdb(grid_main_sum.sparse_grid.get());

					nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*nanovdb_handleI);

					// Serialize NanoVDB to binary vector
					grid_main_final.vector_grid.resize(grid_handle_final.bufferSize());
					memcpy(grid_main_final.vector_grid.data(), grid_handle_final.data(), grid_handle_final.bufferSize());

					//if (grid_main_sum.nano_grid) {
					//	nanovdb::NanoGrid<float>* nanogrid = (nanovdb::NanoGrid<float>*) grid_handle_final.data();

					// Extract min/max from NanoVDB tree

					//nanogrid->tree().extrema(space_data.min_value_reduced, space_data.max_value_reduced);
				}
				// Convert OpenVDB to binary format
				else {

//#ifdef WITH_OPENVDB
//					convert_vdb_base->openvdb_to_vector(grid_main_sum.vdb_grid, grid_main_final.vector_grid);
//					eval_min_max(grid_main_sum.vdb_grid, space_data.min_value_reduced, space_data.max_value_reduced);
//					std::cout << "Active voxels in input grid: " << grid_main_sum.vdb_grid->activeVoxelCount() << std::endl;
//#endif

#ifdef WITH_OPENVDB
					auto openvdb_handleI = convert_vdb_base->sparse_to_openvdb(grid_main_sum.sparse_grid.get());

					convert_vdb_base->openvdb_to_vector(openvdb_handleI, grid_main_final.vector_grid);
					printf("rank #%d: Active voxels in input grid: %lld\n", from_cl.world_rank, openvdb_handleI->activeVoxelCount());

					eval_min_max(openvdb_handleI, space_data.min_value_reduced, space_data.max_value_reduced);
#endif
				}

				// CALL_MPI_BARRIER;

				printf("rank #%d: minI: %e, maxI: %e, reduced: minI: %e, maxI: %e\n", from_cl.world_rank, space_data.min_value, space_data.max_value, space_data.min_value_reduced, space_data.max_value_reduced);
				printf("rank #%d: grid_handle_merged_size: %lld\n", from_cl.world_rank, grid_main_final.vector_grid.size());
				// printf("rank #%d: final grid time: %f\n", from_cl.world_rank, omp_get_wtime() - t);
			}

			//printf("rank #%d: normalize_values: %f, voxels_count: %lld, voxels_count_zero: %lld\n", from_cl.world_rank, omp_get_wtime() - t, voxels_count, voxels_count_zero);

			if (space_data.extracted_type == common::SpaceData::ExtractedType::eDense) {
				grid_main_sum.dense_grid->clear();
			}
		}
		LOG_MeasureStop("finalize_grid");
	}


	/**
	 * @brief Save VDB grid to file
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main_final Grid to save
	 * @param particle_type Type of particle data in grid
	 * @param only_rank0 If true, only rank 0 saves (default behavior for merged grids)
	 * 
	 * @details Saves in multiple formats based on grid type:
	 *   - .vdb: OpenVDB format
	 *   - .nvdb: NanoVDB format
	 *   - .bin: Binary particle data
	 *   - .part: Raw particle format
	 *   - .raw: Raw volume data
	 */
	void save_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main_final, common::vdb::VDBParticleType particle_type, bool only_rank0)
	{
		LOG_MeasureStart("save_vdb");

		if (!only_rank0 || from_cl.world_rank == 0 || (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) &&
			(space_data.anim_type == common::SpaceData::AnimType::eFrameExtract && space_data.frame == space_data.anim_start + from_cl.world_rank || space_data.anim_type == common::SpaceData::AnimType::eAllPath)) {
			
			// Build output filename
			std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);

			// Add frame/rank suffix for animation sequences
			if (!only_rank0 || space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
				char temp[1024];
				sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
				full_filepath = full_filepath + "_" + std::string(temp);
			}

			// Determine file extension based on data type

			if (particle_type == common::vdb::VDBParticleType::eRawParticles) {
				full_filepath = full_filepath + std::string(".part");
			}
			else if (particle_type == common::vdb::VDBParticleType::eVector) {
				if (space_data.extracted_type == common::SpaceData::ExtractedType::eParticle) {
					full_filepath = full_filepath + std::string(".bin");
				}
				else if (from_cl.use_nanovdb) {
					full_filepath = full_filepath + std::string(".nvdb");
				}
				else {
					full_filepath = full_filepath + std::string(".vdb");
				}				
			}
			else if (particle_type == common::vdb::VDBParticleType::eNanoVDB) {
				full_filepath = full_filepath + std::string(".nvdb");
			}
			else {
				full_filepath = full_filepath + std::string(".vdb");
			}
			space_data.full_filepath = full_filepath;

			// Save based on particle type
			if (particle_type == common::vdb::VDBParticleType::eOpenVDB) {
#ifdef WITH_OPENVDB
				//openvdb::io::File(full_filepath).write({ grid_main_final.vdb_grid });
				auto openvdb_handleI = convert_vdb_base->sparse_to_openvdb(grid_main_final.sparse_grid.get());
				openvdb::io::File(full_filepath).write({ openvdb_handleI });
#endif

				printf("rank #%d: finished: %s\n", from_cl.world_rank, full_filepath.c_str());
			}
			// Save serialized VDB data (NanoVDB or binary particle data)
			else if (particle_type == common::vdb::VDBParticleType::eVector) {

				// Open a file in binary mode
				std::ofstream output_file(full_filepath, std::ios::binary);
				if (!output_file) {
					printf("Unable to open file for writing: %s\n", full_filepath.c_str());
					return;
				}

				// Write the vector data to file
				output_file.write((char*)grid_main_final.vector_grid.data(), grid_main_final.vector_grid.size());

				// Close the file
				output_file.close();

				printf("rank #%d: finished: %s\n", from_cl.world_rank, full_filepath.c_str());

			}
			else if (particle_type == common::vdb::VDBParticleType::eNanoVDB) {
//				nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*grid_main_final.nano_grid);

				auto nanovdb_handleI = convert_vdb_base->sparse_to_nanovdb(grid_main_final.sparse_grid.get());

				nanovdb::GridHandle<nanovdb::HostBuffer> grid_handle_final = nanovdb::tools::createNanoGrid(*nanovdb_handleI);

//				if (grid_main_final.nano_grid) {
//#if OPENVDB_VERSION == 11
//					grid_handle_final = nanovdb::createNanoGrid(*grid_main_final.nano_grid);
//#else
//					grid_handle_final = nanovdb::tools::createNanoGrid(*grid_main_final.nano_grid);
//#endif				
//				}
//				else if (grid_main_final.nano_grid3) {
//#if OPENVDB_VERSION == 11
//					grid_handle_final = nanovdb::createNanoGrid(*grid_main_final.nano_grid3);
//#else
//					grid_handle_final = nanovdb::tools::createNanoGrid(*grid_main_final.nano_grid3);
//#endif
//				}

				nanovdb::io::writeGrid(space_data.full_filepath, grid_handle_final);
				printf("rank #%d: finished: %s\n", from_cl.world_rank, full_filepath.c_str());
			}
			// Save raw particle format
			else if (particle_type == common::vdb::VDBParticleType::eRawParticles) {
				std::vector<uint8_t> grid_handle_main;
				grid_main_final.raw_particles.serialize(grid_handle_main);
				
				// Save serialized particle data
				std::ofstream output_file(full_filepath, std::ios::binary);
				if (!output_file) {
					printf("Unable to open file for writing: %s\n", full_filepath.c_str());
					return;
				}

				// Write the content of the vector to the file
				output_file.write((char*)grid_handle_main.data(), grid_handle_main.size());

				// Close the file
				output_file.close();

				printf("rank #%d: finished: %s\n", from_cl.world_rank, full_filepath.c_str());
			}
			// Save dense volume format
			else if (particle_type == common::vdb::VDBParticleType::eDense) {
				if (space_data.extracted_dense_type == common::SpaceData::ExtractedDenseType::eRAW)
					save_raw_volume(convert_vdb_base, from_cl, space_data, grid_main_final, only_rank0);
				else if (space_data.extracted_dense_type == common::SpaceData::ExtractedDenseType::eVTI)
					save_vti_volume(convert_vdb_base, from_cl, space_data, grid_main_final, only_rank0);				
			}
			else {
				printf("Unknown Type for Saving\n");
			}
		}

		LOG_MeasureStop("save_vdb");
	}

	/**
	 * @brief Save raw volume data to binary file
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main Grid containing dense volume data
	 * @param only_rank0 If true, only rank 0 saves
	 * 
	 * @details Saves dense voxel data as raw binary float array.
	 *          Filename includes grid dimensions (e.g., _512_512_512_float.raw)
	 */
	void save_raw_volume(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, bool only_rank0)
	{
		LOG_MeasureStart("save_raw_volume");

		if (grid_main.dense_grid->data_density.size() > 0) {
			if (!only_rank0 || from_cl.world_rank == 0) {
				// Build output filename with dimensions
				std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);
				if (!only_rank0 || space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
					char temp[1024];
					sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
					full_filepath = full_filepath + "_" + std::string(temp);
				}

				full_filepath = full_filepath
					+ std::string("_") + std::to_string(grid_main.dense_grid->x())
					+ std::string("_") + std::to_string(grid_main.dense_grid->y())
					+ std::string("_") + std::to_string(grid_main.dense_grid->z())
					+ std::string("_float.raw");

				space_data.full_filepath = full_filepath;

				// Write raw float data to file
				std::ofstream out(full_filepath.c_str(), std::ios::binary);
				out.write((const char*)grid_main.dense_grid->data_density.data(), grid_main.dense_grid->data_density.size() * sizeof(grid_main.dense_grid->data_density[0]));

				printf("rank #%d: finished: %s\n", from_cl.world_rank, full_filepath.c_str());
			}
		}

		LOG_MeasureStop("save_raw_volume");
	}

	/**
	 * @brief Save dense volume data to VTK ImageData (.vti) file
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main Grid containing dense volume data
	 * @param only_rank0 If true, only rank 0 saves
	 */
	void save_vti_volume(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, bool only_rank0)
	{
		LOG_MeasureStart("save_vti_volume");

		if (grid_main.dense_grid->data_density.size() > 0) {
			if (!only_rank0 || from_cl.world_rank == 0) {
				std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);
				if (!only_rank0 || space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
					char temp[1024];
					sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
					full_filepath = full_filepath + "_" + std::string(temp);
				}

				full_filepath = full_filepath
					+ std::string("_") + std::to_string(grid_main.dense_grid->x())
					+ std::string("_") + std::to_string(grid_main.dense_grid->y())
					+ std::string("_") + std::to_string(grid_main.dense_grid->z())
					+ std::string("_float.vti");

				space_data.full_filepath = full_filepath;

#ifdef WITH_VTK
				const size_t dim_x = grid_main.dense_grid->x();
				const size_t dim_y = grid_main.dense_grid->y();
				const size_t dim_z = grid_main.dense_grid->z();
				const vtkIdType num_points = static_cast<vtkIdType>(grid_main.dense_grid->data_density.size());

				auto image_data = vtkSmartPointer<vtkImageData>::New();
				image_data->SetDimensions(static_cast<int>(dim_x), static_cast<int>(dim_y), static_cast<int>(dim_z));
				const double spacing = static_cast<double>(space_data.transform_scale);
				image_data->SetOrigin(
					static_cast<double>(grid_main.dense_grid->offset[0]) * spacing,
					static_cast<double>(grid_main.dense_grid->offset[1]) * spacing,
					static_cast<double>(grid_main.dense_grid->offset[2]) * spacing);
				image_data->SetSpacing(spacing, spacing, spacing);

				auto density_array = vtkSmartPointer<vtkFloatArray>::New();
				density_array->SetName("density");
				density_array->SetNumberOfComponents(1);
				density_array->SetNumberOfTuples(num_points);

				for (vtkIdType i = 0; i < num_points; i++) {
					density_array->SetValue(i, grid_main.dense_grid->data_density[static_cast<size_t>(i)]);
				}

				image_data->GetPointData()->SetScalars(density_array);

				auto writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
				writer->SetFileName(full_filepath.c_str());
				writer->SetInputData(image_data);
				writer->SetDataModeToBinary();

				if (!writer->Write()) {
					printf("Unable to write file: %s\n", full_filepath.c_str());
					return;
				}

				printf("rank #%d: finished: %s\n", from_cl.world_rank, full_filepath.c_str());
#else
				printf("VTK support is not enabled, cannot write file: %s\n", full_filepath.c_str());
#endif
			}
		}

		LOG_MeasureStop("save_vti_volume");
	}

	/**
	 * @brief Save raw particle data as OpenVDB point cloud
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main Grid containing raw particle data
	 * 
	 * @details Converts raw particle attributes to OpenVDB points format with
	 *          support for multiple attributes (position, density, velocity, etc.)
	 *          Only rank 0 performs the conversion and save.
	 */
	void save_raw_particles_to_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main)
	{
		LOG_MeasureStart("save_raw_particles_to_vdb");

#ifdef WITH_OPENVDB
		if (/*from_cl.use_rawpart2vdb &&*/ grid_main.raw_particles.data.size() > 0) {
			if (from_cl.world_rank == 0) {
				std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);

				if (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge) {
					char temp[1024];
					sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
					full_filepath = full_filepath + "_" + std::string(temp);
					//space_data.anim_task_counter++;
				}

				full_filepath = full_filepath + std::string(".part.vdb");
				space_data.full_filepath = full_filepath;

				// positions
				std::vector<openvdb::Vec3f> points;
				for (size_t type = 0; type < grid_main.raw_particles.data.size(); type++) {
					const std::string attrib_name = grid_main.raw_particles.data[type].name;
					if (attrib_name == "position") {
						const size_t num_points = grid_main.raw_particles.data[type].values.size() / 3;
						points.resize(num_points);
#pragma omp parallel for
						for (size_t i = 0; i < num_points; i++) {
							points[i] = openvdb::Vec3f(grid_main.raw_particles.data[type].values[i * 3 + 0], grid_main.raw_particles.data[type].values[i * 3 + 1], grid_main.raw_particles.data[type].values[i * 3 + 2]);
						}
						break;
					}
				}

				// Compute voxel size and create the transform				
				openvdb::points::PointAttributeVector positions_wrapper(points);
				float vSize = openvdb::points::computeVoxelSize(positions_wrapper, 8);
				openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(vSize);

				// Create a PointIndexGrid. This can be done automatically on creation of the grid, however as this index grid is
				// required for the position and radius attributes, we create one we can use for both attribute creations.
				openvdb::tools::PointIndexGrid::Ptr pt_index = openvdb::tools::createPointIndexGrid<openvdb::tools::PointIndexGrid>(positions_wrapper, *transform);

				// Create a PointDataGrid containing these four points and using the point index grid. This requires the positions
				// wrapper.
				openvdb::points::PointDataGrid::Ptr pt_grid = openvdb::points::createPointDataGrid<openvdb::points::NullCodec, openvdb::points::PointDataGrid>(
					*pt_index, positions_wrapper, *transform);

				pt_grid->setName("points");

				// Attribute definition
				const openvdb::NamePair vdb_var = openvdb::points::TypedAttributeArray<float>::attributeType();
				const openvdb::NamePair vdb_var3 = openvdb::points::TypedAttributeArray<openvdb::Vec3f>::attributeType();

				for (size_t type = 0; type < grid_main.raw_particles.data.size(); type++) {
					const std::string attrib_name = grid_main.raw_particles.data[type].name;
					if (attrib_name == "position") {
						continue;
					}

					if (grid_main.raw_particles.data[type].num_comp == 1) {
						// Add the density attribute to the grid
						openvdb::points::appendAttribute(pt_grid->tree(), attrib_name, vdb_var);
						openvdb::points::populateAttribute(pt_grid->tree(), pt_index->tree(), attrib_name, openvdb::points::PointAttributeVector(grid_main.raw_particles.data[type].values));
					}
					else if (grid_main.raw_particles.data[type].num_comp == 3) {
						const size_t num_points = grid_main.raw_particles.data[type].values.size() / 3;
						std::vector<openvdb::Vec3f> values_vector(num_points);

#pragma omp parallel for
						for (size_t i = 0; i < num_points; i++) {
							values_vector[i] = openvdb::Vec3f(grid_main.raw_particles.data[type].values[i * 3 + 0], grid_main.raw_particles.data[type].values[i * 3 + 1], grid_main.raw_particles.data[type].values[i * 3 + 2]);
						}

						// Add the density attribute to the grid
						openvdb::points::appendAttribute(pt_grid->tree(), attrib_name, vdb_var3);
						openvdb::points::populateAttribute(pt_grid->tree(), pt_index->tree(), attrib_name, openvdb::points::PointAttributeVector(values_vector));
					}
				}

				openvdb::io::File file(full_filepath);
				file.setCompression(openvdb::io::COMPRESS_BLOSC | openvdb::io::COMPRESS_ACTIVE_MASK);
				file.write({ pt_grid });
				file.close();

				printf("finished: %s\n", full_filepath.c_str());
			}
		}
#endif
		LOG_MeasureStop("save_raw_particles_to_vdb");
	}

	/**
	 * @brief Save raw particle data as VTK PolyData (.vtp)
	 * @param convert_vdb_base Converter instance
	 * @param from_cl Command line parameters
	 * @param space_data Spatial data configuration
	 * @param grid_main Grid containing raw particle data
	 *
	 * @details Converts raw particle attributes to VTK point data arrays.
	 *          Supports scalar (1 component) and vector (3 component) attributes.
	 *          Only rank 0 performs the conversion and save.
	 */
	void save_raw_particles_to_vtp(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& space_data, common::vdb::VDBParticles& grid_main, bool only_rank0)
	{
		LOG_MeasureStart("save_raw_particles_to_vtp");

#ifdef WITH_VTK
		if (/*from_cl.use_rawpart2vdb &&*/ grid_main.raw_particles.data.size() > 0) {
			if (!only_rank0 || from_cl.world_rank == 0) {
				std::string full_filepath = from_cl.output_path + "/" + convert_vdb_base->get_type_name(space_data.particle_type) + "_" + convert_vdb_base->get_dataset_name(space_data.block_name_id);

				if (!only_rank0 || (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge)) {
					char temp[1024];
					sprintf(temp, "%d_%05d", space_data.anim_task_counter, from_cl.world_rank);
					full_filepath = full_filepath + "_" + std::string(temp);
				}

				full_filepath = full_filepath + std::string(".part.vtp");
				space_data.full_filepath = full_filepath;

				const common::vdb::RawParticles::ParticleData* position_data = nullptr;
				for (size_t type = 0; type < grid_main.raw_particles.data.size(); type++) {
					if (grid_main.raw_particles.data[type].name == "position") {
						position_data = &grid_main.raw_particles.data[type];
						break;
					}
				}

				if (!position_data || position_data->values.size() < 3 || position_data->num_comp != 3) {
					printf("Unable to export VTP: missing valid 'position' attribute\n");
					return;
				}

				const size_t num_points = position_data->values.size() / 3;

				auto points = vtkSmartPointer<vtkPoints>::New();
				points->SetDataTypeToFloat();
				points->SetNumberOfPoints(num_points);

#pragma omp parallel for
				for (size_t i = 0; i < num_points; i++) {
					points->SetPoint(
						static_cast<vtkIdType>(i),
						position_data->values[i * 3 + 0],
						position_data->values[i * 3 + 1],
						position_data->values[i * 3 + 2]);
				}

				auto vertices = vtkSmartPointer<vtkCellArray>::New();
				for (size_t i = 0; i < num_points; i++) {
					vertices->InsertNextCell(1);
					vertices->InsertCellPoint(static_cast<vtkIdType>(i));
				}

				auto poly_data = vtkSmartPointer<vtkPolyData>::New();
				poly_data->SetPoints(points);
				poly_data->SetVerts(vertices);

				for (size_t type = 0; type < grid_main.raw_particles.data.size(); type++) {
					const std::string attrib_name = grid_main.raw_particles.data[type].name;
					if (attrib_name == "position") {
						continue;
					}

					if (grid_main.raw_particles.data[type].num_comp == 1) {
						auto values = vtkSmartPointer<vtkFloatArray>::New();
						values->SetName(attrib_name.c_str());
						values->SetNumberOfComponents(1);
						values->SetNumberOfTuples(num_points);

						const size_t max_points = std::min(num_points, grid_main.raw_particles.data[type].values.size());
						for (size_t i = 0; i < max_points; i++) {
							values->SetValue(static_cast<vtkIdType>(i), grid_main.raw_particles.data[type].values[i]);
						}

						poly_data->GetPointData()->AddArray(values);
					}
					else if (grid_main.raw_particles.data[type].num_comp == 3) {
						auto vectors = vtkSmartPointer<vtkFloatArray>::New();
						vectors->SetName(attrib_name.c_str());
						vectors->SetNumberOfComponents(3);
						vectors->SetNumberOfTuples(num_points);

						const size_t point_count = std::min(num_points, grid_main.raw_particles.data[type].values.size() / 3);
						for (size_t i = 0; i < point_count; i++) {
							vectors->SetTuple3(
								static_cast<vtkIdType>(i),
								grid_main.raw_particles.data[type].values[i * 3 + 0],
								grid_main.raw_particles.data[type].values[i * 3 + 1],
								grid_main.raw_particles.data[type].values[i * 3 + 2]);
						}

						poly_data->GetPointData()->AddArray(vectors);
					}
				}

				auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
				writer->SetFileName(full_filepath.c_str());
				writer->SetInputData(poly_data);
				writer->SetDataModeToBinary();

				if (!writer->Write()) {
					printf("Unable to write file: %s\n", full_filepath.c_str());
					return;
				}

				printf("finished: %s\n", full_filepath.c_str());

				if (!only_rank0 || (from_cl.world_rank == 0 && (space_data.anim_type != common::SpaceData::AnimType::eNone && space_data.anim_type != common::SpaceData::AnimType::eAllMerge))) {
					std::string pfull_filepath = full_filepath;
					const std::string part_ext = ".part.vtp";
					if (pfull_filepath.size() >= part_ext.size() && pfull_filepath.compare(pfull_filepath.size() - part_ext.size(), part_ext.size(), part_ext) == 0) {
						pfull_filepath.replace(pfull_filepath.size() - part_ext.size(), part_ext.size(), ".part.pvtp");
					}
					else {
						pfull_filepath += ".pvtp";
					}

					auto pwriter = vtkSmartPointer<vtkXMLPPolyDataWriter>::New();
					pwriter->SetFileName(pfull_filepath.c_str());
					pwriter->SetNumberOfPieces(from_cl.world_size);
					pwriter->SetStartPiece(0);
					pwriter->SetEndPiece(from_cl.world_size - 1);
					pwriter->SetInputData(poly_data);

					if (!pwriter->Write()) {
						printf("Unable to write file: %s\n", pfull_filepath.c_str());
						return;
					}

					printf("finished: %s\n", pfull_filepath.c_str());
				}
			}
		}
#endif
		LOG_MeasureStop("save_raw_particles_to_vtp");
	}
}