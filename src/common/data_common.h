#pragma once

#include <string>
#include <float.h>
#include <cstring>
#include <cstdint>
#include <vector>
#include <memory>
#include <iostream>

#include "raw_common.h"
// Note: dense_common.h includes data_common.h (not the other way around to avoid circular deps).
// Include dense_common.h explicitly where VoxelCPUDenseManager / VoxelGPUDenseManager are needed.

#define FDATA_EPSILON 1e-8f  // Small epsilon value for floating-point comparisons

/**
 * Namespace common provides functionality for handling spatial data conversions
 */
namespace common {

	namespace vdb {

		/**
		 * @brief Type of VDB particle representation
		 */
		enum VDBParticleType
		{
			eDense,        ///< Dense regular grid
			eVector,       ///< Serialized binary format
			eNanoVDB,      ///< NanoVDB sparse grid
			eOpenVDB,      ///< OpenVDB sparse grid
			eRawParticles  ///< Raw point cloud
		};

		/**
		 * @brief Base class for dense regular grid representation of particle data.
		 *
		 * Holds the host-side data buffers and grid metadata. Concrete subclasses
		 * (VoxelCPUDenseManager, VoxelGPUDenseManager) override clear() and create()
		 * for CPU-only and CPU+GPU paths respectively.
		 *
		 * Also provides backward-compatible virtual stubs for the sparse-grid interface
		 * (set_transform, init, update, serialize, …) so that existing sparse-grid
		 * implementations that inherit from this class continue to compile unchanged.
		 */
		class VoxelDenseManager {
		public:
			virtual ~VoxelDenseManager() = default;

			// ── Host data ──────────────────────────────────────────────────────────
			std::vector<float> data_density;   ///< Primary density/value data
			std::vector<float> data_temp;      ///< Temporary accumulation buffer for normalisation (only allocated when needed)
			size_t dims[3]   = { 0, 0, 0 };   ///< Grid dimensions [x, y, z]
			size_t offset[3] = { 0, 0, 0 };   ///< Grid offset in global coordinate space

			// ── Lifecycle ──────────────────────────────────────────────────────────

			/** @brief Clear all grid data and reset dimensions. */
			virtual void clear() {
				data_density.clear();
				data_temp.clear();
				memset(dims,   0, 3 * sizeof(size_t));
				memset(offset, 0, 3 * sizeof(size_t));
			}

			/**
			 * @brief Allocate and zero-initialise grid buffers.
			 * @param x  Width  of the grid.
			 * @param y  Height of the grid.
			 * @param z  Depth  of the grid.
			 * @param allocate_data_temp  If true, allocate data_temp buffer for normalization.
			 */
			virtual void create(size_t x, size_t y, size_t z, bool allocate_data_temp = false) {
				dims[0] = x;  dims[1] = y;  dims[2] = z;
				data_density.resize(size());
				memset(data_density.data(), 0, memsize());
				if (allocate_data_temp) {
					data_temp.resize(size());
					memset(data_temp.data(), 0, memsize());
				}
			}

			// ── Dimension accessors ────────────────────────────────────────────────

			/** @brief Grid width. */
			size_t x() const { return dims[0]; }
			/** @brief Grid height. */
			size_t y() const { return dims[1]; }
			/** @brief Grid depth. */
			size_t z() const { return dims[2]; }
			/** @brief Total number of voxels. */
			size_t size()    const { return dims[0] * dims[1] * dims[2]; }
			/** @brief Total memory size in bytes. */
			size_t memsize() const { return dims[0] * dims[1] * dims[2] * sizeof(float); }
			/**
			 * @brief Convert 3-D coordinates to a linear index.
			 * @param x X coordinate  @param y Y coordinate  @param z Z coordinate
			 * @return Linear array index.
			 */
			size_t get_index(size_t x, size_t y, size_t z) const {
				return x + y * dims[0] + z * dims[0] * dims[1];
			}
		};

		class VoxelSparseManager {
		public:
			// ── Sparse-grid interface (backward-compat stubs) ──────────────────────

			virtual void set_transform(double ts) { transform_scale = ts; }
			virtual void init(unsigned int /*expected_voxels*/) {}
			// virtual int  update() { return 0; }
			virtual void insertOrUpdatePackedSequential(uint64_t /*key*/, float /*value*/) {}
			virtual void serialize(std::vector<uint8_t> & /*bin_data*/) {}
			virtual void deserialize(std::vector<uint8_t> & /*bin_data*/) {}
			virtual void merge(VoxelSparseManager* /*other*/) {}
			virtual void merge(std::vector<uint8_t>& /*bin_data*/) {}
			virtual size_t mem_size() const { return 0; }

			virtual void insertOrUpdate(void* /*h_voxels*/, size_t /*num_voxels*/) {}
			//virtual int extractAll(Voxel** h_output_voxels) {}
			virtual void clear() {}

		public:
			double transform_scale = 0.0;
		};		

		/**
		 * @brief Container for different VDB grid representations
		 *
		 * Can hold particle data in various formats: dense grid, sparse NanoVDB,
		 * OpenVDB, serialized binary, or raw particle data.
		 */
		class VDBParticles
		{
		public:
			std::shared_ptr<VoxelDenseManager> dense_grid;  ///< Dense regular grid representation
			std::vector<uint8_t> vector_grid;  ///< Serialized binary grid data (for MPI transfer)
			std::shared_ptr<VoxelSparseManager> sparse_grid;  ///< Simple sparse grid representation (i, j, k, value)
			RawParticles raw_particles;  ///< Raw particle point cloud data

			VDBParticleType type;  ///< Current representation type
		};
	}

	/**
	 * SpaceData class stores and manages configuration and metadata for spatial data conversions
	 * Contains all parameters needed for particle-to-VDB conversion operations
	 */
	class SpaceData {
	public:
		/**
		 * Message types for inter-process communication
		 */
		enum class MessageType {
			eExit = -1,   // Signal to exit the process
			eEmpty = 0,   // Empty message
			eInfo = 1,    // Informational message
			eData = 2,    // Data processing message
			eBBOX = 3     // Bounding box message
		};

		/**
		 * Density kernel types for SPH (Smoothed Particle Hydrodynamics) calculations
		 */
		enum class DenseType {
			eNone = 0,        // No density computation
			eCubic = 1,       // Cubic spline kernel
			eQuintic = 2,     // Quintic spline kernel
			eWendlandC2 = 3,  // Wendland C2 kernel
			eWendlandC4 = 4,  // Wendland C4 kernel
			eWendlandC6 = 5,  // Wendland C6 kernel (recommended)
			eWendlandC8 = 6,  // Wendland C8 kernel
		};

		/**
		 * Normalization methods for density calculations
		 */
		enum class DenseNorm {
			eNone = 0,              // No normalization
			eCount = 1,             // Normalize by particle count
			eSPHInterpolation = 2,  // SPH interpolation normalization
		};

		/**
		 * Output data format types
		 */
		enum class ExtractedType {
			eSparse = 0,    // Sparse grid (default)
			eDense = 1,     // Dense volumetric grid
			eParticle = 2   // Raw particle data
		};

		enum class ExtractedParticleType {
			eNone = 0,   // default
			eRAW = 1,    // Raw particle data
			eVDB = 2,    // OpenVDB point cloud format
			eVTP = 3,    // VTK PolyData format
		};
		
		enum class ExtractedDenseType {
			eNone = 0,   // default
			eRAW = 1,    // Raw particle data
			eVTI = 2,    // VTK Image Data format
		};		

		/**
		 * Animation sequence processing modes
		 */
		enum class AnimType {
			eNone = 0,          // No animation (default)
			eAllPath = 1,       // Process all frames separately
			eAllMerge = 2,      // Merge all frames
			eFrameExtract = 3,  // Extract specific frame
		};

	public:
		// ========== Communication and Processing ==========
		MessageType message_type = MessageType::eEmpty;  // Type of message for inter-process communication
		int particle_type = 0;        // Particle type identifier
		int block_name_id = 0;        // Data block identifier

		// ========== Grid Configuration ==========
		float grid_transform = 1.0f;  // Grid transformation scale
		ExtractedType extracted_type = ExtractedType::eSparse;  // Output format type
		ExtractedParticleType extracted_particle_type = ExtractedParticleType::eNone;  // Particle output format type
		DenseType dense_type = DenseType::eNone;                // Density kernel type
		DenseNorm dense_norm = DenseNorm::eNone;                // Density normalization method
		ExtractedDenseType extracted_dense_type = ExtractedDenseType::eNone;  // Dense output format type
		float object_size = 1000.0f;                            // Object size in world units

		// ========== Value Ranges ==========
		float min_value = 0.0f;          // Minimum attribute value
		float max_value = 1.0f;          // Maximum attribute value
		float min_rho = 0.0f;            // Minimum density value
		float max_rho = 1.0f;            // Maximum density value
		float min_value_reduced = 0.0f;  // Reduced minimum after MPI reduction
		float max_value_reduced = 1.0f;  // Reduced maximum after MPI reduction

		// ========== Particle Parameters ==========
		float particle_radius_multiplier = 0.0f;  // Multiplier for particle radius
		double particle_radius_const = 0.0; // Fixed particle radius (0.0 = use smoothing length)
		float filter_min = -FLT_MAX;     // Minimum filter threshold
		float filter_max = FLT_MAX;      // Maximum filter threshold

		// ========== Counts and Statistics ==========
		size_t particles_count = 0;       // Total particle count
		size_t voxels_count = 0;          // Total voxel count
		size_t particles_count_local = 0; // Local particle count (for MPI)

		// ========== Bounding Box Configuration ==========
		float bbox_min[3] = { 0.0f, 0.0f, 0.0f };          // Min bounding box (world coordinates)
		float bbox_max[3] = { 1000.0f, 1000.0f, 1000.0f }; // Max bounding box (world coordinates)
		int bbox_dim = 100;                                 // Grid resolution
		int bbox_min_orig_local[3] = { 0,0,0 };             // Local min bbox (grid coordinates)
		int bbox_max_orig_local[3] = { 0,0,0 };             // Local max bbox (grid coordinates)
		int bbox_min_orig[3] = { 0,0,0 };                   // Global min bbox (grid coordinates)
		int bbox_max_orig[3] = { 0,0,0 };                   // Global max bbox (grid coordinates)
		double bbox_size_orig = 0.0;                        // Original bounding box size

		// ========== Local Values (MPI) ==========
		float min_value_local = 0.0f;  // Local minimum value per MPI rank
		float max_value_local = 0.0f;  // Local maximum value per MPI rank
		double transform_scale = 0.0;  // Computed transform scale

		// ========== Animation Parameters ==========
		int frame = 0;                      // Current frame number
		AnimType anim_type = AnimType::eNone;  // Animation processing mode
		std::string full_filepath;          // Output VDB file path
		int anim_task_counter = 0;          // Animation task counter
		int anim_start = -1;                // Start frame (-1 = not set)
		int anim_end = -1;                  // End frame (-1 = not set)
		int anim_step = -1;                 // Step frame (-1 = not set)

		// ========== SPH Radius Calculation ==========
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
		int calc_radius_neigh = -1;  // Number of neighbors for radius calculation
		DenseType calc_radius_neigh_rho_kernel = DenseType::eWendlandC6;  // Kernel for radius calculation
#endif

#ifdef WITH_CUDAKDTREE
		// ========== Voxel-Centric Dense Loop ==========
		/// When true: dense grid conversion loops over voxels and queries nearest particles
		/// via KD-tree instead of looping over particles and splatting. Requires calc_radius_neigh > 0.
		bool use_dense_loop_over_voxels = false;
#endif
		std::string calc_radius_neigh_file = "";  // Pre-computed radius file

		// ========== Additional Filtering ==========
		bool use_bbox_sphere = false;           // Use spherical bounding region
		float bbox_sphere_pos[3] = { 0.0f,0.0f,0.0f };  // Sphere center position
		float bbox_sphere_r = 0.0f;             // Sphere radius
		bool use_simple_density = false;        // Use simplified density calculation
		bool use_norm_value = true;		// Use normalized values
		float offset_position[3] = { 0.0f,0.0f,0.0f };  // Position offset

		/**
		 * @brief Print all current SpaceData configuration values.
		 * 
		 * This function outputs all member variables of the SpaceData structure
		 * to standard output, useful for debugging and verification of spatial
		 * data conversion parameters.
		 */
		void print_info() const
		{
			std::cout << "\n=== SpaceData Configuration ===" << std::endl;
			
			// Communication and Processing
			std::cout << "message_type: " << static_cast<int>(message_type) << std::endl;
			std::cout << "particle_type: " << particle_type << std::endl;
			std::cout << "block_name_id: " << block_name_id << std::endl;
			
			// Grid Configuration
			std::cout << "grid_transform: " << grid_transform << std::endl;
			std::cout << "extracted_type: " << static_cast<int>(extracted_type) << std::endl;
			std::cout << "extracted_particle_type: " << static_cast<int>(extracted_particle_type) << std::endl;
			std::cout << "dense_type: " << static_cast<int>(dense_type) << std::endl;
			std::cout << "dense_norm: " << static_cast<int>(dense_norm) << std::endl;
			std::cout << "extracted_dense_type: " << static_cast<int>(extracted_dense_type) << std::endl;
			std::cout << "object_size: " << object_size << std::endl;
			
			// Value Ranges
			std::cout << "min_value: " << min_value << std::endl;
			std::cout << "max_value: " << max_value << std::endl;
			std::cout << "min_rho: " << min_rho << std::endl;
			std::cout << "max_rho: " << max_rho << std::endl;
			std::cout << "min_value_reduced: " << min_value_reduced << std::endl;
			std::cout << "max_value_reduced: " << max_value_reduced << std::endl;
			
			// Particle Parameters
			std::cout << "particle_radius_multiplier: " << particle_radius_multiplier << std::endl;
			std::cout << "particle_radius_const: " << particle_radius_const << std::endl;
			std::cout << "filter_min: " << filter_min << std::endl;
			std::cout << "filter_max: " << filter_max << std::endl;
			
			// Counts and Statistics
			std::cout << "particles_count: " << particles_count << std::endl;
			std::cout << "voxels_count: " << voxels_count << std::endl;
			std::cout << "particles_count_local: " << particles_count_local << std::endl;
			
			// Bounding Box Configuration
			std::cout << "bbox_min: [" << bbox_min[0] << ", " << bbox_min[1] << ", " << bbox_min[2] << "]" << std::endl;
			std::cout << "bbox_max: [" << bbox_max[0] << ", " << bbox_max[1] << ", " << bbox_max[2] << "]" << std::endl;
			std::cout << "bbox_dim: " << bbox_dim << std::endl;
			std::cout << "bbox_min_orig_local: [" << bbox_min_orig_local[0] << ", " << bbox_min_orig_local[1] << ", " << bbox_min_orig_local[2] << "]" << std::endl;
			std::cout << "bbox_max_orig_local: [" << bbox_max_orig_local[0] << ", " << bbox_max_orig_local[1] << ", " << bbox_max_orig_local[2] << "]" << std::endl;
			std::cout << "bbox_min_orig: [" << bbox_min_orig[0] << ", " << bbox_min_orig[1] << ", " << bbox_min_orig[2] << "]" << std::endl;
			std::cout << "bbox_max_orig: [" << bbox_max_orig[0] << ", " << bbox_max_orig[1] << ", " << bbox_max_orig[2] << "]" << std::endl;
			std::cout << "bbox_size_orig: " << bbox_size_orig << std::endl;
			
			// Local Values (MPI)
			std::cout << "min_value_local: " << min_value_local << std::endl;
			std::cout << "max_value_local: " << max_value_local << std::endl;
			std::cout << "transform_scale: " << transform_scale << std::endl;
			
			// Animation Parameters
			std::cout << "frame: " << frame << std::endl;
			std::cout << "anim_type: " << static_cast<int>(anim_type) << std::endl;
			std::cout << "full_filepath: " << full_filepath << std::endl;
			std::cout << "anim_task_counter: " << anim_task_counter << std::endl;
			std::cout << "anim_start: " << anim_start << std::endl;
			std::cout << "anim_end: " << anim_end << std::endl;
			std::cout << "anim_step: " << anim_step << std::endl;
			
			// SPH Radius Calculation
#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
			std::cout << "calc_radius_neigh: " << calc_radius_neigh << std::endl;
			std::cout << "calc_radius_neigh_rho_kernel: " << static_cast<int>(calc_radius_neigh_rho_kernel) << std::endl;
#endif
			std::cout << "calc_radius_neigh_file: " << calc_radius_neigh_file << std::endl;
#ifdef WITH_CUDAKDTREE
			std::cout << "use_dense_loop_over_voxels: " << (use_dense_loop_over_voxels ? "true" : "false") << std::endl;
#endif
			
			// Additional Filtering
			std::cout << "use_bbox_sphere: " << (use_bbox_sphere ? "true" : "false") << std::endl;
			std::cout << "bbox_sphere_pos: [" << bbox_sphere_pos[0] << ", " << bbox_sphere_pos[1] << ", " << bbox_sphere_pos[2] << "]" << std::endl;
			std::cout << "bbox_sphere_r: " << bbox_sphere_r << std::endl;
			std::cout << "use_simple_density: " << (use_simple_density ? "true" : "false") << std::endl;
			std::cout << "use_norm_value: " << (use_norm_value ? "true" : "false") << std::endl;
			std::cout << "offset_position: [" << offset_position[0] << ", " << offset_position[1] << ", " << offset_position[2] << "]" << std::endl;
			
			std::cout << "================================\n" << std::endl;
		}
	};

} // namespace space_converter
