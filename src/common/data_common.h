#pragma once

#include <string>
#include <float.h>
#include <cstring>

#define FDATA_EPSILON 1e-8f // Define a small epsilon value for floating-point comparisons.

// Namespace common provides functionality for handling spatial data conversions.
namespace common {

	// SpaceData class stores and manages configuration and metadata for spatial data conversions.
	class SpaceData {
	public:
		// Enumeration for message types used in communication or processing.
		enum class MessageType {
			eExit = -1,  // Signal to exit the process.
			eEmpty = 0,
			eInfo = 1,   // Informational message.
			eData = 2,    // Data processing message.
			eBBOX = 3 // Bounding box message.
		};

		//// Enumeration for value conversion types.
		//enum class ValueConvertType {
		//    eDefault = 0 // Default conversion type.
		//};

		// Enumeration for density calculation methods.
		enum class DenseType {
			eNone = 0,     // No density computation.
			eCubic = 1,
			eQuintic = 2,
			eWendlandC2 = 3,
			eWendlandC4 = 4,
			eWendlandC6 = 5,
			eWendlandC8 = 6,
		};

		enum class DenseNorm {
			eNone = 0,     // No normalization.
			eCount = 1,
			eSPHInterpolation = 2,
		};

		enum class ExtractedType {
			eSparse = 0,     // Default - sparse.
			eDense = 1,
			eParticle = 2,
		};

		enum class AnimType {
			eNone = 0,
			eAllPath = 1,
			eAllMerge = 2,
			eFrameExtract = 3,
			//eFrameCache = 4,
		};

	public:
		// Public member variables to store configuration and metadata.
		MessageType message_type = MessageType::eEmpty;// Type of the message (see MessageType).
		int particle_type = 0;       // Type of particle being processed.
		int block_name_id = 0;       // Identifier for the block name.

		float grid_transform = 1.0f;    // Transformation applied to the grid.

		//ValueConvertType value_convert;// Type of value conversion (see ValueConvertType).
		ExtractedType extracted_type = ExtractedType::eSparse; // Type of extraction (see ExtractedType).
		DenseType dense_type = DenseType::eNone;    // Type of density calculation (see DenseType).
		DenseNorm dense_norm = DenseNorm::eNone;    // Type of density calculation (see DenseNorm).

		float object_size = 1000.0f;       // Size of the object in the grid.

		float min_value = 0.0f;         // Minimum value in the data.
		float max_value = 1.0f;         // Maximum value in the data.

		float min_rho = 0.0f;           // Minimum density value.
		float max_rho = 1.0f;           // Maximum density value.

		float min_value_reduced = 0.0f; // Reduced minimum value after processing.
		float max_value_reduced = 1.0f; // Reduced maximum value after processing.

		float particle_fix_size = 0.0f; // Flag indicating whether particle length is used.
		float filter_min = -FLT_MAX;        // Minimum filter value.
		float filter_max = FLT_MAX;        // Maximum filter value.

		size_t particles_count = 0;  // Total number of particles.
		size_t voxels_count = 0;     // Total number of voxels.

		float bbox_min[3] = { 0.0f, 0.0f, 0.0f };       // Minimum bounding box coordinates (x, y, z).
		float bbox_max[3] = { 1000.0f, 1000.0f, 1000.0f };       // Maximum bounding box coordinates (x, y, z).
		int bbox_dim = 100;            // Dimensions of the bounding box.

		int bbox_min_orig_local[3] = { 0,0,0 }; // Local original minimum bounding box coordinates.
		int bbox_max_orig_local[3] = { 0,0,0 }; // Local original maximum bounding box coordinates.

		int bbox_min_orig[3] = { 0,0,0 };    // Original minimum bounding box coordinates.
		int bbox_max_orig[3] = { 0,0,0 };    // Original maximum bounding box coordinates.

		//int bbox_min_orig_per_type[3];    // Original minimum bounding box coordinates.
		//int bbox_max_orig_per_type[3];    // Original maximum bounding box coordinates.

		double bbox_size_orig = 0.0;         // Size of the bounding box.

		float min_value_local = 0.0f;   // Local minimum value.
		float max_value_local = 0.0f;   // Local maximum value.
		size_t particles_count_local = 0; // Local particle count.

		double transform_scale = 0.0;  // Scale transformation applied to the data.

		int frame = 0;
		AnimType anim_type = AnimType::eNone;
		std::string full_filepath; //path to saved vdb
		int anim_task_counter = 0;

		int anim_start = -1;               // Start frame (-1 indicates not set).
		int anim_end = -1;                 // End frame (-1 indicates not set).

#if defined(WITH_CUDAKDTREE) || defined(WITH_NANOFLANN)
		int calc_radius_neigh = -1;
		DenseType calc_radius_neigh_rho_kernel = DenseType::eWendlandC6;
#endif
		std::string calc_radius_neigh_file = "";

		bool use_bbox_sphere = false;
		float bbox_sphere_pos[3] = { 0.0f,0.0f,0.0f };
		float bbox_sphere_r = 0.0f;
		bool use_simple_density = false;
		float offset_position[3] = { 0.0f,0.0f,0.0f };

	//public:
		// Constructor to initialize SpaceData with default or provided values.
		// @param export_type: Type of export (e.g., particle type).
		// @param export_dataset: Dataset identifier for the export.
		// @param grid_dim: Dimension of the computational grid.
   //     SpaceData(
   //         //int export_type, 
   //         //int export_dataset,
   //         //int export_extracted_type,
   //         //int export_dense_type, 
   //         //int export_dense_norm,
   //         //int grid_dim, 
   //         //bool use_bbox, 
   //         //float bbox_pos[6]
   //     ) {

   //         //message_type = MessageType::eEmpty;
   //         //particle_type = export_type;
   //         //block_name_id = export_dataset;

   //         //grid_transform = 1.0f;
   //         //value_convert = ValueConvertType::eDefault;
   ////         extracted_type = (ExtractedType)export_extracted_type;
   ////         dense_type = (DenseType)export_dense_type;
			////dense_norm = (DenseNorm)export_dense_norm;

   //         //object_size = 1000.0f;

   //         //min_value = 0.0f;
   //         //max_value = 1.0f;

   //         //min_rho = 0.0f;
   //         //max_rho = 1.0f;

   //         //min_value_reduced = 0.0f;
   //         //max_value_reduced = 1.0f;

   //         //particle_fix_size = 1.2f;
   //         //filter_min = -FLT_MAX;
   //         //filter_max = FLT_MAX;

   //         //particles_count = 0;
   //         //voxels_count = 0;

   //         //if (use_bbox) {
   //         //    bbox_min[0] = bbox_pos[0];
   //         //    bbox_min[1] = bbox_pos[1];
   //         //    bbox_min[2] = bbox_pos[2];

   //         //    bbox_max[0] = bbox_pos[3];
   //         //    bbox_max[1] = bbox_pos[4];
   //         //    bbox_max[2] = bbox_pos[5];
   //         //}
   //         //else {
   //             bbox_min[0] = 0.0f;
   //             bbox_min[1] = 0.0f;
   //             bbox_min[2] = 0.0f;

   //             bbox_max[0] = object_size;
   //             bbox_max[1] = object_size;
   //             bbox_max[2] = object_size;
   //         //}

   //         //bbox_dim = grid_dim;

   //         //memset(bbox_min_orig_local, 0, sizeof(bbox_min_orig_local));
   //         //memset(bbox_max_orig_local, 0, sizeof(bbox_max_orig_local));

   //         //memset(bbox_min_orig, 0, sizeof(bbox_min_orig));
   //         //memset(bbox_max_orig, 0, sizeof(bbox_max_orig));

   //         //memset(bbox_min_orig_per_type, 0, sizeof(bbox_min_orig_per_type));
   //         //memset(bbox_max_orig_per_type, 0, sizeof(bbox_max_orig_per_type));

   //         //bbox_size_orig = 0.0;

   //         //min_value_local = 0.0f;
   //         //max_value_local = 0.0f;
   //         //particles_count_local = 0;

   //         //transform_scale = 0.0;
   //         //frame = 0;
   //         //anim_type = AnimType::eNone;
   //         //anim_task_counter = 0;

   //         //anim_start = -1;               // Start frame (-1 indicates not set).
   //         //anim_end = -1;                 // End frame (-1 indicates not set).
   //     }
	};

} // namespace space_converter
