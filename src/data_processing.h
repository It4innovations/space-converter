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

#pragma once

#include "args_processing.h"
#include "convert_vdb.h"
#include "data_common.h"

 // Namespace space_converter provides functions for managing VDB conversions and related operations.
namespace space_converter {

    // Initialize the VDB converter with command-line arguments and configuration.
    // @param argc: Number of command-line arguments.
    // @param argv: Array of command-line argument strings.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @return A pointer to an initialized ConvertVDBBase object.
    common::vdb::ConvertVDBBase* init_converter(int argc, char** argv, space_converter::FromCL& from_cl, common::SpaceData& space_data);

    // Deinitialize and clean up the VDB converter.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object to deinitialize.
    void deinit_converter(common::vdb::ConvertVDBBase* convert_vdb_base);

    // Print information about the converter and global particle types and blocks.
    // @param convert_vdb_base: Pointer to the initialized ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param types_and_blocks_global: Vector containing global types and blocks information.
    void print_info(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, std::vector<int>& types_and_blocks_global);

    // Calculate the bounding box for the particle data.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object to store bounding box information.
    void find_bbox(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& spaceData, int particle_type = -1);

    // Create the main VDB grid based on particle data and configuration.
    // @param grid_main: Reference to the VDBParticles object representing the main grid.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object with metadata for grid creation.
    void create_grid(common::vdb::VDBParticles& grid_main, space_converter::FromCL& from_cl, common::SpaceData& spaceData);

    // Convert particle data to the VDB grid.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object containing data for conversion.
    // @param grid_main: Reference to the VDBParticles object to populate with converted data.
    void convert_to_grid(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_main);

    // Find the minimum and maximum values in the data for scaling or filtering.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object to store min/max values.
    void find_minmax_value(space_converter::FromCL& from_cl, common::SpaceData& spaceData);
    void find_minmax_reduced_value(space_converter::FromCL& from_cl, common::SpaceData& spaceData);

    // Perform a reduction operation to combine data from multiple grids or processes.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object for metadata.
    // @param grid_main: Reference to the VDBParticles object for the main grid.
    // @param grid_main_sum: Reference to the VDBParticles object to store reduced data.
    void reduction(common::vdb::ConvertVDBBase* convert_vdb_base, space_converter::FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_main, common::vdb::VDBParticles& grid_main_sum);

    // Finalize the VDB grid by applying transformations and optimizations.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object with metadata.
    // @param grid_main_sum: Reference to the VDBParticles object containing reduced data.
    // @param grid_main_final: Reference to the VDBParticles object for the finalized grid.
    void finalize_grid(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_main_sum, common::vdb::VDBParticles& grid_main_final);

    // Save the finalized VDB grid to a file.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object with metadata.
    // @param grid_main_final: Reference to the VDBParticles object for the finalized grid to save.
    void save_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_main_final, common::vdb::VDBParticles::VDBParticleType particle_type, bool only_rank0 = true);

    // Save the finalized RAW grid to a file.
    // @param convert_vdb_base: Pointer to the ConvertVDBBase object.
    // @param from_cl: Reference to the FromCL struct for configuration.
    // @param spaceData: Reference to the SpaceData object with metadata.
    // @param grid_main: Reference to the VDBParticles object for the finalized grid to save.
    void save_raw_volume(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_main, bool only_rank0 = true);

    void save_raw_particles_to_vdb(common::vdb::ConvertVDBBase* convert_vdb_base, FromCL& from_cl, common::SpaceData& spaceData, common::vdb::VDBParticles& grid_main);

    void test_converter(int argc, char** argv, space_converter::FromCL& from_cl);

} // namespace space_converter
