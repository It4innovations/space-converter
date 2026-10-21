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

#include "pluto_vtk_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>

#include <fstream>
#include <sstream>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include "convert_common.h"

#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkRectilinearGridReader.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#define RETURN_NORM_EMPTY return 0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)space_converter::common::calculate_dmagnitude3(v[0], v[1], v[2]);

#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;

#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)

namespace plutovtk {
	namespace io {

		// VTK data storage
		vtkSmartPointer<vtkRectilinearGrid> vtk_grid = nullptr;
		std::vector<std::string> scalar_names;
		
		// Grid dimensions
		int grid_dims[3] = {0, 0, 0};
		int64_t total_cells = 0;
		int64_t local_cells = 0;
		
		// Coordinate arrays
		std::vector<float> x_coords;
		std::vector<float> y_coords;
		std::vector<float> z_coords;
		
		// Cell size (for HSML calculation)
		double cell_size = 0.0;
		
		double steps_time[2];

		void init_lib(
			std::string vtk_file, 
			int world_rank, 
			int world_size,
			std::vector<std::string>& scalar_names_vec
		) {
#ifdef WITH_OPENMP
			steps_time[0] = omp_get_wtime();
#endif

			scalar_names = scalar_names_vec;

			// Read VTK file
			auto reader = vtkSmartPointer<vtkRectilinearGridReader>::New();
			reader->SetFileName(vtk_file.c_str());
			reader->Update();

			vtk_grid = reader->GetOutput();
			
			if (!vtk_grid) {
				throw std::runtime_error("Failed to read VTK file: " + vtk_file);
			}

			// Get grid dimensions
			vtk_grid->GetDimensions(grid_dims);
			
			// Calculate number of cells (dimensions - 1 for cells)
			//int cell_dims[3] = {
			//	grid_dims[0] - 1,
			//	grid_dims[1] - 1,
			//	grid_dims[2] - 1
			//};

			int cell_dims[3] = {
				grid_dims[0],
				grid_dims[1],
				grid_dims[2]
			};
			
			total_cells = (int64_t)cell_dims[0] * cell_dims[1] * cell_dims[2];
			
			// Distribute cells across MPI ranks
			int64_t cells_per_rank = total_cells / world_size;
			int64_t start_cell = world_rank * cells_per_rank;
			local_cells = (world_rank == world_size - 1) ? total_cells - start_cell : cells_per_rank;
			
			// Extract coordinate arrays
			vtkFloatArray* xCoords = vtkFloatArray::SafeDownCast(vtk_grid->GetXCoordinates());
			vtkFloatArray* yCoords = vtkFloatArray::SafeDownCast(vtk_grid->GetYCoordinates());
			vtkFloatArray* zCoords = vtkFloatArray::SafeDownCast(vtk_grid->GetZCoordinates());
			
			if (xCoords && yCoords && zCoords) {
				x_coords.resize(grid_dims[0]);
				y_coords.resize(grid_dims[1]);
				z_coords.resize(grid_dims[2]);
				
				for (int i = 0; i < grid_dims[0]; i++) {
					x_coords[i] = xCoords->GetValue(i);
				}
				for (int i = 0; i < grid_dims[1]; i++) {
					y_coords[i] = yCoords->GetValue(i);
				}
				for (int i = 0; i < grid_dims[2]; i++) {
					z_coords[i] = zCoords->GetValue(i);
				}
				
				// Calculate approximate cell size
				if (grid_dims[0] > 1) {
					cell_size = (x_coords[1] - x_coords[0]);
				}
			}

			if (world_rank == 0) {
				printf("VTK Grid Dimensions: %d x %d x %d\n", grid_dims[0], grid_dims[1], grid_dims[2]);
				printf("Total cells: %lld\n", total_cells);
				printf("Cell size: %f\n", cell_size);
			}
			
			printf("Rank %d: Local cells: %lld\n", world_rank, local_cells);

#ifdef WITH_OPENMP
			steps_time[1] = omp_get_wtime();
#endif
		}

		void finish_lib() {
			vtk_grid = nullptr;
			x_coords.clear();
			y_coords.clear();
			z_coords.clear();
		}

		uint64_t get_particle_type_offset(PlutoVTKParticleType pt) {
			switch (pt) {
			case PlutoVTKParticleType::PLUTO:
				return 0;
			}
			return 0;
		}

		void print_CPU_steps() {
			//printf("init_lib time: %f\n", steps_time[1] - steps_time[0]);
		}

		float get_particle_norm_value(int blocknr, uint64_t id) {
			PlutoVTKParticleType pt = (PlutoVTKParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			if (blocknr >= PlutoVTKBlockType::BTMax && vtk_grid) {
				// Access scalar data from cell data
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;
				
				if (scalar_idx < vtk_grid->GetCellData()->GetNumberOfArrays()) {
					vtkDataArray* array = vtk_grid->GetCellData()->GetArray(scalar_idx);
					if (array && id - offset < array->GetNumberOfTuples()) {
						if (array->GetNumberOfComponents() == 1) {
							RETURN_NORM_VALUE(array->GetComponent(id - offset, 0));
						}
						else if (array->GetNumberOfComponents() == 3) {
							float v[3];
							v[0] = array->GetComponent(id - offset, 0);
							v[1] = array->GetComponent(id - offset, 1);
							v[2] = array->GetComponent(id - offset, 2);
							RETURN_NORM_VECTOR3(v);
						}
					}
				}
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value) {
			PlutoVTKParticleType pt = (PlutoVTKParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			if (blocknr >= PlutoVTKBlockType::BTMax && vtk_grid) {
				// Access scalar data from cell data
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;
				
				if (scalar_idx < vtk_grid->GetCellData()->GetNumberOfArrays()) {
					vtkDataArray* array = vtk_grid->GetCellData()->GetArray(scalar_idx);
					if (array && id - offset < array->GetNumberOfTuples()) {
						if (array->GetNumberOfComponents() == 1) {
							RETURN_ORIG_VALUE(array->GetComponent(id - offset, 0));
						}
						else if (array->GetNumberOfComponents() == 3) {
							float v[3];
							v[0] = array->GetComponent(id - offset, 0);
							v[1] = array->GetComponent(id - offset, 1);
							v[2] = array->GetComponent(id - offset, 2);
							RETURN_ORIG_VECTOR3(v);
						}
					}
				}
			}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id) {
			PlutoVTKParticleType pt = (PlutoVTKParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			if (blocknr >= PlutoVTKBlockType::BTMax && vtk_grid) {
				// Access scalar data from cell data
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;
				
				if (scalar_idx < vtk_grid->GetCellData()->GetNumberOfArrays()) {
					vtkDataArray* array = vtk_grid->GetCellData()->GetArray(scalar_idx);
					if (array && id - offset < array->GetNumberOfTuples()) {
						return array->GetNumberOfComponents();
					}
				}
			}

			RETURN_COMP_EMPTY;
		}

		int get_particle_type(uint64_t id) {
			// All cells are of type PLUTO
			return PlutoVTKParticleType::PLUTO;
		}

		void get_particle_position(uint64_t id, double* pos) {
			// Convert linear cell ID to (i, j, k) indices
			//int cell_dims[3] = {
			//	grid_dims[0] - 1,
			//	grid_dims[1] - 1,
			//	grid_dims[2] - 1
			//};
			int cell_dims[3] = {
				grid_dims[0],
				grid_dims[1],
				grid_dims[2]
			};
			
			int k = id / (cell_dims[0] * cell_dims[1]);
			int j = (id % (cell_dims[0] * cell_dims[1])) / cell_dims[0];
			int i = id % cell_dims[0];
			
			// Get cell center position
			if (i < x_coords.size() /* - 1*/ && j < y_coords.size() /* - 1 */ && k < z_coords.size() /* - 1 */) {
				pos[0] = x_coords[i]; // (x_coords[i] + x_coords[i + 1]) * 0.5;
				pos[1] = y_coords[j]; // (y_coords[j] + y_coords[j + 1]) * 0.5;
				pos[2] = z_coords[k]; // (z_coords[k] + z_coords[k + 1]) * 0.5;
			}
			else {
				pos[0] = pos[1] = pos[2] = 0.0;
			}
		}

		size_t get_local_num_particles() {
			return local_cells;
		}

		size_t get_global_num_particles() {
			return total_cells;
		}

		double get_particle_hsml(uint64_t id) {
			// Return cell size as smoothing length
			return cell_size;
		}

		double get_particle_mass(uint64_t id) {
			// Calculate mass from density and cell volume
			double rho = get_particle_rho(id);
			double volume = cell_size * cell_size * cell_size;
			return rho * volume;
		}

		double get_particle_rho(uint64_t id) {
			// Get density from first scalar field (usually "rho")
			if (vtk_grid && vtk_grid->GetCellData()->GetNumberOfArrays() > 0) {
				vtkDataArray* array = vtk_grid->GetCellData()->GetArray(0);
				if (array && id < array->GetNumberOfTuples()) {
					return array->GetComponent(id, 0);
				}
			}
			return 1.0; // Default density
		}

		int get_particle_rho_blocknr() {
			// First scalar field is rho
			return 0;// PlutoVTKBlockType::BTMax;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			int num_scalars = vtk_grid ? vtk_grid->GetCellData()->GetNumberOfArrays() : 0;
			
			types_and_blocks.resize(PlutoVTKParticleType::PTMax * num_scalars, 0);
			
			// Mark all scalar fields as available for PLUTO type
			for (int i = 0; i < num_scalars; i++) {
				types_and_blocks[PlutoVTKParticleType::PTMax * i + PlutoVTKParticleType::PLUTO] = 1;
			}
		}

		void print_types_and_blocks_local() {
			std::vector<int> types_and_blocks;
			get_types_and_blocks(types_and_blocks);

			printf("\n");
			printf("Type: PLUTO (0)\n");
			
			if (vtk_grid) {
				int num_scalars = vtk_grid->GetCellData()->GetNumberOfArrays();
				for (int i = 0; i < num_scalars; i++) {
					const char* name = vtk_grid->GetCellData()->GetArrayName(i);
					printf("\t%s (%d)\n", name, i);
				}
			}
		}

		void print_types_and_blocks(std::vector<int>& types_and_blocks) {
			printf("\nAll snapshots contain:\n");
			print_types_and_blocks_local();
		}

		std::string get_dataset_name(int blocknr) {
			if (vtk_grid && blocknr >= PlutoVTKBlockType::BTMax) {
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;
				if (scalar_idx < vtk_grid->GetCellData()->GetNumberOfArrays()) {
					const char* name = vtk_grid->GetCellData()->GetArrayName(scalar_idx);
					return name ? name : "unknown";
				}
			}
			return "unknown";
		}

		std::string get_particle_unit(int blocknr) {
			// Units would need to be defined in VTK metadata
			return "";
		}

	} // namespace io
} // namespace plutovtk
