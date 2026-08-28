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

#include "reader_return_macros.h"

// Reader data contract (see docs/SpaceConverter_Code_Analysis_2026-08.md §7):
//   get_particle_mass(id) -> density of the cell times the actual local cell volume
//                            (code units of the PLUTO run), for the single PLUTO type
//   get_particle_rho(id)  -> value of the first selected cell scalar array (usually "rho"),
//                            in the code units stored in the VTK file
//   get_particle_hsml(id) -> maximum of the local grid spacings (dx, dy, dz) of the cell
//   particle id space     -> global cell index space [0, total_cells), decomposed per rank
//                            as [g_start_cell, g_start_cell + local_cells); local id + g_start_cell
//                            gives the global (i, j, k) cell index (x fastest)

namespace plutovtk {
	namespace io {

		// VTK data storage
		vtkSmartPointer<vtkRectilinearGrid> vtk_grid = nullptr;
		std::vector<std::string> scalar_names;

		// VTK cell-array indices of the scalars in use
		// (order given by --scalar-names, or file order when no names were given)
		std::vector<int> selected_scalars;

		// Grid dimensions
		int grid_dims[3] = {0, 0, 0};
		int cell_dims[3] = {1, 1, 1};
		int64_t total_cells = 0;
		int64_t local_cells = 0;
		int64_t g_start_cell = 0;	// first global cell index owned by this rank

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

			// Calculate number of cells (point dimensions - 1 per axis; degenerate axes keep 1)
			for (int i = 0; i < 3; i++) {
				cell_dims[i] = (grid_dims[i] > 1) ? grid_dims[i] - 1 : 1;
			}

			total_cells = (int64_t)cell_dims[0] * cell_dims[1] * cell_dims[2];

			// Distribute cells across MPI ranks
			int64_t cells_per_rank = total_cells / world_size;
			g_start_cell = world_rank * cells_per_rank;
			local_cells = (world_rank == world_size - 1) ? total_cells - g_start_cell : cells_per_rank;

			// Select the scalar cell arrays to expose
			selected_scalars.clear();
			int num_arrays = vtk_grid->GetCellData()->GetNumberOfArrays();
			if (!scalar_names.empty()) {
				// Only the named cell arrays, in the given order
				for (size_t n = 0; n < scalar_names.size(); n++) {
					int found = -1;
					for (int a = 0; a < num_arrays; a++) {
						const char* aname = vtk_grid->GetCellData()->GetArrayName(a);
						if (aname && scalar_names[n] == aname) {
							found = a;
							break;
						}
					}
					if (found >= 0) {
						selected_scalars.push_back(found);
					}
					else if (world_rank == 0) {
						printf("Warning: scalar '%s' not found in VTK file. Available cell arrays:", scalar_names[n].c_str());
						for (int a = 0; a < num_arrays; a++) {
							const char* aname = vtk_grid->GetCellData()->GetArrayName(a);
							printf(" %s", aname ? aname : "(unnamed)");
						}
						printf("\n");
					}
				}
			}
			else {
				for (int a = 0; a < num_arrays; a++) {
					selected_scalars.push_back(a);
				}
			}

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
			selected_scalars.clear();
		}

		// Convert a local particle id to the global cell index of this rank.
		static int64_t get_global_cell_id(uint64_t id) {
			return g_start_cell + (int64_t)id;
		}

		// Decode a global cell index into (i, j, k) cell indices (x fastest).
		static void decode_cell_id(int64_t gid, int& i, int& j, int& k) {
			k = (int)(gid / ((int64_t)cell_dims[0] * cell_dims[1]));
			j = (int)((gid % ((int64_t)cell_dims[0] * cell_dims[1])) / cell_dims[0]);
			i = (int)(gid % cell_dims[0]);
		}

		// Local grid spacing of cell index i along one coordinate axis (0 for degenerate axes).
		static double local_spacing(const std::vector<float>& coords, int i) {
			if (i + 1 < (int)coords.size())
				return (double)coords[i + 1] - (double)coords[i];
			return 0.0;
		}

		// Cell data array for the given index into the selected scalar list (nullptr if out of range).
		static vtkDataArray* get_selected_scalar_array(int scalar_idx) {
			if (!vtk_grid || scalar_idx < 0 || scalar_idx >= (int)selected_scalars.size())
				return nullptr;
			return vtk_grid->GetCellData()->GetArray(selected_scalars[scalar_idx]);
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
				// Access scalar data from cell data (selected scalar list)
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;

				vtkDataArray* array = get_selected_scalar_array(scalar_idx);
				int64_t gid = get_global_cell_id(id - offset);
				if (array && gid < array->GetNumberOfTuples()) {
					if (array->GetNumberOfComponents() == 1) {
						RETURN_NORM_VALUE(array->GetComponent(gid, 0));
					}
					else if (array->GetNumberOfComponents() == 3) {
						float v[3];
						v[0] = array->GetComponent(gid, 0);
						v[1] = array->GetComponent(gid, 1);
						v[2] = array->GetComponent(gid, 2);
						RETURN_NORM_VECTOR3(v);
					}
				}
			}

			RETURN_NORM_EMPTY;
		}

		int get_particle_value(int blocknr, uint64_t id, float* out_value) {
			PlutoVTKParticleType pt = (PlutoVTKParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			if (blocknr >= PlutoVTKBlockType::BTMax && vtk_grid) {
				// Access scalar data from cell data (selected scalar list)
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;

				vtkDataArray* array = get_selected_scalar_array(scalar_idx);
				int64_t gid = get_global_cell_id(id - offset);
				if (array && gid < array->GetNumberOfTuples()) {
					if (array->GetNumberOfComponents() == 1) {
						RETURN_ORIG_VALUE(array->GetComponent(gid, 0));
					}
					else if (array->GetNumberOfComponents() == 3) {
						float v[3];
						v[0] = array->GetComponent(gid, 0);
						v[1] = array->GetComponent(gid, 1);
						v[2] = array->GetComponent(gid, 2);
						RETURN_ORIG_VECTOR3(v);
					}
				}
			}

			RETURN_ORIG_EMPTY;
		}

		int get_particle_value_comp(int blocknr, uint64_t id) {
			PlutoVTKParticleType pt = (PlutoVTKParticleType)get_particle_type(id);
			uint64_t offset = get_particle_type_offset(pt);

			if (blocknr >= PlutoVTKBlockType::BTMax && vtk_grid) {
				// Access scalar data from cell data (selected scalar list)
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;

				vtkDataArray* array = get_selected_scalar_array(scalar_idx);
				int64_t gid = get_global_cell_id(id - offset);
				if (array && gid < array->GetNumberOfTuples()) {
					return array->GetNumberOfComponents();
				}
			}

			RETURN_COMP_EMPTY;
		}

		int get_particle_type(uint64_t id) {
			// All cells are of type PLUTO
			return PlutoVTKParticleType::PLUTO;
		}

		void get_particle_position(uint64_t id, double* pos) {
			// Convert linear global cell ID to (i, j, k) cell indices
			int i, j, k;
			decode_cell_id(get_global_cell_id(id), i, j, k);

			// Get cell center position (node coordinate on degenerate axes)
			if (i < (int)x_coords.size() && j < (int)y_coords.size() && k < (int)z_coords.size()) {
				pos[0] = (i + 1 < (int)x_coords.size()) ? (x_coords[i] + x_coords[i + 1]) * 0.5 : x_coords[i];
				pos[1] = (j + 1 < (int)y_coords.size()) ? (y_coords[j] + y_coords[j + 1]) * 0.5 : y_coords[j];
				pos[2] = (k + 1 < (int)z_coords.size()) ? (z_coords[k] + z_coords[k + 1]) * 0.5 : z_coords[k];
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
			// Return the maximum of the local grid spacings of this cell as smoothing length
			int i, j, k;
			decode_cell_id(get_global_cell_id(id), i, j, k);

			double dx = local_spacing(x_coords, i);
			double dy = local_spacing(y_coords, j);
			double dz = local_spacing(z_coords, k);

			double h = dx;
			if (dy > h) h = dy;
			if (dz > h) h = dz;

			return (h > 0.0) ? h : cell_size;
		}

		double get_particle_mass(uint64_t id) {
			// Calculate mass from density and the actual local cell volume
			int i, j, k;
			decode_cell_id(get_global_cell_id(id), i, j, k);

			double dx = local_spacing(x_coords, i);
			double dy = local_spacing(y_coords, j);
			double dz = local_spacing(z_coords, k);

			// Degenerate axes (a single coordinate plane) contribute a unit extent
			double volume = ((dx > 0.0) ? dx : 1.0) * ((dy > 0.0) ? dy : 1.0) * ((dz > 0.0) ? dz : 1.0);

			double rho = get_particle_rho(id);
			return rho * volume;
		}

		double get_particle_rho(uint64_t id) {
			// Get density from the first selected scalar field (usually "rho")
			vtkDataArray* array = get_selected_scalar_array(0);
			int64_t gid = get_global_cell_id(id);
			if (array && gid < array->GetNumberOfTuples()) {
				return array->GetComponent(gid, 0);
			}
			return 1.0; // Default density
		}

		int get_particle_rho_blocknr() {
			// First scalar field is rho
			return 0;// PlutoVTKBlockType::BTMax;
		}

		void get_types_and_blocks(std::vector<int>& types_and_blocks) {
			int num_scalars = (int)selected_scalars.size();

			// Blocks 0..num_scalars: block 0 is the rho alias (handled through
			// get_particle_rho_blocknr()), block b >= BTMax maps to selected scalar b - BTMax.
			types_and_blocks.resize(PlutoVTKParticleType::PTMax * (num_scalars + PlutoVTKBlockType::BTMax), 0);

			for (int i = 0; i < num_scalars + PlutoVTKBlockType::BTMax; i++) {
				types_and_blocks[PlutoVTKParticleType::PTMax * i + PlutoVTKParticleType::PLUTO] = 1;
			}
		}

		void print_types_and_blocks_local() {
			std::vector<int> types_and_blocks;
			get_types_and_blocks(types_and_blocks);

			printf("\n");
			printf("Type: PLUTO (0)\n");

			int num_blocks = (int)types_and_blocks.size() / PlutoVTKParticleType::PTMax;
			for (int i = 0; i < num_blocks; i++) {
				printf("\t%s (%d)\n", get_dataset_name(i).c_str(), i);
			}
		}

		void print_types_and_blocks(std::vector<int>& types_and_blocks) {
			printf("\nAll snapshots contain:\n");
			print_types_and_blocks_local();
		}

		std::string get_dataset_name(int blocknr) {
			if (blocknr == PlutoVTKBlockType::Rho) {
				// Block 0 is the rho alias (first selected scalar, see get_particle_rho)
				return "Rho";
			}
			if (vtk_grid && blocknr >= PlutoVTKBlockType::BTMax) {
				int scalar_idx = blocknr - PlutoVTKBlockType::BTMax;
				if (scalar_idx < (int)selected_scalars.size()) {
					const char* name = vtk_grid->GetCellData()->GetArrayName(selected_scalars[scalar_idx]);
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
