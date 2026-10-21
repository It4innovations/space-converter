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

#include "ipic3d_hdf5_extract_iolib.h"

#include <iostream>
#include <cstdio>
#include <string.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <thread>

#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef WITH_OPENMP
#include <omp.h>
#endif

#include <hdf5.h>

#include "convert_common.h"

#define RETURN_NORM_EMPTY return 0;
#define RETURN_NORM_VALUE(v) return (float)(v);
#define RETURN_NORM_VECTOR3(v) return (float)space_converter::common::calculate_dmagnitude3(v[0], v[1], v[2]);

#define RETURN_COMP_EMPTY return 0;
#define RETURN_COMP_VALUE(v) return 1;
#define RETURN_COMP_VECTOR3(v) return 3;

#define RETURN_ORIG_EMPTY RETURN_COMP_EMPTY
#define RETURN_ORIG_VALUE(v) out_value[0] = (float)v; RETURN_COMP_VALUE(v)
#define RETURN_ORIG_VECTOR3(v) out_value[0] = (float)v[0]; out_value[1] = (float)v[1]; out_value[2] = (float)v[2]; RETURN_COMP_VECTOR3(v)

namespace ipic3d {
    namespace io {

        // Data structure to hold particle data
        struct ParticleData {
            std::vector<double> x, y, z;    // Position
            std::vector<double> u, v, w;    // Velocity
            std::vector<double> q;          // Charge
            std::vector<uint64_t> id;       // Particle ID
            size_t count;                   // Number of particles
            int species_type;               // Species type (0-3)
        };

        // Data structure to hold structured-grid field/moment data, exposed as one
        // synthetic "particle" per grid point, positioned at the cell center.
        struct GridData {
            std::vector<double> px, py, pz;             // Grid point position
            std::vector<double> ex, ey, ez;             // /fields/E{x,y,z}
            std::vector<double> bx, by, bz;             // /fields/B{x,y,z}
            std::vector<double> rho;                    // /moments/species_N/rho
            std::vector<double> pxx, pxy, pxz, pyy, pyz, pzz; // /moments/species_N/p{XX,XY,XZ,YY,YZ,ZZ}
            size_t count = 0;                           // Number of grid points

            bool has_efield() const { return !ex.empty(); }
            bool has_bfield() const { return !bx.empty(); }
            bool has_rho() const { return !rho.empty(); }
            bool has_pxx() const { return !pxx.empty(); }
            bool has_pxy() const { return !pxy.empty(); }
            bool has_pxz() const { return !pxz.empty(); }
            bool has_pyy() const { return !pyy.empty(); }
            bool has_pyz() const { return !pyz.empty(); }
            bool has_pzz() const { return !pzz.empty(); }
        };

        // Grid geometry read from the companion settings.hdf file.
        struct GridSettings {
            double Dx = 1.0, Dy = 1.0, Dz = 1.0;
            double Lx = 0.0, Ly = 0.0, Lz = 0.0;
            int Nxc = 0, Nyc = 0, Nzc = 0;
            int XLEN = 1, YLEN = 1, ZLEN = 1;
            bool loaded = false;
        };

        // iPIC3D loads every species found in the file at once. Each species
        // contributes two kinds of "particles" sharing one global id space:
        //   real particles (/particles/species_N) and synthetic grid points
        //   (/fields, /moments/species_N), one per structured-grid cell.
        // IdSegment maps a contiguous range of that combined id space back to
        // (species, is_grid, local index within that species' ParticleData/GridData).
        struct IdSegment {
            int species;
            bool is_grid;
            size_t start;
            size_t count;
        };

        std::vector<ParticleData> particle_data_by_species; // indexed by species id, size PTMax
        std::vector<GridData> grid_data_by_species;          // indexed by species id, size PTMax
        std::vector<IdSegment> id_segments;
        size_t g_local_total = 0;
        GridSettings grid_settings;
        int64_t g_total_particles = 0;
        double steps_time[2];

        void print_CPU_steps() {
            // printf("iPIC3D HDF5 reading time: %f seconds\n", steps_time[1] - steps_time[0]);
        }

        // Resolve a combined id into the species/container/local-index it belongs to.
        static bool resolve_id(uint64_t id, int& species, bool& is_grid, size_t& local_idx) {
            for (const auto& seg : id_segments) {
                if (id >= seg.start && (id - seg.start) < seg.count) {
                    species = seg.species;
                    is_grid = seg.is_grid;
                    local_idx = (size_t)(id - seg.start);
                    return true;
                }
            }
            return false;
        }

        float get_particle_norm_value(int blocknr, uint64_t id) {
            int species; bool is_grid; size_t i;
            if (!resolve_id(id, species, is_grid, i)) RETURN_NORM_EMPTY;

            if (!is_grid) {
                const ParticleData& p = particle_data_by_species[species];
                switch (blocknr) {
                case IPIC3DBlockType::Pos: {
                    double pos[3] = { p.x[i], p.y[i], p.z[i] };
                    RETURN_NORM_VECTOR3(pos);
                }
                case IPIC3DBlockType::Vel: {
                    double vel[3] = { p.u[i], p.v[i], p.w[i] };
                    RETURN_NORM_VECTOR3(vel);
                }
                case IPIC3DBlockType::Charge:
                    RETURN_NORM_VALUE(p.q[i]);
                case IPIC3DBlockType::ID:
                    RETURN_NORM_VALUE(p.id[i]);
                default:
                    RETURN_NORM_EMPTY;
                }
            }

            const GridData& g = grid_data_by_species[species];
            switch (blocknr) {
            case IPIC3DBlockType::EField: {
                double v[3] = { g.ex[i], g.ey[i], g.ez[i] };
                RETURN_NORM_VECTOR3(v);
            }
            case IPIC3DBlockType::BField: {
                double v[3] = { g.bx[i], g.by[i], g.bz[i] };
                RETURN_NORM_VECTOR3(v);
            }
            case IPIC3DBlockType::Rho:
                RETURN_NORM_VALUE(g.rho[i]);
            case IPIC3DBlockType::Pxx:
                RETURN_NORM_VALUE(g.pxx[i]);
            case IPIC3DBlockType::Pxy:
                RETURN_NORM_VALUE(g.pxy[i]);
            case IPIC3DBlockType::Pxz:
                RETURN_NORM_VALUE(g.pxz[i]);
            case IPIC3DBlockType::Pyy:
                RETURN_NORM_VALUE(g.pyy[i]);
            case IPIC3DBlockType::Pyz:
                RETURN_NORM_VALUE(g.pyz[i]);
            case IPIC3DBlockType::Pzz:
                RETURN_NORM_VALUE(g.pzz[i]);
            default:
                RETURN_NORM_EMPTY;
            }
        }

        int get_particle_value(int blocknr, uint64_t id, float* out_value) {
            int species; bool is_grid; size_t i;
            if (!resolve_id(id, species, is_grid, i)) RETURN_ORIG_EMPTY;

            if (!is_grid) {
                const ParticleData& p = particle_data_by_species[species];
                switch (blocknr) {
                case IPIC3DBlockType::Pos: {
                    double pos[3] = { p.x[i], p.y[i], p.z[i] };
                    RETURN_ORIG_VECTOR3(pos);
                }
                case IPIC3DBlockType::Vel: {
                    double vel[3] = { p.u[i], p.v[i], p.w[i] };
                    RETURN_ORIG_VECTOR3(vel);
                }
                case IPIC3DBlockType::Charge:
                    RETURN_ORIG_VALUE(p.q[i]);
                case IPIC3DBlockType::ID:
                    RETURN_ORIG_VALUE(p.id[i]);
                default:
                    RETURN_ORIG_EMPTY;
                }
            }

            const GridData& g = grid_data_by_species[species];
            switch (blocknr) {
            case IPIC3DBlockType::EField: {
                double v[3] = { g.ex[i], g.ey[i], g.ez[i] };
                RETURN_ORIG_VECTOR3(v);
            }
            case IPIC3DBlockType::BField: {
                double v[3] = { g.bx[i], g.by[i], g.bz[i] };
                RETURN_ORIG_VECTOR3(v);
            }
            case IPIC3DBlockType::Rho:
                RETURN_ORIG_VALUE(g.rho[i]);
            case IPIC3DBlockType::Pxx:
                RETURN_ORIG_VALUE(g.pxx[i]);
            case IPIC3DBlockType::Pxy:
                RETURN_ORIG_VALUE(g.pxy[i]);
            case IPIC3DBlockType::Pxz:
                RETURN_ORIG_VALUE(g.pxz[i]);
            case IPIC3DBlockType::Pyy:
                RETURN_ORIG_VALUE(g.pyy[i]);
            case IPIC3DBlockType::Pyz:
                RETURN_ORIG_VALUE(g.pyz[i]);
            case IPIC3DBlockType::Pzz:
                RETURN_ORIG_VALUE(g.pzz[i]);
            default:
                RETURN_ORIG_EMPTY;
            }
        }

        int get_particle_value_comp(int blocknr, uint64_t id) {
            int species; bool is_grid; size_t i;
            if (!resolve_id(id, species, is_grid, i)) RETURN_COMP_EMPTY;

            if (!is_grid) {
                switch (blocknr) {
                case IPIC3DBlockType::Pos:
                case IPIC3DBlockType::Vel:
                    RETURN_COMP_VECTOR3(nullptr);
                case IPIC3DBlockType::Charge:
                case IPIC3DBlockType::ID:
                    RETURN_COMP_VALUE(nullptr);
                default:
                    RETURN_COMP_EMPTY;
                }
            }

            switch (blocknr) {
            case IPIC3DBlockType::EField:
            case IPIC3DBlockType::BField:
                RETURN_COMP_VECTOR3(nullptr);
            case IPIC3DBlockType::Rho:
            case IPIC3DBlockType::Pxx:
            case IPIC3DBlockType::Pxy:
            case IPIC3DBlockType::Pxz:
            case IPIC3DBlockType::Pyy:
            case IPIC3DBlockType::Pyz:
            case IPIC3DBlockType::Pzz:
                RETURN_COMP_VALUE(nullptr);
            default:
                RETURN_COMP_EMPTY;
            }
        }

        int get_particle_type(uint64_t id) {
            int species; bool is_grid; size_t i;
            if (!resolve_id(id, species, is_grid, i)) return 0;
            return species;
        }

        void get_particle_position(uint64_t id, double* pos) {
            int species; bool is_grid; size_t i;
            if (!resolve_id(id, species, is_grid, i)) return;

            if (!is_grid) {
                const ParticleData& p = particle_data_by_species[species];
                pos[0] = p.x[i];
                pos[1] = p.y[i];
                pos[2] = p.z[i];
            } else {
                const GridData& g = grid_data_by_species[species];
                pos[0] = g.px[i];
                pos[1] = g.py[i];
                pos[2] = g.pz[i];
            }
        }

        size_t get_local_num_particles() {
            return g_local_total;
        }

        size_t get_global_num_particles() {
            return g_total_particles;
        }

        double get_particle_hsml(uint64_t id) {
            // iPIC3D doesn't have smoothing length, return a default value
            return 1.0;
        }

        double get_particle_mass(uint64_t id) {
            // Use charge as a proxy for mass in iPIC3D
            int species; bool is_grid; size_t i;
            if (resolve_id(id, species, is_grid, i) && !is_grid) {
                return particle_data_by_species[species].q[i];
            }
            return 1.0;
        }

        double get_particle_rho(uint64_t id) {
            int species; bool is_grid; size_t i;
            if (resolve_id(id, species, is_grid, i) && is_grid) {
                return grid_data_by_species[species].rho[i];
            }
            return 1.0;
        }

        int get_particle_rho_blocknr() {
            return IPIC3DBlockType::Rho;
        }

        std::string get_dataset_name(int blocknr) {
            switch (blocknr) {
            case IPIC3DBlockType::Pos:
                return "Position";
            case IPIC3DBlockType::Vel:
                return "Velocity";
            case IPIC3DBlockType::Charge:
                return "Charge";
            case IPIC3DBlockType::ID:
                return "ID";
            case IPIC3DBlockType::EField:
                return "EField";
            case IPIC3DBlockType::BField:
                return "BField";
            case IPIC3DBlockType::Rho:
                return "Rho";
            case IPIC3DBlockType::Pxx:
                return "Pxx";
            case IPIC3DBlockType::Pxy:
                return "Pxy";
            case IPIC3DBlockType::Pxz:
                return "Pxz";
            case IPIC3DBlockType::Pyy:
                return "Pyy";
            case IPIC3DBlockType::Pyz:
                return "Pyz";
            case IPIC3DBlockType::Pzz:
                return "Pzz";
            default:
                return "Unknown";
            }
        }

        void get_types_and_blocks(std::vector<int>& types_and_blocks) {
            types_and_blocks.resize(IPIC3DParticleType::PTMax * IPIC3DBlockType::BTMax, 0);

            for (int s = 0; s < IPIC3DParticleType::PTMax && s < (int)particle_data_by_species.size(); s++) {
                const ParticleData& p = particle_data_by_species[s];
                const GridData& g = grid_data_by_species[s];

                if (p.count > 0) {
                    types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pos] = 1;
                    types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Vel] = 1;
                    types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Charge] = 1;
                    types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::ID] = 1;
                }

                if (g.has_efield()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::EField] = 1;
                if (g.has_bfield()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::BField] = 1;
                if (g.has_rho()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Rho] = 1;
                if (g.has_pxx()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pxx] = 1;
                if (g.has_pxy()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pxy] = 1;
                if (g.has_pxz()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pxz] = 1;
                if (g.has_pyy()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pyy] = 1;
                if (g.has_pyz()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pyz] = 1;
                if (g.has_pzz()) types_and_blocks[IPIC3DBlockType::BTMax * s + IPIC3DBlockType::Pzz] = 1;
            }
        }

        void print_types_and_blocks_local() {
            std::vector<int> types_and_blocks;
            get_types_and_blocks(types_and_blocks);

            printf("\n");
            for (int type = 0; type < IPIC3DParticleType::PTMax; type++) {
                printf("Type: Species_%d (%d)\n", type, type);
                for (int blocknr = 0; blocknr < IPIC3DBlockType::BTMax; blocknr++) {
                    if (types_and_blocks[IPIC3DBlockType::BTMax * type + blocknr] > 0) {
                        printf("\t%s (%d)\n", get_dataset_name(blocknr).c_str(), blocknr);
                    }
                }
            }
        }

        void print_types_and_blocks(std::vector<int>& types_and_blocks) {
            printf("\nAll snapshots contain:\n");
            for (int type = 0; type < IPIC3DParticleType::PTMax; type++) {
                printf("Type: Species_%d (%d)\n", type, type);
                for (int blocknr = 0; blocknr < IPIC3DBlockType::BTMax; blocknr++) {
                    if (types_and_blocks[IPIC3DBlockType::BTMax * type + blocknr] > 0) {
                        printf("\t%s (%d)\n", get_dataset_name(blocknr).c_str(), blocknr);
                    }
                }
            }
        }

        // Replace "{}" placeholders in pattern with the given file number.
        std::string format_filename(const std::string& pattern, int number) {
            std::string result = pattern;
            const std::string placeholder = "{}";
            size_t pos = result.find(placeholder);
            while (pos != std::string::npos) {
                std::string number_str = std::to_string(number);
                result.replace(pos, placeholder.length(), number_str);
                pos = result.find(placeholder, pos + number_str.length());
            }
            return result;
        }

        // H5Fopen with retry/backoff: with hundreds of MPI ranks each opening a
        // multi-GB file on a shared filesystem at the same instant, transient
        // "unable to synchronously open file" errors are expected (metadata
        // server contention), not a real missing/invalid file. Retrying a few
        // times with growing delay avoids aborting the whole job on a hiccup.
        hid_t h5fopen_retry(const std::string& file_path, int max_attempts = 8, int initial_delay_ms = 250) {
            H5E_auto2_t old_func;
            void* old_data;
            H5Eget_auto(H5E_DEFAULT, &old_func, &old_data);

            hid_t file_id = -1;
            int delay_ms = initial_delay_ms;
            for (int attempt = 0; attempt < max_attempts; attempt++) {
                H5Eset_auto(H5E_DEFAULT, nullptr, nullptr); // silence HDF5-DIAG spam while probing
                file_id = H5Fopen(file_path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
                H5Eset_auto(H5E_DEFAULT, old_func, old_data);

                if (file_id >= 0) return file_id;
                if (attempt + 1 < max_attempts) {
                    std::this_thread::sleep_for(std::chrono::milliseconds(delay_ms));
                    delay_ms = std::min(delay_ms * 2, 5000);
                }
            }
            return file_id;
        }

        // Split `total` items into `parts` contiguous chunks, returning the start
        // index and size of chunk `idx` (remainder distributed to the first chunks).
        void partition_range(size_t total, size_t parts, size_t idx, size_t& start, size_t& count) {
            size_t base = total / parts;
            size_t rem = total % parts;
            count = base + (idx < rem ? 1 : 0);
            start = idx * base + std::min(idx, rem);
        }

        // Append a [start, start + count) slice of src's fields onto dst's fields.
        void append_particle_range(const ParticleData& src, ParticleData& dst, size_t start, size_t count) {
            dst.x.insert(dst.x.end(), src.x.begin() + start, src.x.begin() + start + count);
            dst.y.insert(dst.y.end(), src.y.begin() + start, src.y.begin() + start + count);
            dst.z.insert(dst.z.end(), src.z.begin() + start, src.z.begin() + start + count);
            dst.u.insert(dst.u.end(), src.u.begin() + start, src.u.begin() + start + count);
            dst.v.insert(dst.v.end(), src.v.begin() + start, src.v.begin() + start + count);
            dst.w.insert(dst.w.end(), src.w.begin() + start, src.w.begin() + start + count);
            dst.q.insert(dst.q.end(), src.q.begin() + start, src.q.begin() + start + count);
            dst.id.insert(dst.id.end(), src.id.begin() + start, src.id.begin() + start + count);
        }

        // Helper function to read a dataset from HDF5
        template<typename T>
        void read_hdf5_dataset(hid_t file_id, const std::string& dataset_path, std::vector<T>& data) {
            hid_t dataset_id = H5Dopen(file_id, dataset_path.c_str(), H5P_DEFAULT);
            if (dataset_id < 0) {
                throw std::runtime_error("Failed to open dataset: " + dataset_path);
            }

            hid_t dataspace_id = H5Dget_space(dataset_id);
            hssize_t num_elements = H5Sget_simple_extent_npoints(dataspace_id);
            
            data.resize(num_elements);
            
            hid_t memtype;
            if (std::is_same<T, double>::value) {
                memtype = H5T_NATIVE_DOUBLE;
            } else if (std::is_same<T, float>::value) {
                memtype = H5T_NATIVE_FLOAT;
            } else if (std::is_same<T, uint64_t>::value) {
                memtype = H5T_NATIVE_UINT64;
            } else if (std::is_same<T, int64_t>::value) {
                memtype = H5T_NATIVE_INT64;
            } else {
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                throw std::runtime_error("Unsupported data type for HDF5 reading");
            }
            
            herr_t status = H5Dread(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
            
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            
            if (status < 0) {
                throw std::runtime_error("Failed to read dataset: " + dataset_path);
            }
        }

        // Open one HDF5 file (file_index substituted into the "{}" placeholder of
        // hdf5_file_pattern) and read all particle datasets for species/cycle into out.
        void read_particle_file(const std::string& hdf5_file_pattern, int file_index,
                                 const std::string& species_path, const std::string& cycle_suffix,
                                 int species_id, ParticleData& out) {
            std::string file_path = format_filename(hdf5_file_pattern, file_index);
            hid_t file_id = h5fopen_retry(file_path);
            if (file_id < 0) {
                throw std::runtime_error("Failed to open HDF5 file: " + file_path);
            }

            try {
                out.species_type = species_id;
                read_hdf5_dataset(file_id, species_path + "/x" + cycle_suffix, out.x);
                read_hdf5_dataset(file_id, species_path + "/y" + cycle_suffix, out.y);
                read_hdf5_dataset(file_id, species_path + "/z" + cycle_suffix, out.z);
                read_hdf5_dataset(file_id, species_path + "/u" + cycle_suffix, out.u);
                read_hdf5_dataset(file_id, species_path + "/v" + cycle_suffix, out.v);
                read_hdf5_dataset(file_id, species_path + "/w" + cycle_suffix, out.w);
                read_hdf5_dataset(file_id, species_path + "/q" + cycle_suffix, out.q);
                read_hdf5_dataset(file_id, species_path + "/ID" + cycle_suffix, out.id);
                out.count = out.x.size();
            } catch (const std::exception& e) {
                H5Fclose(file_id);
                throw;
            }

            H5Fclose(file_id);
        }

        // Like read_hdf5_dataset, but returns false instead of throwing when the dataset
        // does not exist. Used for /fields and /moments datasets, which are optional.
        template<typename T>
        bool try_read_hdf5_dataset(hid_t file_id, const std::string& dataset_path, std::vector<T>& data) {
            if (H5Lexists(file_id, dataset_path.c_str(), H5P_DEFAULT) <= 0) {
                return false;
            }
            try {
                read_hdf5_dataset(file_id, dataset_path, data);
            } catch (const std::exception&) {
                return false;
            }
            return true;
        }

        // Read a single scalar value from a 0/1-element HDF5 dataset (used for settings.hdf
        // /collective and /topology entries, and per-file /topology/cartesian_coord).
        template<typename T>
        bool try_read_hdf5_scalar(hid_t file_id, const std::string& dataset_path, T& out) {
            std::vector<T> data;
            if (!try_read_hdf5_dataset(file_id, dataset_path, data) || data.empty()) {
                return false;
            }
            out = data[0];
            return true;
        }

        // Read grid geometry from the companion settings.hdf file. Missing file/datasets
        // fall back to the GridSettings defaults (unit spacing, single-rank topology).
        GridSettings read_settings(const std::string& settings_path) {
            GridSettings s;
            hid_t file_id = h5fopen_retry(settings_path);
            if (file_id < 0) {
                printf("Warning: could not open iPIC3D settings file '%s', using unit grid spacing\n", settings_path.c_str());
                return s;
            }

            try_read_hdf5_scalar(file_id, "/collective/Dx", s.Dx);
            try_read_hdf5_scalar(file_id, "/collective/Dy", s.Dy);
            try_read_hdf5_scalar(file_id, "/collective/Dz", s.Dz);
            try_read_hdf5_scalar(file_id, "/collective/Lx", s.Lx);
            try_read_hdf5_scalar(file_id, "/collective/Ly", s.Ly);
            try_read_hdf5_scalar(file_id, "/collective/Lz", s.Lz);

            double tmp;
            if (try_read_hdf5_scalar(file_id, "/collective/Nxc", tmp)) s.Nxc = (int)tmp;
            if (try_read_hdf5_scalar(file_id, "/collective/Nyc", tmp)) s.Nyc = (int)tmp;
            if (try_read_hdf5_scalar(file_id, "/collective/Nzc", tmp)) s.Nzc = (int)tmp;
            if (try_read_hdf5_scalar(file_id, "/topology/XLEN", tmp)) s.XLEN = (int)tmp;
            if (try_read_hdf5_scalar(file_id, "/topology/YLEN", tmp)) s.YLEN = (int)tmp;
            if (try_read_hdf5_scalar(file_id, "/topology/ZLEN", tmp)) s.ZLEN = (int)tmp;

            s.loaded = true;
            H5Fclose(file_id);
            return s;
        }

        // Read this file's MPI-rank coordinate within the iPIC3D cartesian topology
        // (used to place its local field/moment grid tile within the global domain).
        void read_topology_coord(hid_t file_id, int coord[3]) {
            coord[0] = coord[1] = coord[2] = 0;
            std::vector<int> data;
            if (try_read_hdf5_dataset(file_id, "/topology/cartesian_coord", data) && data.size() >= 3) {
                coord[0] = data[0];
                coord[1] = data[1];
                coord[2] = data[2];
            }
        }

        // Read one 3D grid dataset (e.g. /fields/Bx/cycle_N) and report its dimensions.
        bool read_grid_dataset(hid_t file_id, const std::string& dataset_path, std::vector<double>& data, hsize_t dims[3]) {
            if (H5Lexists(file_id, dataset_path.c_str(), H5P_DEFAULT) <= 0) {
                return false;
            }

            hid_t dataset_id = H5Dopen(file_id, dataset_path.c_str(), H5P_DEFAULT);
            if (dataset_id < 0) return false;

            hid_t dataspace_id = H5Dget_space(dataset_id);
            int ndims = H5Sget_simple_extent_ndims(dataspace_id);
            hsize_t file_dims[3] = { 1, 1, 1 };
            if (ndims > 0 && ndims <= 3) {
                H5Sget_simple_extent_dims(dataspace_id, file_dims, nullptr);
                // Datasets with fewer than 3 dims are placed in the trailing axes.
                for (int i = 0; i < ndims; i++) {
                    dims[3 - ndims + i] = file_dims[i];
                }
                for (int i = 0; i < 3 - ndims; i++) {
                    dims[i] = 1;
                }
            } else {
                H5Sclose(dataspace_id);
                H5Dclose(dataset_id);
                return false;
            }

            hssize_t num_elements = H5Sget_simple_extent_npoints(dataspace_id);
            data.resize(num_elements);
            herr_t status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            return status >= 0;
        }

        // Compute the world-space position of every point in a [nx,ny,nz] local grid tile,
        // given this tile's coordinate in the cartesian process topology.
        // NOTE: assumes iPIC3D's standard 1-cell ghost layer on each side of every local
        // field/moment tile; adjust GHOST if a given dataset's layout differs.
        void compute_grid_positions(const hsize_t dims[3], const int coord[3], const GridSettings& s,
                                     std::vector<double>& px, std::vector<double>& py, std::vector<double>& pz) {
            const int64_t GHOST = 1;
            size_t n = (size_t)dims[0] * dims[1] * dims[2];
            px.resize(n); py.resize(n); pz.resize(n);

            int64_t interior_x = (int64_t)dims[0] - 2 * GHOST;
            int64_t interior_y = (int64_t)dims[1] - 2 * GHOST;
            int64_t interior_z = (int64_t)dims[2] - 2 * GHOST;
            if (interior_x < 1) interior_x = (int64_t)dims[0];
            if (interior_y < 1) interior_y = (int64_t)dims[1];
            if (interior_z < 1) interior_z = (int64_t)dims[2];

            size_t idx = 0;
            for (hsize_t i = 0; i < dims[0]; i++) {
                double gx = (coord[0] * interior_x + ((int64_t)i - GHOST)) * s.Dx;
                for (hsize_t j = 0; j < dims[1]; j++) {
                    double gy = (coord[1] * interior_y + ((int64_t)j - GHOST)) * s.Dy;
                    for (hsize_t k = 0; k < dims[2]; k++) {
                        double gz = (coord[2] * interior_z + ((int64_t)k - GHOST)) * s.Dz;
                        px[idx] = gx;
                        py[idx] = gy;
                        pz[idx] = gz;
                        idx++;
                    }
                }
            }
        }

        // Open one HDF5 file and read all /fields and /moments/species_N datasets for the
        // given cycle into a flattened GridData (one entry per grid point). Missing
        // datasets are simply left empty (reported as unavailable by get_types_and_blocks).
        void read_grid_file(const std::string& hdf5_file_pattern, int file_index,
                             int species_id, const std::string& cycle_suffix,
                             const GridSettings& settings, GridData& out) {
            std::string file_path = format_filename(hdf5_file_pattern, file_index);
            hid_t file_id = h5fopen_retry(file_path);
            if (file_id < 0) {
                throw std::runtime_error("Failed to open HDF5 file: " + file_path);
            }

            try {
                int coord[3];
                read_topology_coord(file_id, coord);

                hsize_t dims[3] = { 0, 0, 0 };
                bool have_dims = false;
                std::string moments_path = "/moments/species_" + std::to_string(species_id);

                auto read_field = [&](const char* name, std::vector<double>& dst) {
                    hsize_t local_dims[3];
                    if (read_grid_dataset(file_id, std::string("/fields/") + name + cycle_suffix, dst, local_dims)) {
                        if (!have_dims) { dims[0] = local_dims[0]; dims[1] = local_dims[1]; dims[2] = local_dims[2]; have_dims = true; }
                    }
                };
                read_field("Ex", out.ex); read_field("Ey", out.ey); read_field("Ez", out.ez);
                read_field("Bx", out.bx); read_field("By", out.by); read_field("Bz", out.bz);

                auto read_moment = [&](const char* name, std::vector<double>& dst) {
                    hsize_t local_dims[3];
                    if (read_grid_dataset(file_id, moments_path + "/" + name + cycle_suffix, dst, local_dims)) {
                        if (!have_dims) { dims[0] = local_dims[0]; dims[1] = local_dims[1]; dims[2] = local_dims[2]; have_dims = true; }
                    }
                };
                read_moment("rho", out.rho);
                read_moment("pXX", out.pxx); read_moment("pXY", out.pxy); read_moment("pXZ", out.pxz);
                read_moment("pYY", out.pyy); read_moment("pYZ", out.pyz); read_moment("pZZ", out.pzz);

                if (have_dims) {
                    compute_grid_positions(dims, coord, settings, out.px, out.py, out.pz);
                    out.count = (size_t)dims[0] * dims[1] * dims[2];
                } else {
                    out.count = 0;
                }
            } catch (const std::exception& e) {
                H5Fclose(file_id);
                throw;
            }

            H5Fclose(file_id);
        }

        // Append all of src's grid points onto dst (only when both have data and the
        // same per-point fields available, since grid tiles from different files of the
        // same simulation always carry the same set of fields/moments).
        void append_grid_data(const GridData& src, GridData& dst) {
            if (src.count == 0) return;
            dst.px.insert(dst.px.end(), src.px.begin(), src.px.end());
            dst.py.insert(dst.py.end(), src.py.begin(), src.py.end());
            dst.pz.insert(dst.pz.end(), src.pz.begin(), src.pz.end());
            dst.ex.insert(dst.ex.end(), src.ex.begin(), src.ex.end());
            dst.ey.insert(dst.ey.end(), src.ey.begin(), src.ey.end());
            dst.ez.insert(dst.ez.end(), src.ez.begin(), src.ez.end());
            dst.bx.insert(dst.bx.end(), src.bx.begin(), src.bx.end());
            dst.by.insert(dst.by.end(), src.by.begin(), src.by.end());
            dst.bz.insert(dst.bz.end(), src.bz.begin(), src.bz.end());
            dst.rho.insert(dst.rho.end(), src.rho.begin(), src.rho.end());
            dst.pxx.insert(dst.pxx.end(), src.pxx.begin(), src.pxx.end());
            dst.pxy.insert(dst.pxy.end(), src.pxy.begin(), src.pxy.end());
            dst.pxz.insert(dst.pxz.end(), src.pxz.begin(), src.pxz.end());
            dst.pyy.insert(dst.pyy.end(), src.pyy.begin(), src.pyy.end());
            dst.pyz.insert(dst.pyz.end(), src.pyz.begin(), src.pyz.end());
            dst.pzz.insert(dst.pzz.end(), src.pzz.begin(), src.pzz.end());
            dst.count += src.count;
        }

        // If settings_file is empty, derive "settings.hdf" living next to hdf5_file.
        std::string default_settings_path(const std::string& hdf5_file) {
            size_t slash = hdf5_file.find_last_of("/\\");
            std::string dir = (slash == std::string::npos) ? "" : hdf5_file.substr(0, slash + 1);
            return dir + "settings.hdf";
        }

        // Find every species_N present under /particles or /moments in this file.
        std::vector<int> discover_species(hid_t file_id) {
            std::vector<int> species;
            for (int s = 0; s < IPIC3DParticleType::PTMax; s++) {
                std::string p = "/particles/species_" + std::to_string(s);
                std::string m = "/moments/species_" + std::to_string(s);
                if (H5Lexists(file_id, p.c_str(), H5P_DEFAULT) > 0 || H5Lexists(file_id, m.c_str(), H5P_DEFAULT) > 0) {
                    species.push_back(s);
                }
            }
            return species;
        }

        // Find the highest "cycle_N" dataset name under the given group (e.g.
        // "/particles/species_0/x" or "/fields/Bx"). Returns -1 if the group is missing
        // or has no cycle_* entries.
        int discover_last_cycle(hid_t file_id, const std::string& group_path) {
            if (H5Lexists(file_id, group_path.c_str(), H5P_DEFAULT) <= 0) return -1;

            hid_t group_id = H5Gopen(file_id, group_path.c_str(), H5P_DEFAULT);
            if (group_id < 0) return -1;

            int max_cycle = -1;
            H5G_info_t info;
            if (H5Gget_info(group_id, &info) >= 0) {
                const std::string prefix = "cycle_";
                for (hsize_t i = 0; i < info.nlinks; i++) {
                    char name_buf[256];
                    ssize_t len = H5Lget_name_by_idx(group_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, name_buf, sizeof(name_buf), H5P_DEFAULT);
                    if (len > 0 && (size_t)len < sizeof(name_buf)) {
                        std::string name(name_buf, (size_t)len);
                        if (name.rfind(prefix, 0) == 0) {
                            int c = std::atoi(name.c_str() + prefix.size());
                            if (c > max_cycle) max_cycle = c;
                        }
                    }
                }
            }
            H5Gclose(group_id);
            return max_cycle;
        }

        // Probe one file to find which species are present and the most recent cycle
        // written for them, so the caller can load everything without being told.
        void discover_species_and_cycle(const std::string& hdf5_file, int file_index,
                                         std::vector<int>& species, int& cycle_id) {
            std::string file_path = format_filename(hdf5_file, file_index);
            hid_t file_id = h5fopen_retry(file_path);
            if (file_id < 0) {
                throw std::runtime_error("Failed to open HDF5 file: " + file_path);
            }

            species = discover_species(file_id);

            int cycle = -1;
            for (int s : species) {
                cycle = std::max(cycle, discover_last_cycle(file_id, "/particles/species_" + std::to_string(s) + "/x"));
                cycle = std::max(cycle, discover_last_cycle(file_id, "/moments/species_" + std::to_string(s) + "/rho"));
            }
            cycle = std::max(cycle, discover_last_cycle(file_id, "/fields/Bx"));
            cycle_id = (cycle < 0) ? 0 : cycle;

            H5Fclose(file_id);
        }

        void init_lib(std::string hdf5_file, int world_rank, int world_size,
                      int num_files, std::string settings_file) {
#ifdef WITH_OPENMP
            steps_time[0] = omp_get_wtime();
#endif
            std::string settings_path = settings_file.empty() ? default_settings_path(hdf5_file) : settings_file;
            GridSettings settings = read_settings(settings_path);

            // Discover available species and the latest cycle from the file(s) this
            // rank will read; with a single restart dataset every file shares the
            // same species/cycle layout, so probing this rank's first file suffices.
            size_t probe_file_index = 0;
            if (num_files >= world_size) {
                size_t file_start, file_count;
                partition_range((size_t)num_files, (size_t)world_size, (size_t)world_rank, file_start, file_count);
                probe_file_index = file_start;
            } else {
                for (int f = 0; f < num_files; f++) {
                    size_t rank_start, rank_count;
                    partition_range((size_t)world_size, (size_t)num_files, (size_t)f, rank_start, rank_count);
                    if ((size_t)world_rank >= rank_start && (size_t)world_rank < rank_start + rank_count) {
                        probe_file_index = (size_t)f;
                        break;
                    }
                }
            }

            std::vector<int> species_list;
            int cycle_id;
            discover_species_and_cycle(hdf5_file, (int)probe_file_index, species_list, cycle_id);
            std::string cycle_suffix = "/cycle_" + std::to_string(cycle_id);

            if (world_rank == 0) {
                printf("Reading iPIC3D HDF5 file(s): %s\n", hdf5_file.c_str());
                printf("num_files: %d, world_size: %d, cycle: %d, species: [", num_files, world_size, cycle_id);
                for (size_t i = 0; i < species_list.size(); i++) printf("%s%d", i ? "," : "", species_list[i]);
                printf("]\n");
                printf("Settings file: %s (%s)\n", settings_path.c_str(), settings.loaded ? "loaded" : "not found, using defaults");
            }

            particle_data_by_species.assign(IPIC3DParticleType::PTMax, ParticleData());
            grid_data_by_species.assign(IPIC3DParticleType::PTMax, GridData());
            for (int s = 0; s < IPIC3DParticleType::PTMax; s++) {
                particle_data_by_species[s].species_type = s;
                particle_data_by_species[s].count = 0;
            }

            for (int species_id : species_list) {
                std::string species_path = "/particles/species_" + std::to_string(species_id);
                ParticleData& result = particle_data_by_species[species_id];
                GridData& grid_result = grid_data_by_species[species_id];

                if (num_files >= world_size) {
                    // Each rank owns one or more whole files exclusively; concatenate them.
                    size_t file_start, file_count;
                    partition_range((size_t)num_files, (size_t)world_size, (size_t)world_rank, file_start, file_count);

                    for (size_t i = 0; i < file_count; i++) {
                        ParticleData file_data;
                        read_particle_file(hdf5_file, (int)(file_start + i), species_path, cycle_suffix, species_id, file_data);
                        append_particle_range(file_data, result, 0, file_data.count);

                        GridData file_grid;
                        read_grid_file(hdf5_file, (int)(file_start + i), species_id, cycle_suffix, settings, file_grid);
                        append_grid_data(file_grid, grid_result);
                    }
                } else {
                    // Fewer files than ranks: a group of ranks shares each file and
                    // splits its particles among the group.
                    size_t rank_start = 0, rank_count = 0;
                    int file_index = -1;
                    size_t local_index = 0, group_size = 1;

                    for (int f = 0; f < num_files; f++) {
                        partition_range((size_t)world_size, (size_t)num_files, (size_t)f, rank_start, rank_count);
                        if ((size_t)world_rank >= rank_start && (size_t)world_rank < rank_start + rank_count) {
                            file_index = f;
                            local_index = (size_t)world_rank - rank_start;
                            group_size = rank_count;
                            break;
                        }
                    }

                    ParticleData file_data;
                    read_particle_file(hdf5_file, file_index, species_path, cycle_suffix, species_id, file_data);

                    size_t local_start, local_count;
                    partition_range(file_data.count, group_size, local_index, local_start, local_count);
                    append_particle_range(file_data, result, local_start, local_count);

                    // The grid (field/moment) tile belongs to the whole file, not to a
                    // particle sub-range, so only the first rank in the sharing group
                    // loads it once to avoid duplicating the same grid points.
                    if (local_index == 0) {
                        GridData file_grid;
                        read_grid_file(hdf5_file, file_index, species_id, cycle_suffix, settings, file_grid);
                        append_grid_data(file_grid, grid_result);
                    }
                }

                result.count = result.x.size();
            }

            // Build the combined id space: for each species, its particles followed
            // by its grid points (see IdSegment / resolve_id).
            id_segments.clear();
            size_t offset = 0;
            for (int s = 0; s < IPIC3DParticleType::PTMax; s++) {
                size_t pcount = particle_data_by_species[s].count;
                if (pcount > 0) {
                    id_segments.push_back({ s, false, offset, pcount });
                    offset += pcount;
                }
                size_t gcount = grid_data_by_species[s].count;
                if (gcount > 0) {
                    id_segments.push_back({ s, true, offset, gcount });
                    offset += gcount;
                }
            }
            g_local_total = offset;
            g_total_particles = (int64_t)g_local_total;

            // for (int s : species_list) {
            //     printf("Rank %d species %d: %lld particles, %lld grid points\n", world_rank, s,
            //            (long long)particle_data_by_species[s].count, (long long)grid_data_by_species[s].count);
            // }

#ifdef WITH_OPENMP
            steps_time[1] = omp_get_wtime();
#endif
        }

        void finish_lib() {
            particle_data_by_species.clear();
            grid_data_by_species.clear();
            id_segments.clear();
            g_local_total = 0;
        }

    } // namespace io
} // namespace ipic3d
