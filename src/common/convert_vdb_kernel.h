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

#pragma once

#include <float.h>
#include <vector>
#include <cmath>
#include "data_common.h"
#include "data_cache.h"
#include "../utility/dense_utility.h"

namespace space_converter {
	namespace common {
		namespace vdb {

			class VoxelSparseManager;
			class VDBParticles;
			class VoxelDenseManager;

			namespace kernel {

#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif
					inline bool near_zero(float x)
				{
					return fabsf(x) <= FDATA_EPSILON;
				}

#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif
					inline bool near_zero_norm3(float x, float y, float z)
				{
					return x * x + y * y + z * z <= 1e-12f;
				}

#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif
					inline int ifloor(double v)
				{
					int i = static_cast<int>(v);
					return (static_cast<double>(i) > v) ? i - 1 : i;
				}

#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif
					inline int iceil(double v)
				{
					int i = static_cast<int>(v);
					return (static_cast<double>(i) < v) ? i + 1 : i;
				}

				/*
				 * ── Coordinate spaces used by the conversion kernels ─────────────────────
				 *
				 *  world space         raw particle positions (reader units, offset_position
				 *                      already subtracted where applicable)
				 *  normalized space    p_norm = (pos - bbox_min_orig) / bbox_size_orig,
				 *                      in [0,1] over the cube-symmetrized dataset bbox
				 *  object space        p_norm * object_size (the bspace client works here;
				 *                      bbox_min/bbox_max zoom boxes are in this space)
				 *  voxel space         p_vox = p_norm * bbox_dim / scale_space_diagonal
				 *
				 *  scale_space_diagonal = |bbox_max - bbox_min| / (object_size * sqrt(3))
				 *                       = zoom factor (1 for the full box, < 1 zoomed in)
				 *  len2pix              = bbox_dim / (scale_space_diagonal * bbox_size_orig)
				 *                       = world -> voxel conversion factor
				 *  transform_scale      = scale_space_diagonal * object_size / bbox_dim
				 *                       = voxel edge length in object-space units
				 *
				 *  Dense grids are a bbox_dim^3 window at offset[] = bbox_min_norm *
				 *  bbox_dim / scale_space_diagonal inside the full voxel space. Sparse
				 *  grids store absolute voxel-space coordinates. SPH kernel weights W are
				 *  normalized in VOXEL units (h passed as voxels), so the sum of W over
				 *  voxels approximates 1 per particle.
				 */

#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif
					inline double get_particle_radius(
						uint64_t pid,
						int bbox_dim,
						int* bbox_min_orig,
						double bbox_size_orig,
						double scale_space_diagonal,
						float particle_radius_multiplier,
						const float* radius_particles,
						double particle_radius_const
					) {
					// Fixed radius is specified directly in voxel units
					if (particle_radius_const > 0.0) {
						return particle_radius_const;
					}

					// Conversion factor from world space to voxel space
					double norm_fac = (double)bbox_dim / scale_space_diagonal;

					double radius = 0.0;

					// Per-particle radius (hsml from the reader, or the k-NN estimate)
					if (radius_particles != nullptr) {
						radius = radius_particles[pid];
					}

					// Guard against non-finite radii (the k-NN search stores infinity for
					// particles with fewer than k neighbors); an infinite radius would make
					// fill_voxels walk the whole grid
					if (!(radius >= 0.0) || radius > 1e30) {
						radius = 0.0;
					}

					if (particle_radius_multiplier != 0.0f) {
						radius *= particle_radius_multiplier;
					}

					// Convert radius to voxel units
					double len_to_pix = norm_fac / bbox_size_orig;
					double radiusxyz_max = radius * len_to_pix;

					return radiusxyz_max;
				}

#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif            
					inline void fill_voxels(
						//common::vdb::VoxelDenseManager& grid,
						float* data_density,
						float* data_temp,
						size_t* offset,
						size_t* dims,

						size_t pid,
						float value,
						int bbox_dim,
						int* bbox_min_orig,
						double bbox_size_orig,
						double scale_space_diagonal,
						common::SpaceData::DenseType dense_type,
						common::SpaceData::DenseNorm dense_norm,
						float particle_radius_multiplier,
						const float* radius_particles,
						double particle_radius_const,
						int block_name_id,
						double* pos
					) {

					double norm_fac = (double)bbox_dim / scale_space_diagonal;
					double len2pix = norm_fac / bbox_size_orig;

					if (block_name_id == 0) { //only pos
						value = 1.0f;
					}

					double hsml = get_particle_radius(
						pid,
						bbox_dim,
						bbox_min_orig,
						bbox_size_orig,
						scale_space_diagonal,
						particle_radius_multiplier,
						radius_particles,
						particle_radius_const
					);

					int iradius = static_cast<int>(hsml);

					///////////////////////////////////////////////////////////////////////////////
					double dpx = ((double)pos[0] - (double)bbox_min_orig[0]) * len2pix;
					double dpy = ((double)pos[1] - (double)bbox_min_orig[1]) * len2pix;
					double dpz = ((double)pos[2] - (double)bbox_min_orig[2]) * len2pix;

					int px = ifloor(dpx);
					int py = ifloor(dpy);
					int pz = ifloor(dpz);

					///////////////////////////////////////////////////////////////////////////////

					const size_t dim0 = dims[0];
					const size_t dim1 = dims[1];
					const size_t dim2 = dims[2];
					const size_t slice_size = dim0 * dim1;

					const int off0 = offset[0];
					const int off1 = offset[1];
					const int off2 = offset[2];

					const bool has_data_temp = (data_temp != nullptr);
					const bool add_count = has_data_temp && (dense_norm == common::SpaceData::DenseNorm::eCount);
					const bool add_sph = has_data_temp && (dense_norm == common::SpaceData::DenseNorm::eSPHInterpolation);

					// Fast path: kernel support below one voxel -> deposit the whole value into a single voxel (W = 1)
					if (iradius == 0) {
						int osx = px - off0;
						int osy = py - off1;
						int osz = pz - off2;
						if (osx < 0 || osx >= (int)dim0 || osy < 0 || osy >= (int)dim1 || osz < 0 || osz >= (int)dim2)
							return;

						size_t gindex = (size_t)osz * slice_size + (size_t)osy * dim0 + (size_t)osx;

						if (!near_zero(value)) {
#ifdef SPACE_GPU_DEVICE_COMPILE
							atomicAdd(data_density + gindex, value);
#else
#pragma omp atomic
							data_density[gindex] += value;
#endif
						}

						if (add_count || add_sph) {
							float n = 1.0f;
#ifdef SPACE_GPU_DEVICE_COMPILE
							atomicAdd(data_temp + gindex, n);
#else
#pragma omp atomic
							data_temp[gindex] += n;
#endif
						}
						return;
					}

					const double h = hsml;
					const double h_inv = 1.0 / h;
					const double h2 = h * h;

					// Kernel normalization is constant per particle (contains pow(h_inv, 3)) -> hoist it out of the voxel loops
					const double W_norm = utility::dense::sph_kernel::W_norm(dense_type, h_inv);

					// Iterate only the intersection of the kernel support sphere with the allocated grid window.
					// For zoomed-in extractions (small bbox) the support can span thousands of voxels while the
					// particle sits mostly outside the grid - clamping the ranges avoids walking all of that.
					int sz_lo = ifloor(dpz - h); if (sz_lo < off2) sz_lo = off2;
					int sz_hi = iceil(dpz + h); if (sz_hi > off2 + (int)dim2 - 1) sz_hi = off2 + (int)dim2 - 1;
					int sy_lo = ifloor(dpy - h); if (sy_lo < off1) sy_lo = off1;
					int sy_hi = iceil(dpy + h); if (sy_hi > off1 + (int)dim1 - 1) sy_hi = off1 + (int)dim1 - 1;
					const int sx_min = off0;
					const int sx_max = off0 + (int)dim0 - 1;

					for (int sz = sz_lo; sz <= sz_hi; sz++) {
						double dz = dpz - sz;
						double dz2 = dz * dz;
						if (dz2 >= h2) continue;

						size_t z_offset = (size_t)(sz - off2) * slice_size;

						for (int sy = sy_lo; sy <= sy_hi; sy++) {
							double dy = dpy - sy;
							double dy2 = dy * dy;
							double rem = h2 - dz2 - dy2;
							if (rem <= 0.0) continue; // whole row outside kernel support (W = 0)

							// Limit x to the kernel support: |dpx - sx| < sqrt(rem)
#ifdef SPACE_GPU_DEVICE_COMPILE
							double xr = sqrt(rem);
#else
							double xr = std::sqrt(rem);
#endif
							int sx_lo = iceil(dpx - xr); if (sx_lo < sx_min) sx_lo = sx_min;
							int sx_hi = ifloor(dpx + xr); if (sx_hi > sx_max) sx_hi = sx_max;

							size_t yz_offset = z_offset + (size_t)(sy - off1) * dim0;
							double dyz2 = dy2 + dz2;

							for (int sx = sx_lo; sx <= sx_hi; sx++) {
								double dx = dpx - sx;

#ifdef SPACE_GPU_DEVICE_COMPILE
								double distance_norm = sqrt(dx * dx + dyz2);
#else
								double distance_norm = std::sqrt(dx * dx + dyz2);
#endif

								double q = distance_norm * h_inv;
								double W = W_norm * utility::dense::sph_kernel::W_value(dense_type, q);

								float d = static_cast<float>(W * value);

								size_t gindex = yz_offset + (size_t)(sx - off0);

								if (!near_zero(d)) {
#ifdef SPACE_GPU_DEVICE_COMPILE
									atomicAdd(data_density + gindex, d);
#else
#pragma omp atomic
									data_density[gindex] += d;
#endif
								}

								if (add_count || add_sph) {
									float n = add_count ? 1.0f : static_cast<float>(W);
									if (!near_zero(n)) {
#ifdef SPACE_GPU_DEVICE_COMPILE
										atomicAdd(data_temp + gindex, n);
#else
#pragma omp atomic
										data_temp[gindex] += n;
#endif
									}
								}
							}
						}
					}
				}

				/**
				 * @brief Process a single particle for grid conversion
				 *
				 * This function can be called from both CPU and CUDA kernels.
				 * Returns true if particle should be processed, false if it should be skipped.
				 */
#if defined(__CUDACC__) || defined(__HIPCC__)
				__host__ __device__
#endif
					inline bool process_particle(
						size_t cached_idx,
						const float* pos_particles,
						const float* radius_particles,
						const float* value_particles,
						const float* offset_position,
						int* bbox_min_orig,
						double bbox_size_orig,
						int bbox_dim,
						double scale_space_diagonal,
						float particle_radius_multiplier,
						double bbox_x_min_norm,
						double bbox_x_max_norm,
						double bbox_y_min_norm,
						double bbox_y_max_norm,
						double bbox_z_min_norm,
						double bbox_z_max_norm,
						float filter_min,
						float filter_max,
						const float* bbox_sphere_pos,
						float bbox_sphere_r,
						common::SpaceData::AnimType anim_type,
						int frame_req,
						int frame,
						bool use_simple_density,
						double particle_radius_const,

						bool is_dense_grid,
						// Output parameters
						double* out_Pos,
						double& out_px_norm,
						double& out_py_norm,
						double& out_pz_norm,
						int& out_px,
						int& out_py,
						int& out_pz,
						float& out_v_orig
					) {
					// Filter by animation frame if extracting specific frame
					if (anim_type == common::SpaceData::AnimType::eFrameExtract) {
						if (frame_req != frame)
							return false;
					}

					// Get position from cached data
					out_Pos[0] = pos_particles[cached_idx * 3 + 0];
					out_Pos[1] = pos_particles[cached_idx * 3 + 1];
					out_Pos[2] = pos_particles[cached_idx * 3 + 2];

					// Apply position offset if specified
					if (offset_position[0] != 0.0f || offset_position[1] != 0.0f || offset_position[2] != 0.0f) {
						out_Pos[0] -= offset_position[0];
						out_Pos[1] -= offset_position[1];
						out_Pos[2] -= offset_position[2];
					}

					// Normalize position to [0,1] range based on bounding box
					out_px_norm = ((double)out_Pos[0] - (double)bbox_min_orig[0]) / bbox_size_orig;
					out_py_norm = ((double)out_Pos[1] - (double)bbox_min_orig[1]) / bbox_size_orig;
					out_pz_norm = ((double)out_Pos[2] - (double)bbox_min_orig[2]) / bbox_size_orig;

					// Check if particle is within bounding box
					// Dense grids need extended bounds (particle radius) for proper splatting
					if (is_dense_grid) {
						// Radius in voxel units, converted to NORMALIZED units for this test
						// (a previous version compared the voxel-unit radius against the
						// normalized coordinates directly, inflating the tolerance by a
						// factor of bbox_dim / scale_space_diagonal)
						double radius_voxels = get_particle_radius(
							cached_idx,
							bbox_dim,
							bbox_min_orig,
							bbox_size_orig,
							scale_space_diagonal,
							particle_radius_multiplier,
							radius_particles,
							particle_radius_const
						);
						double radius_norm = radius_voxels * scale_space_diagonal / (double)bbox_dim;

						if (out_px_norm + radius_norm < bbox_x_min_norm || out_px_norm - radius_norm > bbox_x_max_norm)
							return false;

						if (out_py_norm + radius_norm < bbox_y_min_norm || out_py_norm - radius_norm > bbox_y_max_norm)
							return false;

						if (out_pz_norm + radius_norm < bbox_z_min_norm || out_pz_norm - radius_norm > bbox_z_max_norm)
							return false;
					}
					else {
						// Simple point-in-box test for sparse grids
						if (out_px_norm < bbox_x_min_norm || out_px_norm > bbox_x_max_norm)
							return false;

						if (out_py_norm < bbox_y_min_norm || out_py_norm > bbox_y_max_norm)
							return false;

						if (out_pz_norm < bbox_z_min_norm || out_pz_norm > bbox_z_max_norm)
							return false;
					}

					// Convert normalized position to voxel coordinates (floor, not
					// truncation: truncation would fold the two voxels around zero into
					// voxel 0 when the zoom box crosses the origin)
					out_px = ifloor((double)out_px_norm * (double)bbox_dim / scale_space_diagonal);
					out_py = ifloor((double)out_py_norm * (double)bbox_dim / scale_space_diagonal);
					out_pz = ifloor((double)out_pz_norm * (double)bbox_dim / scale_space_diagonal);

					// Get particle value (density, temperature, etc.)
					out_v_orig = value_particles[cached_idx];

					if (use_simple_density) {
						// Simple density mode: uniform weight per particle. The value
						// filter is skipped (values are not deposited), but the SPATIAL
						// sphere filter below still applies — a previous version returned
						// here and silently disabled --bbox-sphere in this mode.
						out_v_orig = 1.0f;
					}
					else {
						// Apply value range filter
						if (out_v_orig < filter_min || out_v_orig > filter_max)
							return false;
					}

					// Apply spherical bounding filter if specified (radius > 0)
					if (bbox_sphere_r > 0.0f)
					{
						// Compute squared distance from particle to sphere center
						// Using squared distance avoids expensive sqrt operation
						double distSquared =
							(out_Pos[0] - bbox_sphere_pos[0]) * (out_Pos[0] - bbox_sphere_pos[0]) +
							(out_Pos[1] - bbox_sphere_pos[1]) * (out_Pos[1] - bbox_sphere_pos[1]) +
							(out_Pos[2] - bbox_sphere_pos[2]) * (out_Pos[2] - bbox_sphere_pos[2]);

						// Skip particles outside the spherical region
						if (distSquared > (bbox_sphere_r * bbox_sphere_r)) {
							return false;
						}
					}

					return true;
				};

				//////////////////

#ifdef WITH_GPU_CUDA

			/**
			 * @brief GPU implementation of bounding box calculation
			 *
			 * Launches CUDA kernel to compute min/max coordinates across all particles.
			 * Implemented in convert_vdb_kernel.cu
			 */
				void find_bbox_gpu(
					const float* d_pos_particles,
					size_t num_particles,
					const float* offset_position,
					float* bbox_min,
					float* bbox_max
				);

				/**
				 * @brief Main GPU dispatcher for converting I/O library data to grids
				 *
				 * Routes to appropriate GPU conversion function based on grid type.
				 * Implemented in convert_vdb_kernel.cu
				 */
				void convert_iolib_to_grid_gpu(
					const float* pos_particles,
					const size_t* particle_ids,
					const float* radius_particles,
					const float* value_particles,
					size_t num_particles,
					float particle_radius_multiplier,
					int bbox_dim,
					int* bbox_min_orig,
					double bbox_size_orig,
					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int block_name_id,
					float object_size,
					float& min_value,
					float& max_value,
					size_t& particles_count,
					VDBParticles& grid,
					double& transform_scale,
					float filter_min,
					float filter_max,
					common::SpaceData::AnimType anim_type,
					int frame_req,
					int frame,
					float* bbox_sphere_pos,
					float bbox_sphere_r,
					bool use_simple_density,
					double particle_radius_const,

					double bbox_x_min_norm,
					double bbox_x_max_norm,
					double bbox_y_min_norm,
					double bbox_y_max_norm,
					double bbox_z_min_norm,
					double bbox_z_max_norm,
					double scale_space_diagonal,

					float* offset_position,

					bool use_dense_loop_over_voxels = false,
					int calc_radius_neigh = 16,
					// Pre-built KD-tree for loop-over-voxels (null = disabled)
					void* d_kdtree_nodes = nullptr,
					int     kdtree_N = 0
				);


#endif // WITH_GPU_CUDA

				/**
				 * @brief CPU implementation of bounding box calculation using OpenMP
				 *
				 * Parallel reduction to compute min/max coordinates across all particles.
				 * Implemented in convert_vdb_kernel.cpp
				 */
				void find_bbox_cpu(
					const float* pos_particles,
					size_t num_particles,
					const float* offset_position,
					float* bbox_min,
					float* bbox_max
				);

				void convert_iolib_to_grid_cpu(
					const float* pos_particles,
					const size_t* particle_ids,
					const float* radius_particles,
					const float* value_particles,
					size_t num_particles,
					float particle_radius_multiplier,
					int bbox_dim,
					int* bbox_min_orig,
					double bbox_size_orig,
					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int block_name_id,
					float object_size,
					float& min_value,
					float& max_value,
					size_t& particles_count,
					VDBParticles& grid,
					double& transform_scale,
					float filter_min,
					float filter_max,
					common::SpaceData::AnimType anim_type,
					int frame_req,
					int frame,
					float* bbox_sphere_pos,
					float bbox_sphere_r,
					bool use_simple_density,
					double particle_radius_const,

					double bbox_x_min_norm,
					double bbox_x_max_norm,
					double bbox_y_min_norm,
					double bbox_y_max_norm,
					double bbox_z_min_norm,
					double bbox_z_max_norm,
					double scale_space_diagonal,

					float* offset_position,

					bool use_dense_loop_over_voxels = false,
					int calc_radius_neigh = 16,
					// Pre-built KD-tree for loop-over-voxels (null = disabled)
					void* kdtree_nodes = nullptr,
					int     kdtree_N = 0
				);



				// ── Voxel-centric dense grid conversion (loop over voxels, KD-tree particle query) ───────────

#ifdef WITH_CUDAKDTREE

			/**
			 * @brief CPU voxel-centric dense grid conversion using a pre-built KD-tree.
			 *
			 * @param tree_nodes  Pre-built float4 KD-tree (from cache_manager.voxel_kdtree_per_ptype).
			 *                    Each node: (x,y,z) = offset-adjusted world position, w = particle idx.
			 * @param tree_N      Number of tree nodes.
			 */
				void convert_to_dense_grid_loop_over_voxels_cpu(
					float4* tree_nodes,
					int           tree_N,
					const float* radius_particles,
					const float* value_particles,
					float  particle_radius_multiplier,
					int* bbox_min_orig,
					double bbox_size_orig,
					int    bbox_dim,

					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int    block_name_id,
					float* offset_position,
					bool   use_simple_density,
					double particle_radius_const,
					double scale_space_diagonal,

					int calc_radius_neigh,

					VoxelDenseManager* grid,
					float& min_value,
					float& max_value,
					size_t& particles_count
				);

#if defined(WITH_GPU_CUDA)
				/**
				 * @brief GPU voxel-centric dense grid conversion using a pre-built device KD-tree.
				 *
				 * @param d_tree  Pre-built device float4 KD-tree (cache_manager.d_voxel_kdtree_per_ptype).
				 * @param tree_N  Number of tree nodes.
				 */
				void convert_to_dense_grid_loop_over_voxels_gpu(
					float4* d_tree,
					int           tree_N,
					const float* d_radius_particles,
					const float* d_value_particles,
					float  particle_radius_multiplier,
					int* bbox_min_orig,
					double bbox_size_orig,
					int    bbox_dim,

					common::SpaceData::DenseType dense_type,
					common::SpaceData::DenseNorm dense_norm,
					int    block_name_id,
					float* offset_position,
					bool   use_simple_density,
					double particle_radius_const,
					double scale_space_diagonal,

					int calc_radius_neigh,

					VoxelDenseManager* grid,
					float& min_value,
					float& max_value,
					size_t& particles_count
				);
#endif // WITH_GPU_CUDA

#endif // WITH_CUDAKDTREE

				} // namespace kernel
		} // namespace vdb
	} // namespace common
}// namespace space_converter