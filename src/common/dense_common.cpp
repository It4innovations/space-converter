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

#include "dense_common.h"

#ifdef WITH_OPENMP
# include <omp.h>
#endif

#include <mpi.h>

#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>

#include "utility/dense_utility.h"

namespace common {
	namespace vdb {
#if 0
		double custom_asinh(double val)
		{
			return log(val + sqrt(1. + val * val));
		}

		//float custom_xexp(float v)
		//{
		//	return std::exp(v + -20.0f);
		//}

		void ConvertVDBBase::fill_voxels_v1(common::vdb::DenseParticles& grid,
			size_t pid, int px, int py, int pz, float value,
			int bbox_dim, int* bbox_min_orig, int* bbox_max_orig, double scale_space_diagonal) {

			grid.type = 1;
			///////////////////////////////////
			double radius = get_particle_radius(pid);

			//if (radius < min_radius_orig) min_radius_orig = radius;
			//if (radius > max_radius_orig) max_radius_orig = radius;

			double hsml = get_particle_hsml(pid);
			double mass = get_particle_mass(pid);
			double rho = get_particle_rho(pid);

			double radiusx_norm = ((double)radius) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double radiusy_norm = ((double)radius) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double radiusz_norm = ((double)radius) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double radiusx = (double)radiusx_norm * (double)bbox_dim / scale_space_diagonal;
			double radiusy = (double)radiusy_norm * (double)bbox_dim / scale_space_diagonal;
			double radiusz = (double)radiusz_norm * (double)bbox_dim / scale_space_diagonal;

			double radius_max_norm = std::max(radiusx, std::max(radiusy, radiusz));
			//if (radius_max_norm < min_radius) min_radius = radius_max_norm;
			//if (radius_max_norm > max_radius) max_radius = radius_max_norm;

			double hsmlx_norm = ((double)hsml) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double hsmly_norm = ((double)hsml) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double hsmlz_norm = ((double)hsml) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double hsmlx = (double)hsmlx_norm * (double)bbox_dim / scale_space_diagonal;
			double hsmly = (double)hsmly_norm * (double)bbox_dim / scale_space_diagonal;
			double hsmlz = (double)hsmlz_norm * (double)bbox_dim / scale_space_diagonal;

			double hsml_max_norm = std::max(hsmlx, std::max(hsmly, hsmlz)); //TODO

			//std::vector<double> densities(radius_max_norm + 1);
			//generate_normalized_gaussian(densities, 3.0);

			double sigma = 1.0 / 3.14159265358979323846;
			double sigma_hsml = sigma / (hsml_max_norm * hsml_max_norm * hsml_max_norm);

			int iradiusx = static_cast<int>(radiusx);
			int iradiusy = static_cast<int>(radiusy);
			int iradiusz = static_cast<int>(radiusz);

			int sx_min = px - iradiusx;
			int sx_max = px + iradiusx + 1;

			//#	pragma omp parallel for
			for (int sx = sx_min; sx < sx_max; sx++) {
				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

						if (sx < 0 || sy < 0 || sz < 0)
							continue;

						if (sx > grid.x() - 1 || sy > grid.y() - 1 || sz > grid.z() - 1)
							continue;

						double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
						double q = distance / hsml_max_norm;

						double W = 0.0;
						if (q >= 0.0 && q < 1.0) {
							W = sigma_hsml * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
						}
						else if (q >= 1.0 && q < 2.0) {
							W = sigma_hsml * 0.25 * (2 - q) * (2 - q) * (2 - q);
						}
						else {
							//W = 0.0;
							continue;
						}

						double w = mass / rho;

						//grid.add_color(sx, sy, sz, w * value * W);
						//grid.add_density(sx, sy, sz, w * W);

						size_t gindex = sx + sy * grid.dims[0] + sz * grid.dims[0] * grid.dims[1];

						//#pragma omp atomic
						//						grid.data_color[gindex] += w * W;

#pragma omp atomic
						grid.data_density[gindex] += w * value * W;

					}
				}
			}
			///////////////////////////////////
		}

		void ConvertVDBBase::fill_voxels_v2(common::vdb::DenseParticles& grid,
			size_t pid, int px, int py, int pz, float value,
			int bbox_dim, int* bbox_min_orig, int* bbox_max_orig, double scale_space_diagonal, float particle_fix_size) {

			grid.type = 2;
			////////////////////READING//////////////////////////
			//config
			float fix_size = particle_fix_size; // 1.0;
			float size_fac = 0.5f;// 1.0;
			float col_fac = 1.0f;
			float smooth_fac = 1.0f;
			float minrad_pix = 1.0f;
			float brightness = 1.f;

			//particle
			particle_sim_v2 p;

			//Reading positions
			double pos[3];
			get_particle_position(pid, pos);
			p.x = pos[0];
			p.y = pos[1];
			p.z = pos[2];

			//TODO: interpol_mode -> Reading velocities + Reading ids

			//Reading smoothing
			double hsml = get_particle_hsml(pid);

			//p.r = (fix_size == 0.0) ? hsml * size_fac : fix_size;
			p.r = (fix_size == 0.0 || hsml == 0.0) ? fix_size : hsml * size_fac;

			//Reading colors == values
			p.e = value * col_fac;

			//Reading intensity
			p.I = get_particle_rho(pid); //1;

			////////////////////NORMALIZING//////////////////////////
			bool log_int = false;
			bool log_col = false;
			bool asinh_col = false;
			bool col_vector = false;

			if (log_int)
			{
				if (p.I > 0)
				{
					p.I = log10(p.I);
				}
				else
					p.I = -38;
			}

			if (log_col)
			{
				if (p.e > 0)
				{
					p.e = log10(p.e);
				}
				else
					p.e = -38;
			}
			else
			{
				if (asinh_col)
					p.e = custom_asinh(p.e);
			}

			//////////////////PROJECTING////////////////////////
			bool classic = true;

			//constants
			double pi = 3.141592653589793238462643383279502884197;
			float h2sigma = 0.5 * pow(pi, -1. / 6.); // 0.413
			float sqrtpi = sqrt(pi); // 1.77245

			float bfak = 0.5 * pow(pi, -5. / 6.); // 0.19261
			float rfac = 0.75;

			if (!classic) {
				rfac = 1.;
			}

			if (classic) {
				p.I *= 0.5f * bfak / p.r;
				p.r *= 2;
			}
			else {
				p.I *= 8.f / (pi * p.r * p.r * p.r); // SPH kernel normalisation
				p.I *= (h2sigma * sqrtpi * p.r); // integral through the center
			}

			p.r *= smooth_fac; // Smoothing factor

			float rcorr = sqrt(p.r * p.r + minrad_pix * minrad_pix) / p.r;
			p.r *= rcorr;
			if (classic) {
				p.I /= rcorr;
			}
			else {
				p.I /= rcorr * rcorr;
			}

			////////////////COLORIZE///////////////////////
			p.e *= p.I * brightness; //TODO: 2x grids

			///////////////////////////////////////////////
			particle_sim_v2 pp = p;
			float rfacr = pp.r * rfac;

			//float posx = pp.x, posy = pp.y, posz = pp.y;
			//int minx = int(posx - rfacr + 1);
			//minx = std::max(minx, bbox_dim);
			//int maxx = int(posx + rfacr + 1);
			//maxx = std::min(maxx, bbox_dim);

			//int miny = int(posy - rfacr + 1);
			//miny = std::max(miny, bbox_dim);
			//int maxy = int(posy + rfacr + 1);
			//maxy = std::min(maxy, bbox_dim);

			//int minz = int(posz - rfacr + 1);
			//minz = std::max(minz, bbox_dim);
			//int maxz = int(posz + rfacr + 1);
			//maxz = std::min(maxz, bbox_dim);

			double radius = rfacr;
			double radiusx_norm = ((double)radius) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double radiusy_norm = ((double)radius) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double radiusz_norm = ((double)radius) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double radiusx = (double)radiusx_norm * (double)bbox_dim / scale_space_diagonal;
			double radiusy = (double)radiusy_norm * (double)bbox_dim / scale_space_diagonal;
			double radiusz = (double)radiusz_norm * (double)bbox_dim / scale_space_diagonal;

			double radius_max_norm = std::max(radiusx, std::max(radiusy, radiusz));

			float radsq = rfacr * rfacr;
			float sigma = h2sigma * pp.r;
			float stp = -1.f / (sigma * sigma);

			//a_eq_e==true
			float a = -pp.e;

			int iradiusx = static_cast<int>(radiusx);
			int iradiusy = static_cast<int>(radiusy);
			int iradiusz = static_cast<int>(radiusz);

			int sx_min = px - iradiusx;
			int sx_max = px + iradiusx + 1;

			exptable<float> xexp(-20.0);

			//#	pragma omp parallel for
			for (int sx = sx_min; sx < sx_max; sx++) {
				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

						if (sx < 0 || sy < 0 || sz < 0)
							continue;

						if (sx > grid.x() - 1 || sy > grid.y() - 1 || sz > grid.z() - 1)
							continue;

						float pre1 = xexp(stp * (sx - px) * (sx - px));
						float pre2 = xexp(stp * (sy - py) * (sy - py));
						float pre3 = xexp(stp * (sz - pz) * (sz - pz));

						//float pre1 = std::exp(stp * (sx - px) * (sx - px));
						//float pre2 = std::exp(stp * (sy - py) * (sy - py));
						//float pre3 = std::exp(stp * (sz - pz) * (sz - pz));

						float att = pre1 * pre2 * pre3;

						//float result = att * a;

						//grid.add(sx, sy, sz, pp.e);
						//grid.add1(sx, sy, sz, att * -pp.I);

						size_t gindex = sx + sy * grid.dims[0] + sz * grid.dims[0] * grid.dims[1];

						//#pragma omp atomic
						//						grid.data_color[gindex] += pp.e;

#pragma omp atomic
						grid.data_density[gindex] += att * a;
					}
				}
			}

			//todo: a_eq_e == false: lpic[x][y].r += xexp.expm1(att*a.r)*(lpic[x][y].r-q.r);
		}

#if 1
		void ConvertVDBBase::fill_voxels_v3(common::vdb::DenseParticles& grid,
			size_t pid, float value,
			int bbox_dim, int* bbox_min_orig, int* bbox_max_orig,
			double scale_space_diagonal,
			int dense_type, float particle_fix_size) {

			grid.type = 3;

			double pos[3];
			get_particle_position(pid, pos);

			double hsml = get_particle_hsml(pid);
			double rho = get_particle_rho(pid);
			double mass = get_particle_mass(pid);

			double radius = 2.0 * hsml;
			double radiusx_norm = ((double)radius) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double radiusy_norm = ((double)radius) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double radiusz_norm = ((double)radius) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double radiusx = (double)radiusx_norm * (double)bbox_dim;// / scale_space_diagonal;
			double radiusy = (double)radiusy_norm * (double)bbox_dim;// / scale_space_diagonal;
			double radiusz = (double)radiusz_norm * (double)bbox_dim;// / scale_space_diagonal;

			if (particle_fix_size != 0.0f) {
				radiusx = particle_fix_size;
				radiusy = particle_fix_size;
				radiusz = particle_fix_size;
			}

			double radiusxyz_max = std::max(radiusx, std::max(radiusy, radiusz));
			// if (radiusxyz_max < 0.5) {
			// 	radiusxyz_max = 0.5;
			// }

			int iradiusx = static_cast<int>(radiusx);
			int iradiusy = static_cast<int>(radiusy);
			int iradiusz = static_cast<int>(radiusz);

			///////////////////////////////////////////////////////////////////////////////
			double px_norm = ((double)pos[0] - (double)bbox_min_orig[0]) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double py_norm = ((double)pos[1] - (double)bbox_min_orig[1]) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double pz_norm = ((double)pos[2] - (double)bbox_min_orig[2]) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double dpx = (double)px_norm * (double)bbox_dim / scale_space_diagonal;
			double dpy = (double)py_norm * (double)bbox_dim / scale_space_diagonal;
			double dpz = (double)pz_norm * (double)bbox_dim / scale_space_diagonal;

			int px = static_cast<int>(dpx);
			int py = static_cast<int>(dpy);
			int pz = static_cast<int>(dpz);
			///////////////////////////////////////////////////////////////////////////////			

			//#	pragma omp parallel for
			for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

						int osx = sx - grid.offset[0];
						int osy = sy - grid.offset[1];
						int osz = sz - grid.offset[2];

						if (osx < 0 || osy < 0 || osz < 0)
							continue;

						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
							continue;

						float density = value;
						float norm = 1.0f;

						//if (iradiusx != 0 || iradiusy != 0 || iradiusz != 0) 
						{
							//int version = 4;
							//Gaussian Kernel v1
							if (dense_type == 1) {
								double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
								double bandwidth = radiusxyz_max;

								double pi = 3.141592653589793238462643383279502884197;
								double factor = 1.0 / (pow(2 * pi, 1.5) * pow(bandwidth, 3));
								double W = factor * exp(-0.5 * (distance * distance) / (bandwidth * bandwidth));

								density = (float)(W * value);
							}

							else if (dense_type == 2) {
								double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
								double q = distance / radiusxyz_max;

								double pi = 3.141592653589793238462643383279502884197;

								double sigma = 1.0 / 3.14159265358979323846;
								double sigma_hsml = sigma / (radiusxyz_max * radiusxyz_max * radiusxyz_max);

								double W = 0.0;
								if (q >= 0.0 && q < 1.0) {
									W = sigma_hsml * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
								}
								else if (q >= 1.0 && q < 2.0) {
									W = sigma_hsml * 0.25 * (2 - q) * (2 - q) * (2 - q);
								}
								else {
									//W = 0.0;
									continue;
								}

								double w = mass / rho;

								density = (float)(w * W * value);
								//norm = (float)(w * W);
							}

							else if (dense_type == 3) {
								double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
								double q = distance / radiusxyz_max;

								double pi = 3.141592653589793238462643383279502884197;

								double sigma_hsml = 8.0 / (pi * radiusxyz_max * radiusxyz_max * radiusxyz_max);

								double W = 0.0;
								if (q >= 0.0 && q < 0.5) {
									W = sigma_hsml * (1.0 - 6.0 * q * q + 6 * q * q * q);
								}
								else if (q >= 0.5 && q < 1.0) {
									W = sigma_hsml * 2.0 * (1.0 - q) * (1.0 - q) * (1.0 - q);
								}
								else {
									//W = 0.0;
									continue;
								}

								density = (float)(W * value);
								//norm = (float)W;
							}

							else if (dense_type == 4) {
								double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
								double q = distance / radiusxyz_max;

								double W = 0.0;
								if (q >= 0.0 && q < 1.0) {
									W = (1.0 - 1.5 * q * q + 0.75 * q * q * q);
								}
								else if (q >= 1.0 && q < 2.0) {
									W = 0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q);
								}
								else {
									//W = 0.0;
									continue;
								}

								density = (float)(W * value);
								//norm = (float)W;
							}

							else if (dense_type == 5) {
#if 1
								if (hsml == 0.0) {
									hsml = 2.0 * pow((mass / rho), 1.0 / 3.0);
								}

								double weight = mass / rho;
								double radkernel = 2.0; // Kernel radius in smoothing lengths(usually 2 for cubic spline)
								double radkernel2 = radkernel * radkernel;

								double hsmlxyz_max = radiusxyz_max / 2.0;

								// Skip particles with non - positive weight or smoothing length
								if (weight <= 0.0 || hsmlxyz_max <= 0.0)
									continue;

								double inv_hi = 1.0 / hsmlxyz_max;
								double inv_hi2 = inv_hi * inv_hi;
								double kernel_radius = radkernel * hsmlxyz_max;
								double kernel_radius2 = kernel_radius * kernel_radius;

								double pi = 3.141592653589793238462643383279502884197;
								double cnormk3D_i = 1.0 / (pi * hsmlxyz_max * hsmlxyz_max * hsmlxyz_max);
								double weight_i = weight * cnormk3D_i;

								double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
								double q = distance / kernel_radius;

								double kernel_weight = 0.0;
								//w_cubic
								if (q >= 0.0 && q < 1.0) {
									kernel_weight = (1.0 - 1.5 * q * q + 0.75 * q * q * q);
								}
								else if (q >= 1.0 && q < 2.0) {
									kernel_weight = 0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q);
								}
								else {
									//kernel_weight = 0.0;
									continue;
								}

								double total_weight = weight_i * kernel_weight;

								density = (float)(total_weight * value);
								//norm = (float)total_weight;							 
#endif
							}
							else if (dense_type == 6) {
								norm = 0.0f;
							}
							else if (dense_type == 7) {
								//do nothing
							}
						}


						size_t gindex = grid.get_index(osx, osy, osz);

#pragma omp atomic
						grid.data_density[gindex] += density;

#ifndef WITH_NO_DATA_TEMP
#pragma omp atomic
						grid.data_temp[gindex] += norm;
#endif						
					}
				}
			}

			//todo: a_eq_e == false: lpic[x][y].r += xexp.expm1(att*a.r)*(lpic[x][y].r-q.r);
		}
#else

		void ConvertVDBBase::fill_voxels_v3(common::vdb::DenseParticles& grid,
			size_t pid, float value,
			int bbox_dim, int* bbox_min_orig, int* bbox_max_orig,
			double scale_space_diagonal,
			int dense_type, float particle_fix_size) {

			grid.type = 3;

			double pos[3];
			get_particle_position(pid, pos);

			double hsml = get_particle_hsml(pid);
			double rho = get_particle_rho(pid);
			double mass = get_particle_mass(pid);

			// if (hsml == 0.0) {
			// 	hsml = 1.2 * pow(mass / rho, 1.0 / 3.0);
			// }

			double radius_orig = 2.0 * hsml;

			////////////////////////////////////////////////////////////////////////////////////
			double radiusx_norm = ((double)radius_orig) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double radiusy_norm = ((double)radius_orig) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double radiusz_norm = ((double)radius_orig) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double radiusx = (double)radiusx_norm * (double)bbox_dim;
			double radiusy = (double)radiusy_norm * (double)bbox_dim;
			double radiusz = (double)radiusz_norm * (double)bbox_dim;

			if (particle_fix_size != 0.0f) {
				radiusx = particle_fix_size;
				radiusy = particle_fix_size;
				radiusz = particle_fix_size;
			}

			// double radius_min = std::min(radiusx, std::min(radiusy, radiusz));

			// // fix min radius
			// if (radius_min < 1.0) {
			// 	radiusx = radiusx * (1.0 / radius_min);
			// 	radiusy = radiusy * (1.0 / radius_min);
			// 	radiusz = radiusz * (1.0 / radius_min);
			// }

			double radius = std::max(radiusx, std::max(radiusy, radiusz));

			int iradiusx = static_cast<int>(radiusx);
			int iradiusy = static_cast<int>(radiusy);
			int iradiusz = static_cast<int>(radiusz);

			///////////////////////////////////////////////////////////////////////////////
			double px_norm = ((double)pos[0] - (double)bbox_min_orig[0]) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double py_norm = ((double)pos[1] - (double)bbox_min_orig[1]) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double pz_norm = ((double)pos[2] - (double)bbox_min_orig[2]) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double px = (double)px_norm * (double)bbox_dim / scale_space_diagonal;
			double py = (double)py_norm * (double)bbox_dim / scale_space_diagonal;
			double pz = (double)pz_norm * (double)bbox_dim / scale_space_diagonal;

			int ipx = static_cast<int>(px);
			int ipy = static_cast<int>(py);
			int ipz = static_cast<int>(pz);
			///////////////////////////////////////////////////////////////////////////////

			int sx_min = ipx - iradiusx;
			int sx_max = ipx + iradiusx + 1;

			int sy_min = ipy - iradiusy;
			int sy_max = ipy + iradiusy + 1;

			int sz_min = ipz - iradiusy;
			int sz_max = ipz + iradiusz + 1;

			//#	pragma omp parallel for
			for (int sx = sx_min; sx < sx_max; sx++) {
				for (int sy = sy_min; sy < sy_max; sy++) {
					for (int sz = sz_min; sz < sz_max; sz++) {

						int osx = sx - grid.offset[0];
						int osy = sy - grid.offset[1];
						int osz = sz - grid.offset[2];

						if (osx < 0 || osy < 0 || osz < 0)
							continue;

						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
							continue;

						float density = value;
						float norm = 1.0f;

						//double sx_orig = ((double)sx / (double)bbox_dim) * ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]) + (double)bbox_min_orig[0];
						//double sy_orig = ((double)sy / (double)bbox_dim) * ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]) + (double)bbox_min_orig[1];
						//double sz_orig = ((double)sz / (double)bbox_dim) * ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]) + (double)bbox_min_orig[2];

						double distance = sqrt(((double)sx - px) * ((double)sx - px) + ((double)sy - py) * ((double)sy - py) + ((double)sz - pz) * ((double)sz - pz));

						//int version = 4;
						//Gaussian Kernel v1
						if (dense_type == 1) {
							//double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
							double bandwidth = radius;

							double pi = 3.141592653589793238462643383279502884197;
							double factor = 1.0 / (pow(2 * pi, 1.5) * pow(bandwidth, 3));
							double W = factor * exp(-0.5 * (distance * distance) / (bandwidth * bandwidth));

							//density = (float)(att * value);
							//double w = mass / rho;

							// density = (float)(w * W * value);
							// norm = (float)(w * W);
							density = (float)(W * value);
						}

						else if (dense_type == 2) {
							//double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
							double q = distance / radius;

							double pi = 3.141592653589793238462643383279502884197;

							double sigma = 1.0 / 3.14159265358979323846;
							double sigma_hsml = sigma / (radius * radius * radius);

							double W = 0.0;
							if (q >= 0.0 && q < 1.0) {
								W = sigma_hsml * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
							}
							else if (q >= 1.0 && q < 2.0) {
								W = sigma_hsml * 0.25 * (2 - q) * (2 - q) * (2 - q);
							}
							else {
								//W = 0.0;
								continue;
							}

							double w = mass / rho;

							// density = (float)(w * W * value);
							// norm = (float)(w * W);
							density = (float)(W * value);
						}

						else if (dense_type == 3) {
							//double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
							double q = distance / radius;

							double pi = 3.141592653589793238462643383279502884197;

							double sigma_hsml = 8.0 / (pi * radius * radius * radius);

							double W = 0.0;
							if (q >= 0.0 && q < 0.5) {
								W = sigma_hsml * (1.0 - 6.0 * q * q + 6 * q * q * q);
							}
							else if (q >= 0.5 && q < 1.0) {
								W = sigma_hsml * 2.0 * (1.0 - q) * (1.0 - q) * (1.0 - q);
							}
							else {
								//W = 0.0;
								continue;
							}

							density = (float)(W * value);
							norm = (float)W;
						}

						else if (dense_type == 4) {
							//double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
							double q = distance / radius;

							double W = 0.0;
							if (q >= 0.0 && q < 1.0) {
								W = (1.0 - 1.5 * q * q + 0.75 * q * q * q);
							}
							else if (q >= 1.0 && q < 2.0) {
								W = 0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q);
							}
							else {
								//W = 0.0;
								continue;
							}

							density = (float)(W * value);
							norm = (float)W;
						}

						else if (dense_type == 5) {
#if 1
							if (hsml == 0.0) {
								hsml = 2.0 * pow((mass / rho), 1.0 / 3.0);
							}

							double weight = mass / rho;
							double radkernel = 2.0; // Kernel radius in smoothing lengths(usually 2 for cubic spline)
							double radkernel2 = radkernel * radkernel;

							double hsmlxyz_max = radius / 2.0;

							// Skip particles with non - positive weight or smoothing length
							if (weight <= 0.0 || hsmlxyz_max <= 0.0)
								continue;

							double inv_hi = 1.0 / hsmlxyz_max;
							double inv_hi2 = inv_hi * inv_hi;
							double kernel_radius = radkernel * hsmlxyz_max;
							double kernel_radius2 = kernel_radius * kernel_radius;

							double pi = 3.141592653589793238462643383279502884197;
							double cnormk3D_i = 1.0 / (pi * hsmlxyz_max * hsmlxyz_max * hsmlxyz_max);
							double weight_i = weight * cnormk3D_i;

							//double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
							double q = distance / kernel_radius;

							double kernel_weight = 0.0;
							//w_cubic
							if (q >= 0.0 && q < 1.0) {
								kernel_weight = (1.0 - 1.5 * q * q + 0.75 * q * q * q);
							}
							else if (q >= 1.0 && q < 2.0) {
								kernel_weight = 0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q);
							}
							else {
								//kernel_weight = 0.0;
								continue;
							}

							double total_weight = weight_i * kernel_weight;

							density = (float)(total_weight * value);
							norm = (float)total_weight;
#endif
						}


						size_t gindex = grid.get_index(osx, osy, osz);

#pragma omp atomic
						grid.data_density[gindex] += density;

#ifndef WITH_NO_DATA_TEMP						
#pragma omp atomic
						grid.data_temp[gindex] += norm;
#endif						
					}
				}
			}

			//todo: a_eq_e == false: lpic[x][y].r += xexp.expm1(att*a.r)*(lpic[x][y].r-q.r);
		}
#endif
		const int maxcoltable = 1000;
		const int radkernel2 = 4;

		std::vector<double> coltable;
		double dq2table = 0.0f, ddq2table = 0.0f;

		double wfunc(double q2) {
			double q;
			double w_cubic;

			if (q2 < 1.0f) {
				q = std::sqrt(q2);
				w_cubic = 1.0f - 1.5f * q2 + 0.75f * q2 * q;
			}
			else if (q2 < 4.0f) {
				q = std::sqrt(q2);
				w_cubic = 0.25f * std::pow(2.0f - q, 3);
			}
			else {
				w_cubic = 0.0f;
			}

			return w_cubic;
		}

		void setup_integratedkernel() {
			if (coltable.size() > 0)
				return;

			coltable.resize(maxcoltable + 1);

			constexpr int npts = 100;

			// Setup the increment steps
			dq2table = (double)radkernel2 / (double)maxcoltable;
			ddq2table = 1.0f / (double)dq2table;

			for (int i = 0; i < maxcoltable; ++i) {
				// Tabulate for cylindrical r**2 between 0 and radkernel**2
				double rxy2 = i * dq2table;

				// Integrate z between 0 and sqrt(radkernel^2 - rxy^2)
				double deltaz = std::sqrt(radkernel2 - rxy2);
				double dz = deltaz / static_cast<double>(npts - 1);
				double coldens = 0.0f;

				if (std::isnan(deltaz)) {
					std::cerr << "WARNING: NaN in kernel table setup" << std::endl;
				}

				for (int j = 1; j <= npts; ++j) {
					double z = (j - 1) * dz;
					double q2 = rxy2 + z * z;
					double wkern = wfunc(q2);

					// Apply trapezoidal rule
					if (j == 1 || j == npts) {
						coldens += 0.5f * wkern * dz;
					}
					else {
						coldens += wkern * dz;
					}
				}
				double cnormk3D = 1.0 / 3.141592653589793238462643383279502884197;
				coltable[i] = 2.0f * coldens * cnormk3D;
			}

			// Final table entry
			coltable[maxcoltable] = 0.0f;
		}

		double wfromtable(double q2) {
			setup_integratedkernel();

			// Find nearest index in the table
			int index = std::max(static_cast<int>(q2 * ddq2table), 0);
			int index1 = std::min(index + 1, maxcoltable);

			// Compute increment along from this index
			double dxx = q2 - index * dq2table;

			// Compute gradient
			double dwdx = (coltable[index1] - coltable[index]) * ddq2table;

			// Compute value of integrated kernel
			return coltable[index] + dwdx * dxx;
		}

		void ConvertVDBBase::fill_voxels_v4(
			common::vdb::DenseParticles& grid,
			size_t pid,
			float value,
			int bbox_dim,
			int* bbox_min_orig,
			int* bbox_max_orig,
			double scale_space_diagonal,
			common::SpaceData::DenseType dense_type,
			float particle_fix_size,
			int particle_type,
			int block_name_id,
			double* pos
		) {

			grid.type = 4;

			double hsml = 0.0;
			double rho = 0.0;
			double mass = 0.0;

			double norm_fac = (double)bbox_dim / scale_space_diagonal;

			if (block_name_id == 0) { //only pos
				value = 1.0f;
			}

			if (particle_type == 0) {  //only sph
				hsml = get_particle_hsml(pid);
				rho = get_particle_rho(pid);
				mass = get_particle_mass(pid);
			}

			if (particle_fix_size != 0.0f) {
				hsml = particle_fix_size * pow((mass / rho), 1.0 / 3.0); // todo
			}

			double radius = 2.0 * hsml;
			if (radius_particles_per_ptype.size() > 0 && radius_particles_per_ptype[particle_type].size() > 0) {
				radius = radius_particles_per_ptype[particle_type][pid - particles_ptype_offset[particle_type]];
			}

			double radiusx_norm = ((double)radius) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double radiusy_norm = ((double)radius) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double radiusz_norm = ((double)radius) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			//double radiusxyz_norm_max = std::max(radiusx_norm, std::max(radiusy_norm, radiusz_norm));
			//if (radiusxyz_norm_max < 1.0 / norm_fac) {
			//	radiusxyz_norm_max = 1.0 / norm_fac;
			//}

			double radiusx = (double)radiusx_norm * norm_fac;
			double radiusy = (double)radiusy_norm * norm_fac;
			double radiusz = (double)radiusz_norm * norm_fac;

			//radiusx = std::max(radiusx, 0.5);
			//radiusy = std::max(radiusy, 0.5);
			//radiusz = std::max(radiusz, 0.5);

			// if (particle_fix_size != 0.0f) {
			// 	radiusx = particle_fix_size;
			// 	radiusy = particle_fix_size;
			// 	radiusz = particle_fix_size;
			// }

			double radiusxyz_max = std::max(radiusx, std::max(radiusy, radiusz));

			int iradiusx = static_cast<int>(radiusx);
			int iradiusy = static_cast<int>(radiusy);
			int iradiusz = static_cast<int>(radiusz);

			///////////////////////////////////////////////////////////////////////////////
			double px_norm = ((double)pos[0] - (double)bbox_min_orig[0]) / ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
			double py_norm = ((double)pos[1] - (double)bbox_min_orig[1]) / ((double)bbox_max_orig[1] - (double)bbox_min_orig[1]);
			double pz_norm = ((double)pos[2] - (double)bbox_min_orig[2]) / ((double)bbox_max_orig[2] - (double)bbox_min_orig[2]);

			double dpx = (double)px_norm * norm_fac;
			double dpy = (double)py_norm * norm_fac;
			double dpz = (double)pz_norm * norm_fac;

			int px = static_cast<int>(dpx);
			int py = static_cast<int>(dpy);
			int pz = static_cast<int>(dpz);

			///////////////////////////////////////////////////////////////////////////////
			//double local_norm = 0.0;

			for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

						int osx = sx - grid.offset[0];
						int osy = sy - grid.offset[1];
						int osz = sz - grid.offset[2];

						if (osx < 0 || osy < 0 || osz < 0)
							continue;

						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
							continue;

						double density = 1.0;
						double norm = 1.0;

						//double sx_norm = ((double)sx / norm_fac);
						//double sy_norm = ((double)sy / norm_fac);
						//double sz_norm = ((double)sz / norm_fac);

						//double radius_min = (0.5 / (double)bbox_dim) * ((double)bbox_max_orig[0] - (double)bbox_min_orig[0]);
						//double hsmooth_norm = hsml; // TODO min HSML

						//double distance_norm = sqrt((sx_norm - px_norm) * (sx_norm - px_norm) + (sy_norm - py_norm) * (sy_norm - py_norm) + (sz_norm - pz_norm) * (sz_norm - pz_norm));

						double dx = dpx - sx;
						double dy = dpy - sy;
						double dz = dpz - sz;
						double distance_norm = std::sqrt(dx * dx + dy * dy + dz * dz);

						if (dense_type == common::SpaceData::DenseType::eType1 /*|| dense_type == common::SpaceData::DenseType::eType2 log10*/) {
							//double sigma = radiusxyz_max;

							// Gaussian kernel formula
							//double W = std::exp(-1.0 * (distance_norm * distance_norm) / (2 * sigma * sigma));												

							if (distance_norm > radiusxyz_max)
								continue;

							double h = radiusxyz_max;

							double pi = 3.141592653589793238462643383279502884197;
							double hpi = h * sqrt(pi);
							double W = (1.0 / (hpi * hpi * hpi)) * std::exp(-1.0 * distance_norm * distance_norm / (h * h));

							double v = (mass / rho) / (radiusx * radiusy * radiusz);
							density = W;
							//norm = W;
							//norm = v * W;
							norm = 0.0;
						}
						else if (dense_type == common::SpaceData::DenseType::eType2) {
							double norm_fac_inverse = radiusxyz_max;
							//double norm_fac_inverse = 1.0 / norm_fac;
							double distance_max = std::sqrt(norm_fac_inverse * norm_fac_inverse + norm_fac_inverse * norm_fac_inverse + norm_fac_inverse * norm_fac_inverse);
							double W = (1.0 - distance_norm / distance_max);
							if (W < 0.0)
								continue;

							if (W > 1.0)
								W = 1.0;

							density = W;
							//norm = W;
						}
						else if (dense_type == common::SpaceData::DenseType::eType3) {
							//double distance = sqrt((sx - px) * (sx - px) + (sy - py) * (sy - py) + (sz - pz) * (sz - pz));
							double bandwidth = radiusxyz_max;

							double pi = 3.141592653589793238462643383279502884197;
							double factor = 1.0 / (pow(2 * pi, 1.5) * pow(bandwidth, 3));
							double W = factor * exp(-0.5 * (distance_norm * distance_norm) / (bandwidth * bandwidth));

							//density = (float)(att * value);
							//double w = mass / rho;

							// density = (float)(w * W * value);
							// norm = (float)(w * W);
							density = W;
							//norm = w * W;
						}
						else if (dense_type == common::SpaceData::DenseType::eType4) {
							//double weight = mass / rho;
							//double radkernel = 2.0; // Kernel radius in smoothing lengths(usually 2 for cubic spline)
							//double radkernel2 = radkernel * radkernel;

							//double hsmlxyz_max = radiusxyz_max / 2.0;

							//// Skip particles with non - positive weight or smoothing length
							//if (weight <= 0.0 || hsmlxyz_max <= 0.0)
							//	continue;

							//double inv_hi = 1.0 / hsmlxyz_max;
							//double inv_hi2 = inv_hi * inv_hi;
							//double kernel_radius = radkernel * hsmlxyz_max;
							//double kernel_radius2 = kernel_radius * kernel_radius;

							//double pi = 3.141592653589793238462643383279502884197;
							//double cnormk3D_i = 1.0 / (pi * hsmlxyz_max * hsmlxyz_max * hsmlxyz_max);
							//double weight_i = weight * cnormk3D_i;

							//double q = distance / kernel_radius;

							//double total_weight = weight_i * kernel_weight;

							double q = distance_norm / radiusxyz_max; //(2.0 * hsmooth_norm); // TODO min HSML

							double h = radiusxyz_max / 2.0;

							double weight = mass / (rho * h * h * h); //m/h**ndim,m/(rho h**ndim)

							double W = 0;
							//w_cubic
							if (q >= 0.0 && q < 1.0) {
								W = (1.0 - 1.5 * q * q + 0.75 * q * q * q);
							}
							else if (q >= 1.0 && q < 2.0) {
								W = 0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q);
							}
							else {
								//W = 0.0;
								continue;
							}

							//density = (float)(weight * hsml * W * value);
							density = weight * h * W;
							//norm = (float)(weight * hsml * W);
							//norm = weight * h * W;
						}
						else if (dense_type == common::SpaceData::DenseType::eType5) {
#if 0							
							double q = distance_norm / radiusxyz_max;
							double q2 = q * q;

							if (q2 > radkernel2)
								continue;

							double h = radiusxyz_max / 2.0;

							double wab = wfromtable(q2);
							double weight = mass / (rho * h * h * h);
							density = weight * h * wab;
#endif
							//double h = radiusxyz_max / 4.0;// / 2.0;
							//float r = distance_norm;
							//float w = 0.0f; //W(r, h);
							//float q = r / h;
							//double PI = 3.141592653589793238462643383279502884197;
							//float sigma = 1.0f / (PI * std::pow(h, 3)); // Normalization factor
							//if (q >= 0 && q <= 1) {
							//	w = sigma * (1 - 1.5f * q * q + 0.75f * q * q * q);
							//} else if (q > 1 && q <= 2) {
							//	w = sigma * 0.25f * std::pow(2 - q, 3);
							//}else{
							//	continue;
							//}					

							//density = mass / rho * w;
							//norm = w;

							//if (distance_norm > radiusxyz_max)
							//	continue;

							//double h = radiusxyz_max;

							//double pi = 3.141592653589793238462643383279502884197;
							//double hpi = h * sqrt(pi);
							//double W = (1.0 / (hpi * hpi * hpi)) * std::exp(-1.0 * distance_norm * distance_norm / (h * h));

							//density = W; //(mass / rho)
							//local_norm += W;
							////norm = W;
							//norm = 0.0;

							//Robert Wissing, 2024
							double h = radiusxyz_max / 2.0;

							double PI = 3.141592653589793238462643383279502884197;

							double q = distance_norm / h;
							double Wkernel = 0.0;
							if (q < 1) {
								Wkernel = (8.0 / (PI * std::pow(h, 3))) * (1 - 1.5 * std::pow(q, 2) + 0.75 * std::pow(q, 3));
							}
							else if (q >= 1 && q < 2) {
								Wkernel = (8.0 / (PI * std::pow(h, 3))) * (std::pow(2 - q, 3) / 4.0);
							}
							else {
								continue;
							}

							double Vpart = mass / rho;

							density = Vpart * Wkernel;
							//norm = 0; // lower density
							//norm = 1; //higher density
							norm = density;

						}
						else if (dense_type == common::SpaceData::DenseType::eType6) {
							////////Monaghan 1992

							//double sigma = radiusxyz_max;

							// Gaussian kernel formula
							//double W = std::exp(-1.0 * (distance_norm * distance_norm) / (2 * sigma * sigma));												

							double h = radiusxyz_max / 2.0;
							//if (distance_norm > b * h)
							//	continue;

							double PI = 3.141592653589793238462643383279502884197;

							double q = distance_norm / radiusxyz_max;
							double W = 0;
							if (q < 1) {
								W = (1.0 / (PI * std::pow(h, 3))) * (1 - 1.5 * q * q + 0.75 * std::pow(q, 3));
							}
							else if (q >= 1 && q < 2) {
								W = (1.0 / (PI * std::pow(h, 3))) * (std::pow(2 - q, 3) / 4.0);
							}
							else {
								continue;
							}

							density = W;
							//norm = 0; // lower density
							//norm = 1; //higher density

						}
						else if (dense_type == common::SpaceData::DenseType::eType7) {
							double pi = 3.141592653589793238462643383279502884197;
							double C = 495.0 / (32.0 * pi);
							double r = distance_norm;
							double h = radiusxyz_max;
							double eta = r / h;
							if (eta >= 2.0) //orig >1
								continue;

							double Wrh = C * pow(1.0 - eta, 6) * (1.0 + 6.0 * eta + (35.0 / 3.0) * eta * eta) / (h * h * h);
							density = Wrh;

							norm = 0.0;
						}

						//size_t gindex = grid.get_index(osx, osy, osz);

						//kernels[iradiusx * iradiusy * (sz - pz + iradiusz) + iradiusx * (sy - py + iradiusy) + (sx - px + iradiusx)] = density;
						//kernels_sum += norm;

// #if 0
// 						kernels[sx + sy * kernel_dim_x + sz * kernel_dim_x * kernel_dim_y] = density;
// 						kernels_sum += density;
// #else
						size_t gindex = grid.get_index(osx, osy, osz);

						float d = density * value;
						float n = norm;

#pragma omp atomic
						grid.data_density[gindex] += d;

						//#pragma omp atomic compare
						//						grid.data_density[gindex] = std::max(grid.data_density[gindex], d);

#ifndef WITH_NO_DATA_TEMP
#pragma omp atomic
						grid.data_temp[gindex] += n;
#endif
						// #endif

												//#pragma omp atomic
												//						grid.data_density[gindex] += density;
												//
												//#pragma omp atomic
												//						grid.data_temp[gindex] += norm;
					}
				}
			}

			// #if 0
			// 			for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
			// 				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
			// 					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

			// 						int osx = sx - grid.offset[0];
			// 						int osy = sy - grid.offset[1];
			// 						int osz = sz - grid.offset[2];

			// 						if (osx < 0 || osy < 0 || osz < 0)
			// 							continue;

			// 						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
			// 							continue;

			// 						size_t gindex = grid.get_index(osx, osy, osz);

			// #pragma omp atomic
			// 						grid.data_density[gindex] += d;
			// 					}
			// 				}
			// 			}
			// #endif

			// #if 0
			// 			double density_check = 0.0;

			// 			//for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
			// 			//	for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
			// 			//		for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {
			// 			for (int sx = 0; sx < kernel_dim_x; ++sx) {
			// 				for (int sy = 0; sy < kernel_dim_y; ++sy) {
			// 					for (int sz = 0; sz < kernel_dim_z; ++sz) {

			// 						int osx = (sx + px - radiusx) - grid.offset[0];
			// 						int osy = (sy + py - radiusy) - grid.offset[1];
			// 						int osz = (sz + pz - radiusz) - grid.offset[2];

			// 						if (osx < 0 || osy < 0 || osz < 0)
			// 							continue;

			// 						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
			// 							continue;

			// 						double W = kernels[sx + sy * kernel_dim_x + sz * kernel_dim_x * kernel_dim_y];

			// 						density_check += W / kernels_sum;
			// 						float density = value * W / kernels_sum;

			// 						float norm = 0.0;

			// 						size_t gindex = grid.get_index(osx, osy, osz);

			// #pragma omp atomic
			// 						grid.data_density[gindex] += density;

			// #pragma omp atomic
			// 						grid.data_temp[gindex] += norm;
			// 					}
			// 				}
			// 			}

			// 			//check
			// 			if (fabs(density_check - 1.0) > FDATA_EPSILON) {
			// 				std::cerr << "error: check sum is not one: " << density_check << std::endl;
			// 			}
			// 			//todo: a_eq_e == false: lpic[x][y].r += xexp.expm1(att*a.r)*(lpic[x][y].r-q.r);
			// #endif
		}
#endif

		void ConvertVDBBase::fill_voxels_v5(
			common::vdb::DenseParticles& grid,
			size_t pid,
			float value,
			int bbox_dim,
			int* bbox_min_orig,
			double bbox_size_orig,
			double scale_space_diagonal,
			common::SpaceData::DenseType dense_type,
			common::SpaceData::DenseNorm dense_norm,
			float particle_fix_size,
			int particle_type,
			int block_name_id,
			double* pos
		) {

			//grid.type = 4;

			//double hsml = 0.0;
			//double rho = 0.0;
			//double mass = 0.0;

			double norm_fac = (double)bbox_dim / scale_space_diagonal;
			double len2pix = norm_fac / bbox_size_orig;

			if (block_name_id == 0) { //only pos
				value = 1.0f;
			}

			//0.11642141634818093
			//0.67700000000000005

			//utility::dense::gadget::GadgetPhysical GU = utility::dense::gadget::create_gadget_physical(redshift, hubble_param);
			// test - reset values
			//GU.x_physical = 1;
			//GU.m_physical = 1;
			//GU.rho_cgs = 1;
			//GU.x_cgs = 1;

			double hsml = get_particle_radius(
				pid,
				bbox_dim,
				bbox_min_orig,
				bbox_size_orig,
				scale_space_diagonal,
				particle_fix_size,
				particle_type
			);

			//double pos_physical[3] = {
			//	((double)pos[0] * GU.x_physical),
			//	((double)pos[1] * GU.x_physical),
			//	((double)pos[2] * GU.x_physical)
			//};

			//double bbox_min_orig_physical[3] = {
			//	((double)bbox_min_orig[0] * GU.x_physical),
			//	((double)bbox_min_orig[1] * GU.x_physical),
			//	((double)bbox_min_orig[2] * GU.x_physical)
			//};

			//for (int p = 0; p < N; ++p) {
			//double bin_q = rho; //Bin_Q == rho
			//if (bin_q == 0.0 && !calc_mean) continue;

#if 0
			double mass = get_particle_mass(pid);
			double rho = get_particle_rho(pid);

			int Npixels[3] = { grid.x(), grid.y(), grid.z() };

			utility::dense::cic_3d::Quantities3D q = utility::dense::cic_3d::get_quantities_3D(pos_physical, 1.0 / len2pix /*Weights[p]*/, hsml, value, mass, Npixels, len2pix);
			//double x = q.pos_pix[0];
			//double y = q.pos_pix[1];
			//double z = q.pos_pix[2];
			double x = ((double)pos_physical[0] - (double)bbox_min_orig_physical[0]) * len2pix;
			double y = ((double)pos_physical[1] - (double)bbox_min_orig_physical[1]) * len2pix;
			double z = ((double)pos_physical[2] - (double)bbox_min_orig_physical[2]) * len2pix;

			//x -= grid.offset[0];
			//y -= grid.offset[1];
			//z -= grid.offset[2];

			utility::dense::cic_3d::MinMax iRange = utility::dense::cic_3d::pix_index_min_max(x, q.hsml_pix, Npixels[0]);
			utility::dense::cic_3d::MinMax jRange = utility::dense::cic_3d::pix_index_min_max(y, q.hsml_pix, Npixels[1]);
			utility::dense::cic_3d::MinMax kRange = utility::dense::cic_3d::pix_index_min_max(z, q.hsml_pix, Npixels[2]);

			utility::dense::sph_kernel::WendlandC6 kernel;
			size_t wk_size = (iRange.max + 1LL - iRange.min) * (jRange.max + 1LL - jRange.min) * (kRange.max + 1LL - kRange.min);

			int tid = omp_get_thread_num();
			int num_threads = omp_get_num_threads();

			// calc temp var
			size_t data_temp_size = grid.data_temp.size() * sizeof(float) / num_threads;
			size_t offset = data_temp_size * tid;
			char* data_temp = (char*)grid.data_temp.data() + offset;
			char* data_temp2 = data_temp + data_temp_size / 2;

			double* wk = (double*)data_temp;
			double* V = (double*)data_temp2;

			utility::dense::cic_3d::CalcWeights result_weights = utility::dense::cic_3d::calculate_weights(
				wk,
				V,
				iRange.min,
				iRange.max,
				jRange.min,
				jRange.max,
				kRange.min,
				kRange.max,
				x, y, z,
				q.hsml_pix,
				q.hsml_inv,
				kernel,
				Npixels[0],
				Npixels[1],
				true);

			double kernel_norm = q.volume / static_cast<double>(result_weights.n_distr_pix);
			double volume_norm = kernel_norm * result_weights.weight_per_pix * q.weight * len2pix;

			for (int i = iRange.min; i <= iRange.max; ++i) {
				for (int j = jRange.min; j <= jRange.max; ++j) {
					for (int k = kRange.min; k <= kRange.max; ++k) {

						int osx = i - grid.offset[0];
						int osy = j - grid.offset[1];
						int osz = k - grid.offset[2];

						if (osx < 0 || osy < 0 || osz < 0)
							continue;

						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
							continue;

						size_t idx = grid.get_index(osx, osy, osz);  //utility::dense::cic_3d::calculate_index(i, j, k, Npixels[0], Npixels[1]);
						size_t idx_block = utility::dense::cic_3d::calculate_index(i - iRange.min, j - jRange.min, k - kRange.min, (jRange.max + 1 - jRange.min), (kRange.max + 1 - kRange.min));

						if (idx_block >= wk_size) {
							printf("ERROR denstiy: idx_block >= wk_size\n", idx_block, wk_size);
						}

						double pix_weight = wk[idx_block] * V[idx_block] * volume_norm;
						//if (pix_weight != 0.0) {
						//	update_image(image, idx, pix_weight, bin_q);
						//}

						float d = pix_weight * bin_q;

#pragma omp atomic
						grid.data_density[idx] += d;
					}
				}
			}
			//}
#endif
			int iradiusx = static_cast<int>(hsml);
			int iradiusy = static_cast<int>(hsml);
			int iradiusz = static_cast<int>(hsml);

			///////////////////////////////////////////////////////////////////////////////
			double dpx = ((double)pos[0] - (double)bbox_min_orig[0]) * len2pix;
			double dpy = ((double)pos[1] - (double)bbox_min_orig[1]) * len2pix;
			double dpz = ((double)pos[2] - (double)bbox_min_orig[2]) * len2pix;

			int px = static_cast<int>(dpx);
			int py = static_cast<int>(dpy);
			int pz = static_cast<int>(dpz);

			///////////////////////////////////////////////////////////////////////////////

			for (int sx = px - iradiusx; sx <= px + iradiusx; sx++) {
				for (int sy = py - iradiusy; sy <= py + iradiusy; sy++) {
					for (int sz = pz - iradiusz; sz <= pz + iradiusz; sz++) {

						int osx = sx - grid.offset[0];
						int osy = sy - grid.offset[1];
						int osz = sz - grid.offset[2];

						if (osx < 0 || osy < 0 || osz < 0)
							continue;

						if (osx > grid.x() - 1 || osy > grid.y() - 1 || osz > grid.z() - 1)
							continue;

						double density = 1.0;

						double dx = dpx - sx;
						double dy = dpy - sy;
						double dz = dpz - sz;
						double distance_norm = std::sqrt(dx * dx + dy * dy + dz * dz);

						double W = 0.0;
						double h = hsml;

						if (iradiusx == 0 && iradiusy == 0 && iradiusz == 0) {
							W = 1.0; // full value
						}
						else {
							double q = distance_norm / h;

							W = utility::dense::sph_kernel::W(dense_type, q, 1.0 / h);
						}

						//if (mass != 0.0 && rho != 0) {
						//	density = (mass / (rho * std::pow(h, 3))) * W;
						//}
						//else {
						//	density = (1.0 / std::pow(h, 3)) * W;
						//}

						density = W;

						//final density
						float d = density * value;

						double norm = 0.0;
						if (dense_norm == common::SpaceData::DenseNorm::eNone) {
							norm = 0.0;
						}
						else if (dense_norm == common::SpaceData::DenseNorm::eCount) {
							norm = 1.0; // count
						}
						else if (dense_norm == common::SpaceData::DenseNorm::eSPHInterpolation) {
							norm = density;
						}

						//final norm
						float n = norm;

						size_t gindex = grid.get_index(osx, osy, osz);

#pragma omp atomic
						grid.data_density[gindex] += d;

#ifndef WITH_NO_DATA_TEMP
#pragma omp atomic
						grid.data_temp[gindex] += n;
#endif
					}
				}
			}
		}

	}// dense
} //common