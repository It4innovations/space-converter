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

#include "convert_common.h"

#include <iostream>
#include <string>
#include <vector>
#include <float.h>
#include <cmath>

const std::string WHITESPACE = " \n\r\t\f\v";

namespace common {

	double calculate_dmagnitude3(double x, double y, double z) {
		return (float)sqrt(x * x + y * y + z * z);
	}

	double calculate_dmagnituden(double* v, int n) {
		double res = 0;

		for (int i = 0; i < n; i++) {
			res += v[i] * v[i];
		}

		return (float)sqrt(res);
	}

	double calculate_fmagnituden(float* v, int n) {
		double res = 0;

		for (int i = 0; i < n; i++) {
			res += v[i] * v[i];
		}

		return (float)sqrt(res);
	}

	std::string ltrim(const std::string& s)
	{
		size_t start = s.find_first_not_of(WHITESPACE);
		return (start == std::string::npos) ? "" : s.substr(start);
	}

	std::string rtrim(const std::string& s)
	{
		size_t end = s.find_last_not_of(WHITESPACE);
		return (end == std::string::npos) ? "" : s.substr(0, end + 1);
	}

	std::string trim(const std::string& s) {
		return rtrim(ltrim(s));
	}

	///////////////////////
	// Gaussian kernel function
	// h is bandwidth for controlling the spread of the Gaussian function
	double gaussian_kernel(double x, double h) {
		return std::exp(-(x * x) / (2 * h * h)) / std::sqrt(2 * 3.14159265358979323846 * h * h);
	}

	void generate_normalized_gaussian(std::vector<double>& densities, double h) {
		// Generate N Gaussian-distributed values and compute their densities
		double total_sum = 0;
		for (int i = 0; i < densities.size(); i++) {
			densities[i] = gaussian_kernel((float)i, h);
			total_sum += densities[i];
		}

		// Scale the densities
		for (int i = 0; i < densities.size(); i++) {
			densities[i] /= total_sum;
		}

		//Check
		//double total_sum_one = std::accumulate(densities.begin(), densities.end(), 0.0);
	}
} //common
