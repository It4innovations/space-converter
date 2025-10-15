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
 * Source: SPHtoGrid: https://github.com/LudwigBoess/SPHtoGrid.jl, SPHKernels: https://github.com/LudwigBoess/SPHKernels.jl
 */

#pragma once

#include <cmath>
#include <stdexcept>
#include <cstdint>

#include <vector>
#include <tuple>
#include <utility>
#include "data_common.h"

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define CUDA_CALLABLE __host__ __device__
#else
#define CUDA_CALLABLE
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace utility {
    namespace dense {
        namespace sph_kernel {

            class AbstractSPHKernel {
            public:
                virtual ~AbstractSPHKernel() = default;

                virtual CUDA_CALLABLE double kernel_norm(double h_inv) const = 0;
                virtual CUDA_CALLABLE double kernel_value(double u) const = 0;
                virtual CUDA_CALLABLE double kernel_deriv(double u) const = 0;
                virtual CUDA_CALLABLE double bias_correction(double density, double m, double h_inv, int n_neighbours) const = 0;

                CUDA_CALLABLE double W(double u, double h_inv) const {
                    return kernel_norm(h_inv) * kernel_value(u);
                }
            };

            // --------------------------------------
            // Cubic B-spline kernel (double only)
            // --------------------------------------

            class Cubic : public AbstractSPHKernel {
            public:
                std::int8_t dim;
                double norm;

                CUDA_CALLABLE Cubic() : dim(3), norm(8.0 / M_PI) {}

                CUDA_CALLABLE Cubic(int dimension) : dim(static_cast<std::int8_t>(dimension)) {
                    if (dimension == 1) {
                        norm = 4.0 / 3.0;
                    }
                    else if (dimension == 2) {
                        norm = 40.0 / (7.0 * M_PI);
                    }
                    else if (dimension == 3) {
                        norm = 8.0 / M_PI;
                    }
                    else {
#ifndef __CUDA_ARCH__
                        throw std::invalid_argument("Cubic not defined for dimension!");
#endif
                    }
                }

                CUDA_CALLABLE Cubic(std::int8_t dimension, double normalization)
                    : dim(dimension), norm(normalization) {
                }

                CUDA_CALLABLE double kernel_norm(double h_inv) const override {
#ifdef __CUDA_ARCH__
                    return norm * powf(h_inv, static_cast<float>(dim));
#else
                    return norm * std::pow(h_inv, static_cast<double>(dim));
#endif
                }

                CUDA_CALLABLE double kernel_value(double u) const override {
                    if (u < 0.5) {
                        return 1.0 + 6.0 * (u - 1.0) * u * u;
                    }
                    else if (u < 1.0) {
                        double t1 = 1.0 - u;
                        return 2.0 * t1 * t1 * t1;
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double kernel_deriv(double u) const override {
                    if (u < 0.5) {
                        return u * (18.0 * u - 12.0);
                    }
                    else if (u < 1.0) {
                        double t1 = 1.0 - u;
                        return -6.0 * t1 * t1;
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double bias_correction(double density, double, double, int) const override {
                    return density;
                }
            };

            // --------------------------------------
            // Quintic B-spline kernel (double only)
            // --------------------------------------

            class Quintic : public AbstractSPHKernel {
            public:
                std::int8_t dim;
                double norm;

                CUDA_CALLABLE Quintic() : dim(3), norm(2187.0 / (40.0 * M_PI)) {}

                CUDA_CALLABLE Quintic(int dimension) : dim(static_cast<std::int8_t>(dimension)) {
                    if (dimension == 1) {
                        norm = 243.0 / 40.0;
                    }
                    else if (dimension == 2) {
                        norm = 15309.0 / (478.0 * M_PI);
                    }
                    else if (dimension == 3) {
                        norm = 2187.0 / (40.0 * M_PI);
                    }
                    else {
#ifndef __CUDA_ARCH__
                        throw std::invalid_argument("Quintic not defined for dimension!");
#endif
                    }
                }

                CUDA_CALLABLE Quintic(std::int8_t dimension, double normalization)
                    : dim(dimension), norm(normalization) {
                }

                CUDA_CALLABLE double kernel_norm(double h_inv) const override {
#ifdef __CUDA_ARCH__
                    return norm * powf(h_inv, static_cast<float>(dim));
#else
                    return norm * std::pow(h_inv, static_cast<double>(dim));
#endif
                }

                CUDA_CALLABLE double kernel_value(double u) const override {
                    if (u < 1.0 / 3.0) {
                        double u_m1 = 1.0 - u;
                        double u_m23 = (2.0 / 3.0) - u;
                        double u_m13 = (1.0 / 3.0) - u;

                        double u_m1_5 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1;
                        double u_m23_5 = u_m23 * u_m23 * u_m23 * u_m23 * u_m23;
                        double u_m13_5 = u_m13 * u_m13 * u_m13 * u_m13 * u_m13;

                        return u_m1_5 - 6.0 * u_m23_5 + 15.0 * u_m13_5;
                    }
                    else if (u < 2.0 / 3.0) {
                        double u_m1 = 1.0 - u;
                        double u_m23 = (2.0 / 3.0) - u;

                        double u_m1_5 = u_m1 * u_m1 * u_m1 * u_m1 * u_m1;
                        double u_m23_5 = u_m23 * u_m23 * u_m23 * u_m23 * u_m23;

                        return u_m1_5 - 6.0 * u_m23_5;
                    }
                    else if (u < 1.0) {
                        double u_m1 = 1.0 - u;
                        return u_m1 * u_m1 * u_m1 * u_m1 * u_m1;
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double kernel_deriv(double u) const override {
                    if (u < 1.0 / 3.0) {
                        double u_m1 = 1.0 - u;
                        double u_m23 = (2.0 / 3.0) - u;
                        double u_m13 = (1.0 / 3.0) - u;

                        double u_m1_4 = u_m1 * u_m1 * u_m1 * u_m1;
                        double u_m23_4 = u_m23 * u_m23 * u_m23 * u_m23;
                        double u_m13_4 = u_m13 * u_m13 * u_m13 * u_m13;

                        return -5.0 * u_m1_4 + 30.0 * u_m23_4 - 75.0 * u_m13_4;
                    }
                    else if (u < 2.0 / 3.0) {
                        double u_m1 = 1.0 - u;
                        double u_m23 = (2.0 / 3.0) - u;

                        double u_m1_4 = u_m1 * u_m1 * u_m1 * u_m1;
                        double u_m23_4 = u_m23 * u_m23 * u_m23 * u_m23;

                        return -5.0 * u_m1_4 + 30.0 * u_m23_4;
                    }
                    else if (u < 1.0) {
                        double u_m1 = 1.0 - u;
                        double u_m1_4 = u_m1 * u_m1 * u_m1 * u_m1;
                        return -5.0 * u_m1_4;
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double bias_correction(double density, double, double, int) const override {
                    return density;
                }
            };

            // --------------------------------------
            // Wendland C2 kernel (double only)
            // --------------------------------------

            class WendlandC2 : public AbstractSPHKernel {
            public:
                std::int8_t dim;
                double norm;

                CUDA_CALLABLE WendlandC2() : dim(3), norm(21.0 / (2.0 * M_PI)) {}

                CUDA_CALLABLE WendlandC2(int dimension) : dim(static_cast<std::int8_t>(dimension)) {
                    if (dimension == 2) {
                        norm = 7.0 / M_PI;
                    }
                    else if (dimension == 3) {
                        norm = 21.0 / (2.0 * M_PI);
                    }
                    else {
#ifndef __CUDA_ARCH__
                        throw std::invalid_argument("WendlandC2 not defined for dimension!");
#endif
                    }
                }

                CUDA_CALLABLE WendlandC2(std::int8_t dimension, double normalization)
                    : dim(dimension), norm(normalization) {
                }

                CUDA_CALLABLE double kernel_norm(double h_inv) const override {
#ifdef __CUDA_ARCH__
                    return norm * powf(h_inv, static_cast<float>(dim));
#else
                    return norm * std::pow(h_inv, static_cast<double>(dim));
#endif
                }

                CUDA_CALLABLE double kernel_value(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        t1 = t1 * t1;  // (1-u)^2
                        t1 = t1 * t1;  // (1-u)^4
                        return t1 * (1.0 + 4.0 * u);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double kernel_deriv(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        return -20.0 * u * t1 * t1 * t1;
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double bias_correction(double density, double m, double h_inv, int n_neighbours) const override {
                    double n = kernel_norm(h_inv);

#ifdef __CUDA_ARCH__
                    double wc_correction = 0.0294 *
                        powf(n_neighbours * 0.01f, -0.977f) *
                        m * n;
#else
                    double wc_correction = 0.0294 *
                        std::pow(n_neighbours * 0.01, -0.977) *
                        m * n;
#endif

                    return density - wc_correction;
                }
            };

            // --------------------------------------
            // Wendland C4 kernel (double only)
            // --------------------------------------

            class WendlandC4 : public AbstractSPHKernel {
            public:
                std::int8_t dim;
                double norm;

                CUDA_CALLABLE WendlandC4() : dim(3), norm(495.0 / (32.0 * M_PI)) {}

                CUDA_CALLABLE WendlandC4(int dimension) : dim(static_cast<std::int8_t>(dimension)) {
                    if (dimension == 2) {
                        norm = 9.0 / M_PI;
                    }
                    else if (dimension == 3) {
                        norm = 495.0 / (32.0 * M_PI);
                    }
                    else {
#ifndef __CUDA_ARCH__
                        throw std::invalid_argument("WendlandC4 not defined for dimension!");
#endif
                    }
                }

                CUDA_CALLABLE WendlandC4(std::int8_t dimension, double normalization)
                    : dim(dimension), norm(normalization) {
                }

                CUDA_CALLABLE double kernel_norm(double h_inv) const override {
#ifdef __CUDA_ARCH__
                    return norm * powf(h_inv, static_cast<float>(dim));
#else
                    return norm * std::pow(h_inv, static_cast<double>(dim));
#endif
                }

                CUDA_CALLABLE double kernel_value(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        t1 = t1 * t1;         // (1-u)^2
                        double t4 = t1 * t1;  // (1-u)^4
                        double u2 = u * u;
                        return t1 * t4 * (1.0 + 6.0 * u + (35.0 / 3.0) * u2);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double kernel_deriv(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        double t5 = t1 * t1;  // (1-u)^2
                        t5 = t5 * t5;         // (1-u)^4
                        t5 = t5 * t1;         // (1-u)^5
                        double u2 = u * u;
                        return (-(288.0 / 3.0) * t5 * u2 - (56.0 / 3.0) * u * t5);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double bias_correction(double density, double m, double h_inv, int n_neighbours) const override {
                    double n = kernel_norm(h_inv);

#ifdef __CUDA_ARCH__
                    double wc_correction = 0.01342 *
                        powf(n_neighbours * 0.01f, -1.579f) *
                        m * n;
#else
                    double wc_correction = 0.01342 *
                        std::pow(n_neighbours * 0.01, -1.579) *
                        m * n;
#endif

                    return density - wc_correction;
                }
            };

            // --------------------------------------
            // Wendland C6 kernel (double only)
            // --------------------------------------

            class WendlandC6 : public AbstractSPHKernel {
            public:
                std::int8_t dim;
                double norm;

                CUDA_CALLABLE WendlandC6() : dim(3), norm(1365.0 / (64.0 * M_PI)) {}

                CUDA_CALLABLE WendlandC6(int dimension) : dim(static_cast<std::int8_t>(dimension)) {
                    if (dimension == 2) {
                        norm = 78.0 / (7.0 * M_PI);
                    }
                    else if (dimension == 3) {
                        norm = 1365.0 / (64.0 * M_PI);
                    }
                    else {
#ifndef __CUDA_ARCH__
                        throw std::invalid_argument("WendlandC6 not defined for dimension!");
#endif
                    }
                }

                CUDA_CALLABLE WendlandC6(std::int8_t dimension, double normalization)
                    : dim(dimension), norm(normalization) {
                }

                CUDA_CALLABLE double kernel_norm(double h_inv) const override {
#ifdef __CUDA_ARCH__
                    return norm * powf(h_inv, static_cast<float>(dim));
#else
                    return norm * std::pow(h_inv, static_cast<double>(dim));
#endif
                }

                CUDA_CALLABLE double kernel_value(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        t1 = t1 * t1;  // (1-u)^2
                        t1 = t1 * t1;  // (1-u)^4
                        t1 = t1 * t1;  // (1-u)^8
                        double u2 = u * u;
                        return t1 * (1.0 + 8.0 * u + 25.0 * u2 + 32.0 * u2 * u);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double kernel_deriv(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        double t7 = t1 * t1;
                        t7 = t7 * t7 * t7 * t1;  // (1-u)^7
                        double u2 = u * u;
                        return -22.0 * t7 * u * (16.0 * u2 + 7.0 * u + 1.0);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double bias_correction(double density, double m, double h_inv, int n_neighbours) const override {
                    double n = kernel_norm(h_inv);

#ifdef __CUDA_ARCH__
                    double wc_correction = 0.0116 *
                        powf(n_neighbours * 0.01f, -2.236f) *
                        m * n;
#else
                    double wc_correction = 0.0116 *
                        std::pow(n_neighbours * 0.01, -2.236) *
                        m * n;
#endif

                    if (wc_correction < 0.2 * density) {
                        density -= wc_correction;
                    }

                    return density;
                }
            };

            // --------------------------------------
            // Wendland C8 kernel (double only)
            // --------------------------------------

            class WendlandC8 : public AbstractSPHKernel {
            public:
                std::int8_t dim;
                double norm;

                CUDA_CALLABLE WendlandC8() : dim(3), norm(357.0 / (64.0 * M_PI)) {}

                CUDA_CALLABLE WendlandC8(int dimension) : dim(static_cast<std::int8_t>(dimension)) {
                    if (dimension == 2) {
                        norm = 8.0 / (3.0 * M_PI);
                    }
                    else if (dimension == 3) {
                        norm = 357.0 / (64.0 * M_PI);
                    }
                    else {
#ifndef __CUDA_ARCH__
                        throw std::invalid_argument("WendlandC8 not defined for dimension!");
#endif
                    }
                }

                CUDA_CALLABLE WendlandC8(std::int8_t dimension, double normalization)
                    : dim(dimension), norm(normalization) {
                }

                CUDA_CALLABLE double kernel_norm(double h_inv) const override {
#ifdef __CUDA_ARCH__
                    return norm * powf(h_inv, static_cast<float>(dim));
#else
                    return norm * std::pow(h_inv, static_cast<double>(dim));
#endif
                }

                CUDA_CALLABLE double kernel_value(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        double t9 = t1 * t1;  // (1-u)^2
                        t9 = t9 * t9;         // (1-u)^4
                        t9 = t9 * t9;         // (1-u)^8
                        t9 = t9 * t1;         // (1-u)^9
                        double u2 = u * u;
                        double u3 = u2 * u;
                        double u4 = u2 * u2;
                        return t1 * t9 * (5.0 + 50.0 * u + 210.0 * u2 + 450.0 * u3 + 429.0 * u4);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double kernel_deriv(double u) const override {
                    if (u < 1.0) {
                        double t1 = 1.0 - u;
                        double t9 = t1 * t1;  // (1-u)^2
                        t9 = t9 * t9;         // (1-u)^4
                        t9 = t9 * t9;         // (1-u)^8
                        t9 = t9 * t1;         // (1-u)^9
                        double u2 = u * u;
                        double u3 = u2 * u;
                        return -26.0 * t9 * u * (231.0 * u3 + 159.0 * u2 + 45.0 * u + 5.0);
                    }
                    else {
                        return 0.0;
                    }
                }

                CUDA_CALLABLE double bias_correction(double density, double m, double h_inv, int n_neighbours) const override {
                    // No bias correction applied for WendlandC8
                    return density;
                }
            };

            inline double part_weight_physical(double pixelSideLength, double x_cgs = 3.085678e21) {
                return pixelSideLength * x_cgs;
            }

            inline CUDA_CALLABLE double W(common::SpaceData::DenseType &dense_type, double u, double h_inv) {
				double W = 0.0;
                if (dense_type == common::SpaceData::DenseType::eCubic) {
                    utility::dense::sph_kernel::Cubic kernel_cubic;
                    W = kernel_cubic.W(u, h_inv);
                }
                else if (dense_type == common::SpaceData::DenseType::eQuintic) {
                    utility::dense::sph_kernel::Quintic kernel_quintic;
                    W = kernel_quintic.W(u, h_inv);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC2) {
                    utility::dense::sph_kernel::WendlandC2 kernel_wendland;
                    W = kernel_wendland.W(u, h_inv);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC4) {
                    utility::dense::sph_kernel::WendlandC4 kernel_wendland;
                    W = kernel_wendland.W(u, h_inv);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC6) {
                    utility::dense::sph_kernel::WendlandC6 kernel_wendland;
                    W = kernel_wendland.W(u, h_inv);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC8) {
                    utility::dense::sph_kernel::WendlandC8 kernel_wendland;
                    W = kernel_wendland.W(u, h_inv);
                }
				return W;
            }

            inline CUDA_CALLABLE double bias_correction(common::SpaceData::DenseType& dense_type, double density, double m, double h_inv, int n_neighbours) {
                double bc = 0.0;
                if (dense_type == common::SpaceData::DenseType::eCubic) {
                    utility::dense::sph_kernel::Cubic kernel_cubic;
                    bc = kernel_cubic.bias_correction(density, m, h_inv, n_neighbours);
                }
                else if (dense_type == common::SpaceData::DenseType::eQuintic) {
                    utility::dense::sph_kernel::Quintic kernel_quintic;
                    bc = kernel_quintic.bias_correction(density, m, h_inv, n_neighbours);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC2) {
                    utility::dense::sph_kernel::WendlandC2 kernel_wendland;
                    bc = kernel_wendland.bias_correction(density, m, h_inv, n_neighbours);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC4) {
                    utility::dense::sph_kernel::WendlandC4 kernel_wendland;
                    bc = kernel_wendland.bias_correction(density, m, h_inv, n_neighbours);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC6) {
                    utility::dense::sph_kernel::WendlandC6 kernel_wendland;
                    bc = kernel_wendland.bias_correction(density, m, h_inv, n_neighbours);
                }
                else if (dense_type == common::SpaceData::DenseType::eWendlandC8) {
                    utility::dense::sph_kernel::WendlandC8 kernel_wendland;
                    bc = kernel_wendland.bias_correction(density, m, h_inv, n_neighbours);
                }
                return bc;
            }

        }  // namespace sph_kernel

    }  // namespace dense
}  // namespace utility

