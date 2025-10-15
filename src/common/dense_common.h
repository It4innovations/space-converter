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

#include "convert_vdb.h"

namespace common {
	namespace vdb {
		struct particle_sim_v2
		{
			float e; //value
			float x, y, z, r, I;
		};

		double custom_asinh(double val);
		//float custom_xexp(float v);

        template<typename T> 
        class exptable
        {
        private:
            T expfac, taylorlimit;
            std::vector<T> tab1, tab2;
            enum {
                nbits = 10,
                dim1 = 1 << nbits,
                mask1 = dim1 - 1,
                dim2 = (1 << nbits) << nbits,
                mask2 = dim2 - 1,
                mask3 = ~mask2
            };

        public:
            exptable(T maxexp)
                : expfac(dim2 / maxexp)
            {
                tab1.resize(dim1);                    
                tab2.resize(dim1);

                using namespace std;
                for (int m = 0; m < dim1; ++m)
                {
                    tab1[m] = exp(m * dim1 / expfac);
                    tab2[m] = exp(m / expfac);
                }
                taylorlimit = sqrt(T(2) * abs(maxexp) / dim2);
            }

            T taylorLimit() const { return taylorlimit; }

            T operator() (T arg) const
            {
                int iarg = int(arg * expfac);
                if (iarg & mask3)
                    return (iarg < 0) ? T(1) : T(0);
                return tab1[iarg >> nbits] * tab2[iarg & mask1];
            }
            T expm1(T arg) const
            {
                if (std::abs(arg) < taylorlimit) 
                    return arg;

                return operator()(arg) - T(1);
            }
        };
	}// vdb
} //common