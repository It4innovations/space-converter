#!/bin/bash

#####################################################################################################################
# Copyright(C) 2023-2025 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
#
# This program is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#
#####################################################################################################################


module use $PWD/easybuild_bar/modules/all/

ml c-blosc/1.21.0-GCC-10.3.0
ml tbb/2020.3-GCCcore-10.2.0
ml CMake/3.24.3-GCCcore-12.2.0
ml Boost/1.81.0-GCC-12.2.0
ml intel/2022b

ROOT_DIR=${PWD}

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_bar
src=${ROOT_DIR}/src

#rm -rf build/space_converter_bar
#-----------space_converter_bar--------------
mkdir ${ROOT_DIR}/build/space_converter_bar
cd ${ROOT_DIR}/build/space_converter_bar

make_d="${src}/space_converter"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 24 install
