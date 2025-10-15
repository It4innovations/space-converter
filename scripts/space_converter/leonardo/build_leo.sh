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

ml purge

ml cmake
ml cuda/11.8
ml openmpi

ROOT_DIR=${PWD}

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter
src=${ROOT_DIR}/src

#rm -rf build/space_converter
#-----------space_converter--------------
mkdir ${ROOT_DIR}/build/space_converter
cd ${ROOT_DIR}/build/space_converter

make_d="${src}/space_converter"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 2 install
