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

############DEPENDENCIES###############
# cd ${lib_dir}
# git-lfs clone -b blender-v4.5-release  https://projects.blender.org/blender/lib-linux_x64.git
#######################################

#rm -rf build/space_converter
#-----------space_converter--------------
mkdir ${ROOT_DIR}/build/space_converter
cd ${ROOT_DIR}/build/space_converter

make_d="${src}/space-converter"

make_d="${make_d} -DTBB_INCLUDE_DIRS=$lib_dir/lib-linux_x64/tbb/include"
make_d="${make_d} -DOPENVDB_INCLUDE_DIRS=$lib_dir/lib-linux_x64/openvdb/include"
make_d="${make_d} -DOPENVDB_LIBRARIES=$lib_dir/lib-linux_x64/openvdb/lib/libopenvdb.so;$lib_dir/lib-linux_x64/tbb/lib/libtbb.so" #$lib_dir/lib-linux_x64/blosc/lib/libblosc.so.1
make_d="${make_d} -DOPENVDB_VERSION=12"

make_d="${make_d} -DWITH_HDF5=OFF"
make_d="${make_d} -DGADGET_READ_ID=OFF"
make_d="${make_d} -DGADGET_MAX_HSML=ON"

make_d="${make_d} -DWITH_OPENVDB=ON"

make_d="${make_d} -DWITH_CUDAKDTREE=ON"

make_d="${make_d} -DWITH_NANOFLANN=OFF"
make_d="${make_d} -DWITH_GENERICIO=ON"
make_d="${make_d} -DWITH_NO_DATA_TEMP=OFF"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 4 install
