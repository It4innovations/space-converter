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


cd ~/he_space

ROOT_DIR=${PWD}

# ml HDF5/1.14.0-gompi-2022b

# ml c-blosc/1.21.0-GCC-10.2.0
# ml Boost/1.79.0-GCC-11.3.0
# ml tbb/2021.10.0-GCCcore-12.2.0
# ml VirtualGL/3.1-GCC-12.3.0

# ml Qt5/5.15.7-GCCcore-12.2.0
ml CMake/3.24.3-GCCcore-12.2.0
ml GCC/12.2.0
ml Mesa/22.2.4-GCCcore-12.2.0
# ml CUDA/12.4.0
# ml OpenMPI/4.1.4-GCC-12.2.0
 ml OpenMPI/4.1.6-NVHPC-24.1-CUDA-12.4.0
 ml GCC/12.2.0

######################################################

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_kar
src=${ROOT_DIR}/src

############DEPENDENCIES###############
# cd ${lib_dir}
# git-lfs clone -b blender-v4.5-release  https://projects.blender.org/blender/lib-linux_x64.git
#######################################

#rm -rf build/space_converter_kar

#-----------space_converter--------------
mkdir ${ROOT_DIR}/build/space_converter_kar
cd ${ROOT_DIR}/build/space_converter_kar

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
make -j 32 install
