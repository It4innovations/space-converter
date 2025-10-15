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


cd /mnt/proj3/open-28-64/blender/he_space

ROOT_DIR=${PWD}

# ml HDF5/1.14.0-gompi-2022b

# ml c-blosc/1.21.0-GCC-10.2.0
# ml Boost/1.79.0-GCC-11.3.0
# ml tbb/2021.10.0-GCCcore-12.2.0
# ml VirtualGL/3.1-GCC-12.3.0

# ml Qt5/5.15.7-GCCcore-12.2.0
# ml CMake/3.24.3-GCCcore-12.2.0
# ml GCC/12.2.0
# ml Mesa/22.2.4-GCCcore-12.2.0
# # ml CUDA/12.4.0
# # ml OpenMPI/4.1.4-GCC-12.2.0
#  ml OpenMPI/4.1.6-NVHPC-24.1-CUDA-12.4.0
#ml CMake/3.23.1-GCCcore-11.3.0
ml OpenMPI/4.1.6-NVHPC-23.11-CUDA-12.2.0
#ml GCC
######################################################

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_cs_p06
src=${ROOT_DIR}/src

export PATH=/mnt/proj3/open-28-64/blender/he_space/install/cmake_p06/cmake-3.28.3/install/bin:$PATH

# export CC='mpicc'
# export CXX='mpic++'

#rm -rf build/space_converter_cs_p06

#-----------space_converter--------------
mkdir ${ROOT_DIR}/build/space_converter_cs_p06
cd ${ROOT_DIR}/build/space_converter_cs_p06

make_d="${src}/space_converter"

# make_d="${make_d} -DTBB_INCLUDE_DIRS=/apps/all/tbb/2021.10.0-GCCcore-12.2.0/include"
# make_d="${make_d} -DOPENVDB_INCLUDE_DIRS=$lib_dir/openvdb_kar/include"
# make_d="${make_d} -DOPENVDB_LIBRARIES=$lib_dir/openvdb_kar/lib64/libopenvdb.so;/apps/all/tbb/2021.10.0-GCCcore-12.2.0/lib64/libtbb.so"
# make_d="${make_d} -DOPENVDB_VERSION=11"

# make_d="${make_d} -DZSTD_LIBRARIES=/apps/all/zstd/1.5.2-GCCcore-12.2.0/lib64/libzstd.so"
# make_d="${make_d} -DBLOSC_LIBRARIES=/apps/all/c-blosc/1.21.0-GCC-10.2.0/lib64/libblosc.so"
# make_d="${make_d} -DZLIB_LIBRARIES=/apps/all/zlib/1.2.12-GCCcore-12.2.0/lib64/libz.so"

make_d="${make_d} -DWITH_HDF5=OFF"
make_d="${make_d} -DGADGET_READ_ID=ON"
make_d="${make_d} -DGADGET_MAX_HSML=ON"

make_d="${make_d} -DWITH_OPENVDB=OFF"
make_d="${make_d} -DWITH_CUDAKDTREE=ON"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 32 install
