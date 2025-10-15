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

#qsub -I -l select=1 -l filesystems=flare -l walltime=1:00:00 -q debug -A XXX

ROOT_DIR=${PWD}

module restore
#ml
ml cmake

######################################################

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_aurora
src=${ROOT_DIR}/src


#export CMAKE_ROOT=$lib_dir/cmake
#export PATH=$lib_dir/cmake/bin:$PATH
#export LD_LIBRARY_PATH=$lib_dir/lib:$lib_dir/blosc/lib:$lib_dir/openvdb/lib64:$lib_dir/tbb/lib64:$lib_dir/boost/lib64:$LD_LIBRARY_PATH
#export CUDAToolkit_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/compilers
#export CUDACXX=/opt/nvidia/hpc_sdk/Linux_x86_64/23.9/compilers/bin/nvcc

#export CC='mpicc'
#export CXX='mpic++'

export CC=gcc
export CXX=g++

#rm -rf build/space_converter_aurora

#-----------space_converter--------------
mkdir ${ROOT_DIR}/build/space_converter_aurora
cd ${ROOT_DIR}/build/space_converter_aurora

make_d="${src}/space_converter"

make_d="${make_d} -DTBB_INCLUDE_DIRS=$lib_dir/tbb/include"
make_d="${make_d} -DOPENVDB_INCLUDE_DIRS=$lib_dir/openvdb/include"
make_d="${make_d} -DOPENVDB_LIBRARIES=$lib_dir/openvdb/lib64/libopenvdb.so;$lib_dir/tbb/lib64/libtbb.so;$lib_dir/blosc/lib/libblosc.so.1"
make_d="${make_d} -DOPENVDB_VERSION=11"

make_d="${make_d} -DBOOST_INCLUDE_DIRS=$lib_dir/boost/include"
make_d="${make_d} -DBOOST_LIBRARIES=$lib_dir/boost/lib64/libboost_iostreams.so"

make_d="${make_d} -DBLOSC_INCLUDE_DIRS=$lib_dir/blosc/include"
make_d="${make_d} -DBLOSC_LIBRARIES=$lib_dir/blosc/lib/libblosc.so"

make_d="${make_d} -DZSTD_LIBRARIES=$lib_dir/lib/libzstd.so"
make_d="${make_d} -DZLIB_LIBRARIES=$lib_dir/lib/libz.so"

# make_d="${make_d} -DWITH_HDF5=OFF"
# make_d="${make_d} -DGADGET_READ_ID=OFF"
# make_d="${make_d} -DGADGET_MAX_HSML=ON"

# #make_d="${make_d} -DMPI_CUSTOM_LIBRARIES=" #/opt/cray/pe/mpich/8.1.28/gtl/lib/libmpi_gtl_cuda.so"

# make_d="${make_d} -DWITH_OPENVDB=OFF"
# make_d="${make_d} -DWITH_CUDAKDTREE=OFF"

make_d="${make_d} -DWITH_HDF5=OFF"
make_d="${make_d} -DGADGET_READ_ID=OFF"
make_d="${make_d} -DGADGET_MAX_HSML=ON"

make_d="${make_d} -DWITH_OPENVDB=ON"

make_d="${make_d} -DWITH_CUDAKDTREE=OFF"
make_d="${make_d} -DWITH_CUDAKDTREE_CPU=OFF"

make_d="${make_d} -DWITH_NANOFLANN=OFF"
make_d="${make_d} -DWITH_GENERICIO=ON"
make_d="${make_d} -DWITH_NO_DATA_TEMP=ON"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 32 install
