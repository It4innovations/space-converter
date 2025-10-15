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

# Prereq:
## /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
## brew install gcc
## brew install open-mpi
## brew install git
## brew install git-lfs
## brew install cmake

source ~/.zprofile

## git clone https://code.it4i.cz/raas/cyclesphi.git
## cd cyclesphi
## make update


ROOT_DIR=${PWD}

CYCLES_DIR_LIB=${PWD}/src/cyclesphi/lib/macos_arm64/
######################################################

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_macos
src=${ROOT_DIR}/src

#rm -rf build/space_converter_macos

#-----------space_converter--------------
mkdir ${ROOT_DIR}/build/space_converter_macos
cd ${ROOT_DIR}/build/space_converter_macos

make_d="${src}/space_converter"

make_d="${make_d} -DTBB_INCLUDE_DIRS=${CYCLES_DIR_LIB}/tbb/include"
make_d="${make_d} -DOPENVDB_INCLUDE_DIRS=${CYCLES_DIR_LIB}/openvdb/include"
make_d="${make_d} -DOPENVDB_LIBRARIES=${CYCLES_DIR_LIB}/openvdb/lib/libopenvdb.dylib;${CYCLES_DIR_LIB}/tbb/lib/libtbb.dylib"
make_d="${make_d} -DOPENVDB_VERSION=11"

make_d="${make_d} -DZSTD_LIBRARIES=${CYCLES_DIR_LIB}/zstd/lib/libzstd.a"
#make_d="${make_d} -DBLOSC_LIBRARIES=${CYCLES_DIR_LIB}/c-blosc/1.21.0-GCC-10.2.0/lib64/libblosc.so"
#make_d="${make_d} -DZLIB_LIBRARIES=${CYCLES_DIR_LIB}/zlib/1.2.12-GCCcore-12.2.0/lib64/libz.so"

make_d="${make_d} -DBOOST_INCLUDE_DIRS=${CYCLES_DIR_LIB}/boost/include"
#make_d="${make_d} -DBOOST_LIBRARIES=${CYCLES_DIR_LIB}/boost/lib/libzstd.a"

make_d="${make_d} -DGADGET_READ_ID=ON"
make_d="${make_d} -DGADGET_MAX_HSML=ON"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake \
    -DCMAKE_C_FLAGS="-Xclang -fopenmp" \
    -DCMAKE_CXX_FLAGS="-Xclang -fopenmp" \
    -DOpenMP_C_FLAGS="-I$lib_dir/openmp/release/usr/local/include" \
    -DOpenMP_CXX_FLAGS="-I$lib_dir/openmp/release/usr/local/include" \
    -DOpenMP_C_LIB_NAMES="omp" \
    -DOpenMP_CXX_LIB_NAMES="omp" \
    -DOpenMP_omp_LIBRARY="$lib_dir/openmp/release/usr/local/lib/libomp.dylib" \
    ${make_d}

#make clean
make -j 1 install
