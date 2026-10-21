#!/bin/bash

#####################################################################################################################
# Copyright(C) 2023-2026 IT4Innovations National Supercomputing Center, VSB - Technical University of Ostrava
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
#
# HIP (AMD MI250X, gfx90a) build of space_converter for LUMI-G.
#

cd /flash/project_465002608/jaromila/space

ml purge

ml LUMI/25.09 partition/D
ml cuda/12.6
ml PrgEnv-gnu

ROOT_DIR=${PWD}

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_lumid
src=${ROOT_DIR}/src

############DEPENDENCIES###############
# cd ${lib_dir}
# git-lfs clone -b blender-v4.5-release  https://projects.blender.org/blender/lib-linux_x64.git
#######################################

# Cray compiler wrappers (PrgEnv-gnu): they add Cray MPICH and, with
# partition/G (craype-accel-amd-gfx90a), also the GTL library needed for
# GPU-aware MPI (WITH_HIP_AWARE_MPI + MPICH_GPU_SUPPORT_ENABLED=1 at runtime).
export CC=cc
export CXX=CC

#rm -rf ${ROOT_DIR}/build/space_converter_lumid

#-----------space_converter_lumid--------------
mkdir -p ${ROOT_DIR}/build/space_converter_lumid
cd ${ROOT_DIR}/build/space_converter_lumid

make_d="${src}/space-converter"

make_d="${make_d} -DTBB_INCLUDE_DIRS=$lib_dir/lib-linux_x64/tbb/include"
make_d="${make_d} -DOPENVDB_INCLUDE_DIRS=$lib_dir/lib-linux_x64/openvdb/include"
make_d="${make_d} -DOPENVDB_LIBRARIES=$lib_dir/lib-linux_x64/openvdb/lib/libopenvdb.so;$lib_dir/lib-linux_x64/tbb/lib/libtbb.so"

make_d="${make_d} -DWITH_HDF5=OFF"
make_d="${make_d} -DGADGET_READ_ID=OFF"
make_d="${make_d} -DGADGET_MAX_HSML=OFF"

make_d="${make_d} -DWITH_OPENVDB=OFF"

make_d="${make_d} -DWITH_CUDAKDTREE=OFF"

make_d="${make_d} -DWITH_NANOFLANN=OFF"
make_d="${make_d} -DWITH_HACC=OFF"

#----------------HIP (AMD GPU)----------------
make_d="${make_d} -DWITH_GPU_HIP=OFF"
make_d="${make_d} -DWITH_HIP_AWARE_MPI=OFF"
make_d="${make_d} -DWITH_HIP_MALLOCMANAGED=OFF"
# NOTE: hipMallocManaged on MI250X maps to slow fine-grained memory - keep OFF
# for performance unless the data does not fit into the 64 GB of a GCD.

# make_d="${make_d} -DCMAKE_HIP_ARCHITECTURES=gfx90a"
# make_d="${make_d} -DCMAKE_HIP_COMPILER=${ROCM_PATH}/bin/amdclang++"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 16 install
