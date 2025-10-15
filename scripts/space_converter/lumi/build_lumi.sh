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

cd /users/jaromila/scratch/space/OpenGadget

ml purge 

ml LUMI/23.09  partition/C
ml PrgEnv-gnu

ml Blosc
ml Boost

ROOT_DIR=${PWD}

lib_dir=${ROOT_DIR}/install
output=${ROOT_DIR}/install/space_converter_lumi
src=${ROOT_DIR}/src

#rm -rf build/space_converter_lumi

#-----------space_converter_lumi--------------
mkdir ${ROOT_DIR}/build/space_converter_lumi
cd ${ROOT_DIR}/build/space_converter_lumi

make_d="${src}/space_converter-public"

make_d="${make_d} -DTBB_INCLUDE_DIRS=${src}/oneTBB/install/include" 
make_d="${make_d} -DTBB_LIBRARIES="

make_d="${make_d} -DBOOST_INCLUDE_DIRS=/appl/lumi/SW/LUMI-23.09/C/EB/Boost/1.82.0-cpeGNU-23.09/include"
make_d="${make_d} -DBOOST_LIBRARIES="

make_d="${make_d} -DNANOVDB_INCLUDE_DIRS="

make_d="${make_d} -DOPENVDB_INCLUDE_DIRS=${lib_dir}/openvdb/include"
make_d="${make_d} -DOPENVDB_LIBRARIES=${lib_dir}/openvdb/lib64/libopenvdb.so;${src}/oneTBB/install/lib64/libtbb.so;/appl/lumi/SW/LUMI-23.09/C/EB/zlib/1.2.13-cpeGNU-23.09/lib/libz.so;/appl/lumi/SW/LUMI-23.09/C/EB/Blosc/1.21.5-cpeGNU-23.09/lib/libblosc.so;/appl/lumi/SW/LUMI-23.09/C/EB/Boost/1.82.0-cpeGNU-23.09/lib/libboost_iostreams-mt-x64.so"

#make_d="${make_d} -DCMAKE_BUILD_TYPE=Debug"
#make_d="${make_d} -DCMAKE_BUILD_TYPE=Release"
make_d="${make_d} -DCMAKE_BUILD_TYPE=RelWithDebInfo"
make_d="${make_d} -DCMAKE_INSTALL_PREFIX=${output}"

cmake ${make_d}
#make clean
make -j 16 install
