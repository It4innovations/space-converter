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

echo "===================RUN"

lib_dir=${ROOT_DIR}/install
install=${ROOT_DIR}/install/space_converter_lumi
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

#${lib_dir}/openvdb/lib64/libopenvdb.so;${src}/oneTBB/install/lib64/libtbb.so;/appl/lumi/SW/LUMI-23.09/C/EB/zlib/1.2.13-cpeGNU-23.09/lib/libz.so;/appl/lumi/SW/LUMI-23.09/C/EB/Blosc/1.21.5-cpeGNU-23.09/lib/libblosc.so;/appl/lumi/SW/LUMI-23.09/C/EB/Boost/1.82.0-cpeGNU-23.09/lib/libboost_iostreams-mt-x64.so

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${lib_dir}/openvdb/lib64/:${src}/oneTBB/install/lib64

#export example_snap="--data-type CHANGA_TIPSY --grid-dim 1000 --output-path ${out} --tipsy-file $HOME/scratch/changa/LOW_512XLARGEMHDVERTDENSWend64FBSB64AB05.03280 --grid-dim 100 --output-path f:\temp\ --port 5000"
export example_snap="--data-type CHANGA_TIPSY --grid-dim 1000 --output-path ${out} --tipsy-file $HOME/scratch/changa/testgalaxyN32-mpi.smp.gcc.ompi-meluxina-p320000-b24-N32-.000050 --grid-dim 100 --output-path f:\temp\ --port 5000"
srun --export=ALL -u -n 512 -c 64 ${install}/bin/space_converter ${example_snap}