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


ROOT_DIR=${PWD}

source ~/.zprofile

#################################################

echo "===================RUN"

lib_dir=${ROOT_DIR}/install

install=${ROOT_DIR}/install/space_converter_macos
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

out=${ROOT_DIR}/out
out_t=${out}/remote
mkdir -p ${out_t}

#################################################
CYCLES_DIR_LIB=${PWD}/src/cyclesphi/lib/macos_arm64/
export DYLD_LIBRARY_PATH=${CYCLES_DIR_LIB}/openvdb/lib:${CYCLES_DIR_LIB}/tbb/lib:${lib_dir}/openmp/release/usr/local/lib:${CYCLES_DIR_LIB}/boost/lib:${DYLD_LIBRARY_PATH}
#################################################

#mpirun -n 4 ${install}/bin/space_converter --data-type GADGET --gadget-file ${data}/very_small_example/snap_081 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --bh-count 1 --grid-dim 100 --output-path ${out} --info 
# mpirun -n 4 
#${install}/bin/space_converter --data-type GADGET --gadget-file ${data}/very_small_example/snap_081 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --bh-count 1 --grid-dim 100 --output-path ${out} --port 5001 --multires

${install}/bin/space_converter --data-type GADGET --gadget-file ${data}/notsosmall_example/snap_091 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 100 --output-path ${out} --port 5001 --bh-count 110 --multires

#mpirun -n 8 ${install}/bin/space_converter --data-type CHANGA_TIPSY --tipsy-file ${data}/changa_tipsy/LOW_64XLARGEMHDVERTDENSWend64FBSB1mergb0.00410  --grid-dim 100 --output-path ${out} --port 8000
