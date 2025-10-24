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

# module use $PWD/easybuild_bar/modules/all/

# ml c-blosc/1.21.0-GCC-10.3.0
# ml tbb/2020.3-GCCcore-10.2.0
ml CMake/3.24.3-GCCcore-12.2.0
# ml Boost/1.81.0-GCC-12.2.0
ml intel/2022b

ROOT_DIR=${PWD}

echo "===================RUN"

lib_dir=${ROOT_DIR}/install
install=${ROOT_DIR}/install/space_converter_bar
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

# export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${lib_dir}/openvdb_bar/lib64
export LD_LIBRARY_PATH=$lib_dir/lib-linux_x64/openvdb/lib:$lib_dir/lib-linux_x64/tbb/lib:$LD_LIBRARY_PATH

# Start the timer
start_time=$(date +%s)

out=${ROOT_DIR}/out
out_t=${out}/t
mkdir -p ${out_t}

snapdir=0
snap=0

snapdir_format=`printf "%03d" $snapdir`
dataset=snapdir_${snapdir_format}/snap_${snapdir_format}.${snap}

out_d=${out_t}/${dataset}
mkdir -p ${out_d}

srun -n 1 ${install}/bin/space_converter --data-type GADGET --gadget-file ${dataset} --max-mem-size 50000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 1000 --output-path ${out_d} --export-data 1 1 ${out_d} >> ${out_d}/log.log 2> ${out_d}/log.err

# End the timer
end_time=$(date +%s)

# Calculate and print the elapsed time
elapsed_time=$((end_time - start_time))
echo "$ Total time: ${elapsed_time} seconds"