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

#################################################

echo "===================RUN"

lib_dir=${ROOT_DIR}/install

install=${ROOT_DIR}/install/space-converter_aurora
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

export LD_LIBRARY_PATH=$lib_dir/lib-linux_x64/openvdb/lib:$lib_dir/lib-linux_x64/tbb/lib:$LD_LIBRARY_PATH

out=${ROOT_DIR}/out
out_t=${out}/remote
mkdir -p ${out_t}

mpirun --np 64 -env OMP_NUM_THREADS=2 --line-buffer ${install}/bin/space_converter --pos-names x y z --vel-names vx vy vz --mass-name mass --rho-name rho --hsml-name hh --data-type GENERICIO --genericio-file HACC/Farpoint/M000p/L1478/HACC000/analysis/Particles/STEP499/m000p-499.select.mpicosmo --output-path ${out_t} --port 5000 # --info #--port 5000 #--calc-radius-neigh 10 #--info
