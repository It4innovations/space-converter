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

#qsub -I -l select=10 -l filesystems=home:eagle -l walltime=1:00:00 -q debug-scaling -A XXX

ROOT_DIR=${PWD}

#ml purge 

#ml PrgEnv-gnu

#################################################

echo "===================RUN"

lib_dir=${ROOT_DIR}/install

install=${ROOT_DIR}/install/space_converter_polaris
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

export LD_LIBRARY_PATH=$lib_dir/lib:$lib_dir/blosc/lib:$lib_dir/openvdb/lib64:$lib_dir/tbb/lib64:$lib_dir/boost/lib64:$LD_LIBRARY_PATH

out=${ROOT_DIR}/out
out_t=${out}/remote
mkdir -p ${out_t}

export MPICH_GPU_SUPPORT_ENABLED=1
mpirun --np 128 -env OMP_NUM_THREADS=2 --line-buffer gpu_bind_kar2.sh ${install}/bin/space_converter --pos-names x y z --vel-names vx vy vz --mass-name mass --rho-name rho --hsml-name hh --data-type GENERICIO --genericio-file HACCdata/m000p.full.mpicosmo.624 --output-path ${out_t} --port 5000 #--calc-radius-neigh 10 #--info
