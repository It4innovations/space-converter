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


cd ${SPACE_ROOT:-${PWD}}   # workspace root with src/, build/, install/ (default: current directory)

ROOT_DIR=${PWD}

#ml Forge 

# ml HDF5/1.14.0-gompi-2022b

# ml c-blosc/1.21.0-GCC-10.2.0
# ml Boost/1.79.0-GCC-11.3.0
# ml tbb/2021.10.0-GCCcore-12.2.0
# ml VirtualGL/3.1-GCC-12.3.0

# ml Qt5/5.15.7-GCCcore-12.2.0
ml CMake/3.24.3-GCCcore-12.2.0
ml GCC/12.2.0
# ml Mesa/22.2.4-GCCcore-12.2.0
# ml CUDA/12.4.0
# ml OpenMPI/4.1.4-GCC-12.2.0

 ml OpenMPI/4.1.6-NVHPC-24.1-CUDA-12.4.0

echo "===================RUN"

lib_dir=${ROOT_DIR}/install

install=${ROOT_DIR}/install/space_converter_kar
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out

export LD_LIBRARY_PATH=$lib_dir/lib-linux_x64/openvdb/lib:$lib_dir/lib-linux_x64/tbb/lib:$LD_LIBRARY_PATH

out_t=${out}/remote
mkdir -p ${out} ${out_t}

# Per-rank GPU binding (CUDA_VISIBLE_DEVICES) for GPU runs
gpu_bind=${src}/space-converter/scripts/space_converter/karolina/gpu.sh

# The examples expect the datasets under ${data}; adjust the file names to your data.

# Dataset info (particle types and data blocks), then exit:
#srun -u -n 4 ${install}/bin/space_converter --data-type GADGET --gadget-file ${data}/very_small_example/snap_081 --max-mem-size 100000 --buffer-size 150.0 --part-alloc-factor 1.5 --grid-dim 100 --output-path ${out} --info

# Batch export of particle type 0, block 1 (dense grid, WendlandC6 kernel, k-NN radius on the GPU):
#srun -u -n 8 --overlap --gres=gpu:8 ${gpu_bind} ${install}/bin/space_converter --data-type CHANGA_TIPSY --tipsy-file ${data}/changa_tipsy/snapshot.00410 --grid-dim 1000 --output-path ${out} --calc-radius-neigh 100 --export-data 0 1 --dense-type 6

# Animation: snapshots 0..1400 ({} is replaced by the frame number), merged into one output:
#srun -u -n 1400 ${install}/bin/space_converter --data-type GADGET --gadget-file ${data}/snapdir_{}/snap_{} --max-mem-size 1000 --buffer-size 150.0 --part-alloc-factor 1.2 --grid-dim 500 --output-path ${out_t} --port 5000 --anim 0 1400 1 --bh-count 1 --anim-merge --raw-particles 1

# HACC GenericIO particles:
#srun -u -n 32 ${install}/bin/space_converter --data-type HACC_GENERICIO --genericio-file ${data}/hacc/m000p-499.mpicosmo --pos-names x y z --vel-names vx vy vz --grid-dim 1000 --output-path ${out} --info

# Server for the BSPACE Blender add-on (listens on --port until killed):
srun --overlap -u -n 128 ${install}/bin/space_converter --data-type CHANGA_TIPSY --tipsy-file ${data}/changa_tipsy/snapshot.00410 --grid-dim 1000 --output-path ${out} --port 5000
