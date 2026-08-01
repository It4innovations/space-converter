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
# Run space_converter (HIP build) on LUMI-G: multi GPU (8 GCDs/node), multi node, MPI.
#

cd /flash/project_465002608/jaromila/space

ml purge

ml LUMI/25.09 partition/G
ml PrgEnv-gnu
ml rocm/6.4.4

ROOT_DIR=${PWD}

echo "===================RUN"

lib_dir=${ROOT_DIR}/install
install=${ROOT_DIR}/install/space_converter_lumi_hip
src=${ROOT_DIR}/src
data=${ROOT_DIR}/data
out=${ROOT_DIR}/out
mkdir -p ${out}

export LD_LIBRARY_PATH=${lib_dir}/lib-linux_x64/openvdb/lib:${lib_dir}/lib-linux_x64/tbb/lib:${ROCM_PATH}/lib:${LD_LIBRARY_PATH}

# GPU-aware MPI (Cray MPICH + GTL) - required by WITH_HIP_AWARE_MPI build
export MPICH_GPU_SUPPORT_ENABLED=1

# Per-rank GPU binding (ROCR_VISIBLE_DEVICES)
gpu_bind=${src}/space-converter/scripts/space_converter/lumi/gpu.sh

# CPU binding for standard-g/dev-g nodes: 7 cores per task, matching the
# closest NUMA domain of each GCD (low-noise cores, core 0 of each L3 reserved)
CPU_BIND="mask_cpu:0xfe000000000000,0xfe00000000000000,0xfe0000,0xfe000000,0xfe,0xfe00,0xfe00000000,0xfe0000000000"

# Small OpenGadget test dataset (already staged on LUMI /scratch).
# Original location on Karolina: /mnt/proj3/open-30-28/OpenGadget/data/verysmall_data/very_small_example
data=/scratch/project_465002608/jaromila/space/data/verysmall_data/very_small_example

# --gpu selects the GPU pipeline at runtime (AMD GPUs in the HIP build).
# Add --cub to write the sort/reduce sparse format (.cub) instead of .vdb.
export example_snap="--data-type GADGET_SIMPLE --gadget-file ${data}/snap_081 --grid-dim 100 --output-path ${out} --gpu --export-data 0 1"

# 2 nodes x 8 ranks = 16 GCDs (1 MPI rank per GCD)
srun --overlap --account=project_465002608 --partition=dev-g --time=00:30:00 \
     --nodes=2 --ntasks-per-node=8 --gpus-per-node=8 \
     --cpu-bind=${CPU_BIND} \
     --export=ALL -u \
     ${gpu_bind} ${install}/bin/space_converter ${example_snap}

# single node, 8 GCDs:
#srun --overlap --account=project_465002608 --partition=dev-g --time=00:30:00 \
#     --nodes=1 --ntasks-per-node=8 --gpus-per-node=8 \
#     --cpu-bind=${CPU_BIND} \
#     --export=ALL -u \
#     ${gpu_bind} ${install}/bin/space_converter ${example_snap}

# inside an existing allocation (salloc / --jobid <id>), e.g.:
#srun --overlap --jobid <JOBID> --ntasks-per-node=8 --gpus-per-node=8 --export=ALL -u \
#     ${gpu_bind} ${install}/bin/space_converter ${example_snap}
