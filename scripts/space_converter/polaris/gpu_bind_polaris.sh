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

# Per-rank GPU binding for Polaris (4x NVIDIA A100 per node).
# Based on the gpu_bind_kar.sh pattern; Polaris uses PBS/PALS (mpiexec),
# so prefer PALS_LOCAL_RANKID and fall back to Slurm variables.

LOCAL_RANK=$PALS_LOCAL_RANKID
if [ -z "$LOCAL_RANK" ]; then
  LOCAL_RANK=$SLURM_LOCALID
fi
if [ -z "$LOCAL_RANK" ]; then
  LOCAL_RANK=0
fi

export CUDA_VISIBLE_DEVICES=$(( LOCAL_RANK % 4 ))
echo "Process ${PALS_RANKID:-$SLURM_PROCID}: Local rank $LOCAL_RANK, Using GPU: $CUDA_VISIBLE_DEVICES"
exec "$@"
