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

#export CUDA_VISIBLE_DEVICES=$(( $PMIX_RANK % 8 ))
# echo "Process $SLURM_PROCID: Local rank $LOCAL_RANK on node $SLURM_NODEID, Using GPU: $CUDA_VISIBLE_DEVICES"
# $@

LOCAL_RANK=$SLURM_LOCALID
if [ -z "$LOCAL_RANK" ]; then
  TASKS_PER_NODE=${SLURM_TASKS_PER_NODE%%(*}
  LOCAL_RANK=$((SLURM_PROCID - SLURM_NODEID * TASKS_PER_NODE))
fi
export CUDA_VISIBLE_DEVICES=$(( $LOCAL_RANK % 8 ))
echo "Process $SLURM_PROCID: Local rank $LOCAL_RANK on node $SLURM_NODEID, Using GPU: $CUDA_VISIBLE_DEVICES"
exec "$@"