#!/bin/bash

#fix_filename="/mnt/proj3/open-30-28/adamiano/Centers_1e5_DFZ19.txt"
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
$@