#!/bin/bash

# Per-rank GPU binding for LUMI-G (8 GCDs of MI250X per node).
# Same purpose as karolina/gpu.sh, but AMD uses ROCR_VISIBLE_DEVICES.

LOCAL_RANK=$SLURM_LOCALID
if [ -z "$LOCAL_RANK" ]; then
  TASKS_PER_NODE=${SLURM_TASKS_PER_NODE%%(*}
  LOCAL_RANK=$((SLURM_PROCID - SLURM_NODEID * TASKS_PER_NODE))
fi
export ROCR_VISIBLE_DEVICES=$(( LOCAL_RANK % 8 ))
echo "Process $SLURM_PROCID: Local rank $LOCAL_RANK on node $SLURM_NODEID, Using GPU: $ROCR_VISIBLE_DEVICES"
exec "$@"
