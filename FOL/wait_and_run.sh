#!/usr/bin/env bash
# Wait until all 8 GPUs are idle, then launch the 50 ns REST2 ladder, one
# replica per GPU (8 replicas, CUDA_VISIBLE_DEVICES=0..7, no MPS).
#
# Detached on purpose: the wait is against another user's training job of
# unknown length, so this must outlive any single tool call.
set -uo pipefail

ROOT=/scratch/ipopa/FOLn-REST2-50ns
RUN="$ROOT/REST2-run1"
LOG="$ROOT/wait_and_run.log"
POLL=60            # seconds between GPU checks
STABLE=2           # consecutive free checks required before launching

exec >>"$LOG" 2>&1

free_gpus() {
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits \
    | awk -F, '{u=$1+0; m=$2+0; if (u<10 && m<3000) c++} END{print c+0}'
}

echo "=== $(date -Is) watcher started; waiting for all 8 GPUs to be idle ==="
echo "    (busy: nicosuter training; will launch REST2 50 ns when free)"

streak=0
while true; do
  n=$(free_gpus)
  if [ "$n" -eq 8 ]; then
    streak=$((streak + 1))
    echo "    $(date -Is) free check $streak/$STABLE (8/8 idle)"
    [ "$streak" -ge "$STABLE" ] && break
  else
    [ "$streak" -ne 0 ] && echo "    $(date -Is) GPUs busy again ($n/8 idle) - resetting"
    streak=0
  fi
  sleep "$POLL"
done

echo "=== $(date -Is) all 8 GPUs free -> launching REST2 50 ns ==="
cd "$RUN" || { echo "FATAL: no run dir $RUN"; exit 1; }

source /scratch/ipopa/folenv.sh
micromamba activate openmm-env
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

start=$(date +%s)
./run.sh
rc=$?
end=$(date +%s)
echo "=== $(date -Is) run.sh exited rc=$rc after $(( (end-start)/60 )) min ==="
exit "$rc"
