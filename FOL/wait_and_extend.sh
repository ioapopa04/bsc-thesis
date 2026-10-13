#!/usr/bin/env bash
# Wait until all 8 GPUs are idle, then EXTEND the completed 50 ns REST2 ladder
# by another 50 ns IN PLACE (same REST2.nc / *_prod1.nc files, appended).
#
# This is NOT run.sh: run.sh only runs prod1, which is already complete. The
# extension goes through the generated REST2.py, whose --extend is the supported
# "lengthen a finished run" path (md-openmm md-run does not expose --extend).
#
# Detached on purpose: the wait is against another user's training job of
# unknown length, so this must outlive any single tool call.
set -uo pipefail

ROOT=/scratch/ipopa/FOLn-REST2-50ns
RUN="$ROOT/REST2-run1"
LOG="$ROOT/wait_and_extend.log"
POLL=60            # seconds between GPU checks
STABLE=2           # consecutive free checks required before launching

# This segment adds 25000 exchanges (x 2 ps) = 50 ns, taking the ladder from
# 25000 exchanges / 50 ns to 50000 exchanges / 100 ns. The absolute step total
# becomes 25000000 + 25000*1000 = 50000000, which the guard below checks.
EXTEND=25000
TARGET_STEP=50000000

exec >>"$LOG" 2>&1

free_gpus() {
  nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader,nounits \
    | awk -F, '{u=$1+0; m=$2+0; if (u<10 && m<3000) c++} END{print c+0}'
}

echo "=== $(date -Is) extension watcher started; waiting for all 8 GPUs to be idle ==="
echo "    (busy: another user's training job; will extend REST2 by +50 ns when free)"

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

echo "=== $(date -Is) all 8 GPUs free -> launching REST2 +50 ns extension ==="
cd "$RUN" || { echo "FATAL: no run dir $RUN"; exit 1; }

source /scratch/ipopa/folenv.sh
micromamba activate openmm-env
export CUDA_DEVICE_ORDER=PCI_BUS_ID
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

# --- guard: never extend twice -------------------------------------------------
step=$(python - <<'PY'
from md_tools.remd import storage
try:
    print(int(storage.ReplicaCheckpoint("REST2_checkpoint.nc").read()["step"]))
except Exception as exc:                       # noqa: BLE001
    print("ERROR:%s" % exc)
PY
)
if [ "${step#ERROR:}" != "$step" ]; then
  echo "FATAL: could not read REST2_checkpoint.nc ($step); refusing to launch"
  exit 1
fi
echo "checkpoint step before launch: $step (target after extension: $TARGET_STEP)"
if [ "$step" -ge "$TARGET_STEP" ]; then
  echo "ABORT: checkpoint is already at/after $TARGET_STEP - this ladder looks already extended."
  echo "       Refusing to launch a second extension. Remove this guard deliberately to override."
  exit 1
fi

start=$(date +%s)
mpirun -n 8 python REST2.py -ng 8 \
  -p ../build/built.pdb \
  --groupfile remd_groupfile.1 \
  -odir . \
  -x REST2.nc \
  --checkpoint REST2_checkpoint.nc \
  -r remd_records/restart_ext1.json \
  -o remd_records/REST2_ext1.out \
  -log remd_records/REST2_ext1.log \
  --extend "$EXTEND"
rc=$?
end=$(date +%s)
echo "=== $(date -Is) extension exited rc=$rc after $(( (end-start)/60 )) min ==="
exit "$rc"
