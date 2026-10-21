#!/usr/bin/env bash
# ============================================================================
# Extend every completed REST2-run1 by EXT_NS ns per state (default 30 -> 60 ns)
# ----------------------------------------------------------------------------
# Uses md-tools' OUT-OF-PLACE extension (--extend-from): the parent REST2-run1/
# is opened read-only and left byte-for-byte unchanged; the new 30 ns are
# written to a fresh sibling directory  <name>/REST2-run1-ext1/  that starts
# from the parent's final checkpoint (positions, velocities, box, state->walker
# mapping, RNG, step count). See MD-tools docs/tutorial/paracetamol/REST2.md s.8.
#
# Needs md-tools >= 0.6.4 (0.6.2 in openmm-env cannot find the manifest of a
# build-md tree). 0.6.4 is installed --no-deps in $OVERLAY and put first on
# PYTHONPATH; openmm-env itself is not modified.
#
# An extension segment is ATOMIC: if it is interrupted it cannot be resumed.
# Delete (or rename) the partial REST2-run1-ext1/ and rerun this script.
#
# USAGE
#   ./extend_all.sh              extend every compound, one after another
#   ./extend_all.sh 5ov8 4ibi    only these
# ENV: GPUS=0,1,2,3,4,5  EXT_NS=30  SEG=ext1  PARENT=REST2-run1
# ============================================================================
set -uo pipefail

BASE=/scratch/ipopa/REST2-conformers
OVERLAY=/scratch/ipopa/software/md-tools-0.6.4-overlay
GPUS="${GPUS:-0,1,2,3,4,5}"
EXT_NS="${EXT_NS:-30}"
SEG="${SEG:-ext1}"
PARENT="${PARENT:-REST2-run1}"
EXCHANGE_INTERVAL=1000; TIMESTEP_FS=2; N_REPLICAS=6
N_EXT=$(( EXT_NS * 1000000 / (EXCHANGE_INTERVAL * TIMESTEP_FS) ))   # 30 ns -> 15000

if [ $# -gt 0 ]; then NAMES=("$@"); else NAMES=(4ibi 4ibk 4j7d 5mwh 5ov8); fi

say() { printf '[%s] %s\n' "$(date '+%F %H:%M:%S')" "$*"; }

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env || { echo "no openmm-env"; exit 1; }
export MD_DATA="${MD_DATA:-/scratch/ipopa/md-data}"
export PYTHONPATH="$OVERLAY${PYTHONPATH:+:$PYTHONPATH}"
say "md_tools $(python -c 'import md_tools.cli.md_openmm as m, importlib.metadata as im; print(im.version("md-tools"))') from $(python -c 'import md_tools,os;print(os.path.dirname(md_tools.__file__))')"
say "extending by $EXT_NS ns = $N_EXT exchanges per state, segment $SEG, GPUs $GPUS"

for name in "${NAMES[@]}"; do
  parent="$BASE/$name/$PARENT"
  seg="$BASE/$name/$PARENT-$SEG"
  log="$BASE/$name/run-$SEG.log"

  if ! grep -q '^run_status: completed' "$parent/remd_records/REST2_prod1.out" 2>/dev/null; then
    say "SKIP $name: parent $parent is not completed"; continue
  fi
  if [ -f "$seg/restart.json" ]; then
    say "SKIP $name: $seg/restart.json exists (already extended)"; continue
  fi
  if [ -e "$seg/REST2.nc" ] || [ -e "$seg/REST2_checkpoint.nc" ]; then
    say "SKIP $name: $seg holds a PARTIAL segment - remove/rename it and rerun"; continue
  fi

  mkdir -p "$seg/eq"
  for f in REST2.py resolved.config remd_groupfile.1 eq/eq_3.xml; do
    cp -p "$parent/$f" "$seg/$f"
  done

  say "RUN  $name: $PARENT -> $PARENT-$SEG  (log $log)"
  (
    cd "$seg" || exit 1
    export CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES="$GPUS"
    mpirun -n "$N_REPLICAS" python REST2.py -ng "$N_REPLICAS" -p ../build/built.pdb \
      --groupfile remd_groupfile.1 \
      --extend "$N_EXT" --extend-from "../$PARENT" \
      --extend-manifest remd_records/restart_prod1.json
  ) >>"$log" 2>&1
  say "RUN  $name exited rc=$?"
done

say "done - status:"
for name in "${NAMES[@]}"; do
  out="$BASE/$name/$PARENT-$SEG/REST2.out"
  printf '  %-6s %-10s %s\n' "$name" \
    "$(grep -oE 'run_status: .*' "$out" 2>/dev/null | tail -1 | cut -d' ' -f2- || true)" \
    "$(grep -oE 'overall +[0-9]+/[0-9]+ +[0-9.]+' "$out" 2>/dev/null | tail -1)"
done
