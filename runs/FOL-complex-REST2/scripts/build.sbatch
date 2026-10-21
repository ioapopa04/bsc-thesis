#!/bin/bash
#SBATCH --job-name=folc-build
#SBATCH --partition=defq
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=/scratch/ipopa/thesis/FOL-complex-REST2/logs/build-%j.out
# CPU-only: build-top (complex), REST2 scaler and build-md for every system named on the command line.
#   sbatch scripts/build.sbatch 4LRH 3QL3
set -uo pipefail
ROOT=/scratch/ipopa/thesis/FOL-complex-REST2
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env
export MD_DATA=/scratch/ipopa/md-data
echo "host $(hostname)  job ${SLURM_JOB_ID:-none}  $(date)"
declare -A IDX=([4LRH]=208 [3QL3]=160)
rc_all=0
for s in "$@"; do
  b=$ROOT/$s/build
  echo "=== $s: build-top"
  ( cd "$b" && md-openmm build-top -i ../prep/complex.pdb --config build-top.config \
      -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) \
    || { echo "$s build-top FAILED"; tail -20 "$b/built.out"; rc_all=1; continue; }
  grep -E 'reused|created|net charge|atoms  |status' "$b/built.out" | head

  echo "=== $s: check residue ${IDX[$s]} is FOL"
  python -I - "$b/built.pdb" "${IDX[$s]}" <<'EOF' || { rc_all=1; continue; }
import sys
from openmm.app import PDBFile
r = list(PDBFile(sys.argv[1]).topology.residues())[int(sys.argv[2]) - 1]
print(f"  residue {sys.argv[2]} = {r.name} chain {r.chain.id} id {r.id}")
sys.exit(0 if r.name == "FOL" else 1)
EOF

  echo "=== $s: scaler"
  cp "$ROOT/packages/FOL-water/molecule.sdf" "$b/FOL.sdf"   # scaler reads bond orders from here
  ( cd "$b" && md-openmm build-top --rest2-scaler -s built.xml -p built.pdb \
      --config scaler.config > scaler.out 2>&1 ) \
    || { echo "$s scaler FAILED"; tail -20 "$b/scaler.out"; rc_all=1; continue; }
  grep -E 'selection|state[0-9]+\.xml' "$b/scaler.out"

  echo "=== $s: build-md"
  ( cd "$ROOT/$s" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
    || { echo "$s build-md FAILED"; tail -20 "$ROOT/$s/build-md.out"; rc_all=1; continue; }
  grep -E 'claimed region|status' "$ROOT/$s/build-md.out"
done
echo "done $(date) rc=$rc_all"
exit $rc_all
