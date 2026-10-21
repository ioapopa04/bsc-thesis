#!/bin/bash
#SBATCH --job-name=rcc-build
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/scratch/ipopa/thesis/REST2-conformers-complex/logs/build-%A_%a.out
# CPU only, one compound per array task: build-top (complex, water-run ligand package reused),
# REST2 scaler (ligand only), build-md for the 40 ns production (REST2-run1) and for a short
# 300-exchange test (test/REST2-test, sharing build/).
#   sbatch scripts/build.sbatch                 all five
#   sbatch --array=2 scripts/build.sbatch       only 4j7d
set -uo pipefail
ROOT=/scratch/ipopa/thesis/REST2-conformers-complex
NAMES=(4ibi 4ibk 4j7d 5mwh 5ov8)
LIGS=(1D8 1DE I31 UWX AXN)
N_REPLICAS=6            # as the water runs
TAU_MAX=0.4             # as the water runs
NS=40                   # 40 ns per state (water: 60 ns)
N_EXCHANGES=$(( NS * 1000000 / (1000 * 2) ))    # 1000 steps x 2 fs per exchange
N_TEST=300

name=${NAMES[$SLURM_ARRAY_TASK_ID]}; lig=${LIGS[$SLURM_ARRAY_TASK_ID]}
b=$ROOT/$name/build
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env
export MD_DATA=/scratch/ipopa/md-data
echo "== $name ($lig) on $(hostname)  job ${SLURM_JOB_ID}  $(date)"

echo "== build-top"
( cd "$b" && md-openmm build-top -i ../prep/complex.pdb --config build-top.config \
    -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) \
  || { echo "build-top FAILED"; grep -v '^$' "$b/built.out" | tail -25; exit 1; }
grep -E 'missing atoms|ligand |reused|net charge|atoms  |solvation|WARNING|status' "$b/built.out" | head -20

echo "== ligand residue in the built topology"
read -r IDX RES < <(python -I - "$b/built.pdb" "$lig" "$ROOT/packages/$name-water/parameter.config" <<'EOF'
import sys, re
from openmm.app import PDBFile
pkg_res = re.search(r"residue_name:\s*(\S+)", open(sys.argv[3]).read()).group(1)
hits = [(i, r.name) for i, r in enumerate(PDBFile(sys.argv[1]).topology.residues(), 1)
        if r.name in (sys.argv[2], pkg_res)]
assert len(hits) == 1, hits
print(*hits[0])
EOF
) || { echo "could not locate the ligand residue"; exit 1; }
echo "   residue :$IDX is $RES"

echo "== scaler"
cp -p "$ROOT/packages/$name-water/molecule.sdf" "$b/$RES.sdf"   # bond orders for the unscaled-torsion classifier
cat > "$b/scaler.config" <<EOF
# Selective REST2: ligand only (residue $IDX = $RES), $N_REPLICAS states, tau 0..$TAU_MAX linear (as the water run).
method: REST2
schedule: {kind: linear, n_states: $N_REPLICAS, tau_min: 0.0, tau_max: $TAU_MAX}
ligand_scaling_dict:
  $RES:
    mask: ":$IDX"
    torsion_exclusions: auto
EOF
( cd "$b" && md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config \
    > scaler.out 2>&1 ) || { echo "scaler FAILED"; tail -20 "$b/scaler.out"; exit 1; }
grep -E 'selection|unscaled torsions|state[0-9]+\.xml' "$b/scaler.out"

write_rest2_config() {   # write_rest2_config <file> <n_exchanges> <comment>
  cat > "$1" <<EOF
# $name complex ($lig): selective REST2, ligand only. $3
# Identical to the water run's REST2.config (REST2-conformers/$name) except the length and the
# ligand-only hot region.
protocol: REST2
solvent: explicit

rest2:
  number_of_replicas: $N_REPLICAS
  tau_max: $TAU_MAX
  exchange_interval_steps: 1000
  number_of_exchanges: $2
  equilibration_steps: 50000
  state_trajectory: true
  rem_log: true
  neighbour_acceptance_report: true
  ligand_scaling_dict:
    $RES: {mask: ":$IDX"}

stages:
  minimization_iterations: 5000
  restrained_nvt_steps: 50000
  restrained_npt_steps: 50000
  unrestrained_npt_steps: 100000

reporting:
  crd_printout_solute: 1000
  crd_printout_whole: 50000
  info_printout: 5000
  checkpoint_printout: 50000
EOF
}

echo "== build-md: production ($N_EXCHANGES exchanges = $NS ns) and test ($N_TEST exchanges)"
write_rest2_config "$ROOT/$name/REST2.config" "$N_EXCHANGES" "$N_EXCHANGES exchanges x 1000 steps x 2 fs = $NS ns per state."
( cd "$ROOT/$name" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
  || { echo "build-md FAILED"; tail -20 "$ROOT/$name/build-md.out"; exit 1; }
mkdir -p "$ROOT/$name/test"; ln -sfn ../build "$ROOT/$name/test/build"
write_rest2_config "$ROOT/$name/test/REST2.config" "$N_TEST" "TEST: $N_TEST exchanges = 0.6 ns per state."
( cd "$ROOT/$name/test" && md-openmm build-md -odir ./REST2-test --config REST2.config > build-md.out 2>&1 ) \
  || { echo "build-md (test) FAILED"; tail -20 "$ROOT/$name/test/build-md.out"; exit 1; }
grep -E 'claimed region|status' "$ROOT/$name/build-md.out" "$ROOT/$name/test/build-md.out"
grep -h 'number_of_exchanges' "$ROOT/$name/input/REST2.in" "$ROOT/$name/test/input/REST2.in"
echo "== $name done  $(date)"
