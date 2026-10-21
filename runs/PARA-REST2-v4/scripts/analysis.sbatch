#!/usr/bin/env bash
#SBATCH --job-name=para4-ana
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ipopa/thesis/PARA-REST2-v4/logs/%x-%j.out
# Paracetamol v4: acceptance + Shu-Yu clustering/convergence, same settings as REST2-PARA
# (default omega/aryl/phenol torsions, stride 1, 0.5 ns windows, no crystal reference).
set -euo pipefail
cd /scratch/ipopa/thesis/PARA-REST2-v4
grep -q '^run_status: completed' REST2-run1/remd_records/REST2_prod1.out || { echo "run not completed"; exit 1; }
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
python3 - <<'PY'
import json
cr = json.load(open("REST2-run1/remd_records/restart_prod1.json"))["completion_report"]
for p in cr["by_neighbouring_pair"]:
    print(f"acceptance {p['state_pair']} tau {p['tau_pair']}: {p['acceptance']:.3f}")
PY
SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
COMMON=(--traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb --sdf build/TYL.sdf
        --stride 1 --no-single-cluster --name paracetamol-v4 --rest2-config REST2.config)
cp -p build/parameter/TYL.sdf build/TYL.sdf
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" --out paracetamol-v4-clustering-shuyu-nosingle-ref \
    > paracetamol-v4-clustering-shuyu-nosingle-ref.log 2>&1
grep -E '^cluster labels|^groups:|^group ' paracetamol-v4-clustering-shuyu-nosingle-ref.log
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 0.5 --plateau 3 \
    --out paracetamol-v4-convergence > paracetamol-v4-convergence.log 2>&1
grep -h converged paracetamol-v4-convergence.log | tail -1
echo "== done $(date)"
