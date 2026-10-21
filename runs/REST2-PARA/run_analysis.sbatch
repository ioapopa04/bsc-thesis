#!/usr/bin/env bash
#SBATCH --job-name=para-ana
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ipopa/thesis/REST2-PARA/logs/%x-%j.out
# Shu-Yu analysis of the full paracetamol run, same settings as bsc-thesis/results/run_fol_para.sbatch
# task 1 (no --torsions-file: the script's default omega/aryl/phenol set; stride 1; 0.5 ns windows;
# no crystal reference), with the current bsc-thesis/shuyu_clustering scripts.
set -euo pipefail
cd /scratch/ipopa/thesis/REST2-PARA
grep -q '^run_status: completed' REST2-run1/remd_records/REST2_prod1.out || { echo "run not completed"; exit 1; }
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
COMMON=(--traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb --sdf build/TYL.sdf
        --stride 1 --no-single-cluster --name paracetamol --rest2-config REST2.config)
echo "== cluster_prune.py  $(date)"
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" --out paracetamol-clustering-shuyu-nosingle-ref \
    > paracetamol-clustering-shuyu-nosingle-ref.log 2>&1
echo "== cluster_convergence.py  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 0.5 --plateau 3 \
    --out paracetamol-convergence > paracetamol-convergence.log 2>&1
R=/scratch/ipopa/thesis/bsc-thesis/results/paracetamol
cp -p paracetamol-clustering-shuyu-nosingle-ref/clusters_overview.png $R/paracetamol_clusters_overview.png
cp -p paracetamol-clustering-shuyu-nosingle-ref/clusters_overview_unmerged.png $R/paracetamol_clusters_overview_unmerged.png
cp -p paracetamol-clustering-shuyu-nosingle-ref/group_torsions.png $R/paracetamol_group_torsions.png
cp -p paracetamol-convergence/convergence.png $R/paracetamol_convergence.png
grep -E '^cluster labels|^group ' paracetamol-clustering-shuyu-nosingle-ref.log
grep -h 'converged' paracetamol-convergence.log | tail -1
echo "== done  $(date)"
