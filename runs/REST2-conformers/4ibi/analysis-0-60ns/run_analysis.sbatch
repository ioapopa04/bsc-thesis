#!/usr/bin/env bash
#SBATCH --job-name=4ibi-ana60
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=slurm-%x-%j.out
# Shu-Yu clustering (with crystal reference + overview figure) and the 5 ns convergence scan
# for 4ibi, 0-60 ns state 0. Submit from analysis-0-60ns/:  sbatch run_analysis.sbatch
set -euo pipefail

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
cd /scratch/ipopa/thesis/REST2-conformers/4ibi/analysis-0-60ns

COMMON=(--traj solute_state0_0-60ns.nc --top ../build/built.solute.pdb --sdf ../build/L01.sdf
        --torsions-file torsions.txt --torsions-file torsions_rigid.txt
        --torsions-file torsions_polarH.txt
        --stride 5 --no-single-cluster --name 4ibi --rest2-config ../REST2.config)

echo "== cluster_prune.py  $(date)"
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" \
    --ref-sdf "$CRYSTAL/4ibi_1D8_A.sdf" --ref-sdf "$CRYSTAL/4ibi_1D8_B.sdf" \
    --out 4ibi-clustering-shuyu-nosingle-60ns-ref \
    > 4ibi-clustering-shuyu-nosingle-60ns-ref.log 2>&1

echo "== cluster_convergence.py  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 5 --plateau 3 \
    --out 4ibi-convergence-60ns \
    > 4ibi-convergence-60ns.log 2>&1

echo "== done  $(date)"
