#!/usr/bin/env bash
#SBATCH --job-name=fol-wvb
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ipopa/thesis/bsc-thesis/results/logs/%x-%j.out
# Folate water vs bound torsion overlay for the v4 runs (FOL-REST2-v4), same script and settings
# as the original-run figure results/folic-acid/folic-acid_water_vs_complexes.png.
set -euo pipefail
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
python /scratch/ipopa/thesis/bsc-thesis/shuyu_clustering/compare_water_complex.py \
    --water /scratch/ipopa/thesis/FOL-REST2-v4/water --complex /scratch/ipopa/thesis/FOL-REST2-v4 \
    --label "v4 (COO-, amide scaled)" --out /scratch/ipopa/thesis/bsc-thesis/results/symfirst-2026-10-15/1_v4-new-rest2/folic-acid-v4_water_vs_bound
