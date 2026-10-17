#!/usr/bin/env bash
#SBATCH --job-name=rep-compare
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/ipopa/thesis/bsc-thesis/replicates/logs/%x-%A_%a.out
# Replicate comparison (replicates/compare_replicates.py), env analysis-env-main: all replicates of
# one molecule pooled and clustered once (symmetry-first t_hdbscan, same torsions/settings as
# run_symfirst.sbatch), then per-replicate populations, timeline, convergence, torsions.
#   sbatch --array=0-1 replicates/compare_replicates.sbatch     0 paracetamol v4   1 folate water v4
# Outputs: replicates/<name>-compare/
set -euo pipefail
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
HERE=/scratch/ipopa/thesis/bsc-thesis/replicates
REPS=/scratch/ipopa/thesis/REST2-replicates
rep() { echo "$1:$2/REST2-run1/solute_state0_prod1.nc:$3/built.solute.pdb"; }
case ${SLURM_ARRAY_TASK_ID} in
0) NAME=paracetamol-v4; R=/scratch/ipopa/thesis/PARA-REST2-v4
   ARGS=(--top "$R/build/built.solute.pdb" --sdf "$R/build/TYL.sdf" --stride 1
         --rep "$(rep rep0 "$R" "$R/build")" --rep "$(rep rep1 "$REPS/PARA/rep1" "$REPS/PARA/rep1/build")") ;;
1) NAME=folic-acid-v4-water; R=/scratch/ipopa/thesis/FOL-REST2-v4/water; T=$R/REST2-run1
   ARGS=(--top "$R/build/built.solute.pdb" --sdf "$R/build/FOL.sdf" --stride 5 --no-single-cluster
         --torsions-file "$T/torsions.txt" --torsions-file "$T/torsions_rigid.txt"
         --torsions-file "$T/torsions_polarH.txt" --rep "$(rep rep0 "$R" "$R/build")")
   for k in 1 2 3; do ARGS+=(--rep "$(rep rep$k "$REPS/FOL/rep$k" "$REPS/FOL/rep$k/build")"); done ;;
esac
echo "== $NAME on $(hostname)  $(date)"
python "$HERE/compare_replicates.py" "${ARGS[@]}" --name "$NAME" --out "$HERE/$NAME-compare"
echo "== done  $(date)"
