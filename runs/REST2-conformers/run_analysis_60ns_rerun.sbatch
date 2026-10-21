#!/usr/bin/env bash
#SBATCH --job-name=ana60-rerun
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-4
#SBATCH --output=/scratch/ipopa/thesis/REST2-conformers/slurm-%x-%A_%a.out
# Re-run of the Shu-Yu analysis on 0-60 ns state 0 for all five compounds with the CURRENT
# bsc-thesis/shuyu_clustering scripts (the old run_analysis_60ns.sbatch and
# 4ibi/analysis-0-60ns/run_analysis.sbatch point at bsc-thesis/shuyu/, since renamed).
# Arguments are identical to those two scripts; only the script path differs.
set -euo pipefail

NAMES=(4ibi 4ibk 4j7d 5mwh 5ov8)
RESNAMES=(L01 L02 L03 L04 L05)
MAXMATCH=(200000 200000 8000000 8000000 8000000)   # automorphism cap, as before
i=${SLURM_ARRAY_TASK_ID}
NAME=${NAMES[$i]}; RES=${RESNAMES[$i]}

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
cd "/scratch/ipopa/thesis/REST2-conformers/$NAME/analysis-0-60ns"

REFS=()
for f in "$CRYSTAL/${NAME}"_*.sdf; do REFS+=(--ref-sdf "$f"); done
COMMON=(--traj solute_state0_0-60ns.nc --top ../build/built.solute.pdb --sdf "../build/$RES.sdf"
        --torsions-file torsions.txt --torsions-file torsions_rigid.txt
        --torsions-file torsions_polarH.txt
        --stride 5 --no-single-cluster --max-matches "${MAXMATCH[$i]}"
        --name "$NAME" --rest2-config ../REST2.config)

OUT="$NAME-clustering-shuyu-nosingle-60ns-ref"
cp -p "$OUT/groups.csv" groups.before-rerun.csv 2>/dev/null || true
echo "== $NAME cluster_prune.py  $(date)"
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" --out "$OUT" > "$OUT.log" 2>&1
cmp -s groups.before-rerun.csv "$OUT/groups.csv" && echo "groups.csv unchanged" || echo "groups.csv CHANGED"

cp -p "$NAME-convergence-60ns/convergence.csv" convergence.before-rerun.csv 2>/dev/null || true
echo "== $NAME cluster_convergence.py  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 5 --plateau 3 \
    --out "$NAME-convergence-60ns" > "$NAME-convergence-60ns.log" 2>&1
cmp -s convergence.before-rerun.csv "$NAME-convergence-60ns/convergence.csv" && echo "convergence.csv unchanged" || echo "convergence.csv CHANGED"
echo "== $NAME done  $(date)"
grep -h "converged" "$NAME-convergence-60ns.log" | tail -1
