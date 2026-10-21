#!/usr/bin/env bash
#SBATCH --job-name=ana60
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH --output=/scratch/ipopa/thesis/REST2-conformers/slurm-%x-%A_%a.out
# Shu-Yu clustering (crystal reference + overview figure) and the 5 ns convergence scan on
# 0-60 ns state 0 (REST2-run1 + REST2-run1-ext1), same settings as 4ibi/analysis-0-60ns.
#   sbatch run_analysis_60ns.sbatch                 all four
#   sbatch --array=2 run_analysis_60ns.sbatch       only 5mwh
set -euo pipefail

NAMES=(4ibk 4j7d 5mwh 5ov8)
RESNAMES=(L02 L03 L04 L05)
# automorphism cap: L03-L05 exceed MD-tools' default 200000 (see shuyu/cluster_prune.py)
MAXMATCH=(200000 8000000 8000000 8000000)
i=${SLURM_ARRAY_TASK_ID}
NAME=${NAMES[$i]}; RES=${RESNAMES[$i]}

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}

SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
ROOT=/scratch/ipopa/thesis/REST2-conformers/$NAME
OUT=$ROOT/analysis-0-60ns
mkdir -p "$OUT"
cd "$OUT"

for f in torsions torsions_rigid torsions_polarH torsions_ring; do
    cp -p "$ROOT/REST2-run1/$f.txt" .
done
if [ ! -s solute_state0_0-60ns.nc ]; then
    echo "== joining run1 + ext1  $(date)"
    python - "$ROOT" <<'EOF'
import sys, mdtraj as md
root = sys.argv[1]
top = f"{root}/build/built.solute.pdb"
a = md.load(f"{root}/REST2-run1/solute_state0_prod1.nc", top=top)
b = md.load(f"{root}/REST2-run1-ext1/solute_state0_prod1.nc", top=top)
assert b.time[0] > a.time[-1], (a.time[-1], b.time[0])
j = a.join(b)
j.save_netcdf("solute_state0_0-60ns.nc")
print(f"{j.n_frames} frames, {j.time[0]:g}-{j.time[-1]:g} ps")
EOF
fi

REFS=()
for f in "$CRYSTAL/${NAME}"_*.sdf; do REFS+=(--ref-sdf "$f"); done

COMMON=(--traj solute_state0_0-60ns.nc --top ../build/built.solute.pdb --sdf "../build/$RES.sdf"
        --torsions-file torsions.txt --torsions-file torsions_rigid.txt
        --torsions-file torsions_polarH.txt
        --stride 5 --no-single-cluster --max-matches "${MAXMATCH[$i]}"
        --name "$NAME" --rest2-config ../REST2.config)

echo "== $NAME cluster_prune.py  $(date)"
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" \
    --out "$NAME-clustering-shuyu-nosingle-60ns-ref" \
    > "$NAME-clustering-shuyu-nosingle-60ns-ref.log" 2>&1

echo "== $NAME cluster_convergence.py  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 5 --plateau 3 \
    --out "$NAME-convergence-60ns" \
    > "$NAME-convergence-60ns.log" 2>&1

echo "== $NAME done  $(date)"
grep -h "converged" "$NAME-convergence-60ns.log" | tail -1
