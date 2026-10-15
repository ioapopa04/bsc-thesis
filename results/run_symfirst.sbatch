#!/usr/bin/env bash
#SBATCH --job-name=symfirst
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/ipopa/bsc-thesis/results/logs/%x-%A_%a.out
# Shu-Yu's SYMMETRY-FIRST clustering (MD-tools main 69c262b, env analysis-env-main) with this
# thesis' figures: shuyu_clustering/cluster_prune.py + cluster_convergence.py. Same inputs as the
# previous analyses (state 0, torsion files, stride, windows, crystal refs), but Shu-Yu's defaults
# for the fit itself (allow_single_cluster=True, no --no-single-cluster). Outputs go to NEW
# folders *-clustering-symfirst / *-convergence-symfirst next to the old ones.
#   sbatch --array=0-7 run_symfirst.sbatch
#   0 paracetamol v4      1 paracetamol (original)
#   2 folate water v4     3 folate 4LRH v4     4 folate 3QL3 v4
#   5 folate water (orig) 6 folate 4LRH (orig) 7 folate 3QL3 (orig)
# NOSINGLE=1 (sbatch --export=ALL,NOSINGLE=1 ...): allow_single_cluster=False, outputs *-symfirst-nosingle
set -euo pipefail
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env-main
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
SHUYU=/scratch/ipopa/bsc-thesis/shuyu_clustering
CRYSTAL=/scratch/ipopa/bsc-thesis/crystal
FOLV4=/scratch/ipopa/FOL-REST2-v4
FOLD=/scratch/ipopa/FOL-complex-REST2
PAIRS=(--pair-2d 5,3 --pair-2d 4,3 --pair-2d 2,3)

fol_args() {   # fol_args <water-root> <torsion-dir>
  ARGS=(--sdf "$1/build/FOL.sdf" --top "$1/build/built.solute.pdb"
        --torsions-file "$2/torsions.txt" --torsions-file "$2/torsions_rigid.txt"
        --torsions-file "$2/torsions_polarH.txt" --stride 5)
}
case ${SLURM_ARRAY_TASK_ID} in
0) NAME=paracetamol-v4; ROOT=/scratch/ipopa/PARA-REST2-v4; TRAJ=REST2-run1/solute_state0_prod1.nc
   ARGS=(--sdf build/TYL.sdf --top build/built.solute.pdb --stride 1); REFS=(); EXTRA=(); STEP=0.5
   CFG=REST2.config ;;
1) NAME=paracetamol; ROOT=/scratch/ipopa/REST2-PARA; TRAJ=REST2-run1/solute_state0_prod1.nc
   ARGS=(--sdf build/TYL.sdf --top build/built.solute.pdb --stride 1); REFS=(); EXTRA=(); STEP=0.5
   CFG=REST2.config ;;
2) NAME=folic-acid-v4-water; ROOT=$FOLV4/water; TRAJ=REST2-run1/solute_state0_prod1.nc
   fol_args "$FOLV4/water" "$FOLV4/water/REST2-run1"; CFG=REST2.config; STEP=5; EXTRA=("${PAIRS[@]}")
   REFS=(--ref-sdf "$CRYSTAL/fol_4LRH_A.sdf" --ref-sdf "$CRYSTAL/fol_3QL3_A.sdf") ;;
3|4) S=$([ "$SLURM_ARRAY_TASK_ID" = 3 ] && echo 4LRH || echo 3QL3)
   NAME=folic-acid-v4-$S-complex; ROOT=$FOLV4/$S/analysis-shuyu; TRAJ=FOL_state0_prod1.nc
   fol_args "$FOLV4/water" "$FOLV4/water/REST2-run1"; CFG=../REST2.config; STEP=5; EXTRA=("${PAIRS[@]}")
   REFS=(--ref-sdf "$CRYSTAL/fol_${S}_A.sdf") ;;
5) NAME=folic-acid; ROOT=/scratch/ipopa/FOLn-REST2-50ns; TRAJ=REST2-run1/solute_state0_prod1.nc
   fol_args "$ROOT" "$ROOT/REST2-run1"; CFG=REST2.config; STEP=5; EXTRA=("${PAIRS[@]}")
   REFS=(--ref-sdf "$CRYSTAL/fol_4LRH_A.sdf" --ref-sdf "$CRYSTAL/fol_3QL3_A.sdf") ;;
6|7) S=$([ "$SLURM_ARRAY_TASK_ID" = 6 ] && echo 4LRH || echo 3QL3)
   NAME=folic-acid-$S-complex; ROOT=$FOLD/$S/analysis-shuyu; TRAJ=FOL_state0_prod1.nc
   ARGS=(--sdf FOL.sdf --top FOL.solute.pdb --torsions-file torsions.txt --torsions-file torsions_rigid.txt
         --torsions-file torsions_polarH.txt --stride 5)
   CFG=../REST2.config; STEP=5; EXTRA=("${PAIRS[@]}"); REFS=(--ref-sdf "$CRYSTAL/fol_${S}_A.sdf") ;;
esac
cd "$ROOT"
echo "== $NAME in $ROOT on $(hostname)  $(date)"
COMMON=(--traj "$TRAJ" "${ARGS[@]}" --name "$NAME" --rest2-config "$CFG")
TAG=symfirst
if [ "${NOSINGLE:-0}" = 1 ]; then COMMON+=(--no-single-cluster); TAG=symfirst-nosingle; fi
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" "${EXTRA[@]}" \
    --out "$NAME-clustering-$TAG" > "$NAME-clustering-$TAG.log" 2>&1 \
  || { echo "cluster_prune FAILED"; tail -20 "$NAME-clustering-$TAG.log"; exit 1; }
grep -E '^n_clusters|^masses|^density noise|symmetry operations|^group |^allow' "$NAME-clustering-$TAG.log"
echo "== $NAME convergence  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns "$STEP" --plateau 3 \
    --out "$NAME-convergence-$TAG" > "$NAME-convergence-$TAG.log" 2>&1 \
  || { echo "convergence FAILED"; tail -20 "$NAME-convergence-$TAG.log"; exit 1; }
grep -h converged "$NAME-convergence-$TAG.log" | tail -1
echo "== $NAME done  $(date)"
