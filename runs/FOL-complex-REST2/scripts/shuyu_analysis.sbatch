#!/usr/bin/env bash
#SBATCH --job-name=ana-folc
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --array=0-1
#SBATCH --output=/scratch/ipopa/thesis/FOL-complex-REST2/logs/%x-%A_%a.out
# Shu-Yu analysis (bsc-thesis/shuyu_clustering: cluster_prune.py + cluster_convergence.py) of the
# bound-folate REST2 runs, with EXACTLY the settings of the water run (bsc-thesis/results/
# run_fol_para.sbatch, task 0): state 0, stride 5, --no-single-cluster, 5 ns windows, plateau 3,
# the water run's torsion files; crystal ref = the complex's own structure only (4LRH A or 3QL3 A).
#   0: 4LRH (FRa)   1: 3QL3 (DHFR + NADP+)
#
# The water torsion files hold ligand-LOCAL indices 0..50 (the water solute was folate alone). Here
# the solute is protein + folate (+ NADP+), so folate is first cut out into a ligand-only trajectory
# and topology; its 51 atoms are in the same package order as in the water run (checked below).
set -euo pipefail

SYSTEMS=(4LRH 3QL3)
SYS=${SYSTEMS[$SLURM_ARRAY_TASK_ID]}
ROOT=/scratch/ipopa/thesis/FOL-complex-REST2/$SYS
WATER=/scratch/ipopa/thesis/FOLn-REST2-50ns
SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
OUT=$ROOT/analysis-shuyu
NAME="folic-acid-$SYS-complex"

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
echo "== $SYS on $(hostname)  $(date)"

mkdir -p "$OUT"
cd "$OUT"
for f in torsions torsions_rigid torsions_polarH torsions_ring; do
    cp -p "$WATER/REST2-run1/$f.txt" .
done
cp -p "$WATER/build/FOL.sdf" FOL.sdf
# Topology: the water run's own folate-only solute PDB (same 51 atoms, same order -- asserted below,
# and it carries the CONECT records cluster_prune.py needs; an mdtraj atom_slice + save_pdb did not).
cp -p "$WATER/build/built.solute.pdb" FOL.solute.pdb

if [ ! -s FOL_state0_prod1.nc ]; then
    echo "== cutting folate out of the state-0 solute trajectory  $(date)"
    python - "$ROOT" "$WATER" <<'EOF'
import sys, mdtraj as md
root, water = sys.argv[1], sys.argv[2]
top = md.load(f"{root}/build/built.solute.pdb")
lig = top.topology.select("resname FOL")
wat = md.load(f"{water}/build/built.solute.pdb").topology
names = [top.topology.atom(i).name for i in lig]
assert names == [a.name for a in wat.atoms], "folate atom order differs from the water run"
t = md.load(f"{root}/REST2-run1/solute_state0_prod1.nc", top=f"{root}/build/built.solute.pdb",
            atom_indices=lig)
sub = top.atom_slice(lig)
t.save_netcdf("FOL_state0_prod1.nc")
print(f"{t.n_frames} frames x {t.n_atoms} atoms, {t.time[0]:g}-{t.time[-1]:g} ps; atom order == water run")
EOF
fi

COMMON=(--traj FOL_state0_prod1.nc --top FOL.solute.pdb --sdf FOL.sdf
        --torsions-file torsions.txt --torsions-file torsions_rigid.txt
        --torsions-file torsions_polarH.txt --stride 5
        --no-single-cluster --name "$NAME" --rest2-config ../REST2.config)
# Each complex is compared with ITS OWN crystal only (the water run is the one that shows both).
REFS=(--ref-sdf "$CRYSTAL/fol_${SYS}_A.sdf")

echo "== $NAME cluster_prune.py  $(date)"
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" \
    --out "$NAME-clustering-shuyu-nosingle-ref" > "$NAME-clustering-shuyu-nosingle-ref.log" 2>&1

echo "== $NAME cluster_convergence.py  $(date)"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 5 --plateau 3 \
    --out "$NAME-convergence" > "$NAME-convergence.log" 2>&1

echo "== $NAME done  $(date)"
grep -h "converged" "$NAME-convergence.log" | tail -1
