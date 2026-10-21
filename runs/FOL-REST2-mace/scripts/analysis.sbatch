#!/usr/bin/env bash
#SBATCH --job-name=folm-ana
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=/scratch/ipopa/thesis/FOL-REST2-mace/logs/%x-%j.out
# Analysis of the MACE folate runs (Shu-Yu's MACE-OMOL bespoke folate, 3H tautomer, -2, REST2 v4),
# same protocol as FOL-REST2-v4 but with the symmetry-first clustering (env analysis-env-main):
#   per system:  check_run.py (acceptance + binding), fragments.py (complex: per-fragment drift and
#                Arg57 / Asp27 / NADP+ contacts for 3QL3), Shu-Yu cluster_prune.py + cluster_convergence.py
#   with "rmsd": RMSD-to-crystal over time + histogram for all three runs.
#   sbatch scripts/analysis.sbatch water 4LRH            sbatch scripts/analysis.sbatch 3QL3 rmsd
# Torsions: regenerated for the MACE molecule (49 atoms, its own atom order)
# with shuyu_clustering/select_torsions.py exactly as for the old water run; the complexes reuse
# the water files on a folate-only trajectory (same atom order, asserted).
set -uo pipefail
ROOT=/scratch/ipopa/thesis/FOL-REST2-mace
SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
cd "$ROOT"
echo "== $* on $(hostname)  $(date)"

micromamba activate openmm-env           # select_torsions.py needs networkx (not in analysis-env)
W=$ROOT/water
if [ ! -s "$W/REST2-run1/torsions.txt" ]; then
  echo "== torsions for the v4 folate"
  ( cd "$W" && python "$SHUYU/select_torsions.py" --sdf build/FOL.sdf --top build/built.pdb \
      --select "resname FOL" --ligand-indices --polar-h --out REST2-run1/torsions ) || { echo "select_torsions FAILED"; exit 1; }
fi
for f in torsions torsions_rigid torsions_polarH; do echo "-- $f.txt"; grep -v '^#' "$W/REST2-run1/$f.txt"; done

for s in "$@"; do
  [ "$s" = rmsd ] && continue
  run=$ROOT/$s/REST2-run1
  grep -q '^run_status: completed' "$run/remd_records/REST2_prod1.out" || { echo "$s: run not completed, skipped"; continue; }
  micromamba activate openmm-env
  echo; echo "=================== $s: acceptance + binding"
  if [ "$s" = water ]; then
    python -I scripts/check_run.py "$run" --states 0 --ligand FOL 2>&1 | sed -n '1,/^binding/p' | grep -v '^binding'
  else
    python -I scripts/check_run.py "$run" --states 0,5,9 --ligand FOL 2>&1
    mkdir -p "$ROOT/$s/analysis"
    extra=""; [ "$s" = 3QL3 ] && extra=contacts
    python -I scripts/fragments.py "$run" "$ROOT/$s/build/built.solute.pdb" $extra 2>&1 | tee "$ROOT/$s/analysis/fragments.txt"
    [ "$s" = 3QL3 ] && python -I scripts/check_run.py "$run" --states 0 --ligand NAP 2>&1 | tail -4
    # folate-only trajectory for the clustering (same 49 atoms, same order as the water topology)
    mkdir -p "$ROOT/$s/analysis-shuyu"
    python -I - "$ROOT/$s" "$W" <<'EOF' || { echo "$s: folate cut FAILED"; continue; }
import sys, mdtraj as md
root, water = sys.argv[1], sys.argv[2]
top = md.load(f"{root}/build/built.solute.pdb")
lig = top.topology.select("resname FOL")
names = [top.topology.atom(i).name for i in lig]
assert names == [a.name for a in md.load(f"{water}/build/built.solute.pdb").topology.atoms], "atom order differs"
t = md.load(f"{root}/REST2-run1/solute_state0_prod1.nc", top=f"{root}/build/built.solute.pdb", atom_indices=lig)
t.save_netcdf(f"{root}/analysis-shuyu/FOL_state0_prod1.nc")
print(f"folate cut: {t.n_frames} frames x {t.n_atoms} atoms")
EOF
  fi

  micromamba activate analysis-env-main
  echo; echo "=================== $s: Shu-Yu clustering"
  if [ "$s" = water ]; then
    cd "$W"; TRAJ=REST2-run1/solute_state0_prod1.nc; NAME=folic-acid-mace-water
    REFS=(--ref-sdf "$CRYSTAL/fol_4LRH_A.sdf" --ref-sdf "$CRYSTAL/fol_3QL3_A.sdf")
  else
    cd "$ROOT/$s/analysis-shuyu"; TRAJ=FOL_state0_prod1.nc; NAME=folic-acid-mace-$s-complex
    REFS=(--ref-sdf "$CRYSTAL/fol_${s}_A.sdf")
  fi
  COMMON=(--traj "$TRAJ" --top "$W/build/built.solute.pdb" --sdf "$W/build/FOL.sdf"
          --torsions-file "$W/REST2-run1/torsions.txt" --torsions-file "$W/REST2-run1/torsions_rigid.txt"
          --torsions-file "$W/REST2-run1/torsions_polarH.txt" --stride 5
          --name "$NAME" --rest2-config "$ROOT/$s/REST2.config")
  for TAG in symfirst-nosingle symfirst; do
    X=(); [ "$TAG" = symfirst-nosingle ] && X=(--no-single-cluster)
    echo "-- $TAG"
    python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${X[@]}" "${REFS[@]}" --pair-2d 5,3 --pair-2d 4,3 --pair-2d 2,3 \
        --out "$NAME-clustering-$TAG" > "$NAME-clustering-$TAG.log" 2>&1 \
      && grep -E '^n_clusters|^masses|^density noise|symmetry operations' "$NAME-clustering-$TAG.log" \
      || { echo "cluster_prune FAILED"; tail -15 "$NAME-clustering-$TAG.log"; }
    python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" "${X[@]}" --step-ns 5 --plateau 3 \
        --out "$NAME-convergence-$TAG" > "$NAME-convergence-$TAG.log" 2>&1 \
      && grep -h converged "$NAME-convergence-$TAG.log" | tail -1 || { echo "convergence FAILED"; tail -15 "$NAME-convergence-$TAG.log"; }
  done
  cd "$ROOT"
done

if [[ " $* " == *" rmsd "* ]]; then
  micromamba activate openmm-env
  echo; echo "=================== RMSD to crystal (water vs 4LRH/3QL3, complexes vs own crystal)"
  mkdir -p analysis
  python -I scripts/rmsd_vs_crystal.py analysis/folic-acid-mace_rmsd_vs_crystal 2>&1 | grep -iv warning
  python -I scripts/rmsd_histogram.py analysis/folic-acid-mace_rmsd_vs_crystal 2>&1 | grep -iv warning
  echo; echo "=================== water vs bound torsions"
  micromamba activate analysis-env-main
  python "$SHUYU/compare_water_complex.py" --water "$ROOT/water" --complex "$ROOT" \
      --label "MACE-OMOL bespoke (3H tautomer, COO-)" --out analysis/folic-acid-mace_water_vs_bound 2>&1 | tail -2
fi
echo "== done  $(date)"
