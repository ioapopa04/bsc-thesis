#!/usr/bin/env bash
#SBATCH --job-name=small-symfirst
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/ipopa/thesis/bsc-thesis/results/logs/%x-%A_%a.out
# Shu-Yu's symmetry-first clustering (cluster_prune.py + cluster_convergence.py, env
# analysis-env-main) for the five small molecules, water and complex:
#   array 0-4  water    4ibi 4ibk 4j7d 5mwh 5ov8
#   array 5-9  complex  4ibi 4ibk 4j7d 5mwh 5ov8
# SET=orig (default): the original runs (REST2-conformers 0-60 ns, REST2-conformers-complex 40 ns)
# SET=v4            : REST2-conformers-v4 (Open Babel pH 7.4 protomers; 4ibi 4ibk 4j7d only)
# Settings as the earlier analyses: state 0, stride 5, 5 ns windows, crystal refs, the water run's
# torsion files (complexes: ligand cut out of the solute trajectory, same atoms and order).
# NOSINGLE=1 -> allow_single_cluster=False (outputs *-symfirst-nosingle), else Shu-Yu's default.
#   sbatch --array=0-9 run_small_symfirst.sbatch
#   sbatch --export=ALL,SET=v4,NOSINGLE=1 --array=0-2,5-7 run_small_symfirst.sbatch
set -euo pipefail
NAMES=(4ibi 4ibk 4j7d 5mwh 5ov8); RESS=(L01 L02 L03 L04 L05); LIGS=(1D8 1DE I31 UWX AXN)
MAXMATCH=(200000 200000 8000000 8000000 8000000)   # automorphism cap, as in the earlier analyses
i=$(( SLURM_ARRAY_TASK_ID % 5 )); KIND=$([ "$SLURM_ARRAY_TASK_ID" -lt 5 ] && echo water || echo complex)
n=${NAMES[$i]}; res=${RESS[$i]}; lig=${LIGS[$i]}; SET=${SET:-orig}
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
T=/scratch/ipopa/thesis

if [ "$SET" = orig ]; then
  W=$T/REST2-conformers/$n; WTRAJ=$W/analysis-0-60ns/solute_state0_0-60ns.nc; TORS=$W/analysis-0-60ns
  C=$T/REST2-conformers-complex/$n; TAGN=$n
else
  W=$T/REST2-conformers-v4/$n-water; WTRAJ=$W/REST2-run1/solute_state0_prod1.nc; TORS=$W/REST2-run1
  C=$T/REST2-conformers-v4/$n-complex; TAGN=$n-v4
  if [ ! -s "$TORS/torsions.txt" ]; then        # torsions for the re-protonated molecule
    micromamba activate openmm-env               # select_torsions.py needs networkx
    ( cd "$W" && python "$SHUYU/select_torsions.py" --sdf "build/$res.sdf" --top build/built.pdb \
        --select "resname $res" --ligand-indices --polar-h --out REST2-run1/torsions ) || exit 1
  fi
fi
micromamba activate analysis-env-main

if [ "$KIND" = water ]; then
  DIR=$W/analysis-symfirst; TRAJ=$WTRAJ; NAME=$TAGN-water
  REFS=(); for f in "$CRYSTAL/${n}"_*.sdf; do REFS+=(--ref-sdf "$f"); done
else
  DIR=$C/analysis-symfirst; TRAJ=$DIR/${res}_state0_prod1.nc; NAME=$TAGN-complex
  REFS=(); for f in "$CRYSTAL/${n}"_*_A.sdf; do REFS+=(--ref-sdf "$f"); done
  [ ${#REFS[@]} -eq 0 ] && for f in "$CRYSTAL/${n}"_*.sdf; do REFS+=(--ref-sdf "$f"); done
fi
mkdir -p "$DIR"; cd "$DIR"
echo "== $NAME ($SET) on $(hostname)  $(date)"
if [ "$KIND" = complex ] && [ ! -s "$TRAJ" ]; then
  python -I - "$C" "$W" "$res" "$TRAJ" "$lig" <<'EOF'
import sys, mdtraj as md
c, w, res, out, lig = sys.argv[1:]
top = md.load(f"{c}/build/built.solute.pdb")
lig = [a.index for a in top.topology.atoms if a.residue.name in (res, lig)]  # crystal name kept in the complex
names = [top.topology.atom(i).name for i in lig]
assert names == [a.name for a in md.load(f"{w}/build/built.solute.pdb").topology.atoms], "ligand atom order differs"
t = md.load(f"{c}/REST2-run1/solute_state0_prod1.nc", top=f"{c}/build/built.solute.pdb", atom_indices=lig)
t.save_netcdf(out)
print(f"ligand cut: {t.n_frames} frames x {t.n_atoms} atoms; atom order == water run")
EOF
fi
COMMON=(--traj "$TRAJ" --top "$W/build/built.solute.pdb" --sdf "$W/build/$res.sdf"
        --torsions-file "$TORS/torsions.txt" --torsions-file "$TORS/torsions_rigid.txt"
        --torsions-file "$TORS/torsions_polarH.txt" --stride 5 --max-matches "${MAXMATCH[$i]}"
        --name "$NAME" --rest2-config "$([ "$KIND" = water ] && echo "$W" || echo "$C")/REST2.config")
TAG=symfirst
if [ "${NOSINGLE:-0}" = 1 ]; then COMMON+=(--no-single-cluster); TAG=symfirst-nosingle; fi
python "$SHUYU/cluster_prune.py" "${COMMON[@]}" "${REFS[@]}" --out "$NAME-clustering-$TAG" \
    > "$NAME-clustering-$TAG.log" 2>&1 || { echo "cluster_prune FAILED"; tail -20 "$NAME-clustering-$TAG.log"; exit 1; }
grep -E '^n_clusters|^masses|^density noise|symmetry operations' "$NAME-clustering-$TAG.log"
python "$SHUYU/cluster_convergence.py" "${COMMON[@]}" --step-ns 5 --plateau 3 --out "$NAME-convergence-$TAG" \
    > "$NAME-convergence-$TAG.log" 2>&1 || { echo "convergence FAILED"; tail -20 "$NAME-convergence-$TAG.log"; exit 1; }
grep -h converged "$NAME-convergence-$TAG.log" | tail -1
echo "== $NAME done  $(date)"
