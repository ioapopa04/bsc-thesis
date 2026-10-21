#!/usr/bin/env bash
# Re-run bsc-thesis/shuyu_clustering/cluster_prune.py (new horizontal clusters_overview layout) for
# folic acid in water, paracetamol, and the 4LRH / 3QL3 complexes, with the SAME arguments as their
# previous runs. Run through Slurm:  srun -p defq -x sbrinz,sbrinz1 -c 8 --mem=32G bash rerun_overviews.sh
set -uo pipefail
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate analysis-env
export OMP_NUM_THREADS=8
SHUYU=/scratch/ipopa/thesis/bsc-thesis/shuyu_clustering
CRYSTAL=/scratch/ipopa/thesis/bsc-thesis/crystal
run() {  # run <dir> <out> <args...>
  local dir=$1 out=$2; shift 2
  cd "$dir" || return 1
  cp -p "$out/groups.csv" "/tmp/groups-before-$(basename "$out").csv" 2>/dev/null
  python "$SHUYU/cluster_prune.py" "$@" --out "$out" > "$out.log" 2>&1
  local rc=$?
  if cmp -s "/tmp/groups-before-$(basename "$out").csv" "$out/groups.csv"; then same="groups.csv unchanged"; else same="groups.csv CHANGED"; fi
  echo "$(basename "$dir")/$out rc=$rc  $same"
}
# folic acid in water (bsc-thesis/results/run_fol_para.sbatch task 0)
run /scratch/ipopa/thesis/FOLn-REST2-50ns folic-acid-clustering-shuyu-nosingle-ref \
  --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb --sdf build/FOL.sdf \
  --torsions-file REST2-run1/torsions.txt --torsions-file REST2-run1/torsions_rigid.txt \
  --torsions-file REST2-run1/torsions_polarH.txt --stride 5 --no-single-cluster \
  --name folic-acid --rest2-config REST2.config \
  --ref-sdf $CRYSTAL/fol_4LRH_A.sdf --ref-sdf $CRYSTAL/fol_3QL3_A.sdf
# paracetamol (task 1)
run /scratch/ipopa/thesis/REST2-PARA paracetamol-clustering-shuyu-nosingle-ref \
  --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb --sdf build/TYL.sdf \
  --stride 1 --no-single-cluster --name paracetamol --rest2-config REST2.config
# complexes (FOL-complex-REST2/scripts/shuyu_analysis.sbatch), own crystal only
for s in 4LRH 3QL3; do
  run /scratch/ipopa/thesis/FOL-complex-REST2/$s/analysis-shuyu folic-acid-$s-complex-clustering-shuyu-nosingle-ref \
    --traj FOL_state0_prod1.nc --top FOL.solute.pdb --sdf FOL.sdf \
    --torsions-file torsions.txt --torsions-file torsions_rigid.txt --torsions-file torsions_polarH.txt \
    --stride 5 --no-single-cluster --name folic-acid-$s-complex --rest2-config ../REST2.config \
    --ref-sdf $CRYSTAL/fol_${s}_A.sdf
done
