#!/usr/bin/env bash
# Re-run all cumulative windows with the population added to the group_torsions
# legend (cluster_prune.py label change). Same command/settings as the series.
set -uo pipefail

ROOT=/scratch/ipopa/thesis/FOLn-REST2-50ns
cd "$ROOT" || exit 1

source /scratch/ipopa/thesis/folenv.sh
micromamba activate analysis-env

SC=/scratch/ipopa/thesis/bsc-thesis/shuyu/cluster_prune.py
OUT=fol-clustering-shuyu-nosingle

for ns in 10 20 30 40 50 60 70 80 90 100; do
  if [ "$ns" = 100 ]; then
    traj="REST2-run1/solute_state0_prod1.nc"
  else
    traj="_window_trajs/solute_state0_${ns}ns.nc"
  fi
  echo "=== $(date -Is) ${ns} ns ==="
  start=$(date +%s)
  python "$SC" \
    --traj "$traj" --top build/built.solute.pdb --sdf build/FOL.sdf \
    --torsions-file REST2-run1/torsions.txt \
    --torsions-file REST2-run1/torsions_rigid.txt \
    --torsions-file REST2-run1/torsions_polarH.txt \
    --stride 5 --no-single-cluster --out "${OUT}/${ns}ns" \
    > "${OUT}/${ns}ns.log" 2>&1
  echo "    ${ns} ns: rc=$? , $(date -Is), $(( $(date +%s) - start )) s"
done
echo "=== $(date -Is) all windows re-run ==="
