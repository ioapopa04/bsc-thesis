#!/usr/bin/env bash
# Cumulative-window clustering for 60/70/80/90/100 ns, continuing the existing
# 10ns..50ns series in fol-clustering-shuyu-nosingle. Same command as before:
# cumulative window = first N*500 solute frames, --stride 5, --no-single-cluster.
set -uo pipefail

ROOT=/scratch/ipopa/FOLn-REST2-50ns
cd "$ROOT" || exit 1

source /scratch/ipopa/folenv.sh
micromamba activate analysis-env

SC=/scratch/ipopa/bsc-thesis/shuyu/cluster_prune.py
OUT=fol-clustering-shuyu-nosingle

for ns in 60 70 80 90 100; do
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
  rc=$?
  end=$(date +%s)
  echo "    ${ns} ns: rc=$rc, $(( end - start )) s"
done
echo "=== $(date -Is) all windows done ==="
