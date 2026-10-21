---
name: fol-rest2-workflow
description: "Folic acid REST2 dataset layout, the traj/topology convention for its analysis scripts, and the distance-scale caveat"
metadata:
  node_type: memory
  type: project
  originSessionId: f50aeccb-185c-4cda-9137-c0df5613115a
  modified: 2026-09-29T08:49:23.046Z
---

Folic acid REST2 dataset: `/scratch/ipopa/thesis/FOLn-REST2` (built 2026-09-29), analysis scripts cloned
to `/scratch/ipopa/thesis/bsc-thesis` (not `~/bsc-thesis` — `$HOME=/home/ipopa` is a separate filesystem
from `/scratch/ipopa`). Env: [[md-tools-openmm-env]].

Run: `REST2-run1`, 8 replicas, tau 0.0..0.4, 5M steps/replica, 28 min on 8x H200 (one rank per
GPU, `CUDA_VISIBLE_DEVICES=0..7`, no MPS). All 7 neighbouring-pair acceptances 0.418-0.454.

## The traj/topology convention (non-obvious)

`select_torsions.py --ligand-indices` writes ligand-local indices 0..50, which are valid against
**both** topologies only because MD-tools happens to build the solute as atoms 0..50 in
`built.pdb`. The analysis steps then need `--top ../build/built.solute.pdb` together with
`solute_state0_prod1.nc` (5000 frames), NOT `built.pdb` + `whole_state0_prod1.nc` (100 frames) —
`md.load()` rejects an atom-count mismatch. Prefer the solute pair: 50x the frames.

## Distance thresholds do not scale with the torsion count

`cluster.py`'s `--assign-radius` (0.2) and Daura's `--cutoff` (0.5) were calibrated on
paracetamol's 2-torsion feature space. Features are scaled 1/sqrt(n), so with 9 torsions the
whole distance scale shifts up — the 1st-percentile pair distance is already 0.571, making both
thresholds effectively no-ops (Daura reported 80% noise for that reason alone). HDBSCAN itself is
scale-adaptive and can be trusted; the fixed-radius methods cannot, unless rescaled.

Folic acid result: 9 rotatable torsions, 1 rigid (amide omega, stays trans, circSD 12.7 deg).
HDBSCAN gives one basin (52%) + 48% noise. That is real, not a tuning failure — 5 torsions are
well-defined (aryl-amide and aryl-NH both planar at ~21-22 deg folded) and 4 are near-free
(2 carboxyls + the NH-CH2-pteridine methylene linkers), so the ensemble is genuinely diffuse.
PCA is correspondingly flat (PC1 14%, PC1-5 53%). Low silhouette (~0.12) confirms no real
sub-clusters; the 3 streaks in `clusters_pc.png` are an artifact of the symmetry-alignment step.
