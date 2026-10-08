# Shu-Yu's torsion-clustering workflow for REST2

Scripts live in `~/bsc-thesis/shuyu`. Run them from the REST2 run directory so input paths
stay relative to that simulation. This protocol clusters the ligand torsions with MD-tools'
`t_hdbscan` and prunes symmetry-equivalent clusters with MD-tools' `TorsionSymmetry`, then
summarises per-group masses, torsion statistics and representative structures.

Requires the MD-tools version that contains `md_tools.analysis` (available in the
`analysis-env` micromamba environment), plus scikit-learn, rdkit and pandas.

## 0. Setup

```bash
source /scratch/ipopa/folenv.sh          # or: micromamba activate analysis-env
cd ~/FOLn-REST2-50ns                     # run directory (contains build/ and REST2-run1/)

SHUYU="$HOME/bsc-thesis/shuyu"
RES=FOL
```

## 1. Select the torsions

```bash
python "$SHUYU/select_torsions.py" \
    --sdf build/FOL.sdf --top build/built.pdb --select "resname FOL" \
    --ligand-indices --polar-h --out REST2-run1/torsions
```

Writes `REST2-run1/torsions.txt` (rotatable), `torsions_rigid.txt`, `torsions_polarH.txt`
(polar O–H/N–H torsions) and `torsions_ring.txt`. `--ligand-indices` numbers the torsions
`0..n-1` within the ligand, matching a solute-only trajectory.

## 2. Cluster and prune by symmetry

```bash
python "$SHUYU/cluster_prune.py" \
    --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb --sdf build/FOL.sdf \
    --torsions-file REST2-run1/torsions.txt \
    --torsions-file REST2-run1/torsions_rigid.txt \
    --torsions-file REST2-run1/torsions_polarH.txt \
    --stride 5 --no-single-cluster --out fol-clustering-shuyu-nosingle
```

Outputs in `--out/`: `torsion_cv.json`, `merged_groups.json`, `groups.csv`,
`representative_<A/B/…>.pdb`, `group_<i>.png`, `group_torsions.png`,
`representatives.png`, `torsions_<i>_<j>_2d.png`.

`--no-single-cluster` runs `t_hdbscan` with `allow_single_cluster=False` (a unimodal
ensemble is reported as noise rather than forced into one cluster). `--stride 5` subsamples
the dense solute trajectory; pass the same stride if you compare runs.

Optional labelling / reference flags:

- `--name 4ibi` — compound name in every figure title.
- `--rest2-config REST2.config` — replicas and τ_max for the overview title; the simulated
  time comes from the trajectory time stamps (override with `--duration-ns`).
- `--ref-sdf ../crystal/4ibi_1D8_A.sdf` (repeatable) — crystal (bioactive) ligand; its
  torsions are drawn as black vertical lines (solid, dashed, … per file) in
  `group_torsions.png` and `clusters_overview.png`. Matched to `--sdf` on heavy-atom element +
  connectivity, so polar-H torsions get no line (no H in crystal structures), and a torsion
  through atoms the crystal cannot tell apart (the two carboxyl O) gets a line per assignment
  (e.g. 0° and 180°). Crystal ligands live in `../crystal/` (from RCSB, one SDF per chain).

`clusters_overview.png` combines the `group_<i>.png` files: molecule on the left (each torsion's
bonds in its own colour, `t<k>` on the central bond), one row per torsion — the full
distribution first, then one panel per merged group — titled with compound, simulation
length, replicas and τ_max.

## Convergence

`cluster_convergence.py` takes the same torsion/clustering flags as `cluster_prune.py`, re-runs
the clustering + symmetry pruning on cumulative windows 0–5, 0–10, … ns (`--step-ns`), and
writes `convergence.csv` / `convergence.png` (merged groups vs time, raw HDBSCAN count and noise
fraction alongside). Converged = `--plateau` (default 3) consecutive windows with the same
number of groups. Several `--traj` segments (run1, ext1, …) are concatenated in order.

## Cumulative windows

`run_windows.sh` / `run_windows_labels.sh` in the run directory re-run the series of
cumulative windows (10…100 ns) with the same settings, into subfolders `10ns/ … 100ns/`.
They call this script as `/scratch/ipopa/bsc-thesis/shuyu/cluster_prune.py`.

The 100 ns state-0 window on the FOL run gives **2 clusters (A, B) and ~87.5 % noise**
(≈13 % of conformations clustered). Compare against Ioana's protocol (`../ioana/`) for a
clustering that assigns the lower-density frames to basins and reports block-bootstrap
population errors and cluster stability.
