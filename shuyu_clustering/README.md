# Shu-Yu's torsion-clustering workflow for REST2

Scripts live in `~/bsc-thesis/shuyu`. Run them from the REST2 run directory so input paths
stay relative to that simulation. This protocol clusters the ligand torsions with Shu-Yu's
**symmetry-first** `t_hdbscan` (MD-tools `main`, 69c262b, 2026-10-15): the molecular symmetry is
enumerated first and built into the distance (minimum over atom relabellings), HDBSCAN runs on
that distance, and every frame is classified by an 18-of-20 neighbour vote (`-2` density noise,
`-1` ambiguous). There is **no merge step** — each cluster is already one conformer up to
symmetry. For the figures every frame is written in one labelling with
`align_all(per_cluster=True)` (each cluster aligned to its own density peak, each peak to the
density peak of cluster A); representatives come from `fit.representatives()`. The previous
cluster-then-merge version (`TorsionSymmetry`) is in git history before this change.

Requires the `analysis-env-main` micromamba environment (MD-tools 0.6.5 from
`/scratch/ipopa/software/MD-tools-main`), plus scikit-learn, rdkit and pandas. The dense
symmetry-aware distance matrix limits a fit to 10 000 frames (`--max-frames` raises it).

## 0. Setup

```bash
micromamba activate analysis-env-main
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

## 2. Cluster (symmetry-first)

```bash
python "$SHUYU/cluster_prune.py" \
    --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb --sdf build/FOL.sdf \
    --torsions-file REST2-run1/torsions.txt \
    --torsions-file REST2-run1/torsions_rigid.txt \
    --torsions-file REST2-run1/torsions_polarH.txt \
    --stride 5 --out fol-clustering-symfirst
```

Outputs in `--out/`: `torsion_cv.json`, `summary.json` (MD-tools' full fit record),
`merged_groups.json` (clusters, masses, noise/ambiguous, representatives, density peaks),
`groups.csv`, `representative_<A/B/…>.pdb`, `group_torsions.png`, `clusters_overview.png`,
`clusters_overview_unaligned.png`, `representatives.png`, `torsions_<i>_<j>_2d.png`.

`allow_single_cluster` is True by default, as in Shu-Yu's tutorial; `--no-single-cluster`
sets it False (the setting used before 2026-10-15). `--stride 5` subsamples
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

`clusters_overview.png`: molecule on the left (each torsion's
bonds in its own colour, `t<k>` on the central bond), one column per torsion — the full
distribution in the top row, then one row per cluster (aligned torsions); with many torsions they wrap into 2–3
bands, with many groups the groups wrap into 2–3 side-by-side blocks, so the figure stays near 16:9 — titled with compound, simulation
length, replicas and τ_max. `clusters_overview_unaligned.png` is the same with each cluster's
torsions as stored (not relabelled), where a conformer held in two labellings shows both images.

## Convergence

`cluster_convergence.py` takes the same torsion/clustering flags as `cluster_prune.py`, re-runs
the symmetry-first clustering on cumulative windows 0–5, 0–10, … ns (`--step-ns`), and
writes `convergence.csv` / `convergence.png` (clusters vs time, unassigned fraction alongside). Converged = `--plateau` (default 3) consecutive windows with the same
number of groups. Several `--traj` segments (run1, ext1, …) are concatenated in order.

## Cumulative windows

`run_windows.sh` / `run_windows_labels.sh` in the run directory re-run the series of
cumulative windows (10…100 ns) with the same settings, into subfolders `10ns/ … 100ns/`.
They call this script as `/scratch/ipopa/bsc-thesis/shuyu/cluster_prune.py`.

With the earlier cluster-then-merge version, the 100 ns state-0 window on the FOL run gave **2 clusters (A, B) and ~87.5 % noise**
(≈13 % of conformations clustered). Compare against Ioana's protocol (`../ioana/`) for a
clustering that assigns the lower-density frames to basins and reports block-bootstrap
population errors and cluster stability.
