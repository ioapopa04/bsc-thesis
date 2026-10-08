# Ioana's torsion-analysis workflow for REST2

The reusable scripts and their shared symmetry helper (`symtools.py`) live in
`~/bsc-thesis/ioana`. Run them from the relevant REST2 run directory so input paths stay
relative to that simulation and generated files (torsion lists, plots, tables, structures)
stay with its outputs.

## 0. Setup (once per session)

```bash
micromamba activate openmm-env
cd ~/PARA-REST2-long/REST2-run1

# Scripts are maintained separately from simulation outputs.
IOANA="$HOME/bsc-thesis/ioana"

# Set these for the ligand and run:
RES=PAR                                # ligand residue name in built.pdb
TOP=../build/built.pdb                 # topology written by MD-tools
SDF=../build/PAR.sdf                   # ligand SDF with bond orders
TRAJ=whole_state0_prod1.nc             # state 0 = unscaled, physical ensemble
```

**Pick `whole_` vs `solute_` deliberately.** They are written at different intervals:
`whole_state*_prod1.nc` holds the solvent box and is often very short (e.g. 100 frames),
while `solute_state*_prod1.nc` is ligand-only and dense (e.g. 5000 frames). For torsion
analysis/clustering you normally want the dense solute trajectory with the matching solute
topology (`built.solute.pdb`), not `built.pdb` — the two topologies differ in atom count and
`md.load` will refuse the mismatch. Check both with
`python -c "import mdtraj as md; print(md.open('solute_state0_prod1.nc').n_frames)"`.

Check the actual SDF name with `ls ../build/*.sdf`. To find the residue name for a new ligand:

```bash
python -c "import mdtraj as md; print({r.name for r in md.load_topology('$TOP').residues})"
```

## 1. Select and classify torsions

```bash
python "$IOANA/select_torsions.py" --sdf "$SDF" --top "$TOP" \
    --select "resname $RES" --out torsions
```

This checks that SDF and topology atom orders match, then writes:

- `torsions.txt`: rotatable torsions and an informational per-bond symmetry number
- `torsions_symops.txt`: whole-molecule symmetry operations acting jointly on all selected torsions; this is the symmetry input used by clustering and PCA
- `torsions_rigid.txt`: amides, double bonds, esters, and triple-bond neighbours
- `torsions_ring.txt`: non-aromatic ring bonds (not analysed by this workflow yet)

Trivial rotors such as CF₃ are dropped. Check the printed classification table for chemical
sense before continuing. Use the matching `torsions.txt` and `torsions_symops.txt` from this
same run together.

## 2. Torsion distributions and diagnostics

```bash
python "$IOANA/torsion_analysis.py" analyse \
    --traj "$TRAJ" --top "$TOP" --torsions torsions.txt --out results
```

Outputs in `results/`:

| Output | What to inspect |
|---|---|
| `torsions_hist.png` | Marginal distribution and free-energy profile for each torsion |
| `jsd_halves.csv` | Convergence: first-half vs second-half Jensen–Shannon divergence; smaller is better |
| `mi_matrix.csv` / `mi_shuffled.csv` | Torsion coupling and shuffled finite-sample baseline |
| `torsion_pairs.png` | 2D scatter of each torsion pair |
| `top_mi_pair_2d.png` | 2D histogram of the most correlated pair |
| `dihedrals.npy` | Raw angles used by `dpca.py` and `cluster.py` |

Optional sanity check: rigid torsions should generally remain in one state (for example, a trans amide).

```bash
python "$IOANA/torsion_analysis.py" analyse \
    --traj "$TRAJ" --top "$TOP" --torsions torsions_rigid.txt --out results_rigid
```

Optional REST2 comparison: inspect whether torsion distributions broaden in higher-temperature/scaled states. State 0 remains the physical reference ensemble.

```bash
python "$IOANA/torsion_analysis.py" analyse \
    --traj "$TRAJ" --top "$TOP" --torsions torsions.txt \
    --compare whole_state1_prod1.nc whole_state2_prod1.nc whole_state3_prod1.nc \
    --out results_compare
```

## 3. Cluster conformations

```bash
python "$IOANA/cluster.py" --traj "$TRAJ" --top "$TOP" \
    --torsions torsions.txt --symops torsions_symops.txt \
    --select "resname $RES" --out clusters
```

Outputs in `clusters/` include:

- `clusters.csv`: frame counts, populations **with a block-bootstrap error**, a **stability** score,
  relative ΔG, the torsion indices used (`torsions` column), mean torsions, and representative frames
- `clusters_pc.png`: PC1 vs PC2, coloured by cluster (noise is grey)
- `clusters_torsions.png`: torsion distributions per cluster
- `cluster_<n>.pdb`: representative structure for each cluster
- `labels.npy`: cluster label for every frame
- `aligned_dihedrals.npy`: symmetry-aligned torsions per frame

### How the distance is built (this decides whether you get clusters at all)

Each frame is compared under every symmetry operation and the smallest combined distance wins:
`d(A,B) = min_g combine_j chord(phi_j(A), phi_j(g B))`, with `chord = 2|sin(dphi/2)|` in `[0, 2]`
(one torsion: `d = 0.5` ↔ 29°, `d = 1.0` ↔ 60°, `d = 2.0` ↔ 180°).

- `--metric max` (**default**): `combine` = the most-changed torsion. A real single-torsion flip is the
  whole distance, so it is not diluted by thermal jitter elsewhere.
- `--metric rms`: `combine` = sqrt(sum of squares). Each torsion contributes `1/n` of the squared
  distance, so one 180° flip is worth about as much as uniform small jitter in *all* torsions. With this
  metric a flexible molecule has no density gap and a density method returns a single blob — which is why
  it is no longer the default.

Torsion selection (the other half of the same problem):

- `--drop-free` (**default on**): remove *near-free rotors*. A torsion whose 1D free-energy profile never
  rises `--min-barrier` kcal/mol (default 2.0) above its minimum defines no metastable state; it only
  adds a dimension and thermal smear. The run prints a `Gmax` table showing which torsions were kept and
  why, so the choice is visible and can be overridden.
- `--torsions-keep 2,3,4` / `--torsions-drop 0,8`: explicit control, e.g. to cluster on a chemically
  chosen subspace (the PABA/amide torsions) instead of everything the classifier kept.

Useful options to test cluster stability:

- `--method daura --cutoff 0.5` for Zivanovic-style Daura clustering
- `--min-cluster-size`, `--min-samples`, and `--merge-gap` for HDBSCAN
- `--allow-single-cluster` restores the old HDBSCAN behaviour. It is **off by default** because with it
  on, a molecule that has no density contrast is silently reported as one cluster rather than as noise.
- `--blocks` sets the number of contiguous blocks for the population error (0 = off).
- `--stability N` (default 5) re-clusters `N` windows of `--stability-frac` (default 50%) of the
  trajectory with the same settings and reports, per cluster, the fraction of a window's frames that
  still land in one cluster. This is the honest test of whether a cluster is a *state*: a cluster that
  reappears in every window is real, one that only exists over the whole trajectory at once is not.
  Set `--stability 0` to skip it (it re-runs the clustering `N` times).

**Read `population_err`, not just `population`.** The block-bootstrap error resamples contiguous chunks of
the trajectory, so it does not pretend that 2 ps frames are independent. If a cluster's population is not
larger than its error, its abundance is not resolved by the trajectory. A small `jsd_halves.csv` value is
*not* sufficient evidence of convergence: JSD on fine bins is dominated by the shape inside a well and
stays small even when probability moves between wells.

## 3b. Shu-Yu-compatible report (same output files as `shuyu/cluster_prune.py`)

Writes the same file set as Shu-Yu's protocol, but from **Ioana's cluster labels**
(`clusters/labels.npy`), so the two protocols can be compared file-for-file. MD-tools'
`TorsionSymmetry` does the symmetry pruning; Ioana's native `clusters/` outputs are kept
alongside, and the dPCA (`dpca.py`) is unchanged.

```bash
python "$IOANA/group_report.py" --traj "$TRAJ" --top "$TOP" --sdf "$SDF" \
    --torsions torsions.txt --labels clusters/labels.npy --stride 1 --out groups
```

Outputs in `groups/`: `torsion_cv.json`, `merged_groups.json`, `groups.csv`,
`representative_<A/B/…>.pdb`, `group_<i>.png`, `group_torsions.png`, `representatives.png`,
`torsions_<i>_<j>_2d.png`. Pass `--stride`/`--skip` and the torsion-selection flags exactly as
given to `cluster.py`, or the labels will not line up with the frames.

## 4. Dihedral PCA

Run PCA after clustering so each frame can be aligned consistently to its cluster's representative. Use the same trajectory, topology, torsions, symmetry operations, frame filtering, **and torsion selection** (`--drop-free`/`--min-barrier`/`--torsions-keep`/`--torsions-drop`) as the clustering command — the defaults match `cluster.py`, so only copy them across if you changed them there:

```bash
python "$IOANA/dpca.py" --traj "$TRAJ" --top "$TOP" \
    --torsions torsions.txt --symops torsions_symops.txt \
    --labels clusters/labels.npy --out pca
```

Outputs in `pca/` include `pc1_pc2.png`, `pc1_pc2_free_energy.png`, `explained_variance.png`, `loadings.csv`, `pcs.npy`, and `aligned_dihedrals.npy`. The PCA uses the `clusters/labels.npy` labels to align symmetry-equivalent conformations within each cluster. If using non-default `--skip` or `--stride`, pass the same values to both `cluster.py` and `dpca.py` so the labels correspond to the same frames.

## 5. Compare against Shu-Yu's clustering (frame-level cross-tabulation)

`compare_clusters.py` recomputes Shu-Yu's labels deterministically (MD-tools `t_hdbscan`
on `torsions.txt` + `torsions_rigid.txt` + `torsions_polarH.txt`, `allow_single_cluster=False`,
same `--stride` as `cluster_prune.py`) and cross-tabulates them against Ioana's labels.

```bash
python "$IOANA/compare_clusters.py" --traj "$TRAJ" --top "$TOP" \
    --sdf "$SDF" --stride 5 --labels clusters/labels.npy \
    --torsions-file torsions.txt --torsions-file torsions_rigid.txt \
    --torsions-file torsions_polarH.txt --out crosstab.csv
```
