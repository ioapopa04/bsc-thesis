# Torsion-analysis workflow for REST2

The reusable analysis scripts live in `~/bsc-thesis`. Run them from the relevant REST2 run directory so input paths stay relative to that simulation and generated files (torsion lists, plots, tables, and structures) stay with its outputs.

## 0. Setup (once per session)

```bash
micromamba activate openmm-env
cd ~/PARA-REST2-long/REST2-run1

# Scripts are maintained separately from simulation outputs.
SCRIPTS="$HOME/bsc-thesis"

# Set these for the ligand and run:
RES=PAR                                # ligand residue name in built.pdb
TOP=../build/built.pdb                 # topology written by MD-tools
SDF=../build/PAR.sdf                   # ligand SDF with bond orders
TRAJ=whole_state0_prod1.nc             # state 0 = unscaled, physical ensemble
```

Check the actual SDF name with `ls ../build/*.sdf`. To find the residue name for a new ligand:

```bash
python -c "import mdtraj as md; print({r.name for r in md.load_topology('$TOP').residues})"
```

## 1. Select and classify torsions

```bash
python "$SCRIPTS/select_torsions.py" --sdf "$SDF" --top "$TOP" \
    --select "resname $RES" --out torsions
```

This checks that SDF and topology atom orders match, then writes:

- `torsions.txt`: rotatable torsions and a symmetry number for each
- `torsions_rigid.txt`: amides, double bonds, esters, and triple-bond neighbours
- `torsions_ring.txt`: non-aromatic ring bonds (not analysed by this workflow yet)

Trivial rotors such as CF₃ are dropped. Check the printed classification table for chemical sense before continuing.

## 2. Torsion distributions and diagnostics

```bash
python "$SCRIPTS/torsion_analysis.py" analyse \
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
python "$SCRIPTS/torsion_analysis.py" analyse \
    --traj "$TRAJ" --top "$TOP" --torsions torsions_rigid.txt --out results_rigid
```

Optional REST2 comparison: inspect whether torsion distributions broaden in higher-temperature/scaled states. State 0 remains the physical reference ensemble.

```bash
python "$SCRIPTS/torsion_analysis.py" analyse \
    --traj "$TRAJ" --top "$TOP" --torsions torsions.txt \
    --compare whole_state1_prod1.nc whole_state2_prod1.nc whole_state3_prod1.nc \
    --out results_compare
```

## 3. Dihedral PCA

```bash
python "$SCRIPTS/dpca.py" \
    --dihedrals results/dihedrals.npy --torsions torsions.txt --out pca
```

Outputs in `pca/` include `pc1_pc2.png`, `pc1_pc2_free_energy.png`, `explained_variance.png`, and `loadings.csv` (which shows the torsion features driving each PC).

## 4. Cluster conformations

```bash
python "$SCRIPTS/cluster.py" --traj "$TRAJ" --top "$TOP" \
    --torsions torsions.txt --select "resname $RES" --out clusters
```

Outputs in `clusters/` include:

- `clusters.csv`: frame counts, populations, relative ΔG, mean torsions, and representative frames
- `clusters_pc.png`: PC1 vs PC2, coloured by cluster (noise is grey)
- `clusters_torsions.png`: torsion distributions per cluster
- `cluster_<n>.pdb`: representative structure for each cluster
- `labels.npy`: cluster label for every frame

Useful options to test cluster stability:

- `--method daura --cutoff 0.5` for Zivanovic-style Daura clustering
- `--min-cluster-size`, `--min-samples`, and `--merge-gap` for HDBSCAN

## 5. Review results locally

From your laptop, copy the generated analysis folders from the remote run directory (adjust the host if needed):

```bash
scp -r popa0002@pharma-jay:~/PARA-REST2-long/REST2-run1/{results,pca,clusters} .
```

Open the PNGs and load `clusters/cluster_*.pdb` in PyMOL to verify that each cluster represents a sensible conformation. The scripts remain in `~/bsc-thesis`; only generated data and analysis results belong in the simulation run directory.
