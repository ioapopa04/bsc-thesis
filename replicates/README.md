# REST2 replicates from different starting conformers

Question (Shu-Yu): does a REST2 run started from a different RDKit conformer give the same ensemble?
Replicate 0 is the existing v4 water run (start = the package's lowest-MMFF ETKDG conformer, seed 1).
Run data: `/scratch/ipopa/thesis/REST2-replicates/{PARA,FOL}/rep<k>` (not in git).

| step | script | what it does |
|---|---|---|
| 1 | `pick_starts.py` / `.sbatch` | cluster rep0 (symmetry-first t_hdbscan), embed ETKDGv3 + MMFF94s candidates of the **package** molecule, label them with `fit.predict`, pick starts greedily (first an uncovered cluster, then farthest point; min distance 0.1) |
| 2 | `build.sbatch` | build each start with the same package (reused, charges identical, checked), same scaler and REST2.config; only `dynamics.seed` = 1000+k differs |
| 3 | `rest2.sbatch` | the run (MPS, 2 ranks per GPU) |
| 4 | `compare_replicates.py` / `.sbatch` | pool all replicates, cluster ONCE, per-replicate populations (block SE), timeline, convergence, torsions, Jensen-Shannon |

Starts: paracetamol has only two distinct trans-amide minima (= clusters A/B), so rep0 + rep1;
folate rep1 (cluster A), rep2 and rep3 (unassigned), all 3.7-4.5 from rep0 in the symmetry-aware distance.

Results (2026-10-17):

* **Folate in water** (4 x 100 ns per state): A 10-13 %, B 10-12 %, unassigned 76-78 % in every
  replicate, spread within ~1.5 block SE, JS(labels) <= 0.002 bits; torsion distributions overlap.
  Start-independent. Pooled data give 2 clusters (rep0 alone: 3; the ~2 % cluster C does not survive).
* **Paracetamol** (2 x 10 ns per state): same two clusters, A 46 % (rep0) vs 57 % (rep1), ~3 SE
  apart although both switch A<->B constantly; most likely too short. Extension not done yet.

Folders: `<name>/` = start selection (`starts/rep*.sdf`, `candidates.csv`, `selection.json`, figures);
`<name>-compare/` = the comparison (`populations`, `timeline`, `convergence`, `torsions`, `agreement.json`).
