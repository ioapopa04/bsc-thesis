---
name: rest2-replicates
description: REST2 start-conformer replicates (Shu-Yu's request): pick_starts.py, PARA rep1 + FOL rep1-3 in thesis/REST2-replicates, jobs 12622-12626
metadata:
  type: project
---
Shu-Yu (2026-10-17): run replicates from different RDKit conformers to check REST2 gives the same result.
- Scripts in thesis/bsc-thesis/replicates/: pick_starts.py (clusters rep0 with symfirst, ETKDG+MMFF candidates of the PACKAGE molecule, fit.predict labels, greedy new-cluster then maximin, min distance 0.1), pick_starts.sbatch, build.sbatch (reuses package by path, checks RMSD to start + identical charges), rest2.sbatch.
- rep0 = existing v4 water run (seed 1). Replicates: seed 1000+k. Runs in /scratch/ipopa/thesis/REST2-replicates/{PARA,FOL}/rep<k>.
- Paracetamol has only 2 distinct trans-amide minima (= clusters A/B); user chose only 2 starts (rep0 B, rep1 A), no cis-amide start.
- Folate: rep1 cluster A, rep2/rep3 noise, all far from rep0 (d 3.7-4.5). No candidates in clusters B/C.
- Jobs: build 12622 (ok), rest2 PARA1 12623, FOL1-3 12624-12626. Analysis to do next: pool reps, cluster once, per-rep populations. See [[symfirst-clustering]], [[fol-rest2-v4]].
- 2026-10-17 all 4 replicate runs completed. Comparison: replicates/compare_replicates.py (+ .sbatch, pooled single symfirst fit, resampling=10000 for folate), outputs replicates/<name>-compare/.
  - Folate: pooled 2 clusters; A 10-13 %, B 10-12 %, unassigned 76-78 % in all 4 reps; spreads within ~1.5 block SE; JS(labels) <= 0.002 bits -> start-independent.
  - Paracetamol (10 ns): A 46 % (rep0) vs 57 % (rep1), ~3 combined SE apart, both switch A<->B constantly -> likely too short; extension to e.g. 50 ns suggested, not done.
