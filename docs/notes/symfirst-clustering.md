---
name: symfirst-clustering
description: bsc-thesis shuyu_clustering ported to Shu-Yu's symmetry-first t_hdbscan (MD-tools main 69c262b) in env analysis-env-main; array job run_symfirst.sbatch
metadata:
  type: project
---
2026-10-15: user asked to adopt Shu-Yu's NEW clustering exactly, keeping the thesis figures.
- MD-tools main worktree /scratch/ipopa/software/MD-tools-main (69c262b, 0.6.5), installed non-editable in env analysis-env-main (clone of analysis-env). analysis-env (0.6.4) untouched.
- bsc-thesis/shuyu_clustering/cluster_prune.py + cluster_convergence.py now: t_hdbscan(molecule, torsion_definitions, coordinates, atom_mapping) -> align_all(per_cluster=True) -> representatives(); no TorsionSymmetry merge; allow_single_cluster default True (old runs used --no-single-cluster). Old version in git history. Unmerged overview replaced by clusters_overview_unaligned.png.
- bsc-thesis/results/run_symfirst.sbatch array 0-7 (para v4/orig, folate water/4LRH/3QL3 v4 and orig); outputs *-clustering-symfirst / *-convergence-symfirst. Jobs 11530, 11532.
- Max 10000 frames for the dense matrix (folate water stride 5 = exactly 10000).
Related: [[fol-rest2-v4]], [[bsc-thesis-protocols]]
