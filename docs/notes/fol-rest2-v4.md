---
name: fol-rest2-v4
description: "Folate REST2 rerun with Shu-Yu's md-tools branch (COO-, amide scaled): separate env, dataset FOL-REST2-v4, jobs"
metadata:
  type: project
---

On 2026-10-14 the user (via their advisor Shu-Yu) asked to rerun folate in water, 4LRH and 3QL3 with
md-tools branch `feature/ligand-protonation-and-ligand-amides` at b3603e0. This supersedes the old
comparability rule in [[fol-complex-rest2-plan]]: the new runs form a NEW comparable set
(water + complexes share one package). They are not comparable to the old neutral runs.

- Installed side by side ([[md-tools-protonation-branch]]): worktree
  `/scratch/ipopa/software/MD-tools-protonation`, env **`openmm-env-protonation`** (editable,
  reports 0.6.4, openbabel 3.2.1, openmm 8.6). The original `openmm-env` (0.6.2) and the main checkout
  are untouched. The install job is `software/install_protonation_env.sbatch`.
- Dataset `/scratch/ipopa/thesis/FOL-REST2-v4/`: `packages/FOL` = new folate param_80e441460558 (Open Babel
  pH 7.4, net -2, same 1H tautomer as before, AM1-BCC), reused by all three systems. NADP+ reused from
  FOL-complex-REST2 (NAGL). No atom_map is needed: the carboxylates are symmetric.
- The scaler confirms REST2 v4: "ligand amides 1 bond(s), scaled". Settings match the old runs:
  water 8 replicas tau 0.4, now 50000 exchanges = 100 ns in one run; complexes 10 replicas, 50 ns.
  Residues: 4LRH :208, 3QL3 :160.
- Slurm jobs 11474 (water, gpu:4), 11475 (4LRH, gpu:5), 11476 (3QL3, gpu:5); see
  `logs/production-jobs.txt`.
- Analysis: job 11488 (water + 4LRH, done); 3QL3 + RMSD queued afterok:11482. Script `scripts/analysis.sbatch`. v4 torsions regenerated (9 rotatable, 1 rigid amide, 0 polar-H).
  acceptance once the runs start.

Results (2026-10-15): all v4 runs + analyses done (3QL3 11482, ana 11488/11489; paracetamol v4 in /scratch/ipopa/thesis/PARA-REST2-v4, jobs 11490-11492).
3QL3 v4: Arg57 salt bridge held ~15 ns then lost, folate drifts out (pose RMSD mean 4.9 Å vs old 6.6). 4LRH v4 conf RMSD median 0.59 Å (old 0.81). Paracetamol v4 = old within noise (A 52.6 / B 46.1 %).
