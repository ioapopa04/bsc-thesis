---
name: rest2-conformers-v4
description: Advisor rule "rerun only if protonation changes" applied to the 5 small molecules; v4 runs for 4ibi/4ibk/4j7d, jobs and analysis chain
metadata:
  type: project
---
2026-10-15, advisor: proceed with compounds with fewer rotatable bonds, with the protonation state + new clustering; if protonation is exactly the same as before, no new simulations.
Open Babel pH 7.4 (and 7.0): 4ibi -2, 4ibk -2 (COOH + enol OH), 4j7d +1 (imidazoline N); 5mwh, 5ov8 unchanged (thesis/REST2-conformers-complex/protonation-check/result.txt).
New runs in /scratch/ipopa/thesis/REST2-conformers-v4 (env openmm-env-protonation, same settings as originals: water 6 reps tau0.4 60 ns one run; complex 40 ns ligand-only hot region). Sims 11591-11596 (logs/jobs.txt), v4 analyses 11597-11608 afterok.
Analysis of all original runs (water+complex) with symfirst: bsc-thesis/results/run_small_symfirst.sbatch (SET=orig|v4, NOSINGLE).
sbrinz1 GPUs/MPS work again as of 2026-10-15 (test job 11590).
Related: [[symfirst-clustering]], [[rest2-conformers-complex]]
