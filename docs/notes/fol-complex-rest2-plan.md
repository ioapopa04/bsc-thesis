---
name: fol-complex-rest2-plan
description: "Folic-acid complex REST2 (4LRH FRα, 3QL3 DHFR): water-run FOL package reused, 10 rungs, Slurm jobs, prep gotchas"
metadata:
  type: project
---

Complex REST2 for folate bound to **4LRH** (human FRα, chain A) and **3QL3** (E. coli DHFR =
re-refined 1RX2, chain A, altloc A, with NADP+). Built 2026-10-08 at
`/scratch/ipopa/thesis/FOL-complex-REST2/`. Production **completed**: Slurm jobs 8050 (4LRH, 4.1 h, 296 ns/day/replica) and 8051 (3QL3, 3.4 h).
8024/8025 had failed because `SLURM_NTASKS` made the serial stages refuse; rest2.sbatch now unsets it. See
[[fol-bound-confs-clusters]], [[fol-rest2-workflow]], [[slurm-always]].

**Comparability rule (user, 2026-10-08): do NOT re-parameterise folate.** Results must be
comparable with the water runs, so the complexes use the exact water package
`LOCAL-OVBPIULPVIDEAO/param_a76881c75970` (neutral COOH, 1H-pterin, AM1-BCC, openff-2.2.1),
copied byte-identical to `packages/FOL-water` (the whole dir is needed, including FOL.pdb/xml/sdf).
I had first built a −2 dianion package after the user picked "−2" in a plan question; they then
rejected it. That build is archived in `_superseded-dianion/`.
- The neutral COOH needs an explicit `atom_map` (`scripts/atom_map.py`). 4LRH: the OH goes on the
  longer crystal C–O (1.35 vs 1.22 Å). 3QL3 was refined as carboxylates (1.25–1.26 Å), so it is
  assigned by environment: gamma OH on OE2 (free), alpha OH on O2. Both alpha oxygens bind Arg57
  bidentately, so a neutral COOH loses one of those two contacts.
- Ladder: **10 rungs**, tau 0..0.4 (with 8 rungs the −2 test gave 0.23–0.25, worst pair 0.17).
  Hot region: folate only, 51 atoms and 12 torsion bonds. Residues: 4LRH :208, 3QL3 :160 (NAP :161).
- NADP+ (−3) is charged with NAGL: sqm AM1 SCF does not converge for it.
- The scaler needs `build/FOL.sdf`. Its improper handling differs from the water run (0 vs 57
  unscaled impropers), which matters only at hot rungs, not for the state-0 ensemble.

Slurm: `scripts/build.sbatch` (CPU, sbrinz2) and `scripts/rest2.sbatch` (gpu:5, 10 ranks, MPS
private to the job, 12 h limit; resubmitting resumes). Check with `scripts/check_run.py <run>`.
Test speed: 122–139 ns/day per replica at 4 ranks/GPU.

## Results (checked 2026-10-11; `logs/check-production.txt`, `analysis/fragments.txt`)
- Acceptance: 4LRH 0.542–0.577 (overall 0.556); 3QL3 0.521–0.566 (overall 0.548). 10 rungs is ample.
- **4LRH: stays bound.** State-0 pteridine RMSD 0.5–0.7 Å, pABA ~1 Å, glutamate 1.9–3.2 Å.
- **3QL3: the glutamate tail leaves within the first 5 ns.** Arg57 to α-COOH distance is 10 Å by
  0–5 ns and 17 Å later (crystal 2.8 Å); glutamate RMSD is 6–13 Å. The pteridine stays anchored:
  Asp27 is at 2.8 Å throughout. The pteridine shifts after about 20 ns (RMSD 1 to 3 Å; NADP+ contact
  3.4 to 4.5 Å). The likely cause is the neutral α-COOH (the user's comparability choice) losing
  the Arg57 salt bridge, but this is not proven. The crystal tail is also disordered (altloc A/B).
  The user has not yet decided whether to accept this or run a −2 sensitivity check.

## Water-vs-complex torsion figure (2026-10-11)
`bsc-thesis/shuyu_clustering/compare_water_complex.py` (analysis-env, run via srun, job 8714) →
`bsc-thesis/results/folic-acid/folic-acid_water_vs_complexes.{png,csv}`. Torsions are 0-based t0..t11.
Not at crystal: 3QL3 t1, t5, t6; 4LRH t0, t3 (and t6 only 26 %).
