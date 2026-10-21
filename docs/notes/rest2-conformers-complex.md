---
name: rest2-conformers-complex
description: "Complex REST2 for the five REST2-conformers ligands (4ibi/4ibk VP35, 4j7d MDM2, 5mwh/5ov8 BRPF1): setup, atom-map decisions, Slurm jobs"
metadata:
  type: project
---

Started 2026-10-12 at `/scratch/ipopa/thesis/REST2-conformers-complex/`. It is the complex counterpart of
[[rest2-conformers-dataset]] and follows the [[fol-complex-rest2-plan]] workflow. Everything runs via
Slurm ([[slurm-always]]).

User requirements: **reuse the exact water-run ligand packages** (no re-parameterisation, so the runs
are comparable with water); **40 ns** per state (water ran 60 ns); 6 replicas, tau 0.4, as in water.

- Full PDB entries were downloaded from RCSB into `crystal/`: bsc-thesis only had the ligand SDFs.
  Each complex is chain A plus its chain-A ligand. Additives (DMS, GOL, SO4, NO3), chain B and
  waters are dropped. Every ligand sits in its own chain's pocket (VP35: one chain-B residue
  within 4 Å, a crystal contact).
- Packages are byte-identical copies in `packages/<name>-water/`.
- atom_map decisions (`<name>/prep/atom_map.txt`):
  - 4ibi/4ibk: the COOH was refined with equal C–O bonds. The OH goes on the oxygen NOT bound to
    Lys251 NZ (4ibi: OAD, 4ibk: OAA).
  - 4j7d: the imidazoline skeleton is symmetric without H. Crystal C2–N1 1.388 vs C2–N3 1.282 Å,
    so N1 carries the H (package N2).
  - 5mwh/5ov8 need no map.
- Bug fixed in `scripts/atom_map.py`: it must match on the bare skeleton, because bond orders
  assigned to the crystal pin one carboxyl O as C=O. Added `--hydroxyl`.
- Ligand residues: 4ibi :124, 4ibk :124, 4j7d :87, 5mwh :113, 5ov8 :112.
- Scripts: `scripts/build.sbatch` (CPU array 0-4; writes REST2-run1 with 20000 exchanges and
  test/REST2-test with 300) and `scripts/rest2.sbatch` (gpu:3, 6 ranks, MPS).
- Tests (jobs 8870-8874, 300 exchanges) passed on 2026-10-12. Every acceptance pair is
  0.31-0.51 (5mwh is the lowest, overall 0.34), so 6 replicas are fine. Every ligand stays
  bound (state-0 COM shift <= 2.6 Å). Speed: 394-445 ns/day per replica at 2 ranks per GPU.
  Results in `logs/check-tests.txt`.
- **Production (40 ns) submitted 2026-10-12 as jobs 8928-8932** (4ibi, 4ibk, 4j7d, 5mwh, 5ov8;
  6 h limit), recorded in `logs/production-jobs.txt`.
- `scripts/check_run.py` selects the ligand by residue name in Python: mdtraj's selection
  language cannot parse names that start with a digit (1D8, 1DE).
