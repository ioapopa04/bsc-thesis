# bsc-thesis — instructions for Claude Code

REST2 (MD-tools/OpenMM) of folate, paracetamol and 5 small ligands in water and bound, plus
torsion clustering. Moved from the sbrinz cluster on 2026-10-21.

- Read `SETUP.md` first (software, paths, Slurm, workflow, dataset index), then `runs/README.md`,
  `replicates/README.md` and `docs/notes/` (project memory from the previous sessions: the why
  behind the protocol choices).
- **Every compute job via the scheduler (sbatch/srun). Never on the login node** — no nohup, `&`,
  setsid. Quick inspection commands are fine.
- Envs by exact name: `openmm-env-protonation` (simulation, MD-tools b3603e0), `analysis-env-main`
  (analysis, MD-tools 69c262b).
- Reuse the committed ligand packages; never re-parameterise. Never edit MD-tools source.
- Paths in scripts may still be sbrinz's `/scratch/ipopa/...` until ported (SETUP.md §2).
- Ask before deleting run data or overwriting results; don't touch `REST2-conformers/extend_all.sh`.
- Commits end with the Co-Authored-By line; push to github.com:ioapopa04/bsc-thesis.
