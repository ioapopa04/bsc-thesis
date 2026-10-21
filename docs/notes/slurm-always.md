---
name: slurm-always
description: User rule — always launch jobs on sbrinz via Slurm (module load slurm/slurm/24.11), never directly on the login node
metadata:
  type: feedback
---

Always launch compute jobs through Slurm (`sbatch`/`srun`) on this machine. Load it with
`module load slurm/slurm/24.11` (not on PATH by default). Partition `defq`; GPU nodes
sbrinz/sbrinz1 (gpu:8 each), CPU nodes sbrinz[2-3]. Also recorded in `~/.claude/CLAUDE.md`.

**Why:** user asked explicitly (2026-10-08); earlier sessions wrongly assumed there was no
Slurm and ran MD on the shared login node (see [[rest2-conformers-dataset]]).
**How to apply:** any MD run, build, clustering or other non-trivial compute → write an sbatch
script with `--gres=gpu:N` as needed; don't use nohup/`--detach`/GPU-polling launchers on the login node.
