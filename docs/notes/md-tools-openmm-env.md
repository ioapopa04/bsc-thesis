---
name: md-tools-openmm-env
description: Where the MD-tools/OpenMM stack lives on this machine and how to activate it
metadata:
  node_type: memory
  type: project
  originSessionId: 9c2d1b0e-65ba-4412-8c0e-189c00b622c4
  modified: 2026-09-29T08:11:30.762Z
---

MD-tools + OpenMM installed 2026-09-29 on scratch, not `$HOME` (the upstream docs assume
`~/software/md-stack`; here `$HOME=/home/ipopa` is separate from `/scratch/ipopa`):

- micromamba root + binary: `/scratch/ipopa/software/md-stack` (`bin/micromamba`, v2.9.0)
- env: `/scratch/ipopa/software/md-stack/envs/openmm-env` — activate with `micromamba activate openmm-env`
- checkout (v0.6.2, commit 2e1aa6c): `/scratch/ipopa/software/MD-tools`
- `MD_DATA=/scratch/ipopa/md-data`
- user.config: `~/.config/md-tools/user.config` — person_id/name are still the placeholder `ipopa`,
  safe to change only until the first `data-register` (changing it later splits datasets by person)
- `~/.bashrc` carries the `micromamba shell init` block plus the `MD_DATA` export

`md-openmm` is not on PATH until the env is activated — intended, not a misconfiguration.

Hardware: 8x NVIDIA H200 NVL (143 GB each), driver 595.58.03 (CUDA 13.2) — host `sbrinz`, a
~22-user shared login node, not an exclusive compute node. (An earlier version of this memory said
2 GPUs; `nvidia-smi` on 2026-09-29 listed 8.) The `cuda-version=13.0` pin in
`environment.yml` is deliberate (PTX is backward- but not forward-compatible); resolved
`cuda-nvrtc 13.0.88` compiles and runs on this driver — verified, not assumed.

The docs are in the checkout at `docs/install.md` and `docs/structure/machine-config.md`;
the published site 404s at `/machine-configuration/` (the real path is `/structure/machine-config/`).
