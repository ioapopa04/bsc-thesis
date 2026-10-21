---
name: md-tools-protonation-branch
description: "Pending MD-tools update (branch feature/ligand-protonation-and-ligand-amides): user installs it only after the queue is empty"
metadata:
  type: project
---

On 2026-10-12 the user received Shu-Yu's instructions: git fetch, check out
`feature/ligand-protonation-and-ligand-amides`, `micromamba create -f environment.yml`, `pip install -e .`.
The user chose to **wait until all queued jobs finish** (option 2) before replacing `openmm-env`.
The queued jobs were 5mwh/5ov8 complex REST2 and the paracetamol resume 9486/9487.

Why waiting matters:
- The branch's environment.yml is still named `openmm-env`, which every Slurm script activates.
  Recreating it would hit running jobs.
- The branch changes run-time REMD code (`remd/driver.py`), so a resume across versions is risky.

What the branch changes:
- Ligands are protonated at pH 7 before parameterisation (FOL becomes COO-).
- "REST2 v4" scales ligand amides.
- Both differ from the existing water and complex runs (neutral COOH; amide omega unscaled). Runs
  built with the branch are therefore not Hamiltonian-identical to them. See [[fol-complex-rest2-plan]].

At the time, the installed md-tools was 0.6.2 (a non-editable copy), and checkout main was 0.6.4.
The branch was fetched only; the checkout stayed on main. See [[md-tools-openmm-env]].
