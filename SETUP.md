# Rebuilding the thesis setup on a new machine

Written 2026-10-21 on sbrinz (the original cluster), for whoever (person or Claude Code instance)
sets this up elsewhere. Read it top to bottom once before running anything.

What this repo holds and what it does not:

| in git | not in git (stays on sbrinz or the user's laptop) |
|---|---|
| analysis scripts (`shuyu_clustering/`, `ioana_clustering/`, `replicates/`, `results/*.sbatch`) | trajectories (`*.nc`, ~TB) |
| every run's inputs: configs, parameter packages, prepared complexes, sbatch scripts (`runs/`) | built systems (`built.xml/pdb`, rebuildable from `runs/`) |
| final figures (`results/`) | per-run analysis outputs (CSV, representatives) |
| exact environment specs (`envs/`) | the conda environments themselves |
| crystal ligands (`crystal/`), project notes (`docs/notes/`) | |

The results in `results/` are final. Do not rerun anything just to "have it" on the new machine.
Rerun only what the user asks for (new replicates, extensions, new molecules).

---

## 1. Software

### 1.1 MD-tools (Shu-Yu's package, the MD and analysis engine)

Repo: <https://github.com/csy0000/MD-tools> (private or semi-private; the user or Shu-Yu has access).
Two pinned commits are used:

| role | commit | where it came from | env |
|---|---|---|---|
| **simulation** (build-top, scaler, build-md, REST2 runs) | `b3603e0` (2026-10-13) | branch `feature/ligand-protonation-and-ligand-amides`; now also in `main` | `openmm-env-protonation` |
| **analysis** (symmetry-first `t_hdbscan`, `align_all`, `predict`) | `69c262b` (2026-10-15) | `main` | `analysis-env-main` |

**Keep the environment names exactly as written.** Every sbatch script activates them by name.

All current results ("v4", "mace", replicates) use these two commits.

Older runs used other versions, which are needed only to reproduce the old results exactly:
- **Original REST2:** md-tools 0.6.2, commit `2e1aa6c`, env `openmm-env`.
- **Old analysis:** 0.6.4 in env `analysis-env`.

Do not install the older versions unless asked.

### 1.2 Install (on a compute node, through the scheduler)

Steps 1–3 below follow `docs/install.md` in the MD-tools checkout; that file is authoritative if
anything below disagrees with it. Run every install step in a batch job, never on a login node.
The templates are `envs/install_protonation_env.sbatch` and `envs/install_analysis_main_env.sbatch`.

```bash
SOFT=/path/on/new/cluster/software            # choose; it replaces /scratch/ipopa/software
mkdir -p $SOFT/md-stack && cd $SOFT/md-stack
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj bin/micromamba
export MAMBA_ROOT_PREFIX=$SOFT/md-stack
eval "$($SOFT/md-stack/bin/micromamba shell hook -s bash)"

git clone https://github.com/csy0000/MD-tools.git $SOFT/MD-tools
git -C $SOFT/MD-tools worktree add $SOFT/MD-tools-protonation b3603e0
git -C $SOFT/MD-tools worktree add $SOFT/MD-tools-main        69c262b

# simulation env
micromamba create -y -n openmm-env-protonation -f $SOFT/MD-tools-protonation/environment.yml
micromamba activate openmm-env-protonation
pip install --no-deps -e $SOFT/MD-tools-protonation

# analysis env: the same stack + md-tools main with the analysis extra (+ networkx for ioana_clustering)
micromamba create -y -n analysis-env-main -f $SOFT/MD-tools-main/environment.yml
micromamba activate analysis-env-main
micromamba install -y -c conda-forge scikit-learn networkx matplotlib
pip install --no-deps $SOFT/MD-tools-main
```

**The CUDA pin.** `environment.yml` pins `cuda-version=13.0`. It needs a driver that supports
CUDA ≥ 13.0; `nvidia-smi` shows this in its top right corner. sbrinz had driver 595.58 (CUDA 13.2).
- If the new cluster's driver is older, do not just raise or lower the pin. Read the header of
  `environment.yml` first, then ask the user.
- The failure looks like `CUDA_ERROR_UNSUPPORTED_PTX_VERSION`.

**Exact versions.** The exact versions on sbrinz are in `envs/*.lock.yml` (`micromamba env export`)
and `envs/*.pip-freeze.txt`.
- Use them to check what was resolved, or as a fallback with
  `micromamba create -n <name> -f envs/<name>.lock.yml`.
- Do not trust the `prefix:` line in those files.
- Key versions: python 3.12, openmm 8.6, openbabel 3.2.1, AM1-BCC via AmberTools `sqm`.

### 1.3 Verify, in a short GPU job

```bash
micromamba activate openmm-env-protonation
md-openmm --version; command -v sqm mpiexec nvidia-cuda-mps-control
python -c "from openmm import Platform as P; print([P.getPlatform(i).getName() for i in range(P.getNumPlatforms())])"   # must list CUDA
python -c "from md_tools.openmm.ligand_protonation import openbabel_version as v; print(v())"
micromamba activate analysis-env-main
python -c "from md_tools.analysis import t_hdbscan, THDBSCAN; print(hasattr(THDBSCAN,'align_all'))"   # True
```

### 1.4 MD-tools machine config and data store

`md-openmm data-register --init` writes `~/.config/md-tools/user.config`. sbrinz had:

```yaml
user:    {person_id: ipopa, name: ipopa}
machine: {md_data: /scratch/ipopa/md-data, openmm: {platform: CUDA, precision: mixed, device_policy: local_rank}}
```

- **Keep `person_id: ipopa`.** Changing it splits datasets by person.
- **Set `MD_DATA`** to a new directory and export it in every job.
- **Copy the paracetamol package.** Its only content needed is CHEMBL112/param_e932f4c4f371.
  The original REST2-PARA build-top.config points at it; a copy is in `runs/REST2-PARA/build/parameter/`.

---

## 2. Paths

The scripts have absolute sbrinz paths baked in. Choose the new locations, then rewrite them in
one pass:

| sbrinz | meaning | new cluster |
|---|---|---|
| `/scratch/ipopa/thesis` | all runs + this repo (`thesis/bsc-thesis`) | `$THESIS` |
| `/scratch/ipopa/software` | micromamba root + MD-tools worktrees | `$SOFT` |
| `/scratch/ipopa/md-data` | `MD_DATA` | `$MD_DATA` |

```bash
cd $THESIS/bsc-thesis
grep -rlE '/scratch/ipopa/(thesis|software|md-data)' --include='*.sbatch' --include='*.sh' --include='*.py' --include='*.config' . \
 | xargs sed -i -e "s#/scratch/ipopa/thesis#$THESIS#g" -e "s#/scratch/ipopa/software#$SOFT#g" -e "s#/scratch/ipopa/md-data#$MD_DATA#g"
```

Commit this as one "port paths to <cluster>" commit, so the diff is reviewable. Logs and run
records (`*.out`, `restart.json`) keep the old paths; they are history, not inputs.

Then lay out run folders as on sbrinz: `cp -a runs/<RUN> $THESIS/<RUN>` for each run you
need (see `runs/README.md`). The scripts expect `$THESIS/<RUN>/...`, not `bsc-thesis/runs/...`.

---

## 3. Slurm: rules and what to change in every sbatch header

**Rule from the user: every compute job goes through the scheduler.** This covers builds, MD,
analysis and installs. Never run them with nohup, `&`, setsid or a bare python on a login node.
Only quick inspection commands are fine interactively.

The sbrinz-specific header lines to adapt:

| line | sbrinz meaning | change to |
|---|---|---|
| `module load slurm/slurm/24.11` | Slurm not on PATH by default | whatever the new cluster needs (often nothing) |
| `--partition=defq` | only partition, no time limit | the new GPU / CPU partitions |
| `--exclude=sbrinz,sbrinz1` | "CPU job: keep off the GPU nodes" | the new CPU partition, or drop it |
| `--gres=gpu:N` | sbrinz/sbrinz1: 8× H200 each | the new gres syntax / GPU type |
| `--time=…` | generous, defq had no limit | respect the new limit; `run.sh` resumes from checkpoints, so resubmitting is safe |

**GPU layout and MPS** (see `runs/FOL-REST2-v4/scripts/rest2.sbatch`):
- One REST2 replica = one MPI rank.
- md-tools allows more than one rank per GPU only under NVIDIA MPS. The rest2 sbatch scripts start
  a private MPS daemon inside the job (pipe dir `/tmp/mps-$SLURM_JOB_ID`) and stop it at exit.
- Check that `nvidia-cuda-mps-control` exists and that MPS is allowed on the new cluster. If it is
  not, request 1 GPU per replica and remove the MPS block.
- The sbatch scripts must `unset SLURM_NTASKS SLURM_NPROCS` before `./run.sh`. Otherwise the serial
  stages (min, eq) refuse with "launcher world size N, MPI communicator size 1".

Layouts used on sbrinz (MPS limit: 3 ranks per GPU):

| system | replicas | GPUs | ranks per GPU |
|---|---|---|---|
| folate water | 8 | 4 | 2 |
| folate complexes | 10 | 5 (or 4) | 2, or 3/3/2/2 |
| small-molecule water / complex | 6 | 3 | 2 |
| paracetamol | 4 | 2 | 2 |

**H200 speeds** for sizing time limits:
- Small-molecule complexes: 394–445 ns/day per replica at 2 ranks per GPU.
- Folate complexes: 122–139 ns/day at 4 ranks per GPU, 296 ns/day at 1–2 ranks per GPU.

---

## 4. The workflow (what each step is, in the order it runs)

Each run folder in `runs/` is one dataset. The generic chain, using `md-openmm` from
`openmm-env-protonation`:

1. **Parameterise the ligand once into a package.** `build-top --parameterize` with
   `packages/<X>.parameterize.config`. The settings are Open Babel pH 7.4 protonation, AM1-BCC
   (sqm) and OpenFF Sage 2.2.1, giving `packages/<X>/` (molecule.sdf, parameters.ffxml,
   parameter.config).

   **Reuse the committed packages. Do not re-parameterise.** AM1-BCC is not bit-reproducible
   across machines, and reuse is what makes water and complex runs comparable. A reused package
   prints `reused` in `built.out`.

2. **Build the system.** `build-top -i <input> --config build/build-top.config` with
   `solute.parameters: <package>` produces built.xml/pdb. For complexes the input is
   `prep/complex.pdb` (chain A + ligand, made by `scripts/prep.py`).
3. **Scale the hot region.** `build-top --rest2-scaler --config build/scaler.config` (needs
   `build/<RES>.sdf` for bond orders). The scaler must report "ligand amides … scaled", which
   is REST2 v4.
4. **Write the run.** `build-md -odir REST2-run1 --config REST2.config` produces `REST2-run1/run.sh`
   (min, eq_1..3, production).
5. **Run.** `sbatch scripts/rest2.sbatch <run>` (GPU, MPS); check with `scripts/check_run.py`.
6. **Analyse** in `analysis-env-main`:
   - **Clustering:** `results/run_symfirst.sbatch` and `results/run_small_symfirst.sbatch` call
     `shuyu_clustering/cluster_prune.py` and `cluster_convergence.py`.
   - **Inputs:** state 0 only, using `solute_state0_prod1.nc` + `build/built.solute.pdb`. Do not use
     the whole-system trajectory, which has 50× fewer frames.
   - **Torsions:** `shuyu_clustering/select_torsions.py` writes torsions.txt / torsions_rigid.txt /
     torsions_polarH.txt.
   - **Two settings, both reported.** `--no-single-cluster` gives the `*-nosingle` figures (the ones
     to use). `allow_single_cluster=True` is Shu-Yu's default; for folate it collapses to 1 cluster.
   - **Frame cap:** ≤ 10000 frames for the dense matrix. Folate water uses stride 5.
   - **Crystal comparisons:** `FOL-REST2-*/scripts/rmsd_vs_crystal.py`, the water-vs-bound scripts,
     and crystal ligands from `crystal/`.
7. **Run the replicates** (`replicates/`, see its README): pick_starts → build → rest2 → compare.

Run-specific twists:
- **FOL-REST2-mace.** Shu-Yu's bespoke MACE-OMOL System XML (`folic_acid_mace_omol/`) is spliced
  into built.xml by `scripts/transplant.py` before the scaler; build-top cannot read a System XML.
  It is a 3H-pterin tautomer, unlike every other folate run (N1-H).
- **Complex packages are byte-identical copies of the water packages.** For a neutral COOH an
  explicit `atom_map` is needed (`scripts/atom_map.py`, `prep/atom_map.txt`). NADP+ (−3) in 3QL3
  is charged with NAGL because sqm does not converge for it.

### How the starting conformer is chosen

- **Normal runs.** `initial_structure` runs ETKDGv3: 10 conformers, seed 20260814, MMFF94s, lowest
  energy kept. build-top then uses the package SDF coordinates.
- **Replicates.** They change only the start (`start.sdf` from `pick_starts.py`) and
  `dynamics.seed` (1000+k; rep0 = 1). The checks are RMSD < 0.05 Å and identical charges.

---

## 5. Dataset index (what was run, where it lives on sbrinz)

| run folder | system | md-tools | length / state | replicas, τmax | status |
|---|---|---|---|---|---|
| FOLn-REST2-50ns | folate water, neutral | 0.6.2 | 100 ns | 8, 0.4 | done (old protocol) |
| FOL-complex-REST2 | folate 4LRH, 3QL3, neutral | 0.6.2 | 50 ns | 10, 0.4 | done (old) |
| REST2-PARA | paracetamol water | 0.6.2 | 10 ns | 4, 0.5 | done (old) |
| REST2-conformers / -complex | 4ibi 4ibk 4j7d 5mwh 5ov8, water / bound | 0.6.2 | 60 / 40 ns | 6, 0.4 | done (5mwh, 5ov8 still current) |
| **FOL-REST2-v4** | folate −2 (COO⁻), water / 4LRH / 3QL3 | b3603e0 | 100 / 50 / 50 ns | 8 / 10 / 10 | done |
| **PARA-REST2-v4** | paracetamol | b3603e0 | 10 ns | 4, 0.5 | done |
| **REST2-conformers-v4** | 4ibi −2, 4ibk −2, 4j7d +1, water / complex | b3603e0 | 60 / 40 ns | 6, 0.4 | done |
| **FOL-REST2-mace** | folate, MACE bespoke torsions | b3603e0 | 100 / 50 / 50 ns | 8 / 10 / 10 | done |
| **REST2-replicates** | PARA rep1; FOL water rep1–3 | b3603e0 | 10 / 100 ns | as rep0 | done, compared |

Headline results are in `results/symfirst-2026-10-15/README.txt` and `replicates/README.md`.
Open items as of 2026-10-21:
- Paracetamol replicates disagree (A 46 % vs 57 %), probably because 10 ns is too short. Extending
  both to 50 ns was offered but not started.
- The MACE 3H tautomer differs from the other folate runs. This was flagged to Shu-Yu.
- 3QL3: the folate glutamate tail leaves Arg57 in both the old (neutral) and the v4 runs.

---

## 6. Notes for a Claude Code instance on the new machine

- **Start in `$THESIS/bsc-thesis`** so `CLAUDE.md` is loaded. Read this file, `runs/README.md`
  and `docs/notes/` (the previous session's project memory; it explains why several choices were
  made) before acting.
- **Transfer the user's rules.** The user's global rules on sbrinz were "always Slurm, never the
  login node" plus the module load. Recreate them in `~/.claude/CLAUDE.md` with the new cluster's
  details.
- **Before the first real job:**
  1. Install the envs (§1).
  2. Port the paths (§2).
  3. Adapt the sbatch headers (§3).
  4. Run one short GPU test: copy a run, set `number_of_exchanges` to ~300 in REST2.config and
     check acceptance plus `CUDA`.

  Only then reproduce or extend production runs.
- **Never re-parameterise** a ligand that has a committed package (§4.1).
- **Never edit MD-tools source.** Changes go through Shu-Yu.
- **Do not touch `REST2-conformers/extend_all.sh`** (user's request).
- **Ask before deleting** anything in run folders or overwriting results.
