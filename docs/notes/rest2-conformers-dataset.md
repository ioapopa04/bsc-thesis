---
name: rest2-conformers-dataset
description: "REST2-conformers dataset: 5 selected ligands, the run_all.sh launcher, and its 6-replica design"
metadata:
  node_type: memory
  type: project
  originSessionId: eb4b7995-c5e8-49b5-b59d-6fba29fd1f8e
  modified: 2026-10-06T08:01:56.092Z
---

`/scratch/ipopa/thesis/REST2-conformers/` (created 2026-10-05) runs the [[fol-rest2-workflow]]
md-tools REST2 protocol in explicit water for **5 unrelated ligands** picked for having 6-7
rotatable bonds, not folic-acid analogues. Inputs are the RDKit 3D SDFs in
`/scratch/ipopa/thesis/bsc-thesis/conformers/{4ibi,4ibk,4j7d,5mwh,5ov8}_lowest_energy_conf.sdf`
(they carry Cl/F/S; OpenFF sage-2.2.1 + AM1-BCC types them all cleanly).

One launcher does build + queue + run: `/scratch/ipopa/thesis/REST2-conformers/run_all.sh`
(resnames L01..L05). Per-compound self-contained root `<base>/<name>/` with build/, input/,
min/, REST2.config and REST2-run1/. Key options: `--build-only`, `--run-only`,
`--status`, `--only <name>`, `--detach`.

**Design points that differ from the folic-acid runs:**
- **6 replicas**, tau 0.0..0.4 linear (not 8). The 8-replica ladder gave 0.42-0.45
  neighbour acceptance; 6 is steeper, so **watch the per-pair report - every pair must be
  >= 0.2**, else rerun with 8 replicas or a lower tau_max. **Validated 2026-10-06 on ALL FIVE
  finishes:** 4ibi 0.404 (pairs 0.400-0.408), 4ibk 0.421 (0.405-0.432), 4j7d 0.433 (0.414-0.443),
  5mwh 0.347 (0.326-0.366), 5ov8 0.357 (0.336-0.377) - every neighbouring pair clears 0.2, so the
  6-replica ladder stands (5mwh is the low one; a more flexible ligand may still want 8).
- **30 ns** per state = `number_of_exchanges: 15000` x 1000 steps x 2 fs (no HMR). All five
  runs (4ibi/4ibk/4j7d/5mwh/5ov8) **completed 2026-10-06** at the full 30 ns (15M steps, 15000
  exchanges). 4j7d ran in a separate later batch (`4j7d-queue.log`) after its build gotcha was
  cleared; it finished last at 00:06. Per-pair acceptance lives in
  `REST2-run1/remd_records/restart_prod1.json` -> `completion_report.by_neighbouring_pair`
  (not in run.log - run.log only has the exceptions-per-pair counter).
- **Build gotcha (hit by 4j7d, 2026-10-05):** an interrupted `--parameterize` leaves a
  hidden `.parameterize-<rand>/` scratch dir inside `build/parameter/`; on the next build
  `md-tools` rejects the whole package with `PackageError: ... unexpected files
  ['.parameterize-...']`. The parameterization itself is fine - just `rm -rf` that hidden
  dir and rerun `./run_all.sh --build-only --only <name>` (parameterise is then skipped,
  solvate re-validates the now-clean package).
- Defaults `WAIT_FOR_GPUS=1`: the script polls `nvidia-smi` until the GPU pool is idle
  (<10% util, <3 GiB) before launching, because sbrinz is a shared login node (another
  user's job usually sits on all 8 GPUs at ~88 GiB each). Set `WAIT_FOR_GPUS=0` to force.

`run_all.sh` is NOT git-tracked. It was run on the login node because Slurm was missed (it exists via `module load slurm/slurm/24.11`); per [[slurm-always]], future runs must go through sbatch.

Logs when started detached: build phase -> `detach.out`, run phase -> `run_all.log`,
per compound -> `<name>/run.log`; progress with `./run_all.sh --status`.

## Shu-Yu clustering of the five runs (2026-10-06)

Ran [[bsc-thesis-protocols]] shuyu protocol on all five (state 0 solute traj, stride 5,
`--no-single-cluster`; resname L01..L05, sdf `build/L0X.sdf`). Torsion selection is clean:
5-6 rotatable torsions each, **no rigid bonds**, polar-H only for L01/L02 (2 each); L01/L03
drop one CF3-type trivial rotor. Outputs in `<name>/<name>-clustering-shuyu-nosingle/` (+ `.log`);
symmetry pruning merged nothing (all groups singletons - the operations exist but no cluster
pair is symmetry-equivalent).

| name | resname | clusters | noise | sym ops |
|------|---------|----------|-------|---------|
| 4ibi | L01 | 14 | 14.4 % | 1 |
| 4ibk | L02 | 14 | 10.5 % | 1 |
| 4j7d | L03 | 3 (0.64/0.05/0.05) | 26.6 % | 4 |
| 5mwh | L04 | 23 (flat, all <8 %) | 19.1 % | 16 |
| 5ov8 | L05 | 4 (0.35/0.32/0.05/0.03) | 24.8 % | 16 |

L03 and L05 are the two with a dominant basin; L04 is a broad flat multi-basin ligand (23
clusters, no cluster >8 %); L01/L02 sit in between with 14 comparable clusters. Contrast the
FOL 100 ns run (2 clusters + 87.5 % noise) - these flexible acyclic-ish ligands cluster far
richer, as expected from the 5-6 rotatable bonds they were selected for.

## 4ibi 60 ns re-analysis with crystal reference (2026-10-08, Slurm job 8030)

`4ibi/analysis-0-60ns/run_analysis.sbatch` (analysis-env, CPU node, ~30 s) →
`4ibi-clustering-shuyu-nosingle-60ns-ref/` + `4ibi-convergence-60ns/`. 12 groups, 18.3 % noise;
cumulative 5 ns scan 11,14,15,13,13,14,13,13,12,12,12,12 → converged at 55 ns (45-55 plateau,
60 ns also 12). Closest to crystal 4IBI: groups A (9.6 %) and C (8.3 %) — t2/t3/t4 within
~13°, t1 ~30° off (−157° vs −127/−117°); they differ only in the carboxyl flip (crystal-ambiguous).

## All five at 60 ns (2026-10-08, Slurm array 8047, `REST2-conformers/run_analysis_60ns.sbatch`)

Cluster counts per cumulative 5 ns window (5..60): 4ibk 9 14 19 19 18 14 15 15 15 16 16 14;
4j7d 2 1 3 2 3 3 3 3 3 3 3 3 (stable 3 from 25 ns); 5mwh 14 22 24 23 23 23 23 23 22 19 17 16
(falling — not converged); 5ov8 8 10 5 4 5 4 4 4 4 4 3 3. Only 4ibi and 4j7d hold their plateau.
Closest cluster to crystal (max |Δ| over heavy-atom torsions): 4ibi A/C (~30°, t1); 4ibk only K
(1.8 %, 34°) — groups A-H have t4 ~165° off; 4j7d A (60 %, 32°); 5mwh O/K/E/C (30-37°, 1-8 %);
5ov8 none closer than 61° (B, 33 %). Crystal matching filters mirror-image graph matches by
stereocentre handedness (needed for 4j7d's two identical chlorophenyls).
