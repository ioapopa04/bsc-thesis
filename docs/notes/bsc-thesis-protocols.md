---
name: bsc-thesis-protocols
description: "bsc-thesis holds two REST2 ligand-torsion protocols (shuyu/ vs ioana/), with two incompatible select_torsions.py and a split micromamba env requirement"
metadata:
  node_type: memory
  type: project
  originSessionId: 2783d413-c5f5-4b61-9077-829c6b83315b
  modified: 2026-10-06T08:01:33.321Z
---

`/scratch/ipopa/thesis/bsc-thesis` holds **two independent analysis protocols**, split on 2026-10-05 into
`shuyu/` and `ioana/` (shared early-stage ligand-selection files stay at top level). Both analyse
the same [[fol-rest2-workflow]] ligand runs.

- `shuyu/` — `select_torsions.py` (compact, `--polar-h`; writes `torsions.txt`/`_rigid`/`_polarH`,
  no symmetry) → `cluster_prune.py` (MD-tools `t_hdbscan` + `TorsionSymmetry`).

## cluster_prune.py automorphism cap (added 2026-10-06)

MD-tools' `TorsionSymmetry` enumerates graph automorphisms and caps the search at
`MAX_AUTOMORPHISMS = 200000` (`md_tools/analysis/_t_symmetry.py`); a truncated search raises
`RuntimeError` rather than silently reporting "no further symmetry". Molecules with several
equivalent methyl/CF3 groups blow past it (the group is factorial in equivalent H). **Added a
backwards-compatible `--max-matches` flag** to `shuyu/cluster_prune.py` (default 200000, so the
FOL runs are unchanged) passed straight to `TorsionSymmetry(mol, torsion_idx, max_matches=...)`.
Counts seen on the REST2-conformers ligands: L01 12, L02 2 (below cap, ran unflagged); L03
2.24M (17 s), L04/L05 663k (4 s) — pass `--max-matches 8000000` for those. Raising the cap above
the true count gives the complete, correct operation set. `shuyu/` is **untracked in git**
(whole-folder `?? shuyu/`), so this edit is not in the tracked diff.
- `ioana/` — `select_torsions.py` (**symmetry version**; writes `torsions.txt` + `torsions_symops.txt`)
  → `torsion_analysis.py` → `cluster.py` (symmetry-aware `--metric max`, `--drop-free`,
  assign-to-cores) → `dpca.py`. `group_report.py` re-emits Shu-Yu's exact file set from Ioana's labels;
  `compare_clusters.py` cross-tabulates the two.

**The two `select_torsions.py` are NOT interchangeable** — Shu-Yu's needs `--polar-h`, Ioana's
writes the `torsions_symops.txt` that `cluster.py`/`dpca.py` require. Each folder keeps its own copy.

## Split env requirement (non-obvious)

No single env runs the whole `ioana/` chain:
- `openmm-env` — has `networkx`, mdtraj, rdkit, sklearn, scipy, pandas → `select_torsions.py`,
  `cluster.py`, `dpca.py`, `torsion_analysis.py`.
- `analysis-env` — has `md_tools.analysis` (needed by `group_report.py`, `compare_clusters.py`,
  `cluster_prune.py`) but **lacks `networkx`**.

So `select_torsions.py` cannot run in analysis-env, and the reporting steps cannot run in openmm-env.
See [[md-tools-openmm-env]].

## 100 ns result (FOLn-REST2-50ns, state 0, stride 5)

Shu-Yu: 2 clusters A 6.7 % / B 5.8 %, **87.5 % noise** (≈13 % clustered). Ioana: **4 clusters
25.8 / 21.6 / 12.1 / 10.4 %, 30 % noise**. Cross-tab: Shu-Yu's A/B are not separable states —
A splits ~52/38/9 over Ioana clusters 0/1/2, and ~71-99 % of every Ioana cluster is Shu-Yu noise.

## shuyu/ additions (2026-10-08)

`cluster_prune.py` gained `--name`, `--rest2-config`, `--duration-ns`, `--ref-sdf` (crystal
torsion vlines) and writes `clusters_overview.png` (all groups in one figure, own 18-colour
torsion palette since MD-tools' TORSION_COLOURS has only 5). New `cluster_convergence.py`
(cumulative 5 ns windows, plateau of 3 = converged) imports helpers from `cluster_prune.py`.
Crystal ligands in `bsc-thesis/crystal/` (4IBI ligand 1D8, chains A/B, from RCSB models API).
Crystal carboxyl O/OH are indistinguishable → that torsion gets both 0° and 180° lines.

## results/ folder (2026-10-08)

`bsc-thesis/results/<name>/` holds clusters_overview / group_torsions / convergence PNGs for
4ibi 4ibk 4j7d 5mwh 5ov8 (60 ns), folic-acid (FOLn-REST2-50ns = 100 ns; refs 4LRH chain A +
3QL3) and paracetamol (REST2-PARA — interrupted at 3.5 of 10 ns, stride 1, 0.5 ns windows, no
crystal ref). README.md there lists sources. FOL/PARA re-run via `results/run_fol_para.sbatch`
(Slurm array 8062). FOL 100 ns: 2 groups, 87.5 % noise, count 2 from 65 ns on; PARA: 2 groups.
