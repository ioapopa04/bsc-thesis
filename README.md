# bsc-thesis — REST2 ligand-torsion analysis

Two independent analysis protocols live here, each in its own folder. They share the
same inputs (ligand SDF + simulation topology + a state-0 solute trajectory) but differ
in torsion selection, clustering algorithm and reporting.

| | **Ioana's protocol** (`ioana/`) | **Shu-Yu's protocol** (`shuyu/`) |
|---|---|---|
| Torsions | symmetry version of `select_torsions.py` → `torsions.txt` + `torsions_symops.txt` | compact `select_torsions.py` → `torsions.txt` / `torsions_rigid.txt` / `torsions_polarH.txt` (`--polar-h`) |
| Clustering | `cluster.py` — symmetry-aware, `--metric max`, `--drop-free`, assign-to-cores, stability + block-bootstrap error | `cluster_prune.py` — MD-tools `t_hdbscan` + `TorsionSymmetry` pruning + summary |
| PCA | `dpca.py` (dihedral PCA) | — |
| Torsion diagnostics | `torsion_analysis.py` | — |
| Report format | Shu-Yu-compatible files via `group_report.py` **plus** native outputs | `groups.csv`, `merged_groups.json`, `torsion_cv.json`, `representative_<X>.pdb`, `group_*.png`, `group_torsions.png`, `representatives.png`, `torsions_<i>_<j>_2d.png` |

**The two `select_torsions.py` are not interchangeable** — Shu-Yu's needs `--polar-h` and
writes no symmetry file; Ioana's writes `torsions_symops.txt` that `cluster.py`/`dpca.py`
require. Each folder keeps its own copy.

Shared, protocol-independent early-stage files (ligand candidate filtering and conformer
generation) sit at the top level: `filter.py`, `conformer_generation_from_SMILES.py`,
`conformers/`, `easy_candidates.csv`, `rotbonds_histogram.png`.

**Moving to a new machine: see `SETUP.md`** (software, paths, Slurm, workflow). Run inputs are in
`runs/`, environment specs in `envs/`, project notes in `docs/notes/`, results in `results/`.

See the per-folder READMEs for the exact commands and required environment
(`source folenv.sh` then `micromamba activate analysis-env`).
