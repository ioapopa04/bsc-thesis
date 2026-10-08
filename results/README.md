# Results — Shu-Yu clustering, crystal comparison and convergence

One folder per compound, three figures each:

- `<name>_clusters_overview.png` — molecule + one row per torsion: full distribution, then each
  merged cluster group; black lines = crystal (bioactive) torsions.
- `<name>_group_torsions.png` — all groups overlaid per torsion, with the same crystal lines.
- `<name>_convergence.png` — number of merged groups vs cumulative simulated time.

Scripts: `../shuyu/cluster_prune.py`, `../shuyu/cluster_convergence.py` (analysis-env, state 0,
`--no-single-cluster`). Crystal ligands: `../crystal/`.

| folder | run | length | replicas / τ_max | stride | windows | crystal reference |
|---|---|---|---|---|---|---|
| 4ibi | REST2-conformers/4ibi (run1 + ext1) | 60 ns | 6 / 0.4 | 5 | 5 ns | 4IBI 1D8, chains A, B |
| 4ibk | REST2-conformers/4ibk (run1 + ext1) | 60 ns | 6 / 0.4 | 5 | 5 ns | 4IBK 1DE, chains A, B |
| 4j7d | REST2-conformers/4j7d (run1 + ext1) | 60 ns | 6 / 0.4 | 5 | 5 ns | 4J7D I31, chain A |
| 5mwh | REST2-conformers/5mwh (run1 + ext1) | 60 ns | 6 / 0.4 | 5 | 5 ns | 5MWH UWX, chains A, B |
| 5ov8 | REST2-conformers/5ov8 (run1 + ext1) | 60 ns | 6 / 0.4 | 5 | 5 ns | 5OV8 AXN, chains A, B |
| folic-acid | FOLn-REST2-50ns | 100 ns | 8 / 0.4 | 5 | 5 ns | 4LRH chain A, 3QL3 |
| paracetamol | REST2-PARA (interrupted at 3.5 of 10 ns) | 3.5 ns | 4 / 0.5 | 1 | 0.5 ns | none (no established target) |

Several lines per crystal file mean atoms the crystal cannot tell apart (carboxyl/sulfonyl O,
ring flips, equivalent methyls); O–H torsions have no line (no H in crystal structures).
Full outputs (logs, CSVs, representatives) stay next to each run.
