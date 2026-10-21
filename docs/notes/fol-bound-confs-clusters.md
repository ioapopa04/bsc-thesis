---
name: fol-bound-confs-clusters
description: "Where the two bound folic-acid conformations (4LRH, 3QL3) land relative to the REST2 clusters, and how the crystal->topology mapping is done"
metadata:
  node_type: memory
  type: project
  originSessionId: 6f1ebc7e-0e65-4dbe-b805-9c145a9a5050
  modified: 2026-09-30T15:08:27.214Z
---

Reference bound conformations of folic acid from Zivanovic et al. JCIM 2015
(`pubs.acs.org/jcisd8/article/55/10/2154`): **4LRH** = human folate receptor alpha
(ligand `FOL`, 8 copies chains A–H), **3QL3** = DHFR (ligand `FOL`, altloc A/B).
The paper's point is that these two bound conformations differ.

Assignment script: `/scratch/ipopa/thesis/FOLn-REST2/REST2-run1/bound_check/assign_bound.py`
(+ `control_rmsd.py`). It maps crystal ligand atoms onto the built topology by
enumerating **all** graph isomorphisms (RDKit `AssignBondOrdersFromTemplate` +
`GetSubstructMatches`), computes exactly the torsions clustering kept
(`1 2 3 4 6 7`), applies the same chord/`max` distance with the min over symmetry
ops, and reads the nearest MD frame's label. Mapping validated three ways: 2
automorphisms = the symops count, elements match, mapped bond lengths 1.21–1.57 Å,
and RDKit `GetBestRMS` agrees with a Kabsch RMSD.

**Result (10 ns run, `clusters_new`): neither bound conformation belongs to any
cluster — both are nearest to NOISE frames.** Nearest *clusters* do differ:
4LRH → cluster 0, 3QL3 → cluster 5, but at d ≈ 1.03 / 0.93 (chord/max units),
far outside the ensemble's own scale (ensemble NN distance median 0.205, **max
0.878**). Heavy-atom RMSD to the nearest frame is 5.25 / 4.47 Å, versus an
ensemble nearest-neighbour RMSD of 0.13 Å (random pair 0.61 Å) — i.e. 20–40×
typical. Holds under Daura too (which is 80% noise). The two conformers are
maximally different from each other (t2 and t4 differ by ~180°/~120°).

So with 10 ns one cannot say the two structures occupy different clusters — only
that both are unsampled. The 50 ns run (`FOLn-REST2-50ns`) is the test.
See [[fol-rest2-workflow]].
