Shu-Yu's symmetry-first clustering (MD-tools main 69c262b), 2026-10-15

1_v4-new-rest2/  REST2 v4 (Shu-Yu's protonation branch): Open Babel pH 7.4 protomers, ligand amide
                 heated. Runs: FOL-REST2-v4 (folate -2; water / 4LRH / 3QL3), PARA-REST2-v4,
                 REST2-conformers-v4 (4ibi -2, 4ibk -2, 4j7d +1; water 60 ns + complex 40 ns).
2_old-rest2/     original REST2: neutral folate, amide not heated.
                 Runs: FOLn-REST2-50ns, FOL-complex-REST2 (4LRH / 3QL3), REST2-PARA,
                 REST2-conformers (water 0-60 ns) + REST2-conformers-complex (40 ns) for all 5 small
                 molecules. 5mwh and 5ov8 have the SAME protomer under the new pH model, so their
                 runs here are also the current ones (no v4 rerun needed).
3_mace-bespoke/  folate with Shu-Yu's MACE-OMOL bespoke parameters (3H-pterin tautomer, -2, REST2 v4):
                 FOL-REST2-mace (water done; 4LRH / 3QL3 waiting for GPUs).

*_clusters_overview_symfirst-nosingle.png  allow_single_cluster=False (use these; same setting as earlier analyses)
*_clusters_overview_symfirst.png           allow_single_cluster=True, Shu-Yu's default (folate collapses to 1 cluster)
*_water_vs_bound.png                       folate torsions, water vs 4LRH vs 3QL3
*_rmsd_vs_crystal_hist.png                 folate conformation RMSD to crystal, density
