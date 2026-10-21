---
name: fol-rest2-mace
description: Folate REST2 with Shu-Yu's MACE-OMOL bespoke parameters (System XML transplanted into MD-tools builds); 3H tautomer differs from all earlier folate runs
metadata:
  type: project
---
2026-10-15. Shu-Yu: thesis/folic_acid_mace_omol/ligand_sage_mace_system.xml "is the ligand parameter, already optimized with bespokefit" (Sage 2.2.1 + MACE-OMOL torsions, AM1-BCC, -2).
- Its molecule is the 3H-pterin tautomer (H on N3); ALL earlier folate runs (FOLn, FOL-complex, FOL-REST2-v4) used the N1-H tautomer from the original SDF. Told user to flag to Shu-Yu/advisor.
- build-top can't read a System XML (only .offxml via refinement_offxml), so thesis/FOL-REST2-mace/scripts/transplant.py swaps every folate term into built.xml before the scaler (keeps built.sage.xml; checks in transplant.json). Charges identical to MD-tools AM1-BCC; 54 torsion quadruplets changed; energies match exactly.
- Same settings as v4 (env openmm-env-protonation). Sims 11621-11623, analyses 11625-11627 (symfirst, both allow_single_cluster settings, RMSD, water-vs-bound). fragments.py names translated to MACE package names.
Related: [[fol-rest2-v4]], [[symfirst-clustering]]
- 2026-10-17: 11622/11623 failed at start (exit 53, no log; sbrinz1 /scratch outage). Resubmitted: 4LRH 12614, 3QL3 12615, ana 12616 (4LRH), 12617 (3QL3+rmsd).
- 2026-10-17 23:20: 12614/12615 changed to gres/gpu:4 (10 ranks, 3/3/2/2 per GPU under MPS) to fit free GPUs; 12614 started on sbrinz1.
