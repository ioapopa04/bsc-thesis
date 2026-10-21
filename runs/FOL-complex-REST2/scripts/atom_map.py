"""Heavy-atom map {deposited name: package name} for a ligand whose package has carboxylic ACIDS.

  python -I scripts/atom_map.py <package-dir> <complex.pdb> <RESNAME>

md-tools refuses an ambiguous graph match whose alternatives are not parameter-equivalent: the two
oxygens of a C(=O)OH are, so an explicit `atom_map` is required for every heavy atom. Among all
graph matches, this picks the one that puts each acid's hydroxyl oxygen on the deposited oxygen
with the LONGER C-O bond (single bond ~1.31 A vs C=O ~1.21 A); ring flips are true symmetries and
any of them is fine. Prints a YAML flow mapping ready to paste into build-top.config, plus the
bond lengths it decided on.
"""
import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

pkg, pdb, resname = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
names = [a["name"] for a in json.loads((pkg / "metadata.json").read_text())["atoms"]]
pmol = Chem.MolFromMolFile(str(pkg / "molecule.sdf"), removeHs=False)
heavy = [a.GetIdx() for a in pmol.GetAtoms() if a.GetAtomicNum() > 1]
template = Chem.RemoveHs(pmol)            # heavy-atom indices of template == package order (H last)
assert [names[i] for i in heavy] == [names[i] for i in range(len(heavy))], "H atoms not last"

lines = [l for l in pdb.read_text().splitlines()
         if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == resname]
crystal = Chem.MolFromPDBBlock("\n".join(lines) + "\nEND\n", removeHs=True, proximityBonding=True)
crystal = AllChem.AssignBondOrdersFromTemplate(template, crystal)
dep_names = [l[12:16].strip() for l in lines]
xyz = crystal.GetConformer().GetPositions()

# acid patterns in the package: (C, O-carbonyl, O-hydroxyl) package indices
acids = template.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OX2H1,OX2]"))
acids = [(c, o1, o2) for c, o1, o2 in acids
         if pmol.GetAtomWithIdx(o2).GetTotalNumHs(includeNeighbors=True) == 1]
print(f"# {len(acids)} carboxylic acid(s) in the package: "
      + ", ".join(f"{names[c]}(={names[o1]})-{names[o2]}H" for c, o1, o2 in acids))

best = None
for m in crystal.GetSubstructMatches(template, uniquify=False, useChirality=False, maxMatches=10000):
    # m[i] = crystal index of package heavy atom i
    ok = all(np.linalg.norm(xyz[m[c]] - xyz[m[oh]]) > np.linalg.norm(xyz[m[c]] - xyz[m[od]])
             for c, od, oh in acids)
    if ok:
        best = m
        break
if best is None:
    sys.exit("no graph match puts every hydroxyl on the longer C-O bond")
for c, od, oh in acids:
    d1 = np.linalg.norm(xyz[best[c]] - xyz[best[od]])
    d2 = np.linalg.norm(xyz[best[c]] - xyz[best[oh]])
    print(f"# {names[c]}: C=O {dep_names[best[od]]}->{names[od]} {d1:.3f} A,  "
          f"C-OH {dep_names[best[oh]]}->{names[oh]} {d2:.3f} A")
amap = {dep_names[best[i]]: names[i] for i in range(template.GetNumAtoms())}
print("atom_map: {" + ", ".join(f"{k}: {v}" for k, v in amap.items()) + "}")
