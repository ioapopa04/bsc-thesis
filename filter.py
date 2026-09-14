import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')   

REFINED_SET_DIR = "./refined-set" 

records = []

for pdb_code in os.listdir(REFINED_SET_DIR):
    complex_dir = os.path.join(REFINED_SET_DIR, pdb_code)
    sdf_path = os.path.join(complex_dir, f"{pdb_code}_ligand.sdf")

    if not os.path.isfile(sdf_path):
        continue  # skip anything that doesn't match expected name structure

    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mol = supplier[0] if len(supplier) > 0 else None

    if mol is None:
        # SDF failed to parse -- flag and skip (mol2 parsing errors are common in PDBBind)
        records.append({"pdb_code": pdb_code, "parse_failed": True})
        continue

    records.append({
        "pdb_code": pdb_code,
        "parse_failed": False,
        "mw": Descriptors.MolWt(mol),
        "rot_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "formal_charge": Chem.GetFormalCharge(mol),
        "num_rings": rdMolDescriptors.CalcNumRings(mol),
        "smiles": Chem.MolToSmiles(mol),
    })

df = pd.DataFrame(records)
print(f"Total complexes processed: {len(df)}")
print(f"Parse failures: {df['parse_failed'].sum()}")

df = df[~df["parse_failed"]]

import matplotlib.pyplot as plt

df["rot_bonds"].value_counts().sort_index().plot(kind="bar", figsize=(10,5))
plt.xlabel("Rotatable bonds")
plt.ylabel("Number of ligands")
plt.title("Rotatable bond distribution — PDBBind refined set (valid molecules)")
plt.tight_layout()
plt.savefig("rotbonds_histogram.png")  # saves it so you can view/share it easily

easy_candidates = df[
    (df["rot_bonds"] >= 3) & (df["rot_bonds"] <= 4) &
    (df["mw"] <= 400) & #molecular-weight filter
    (df["formal_charge"] == 0) & # neutral molecules only
    (df["num_heavy_atoms"] <= 30)
]
print(f"Candidates after filtering: {len(easy_candidates)}")

easy_candidates.to_csv("easy_candidates.csv", index=False)
print(f"Saved {len(easy_candidates)} candidates to easy_candidates.csv")