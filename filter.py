import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd

REFINED_SET_DIR = "path/to/refined-set"  # <-- update this

records = []

for pdb_code in os.listdir(REFINED_SET_DIR):
    complex_dir = os.path.join(REFINED_SET_DIR, pdb_code)
    sdf_path = os.path.join(complex_dir, f"{pdb_code}_ligand.sdf")

    if not os.path.isfile(sdf_path):
        continue  # skip anything that doesn't match expected structure

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