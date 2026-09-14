import os
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import pandas as pd
from rdkit import RDLogger
import matplotlib.pyplot as plt

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

df["rot_bonds"].value_counts().sort_index().plot(kind="bar", figsize=(10,5))
plt.xlabel("Rotatable bonds")
plt.ylabel("Number of ligands")
plt.title("Rotatable bond distribution — PDBBind refined set (valid molecules)")
plt.tight_layout()
plt.savefig("rotbonds_histogram.png")  # saves it so you can view/share it easily

easy_candidates = df[
    (df["rot_bonds"] >= 5) & (df["rot_bonds"] <= 6) &
    (df["mw"] >= 450) & (df["mw"] <= 550) &   # around 500
    (df["formal_charge"] == 0) &
    (df["num_heavy_atoms"] <= 40)
]

print(f"Candidates after filtering: {len(easy_candidates)}")

#easy_candidates.to_csv("easy_candidates.csv", index=False)
#print(f"Saved {len(easy_candidates)} candidates to easy_candidates.csv")


def get_protein_size(pdb_code, refined_set_dir=REFINED_SET_DIR):
    """
    Returns (n_residues, n_atoms) for {pdb_code}_protein.pdb.
    n_residues = count of unique (chain_id, resSeq, iCode) combos among CA atoms.
    n_atoms    = total ATOM record count (heavy + H if present) in the protein file.
    Returns (None, None) if the file is missing or unparseable.
    """
    protein_path = os.path.join(refined_set_dir, pdb_code, f"{pdb_code}_protein.pdb")
    if not os.path.isfile(protein_path):
        return None, None

    seen_residues = set()
    n_atoms = 0

    with open(protein_path, "r") as f:
        for line in f:
            if line.startswith("ATOM"):
                n_atoms += 1
                atom_name = line[12:16].strip()
                if atom_name == "CA":
                    chain_id = line[21]
                    res_seq = line[22:26].strip()
                    i_code = line[26]
                    seen_residues.add((chain_id, res_seq, i_code))

    return len(seen_residues), n_atoms

sizes = easy_candidates["pdb_code"].apply(lambda code: pd.Series(
    get_protein_size(code), index=["n_residues", "n_protein_atoms"]
))

easy_candidates = pd.concat([easy_candidates.reset_index(drop=True), sizes.reset_index(drop=True)], axis=1)

missing = easy_candidates[easy_candidates["n_residues"].isna()]
if len(missing) > 0:
    print(f"WARNING: {len(missing)} candidates missing/unparseable protein.pdb:")
    print(missing["pdb_code"].tolist())

valid = easy_candidates.dropna(subset=["n_residues"]).copy()
valid["n_residues"] = valid["n_residues"].astype(int)
valid["n_protein_atoms"] = valid["n_protein_atoms"].astype(int)

# Rank by protein size (smallest protein first) and take top 5
smallest_5 = valid.sort_values("n_residues", ascending=True).head(5)

print(smallest_5[["pdb_code", "mw", "rot_bonds", "n_residues", "n_protein_atoms"]])
smallest_5.to_csv("easy_candidates.csv", index=False)
print(f"Saved {len(smallest_5)} candidates to easy_candidates.csv")