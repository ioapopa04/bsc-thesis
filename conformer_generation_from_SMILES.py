import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

CANDIDATES_CSV = "easy_candidates.csv"
NUM_LIGANDS = 5
N_EMBED_ATTEMPTS = 50  # candidate geometries generated per ligand; only the lowest-energy one is kept
OUTPUT_DIR = "conformers"

os.makedirs(OUTPUT_DIR, exist_ok=True)

df = pd.read_csv(CANDIDATES_CSV)
selected = df.sort_values("n_residues", ascending=True).head(NUM_LIGANDS).copy()
print("Selected ligands:")
print(selected[["pdb_code", "mw", "rot_bonds", "n_residues", "n_protein_atoms"]])

def get_lowest_energy_conformer(smiles, pdb_code, n_confs=N_EMBED_ATTEMPTS, rms_thresh=0.5):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"[{pdb_code}] SMILES parse failed")
        return None

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.pruneRmsThresh = rms_thresh
    params.useRandomCoords = True

    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params))
    if len(cids) == 0:
        print(f"[{pdb_code}] Embedding failed for all attempts")
        return None

    try:
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
    except Exception:
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=2000)

    best_idx = min(range(len(cids)), key=lambda i: results[i][1])
    best_cid = cids[best_idx]
    best_energy = results[best_idx][1]

    print(f"[{pdb_code}] kept lowest-energy conformer: conf_id={best_cid}, "
          f"E={best_energy:.2f} kcal/mol (out of {len(cids)} embedded)")

    return mol, best_cid, best_energy

for _, row in selected.iterrows():
    out = get_lowest_energy_conformer(row["smiles"], row["pdb_code"])
    if out is None:
        continue
    mol, best_cid, best_energy = out

    out_path = os.path.join(OUTPUT_DIR, f"{row['pdb_code']}_lowest_energy_conf.sdf")
    writer = Chem.SDWriter(out_path)
    writer.write(mol, confId=best_cid)
    writer.close()
    print(f"[{row['pdb_code']}] wrote -> {out_path}")