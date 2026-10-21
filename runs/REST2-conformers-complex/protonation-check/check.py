"""Would MD-tools' new ligand pH step (Open Babel, Shu-Yu's protonation branch) change the
protomer of the five small molecules? Applied to each run's own parameter package molecule."""
import sys
from rdkit import Chem
from md_tools.openmm.ligand_protonation import protonate_molecule, protonation_summary

PKG = "/scratch/ipopa/thesis/REST2-conformers-complex/packages"
for name, res in zip(["4ibi", "4ibk", "4j7d", "5mwh", "5ov8"], ["L01", "L02", "L03", "L04", "L05"]):
    mol = Chem.MolFromMolFile(f"{PKG}/{name}-water/molecule.sdf", removeHs=False)
    print(f"{name} ({res}) {Chem.MolToSmiles(Chem.RemoveHs(mol))}  net {Chem.GetFormalCharge(mol):+d}")
    for ph in (7.4, 7.0):
        out, rec = protonate_molecule(mol, ph, route="sdf")
        print(f"   pH {ph}: {protonation_summary(rec)}")
        if rec["changed"]:
            print(f"           -> {Chem.MolToSmiles(Chem.RemoveHs(out))}  changes {rec.get('changes')}")
