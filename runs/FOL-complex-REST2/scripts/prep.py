"""Prepare the inputs of the folate complex REST2 runs from the deposited crystal structures.

  packages/FOL.sdf   folate DIANION (-2), 3H-pteridin-4-one tautomer, bond orders from the
                     template SMILES, heavy-atom pose and stereo from 4LRH chain A
  packages/NAP.sdf   NADP+ (-3: nicotinamide +1, pyrophosphate -2, 2'-phosphate -2), pose and
                     stereo from 3QL3
  4LRH/prep/complex.pdb    FRa chain A + FOL; NAG, waters and the other 7 chains dropped
  3QL3/prep/complex.pdb    DHFR chain A (altloc A) + FOL (altloc A) + NAP; Mn2+ and waters dropped

Ligands are written as heavy atoms only: md-tools places the package hydrogens itself.
"""
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = Path(sys.argv[1]).resolve()
CRYSTAL = ROOT / "crystal"

FOL_SMILES = "Nc1nc2ncc(CNc3ccc(cc3)C(=O)NC(CCC(=O)[O-])C(=O)[O-])nc2c(=O)[nH]1"
NAP_SMILES = ("NC(=O)c1ccc[n+](c1)C1OC(COP(=O)([O-])OP(=O)([O-])OCC2OC(C(OP(=O)([O-])[O-])C2O)"
              "n2cnc3c(N)ncnc23)C(O)C1O")

DROP = {"HOH", "NAG", "MN"}


def records(pdb, chain, keep_het):
    """ATOM/HETATM lines of one chain, altloc ' ' or 'A' (altloc column blanked)."""
    out = []
    for line in (CRYSTAL / pdb).read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")) or line[21] != chain:
            continue
        if line[16] not in " A":
            continue
        resname = line[17:20].strip()
        if line.startswith("HETATM") and (resname in DROP or resname not in keep_het):
            continue
        out.append(line[:16] + " " + line[17:])
    return out


def ligand_mol(lines, smiles, name):
    """Crystal heavy atoms + template bond orders/charges -> explicit-H mol, stereo from 3D."""
    block = "\n".join(lines) + "\nEND\n"
    crystal = Chem.MolFromPDBBlock(block, removeHs=True, sanitize=True, proximityBonding=True)
    template = Chem.MolFromSmiles(smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(template, crystal)
    # AssignBondOrdersFromTemplate copies formal charges; check the net charge explicitly.
    q = Chem.GetFormalCharge(mol)
    qt = Chem.GetFormalCharge(template)
    assert q == qt, f"{name}: charge {q} != template {qt}"
    Chem.AssignStereochemistryFrom3D(mol)
    # A phosphate P carrying resonance-equivalent O-/=O is not a stereocentre.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "P":
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    molH = Chem.AddHs(mol, addCoords=True)
    molH.SetProp("_Name", name)
    centers = Chem.FindMolChiralCenters(molH, includeUnassigned=True, useLegacyImplementation=False)
    print(f"{name}: {molH.GetNumAtoms()} atoms, {mol.GetNumAtoms()} heavy, charge {q}, "
          f"stereocentres {centers}")
    print(f"  {Chem.MolToSmiles(mol)}")
    return molH


def main():
    (ROOT / "packages").mkdir(exist_ok=True)

    # --- 4LRH chain A -------------------------------------------------------
    lines = records("4LRH.pdb", "A", keep_het={"FOL"})
    fol = [l for l in lines if l[17:20] == "FOL"]
    assert len(fol) == 32, len(fol)
    Chem.MolToMolFile(ligand_mol(fol, FOL_SMILES, "FOL"), str(ROOT / "packages/FOL.sdf"))
    write_complex(ROOT / "4LRH/prep/complex.pdb", lines, "4LRH chain A, FOL; NAG/HOH dropped")

    # --- 3QL3 chain A -------------------------------------------------------
    lines = records("3QL3.pdb", "A", keep_het={"FOL", "NAP"})
    fol = [l for l in lines if l[17:20] == "FOL"]
    nap = [l for l in lines if l[17:20] == "NAP"]
    assert len(fol) == 32 and len(nap) == 48, (len(fol), len(nap))
    ligand_mol(fol, FOL_SMILES, "FOL(3QL3 check)")   # same chemistry must be perceivable here too
    Chem.MolToMolFile(ligand_mol(nap, NAP_SMILES, "NAP"), str(ROOT / "packages/NAP.sdf"))
    write_complex(ROOT / "3QL3/prep/complex.pdb", lines,
                  "3QL3 chain A altloc A, FOL (altloc A) + NAP; MN/HOH dropped")


def write_complex(path, lines, remark):
    path.parent.mkdir(parents=True, exist_ok=True)
    prot = [l for l in lines if l.startswith("ATOM")]
    het = [l for l in lines if l.startswith("HETATM")]
    body = [f"REMARK   1 {remark}"] + prot + ["TER"] + het + ["END"]
    path.write_text("\n".join(body) + "\n")
    print(f"wrote {path}: {len(prot)} protein atoms, {len(het)} ligand atoms")


if __name__ == "__main__":
    main()
