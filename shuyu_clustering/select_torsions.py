#!/usr/bin/env python
"""
Select and classify the torsions of one ligand. No symmetry handling here:
cluster_prune.py (md_tools.analysis TorsionSymmetry) detects symmetry itself.

Bond orders come from the ligand SDF; atom indices are written in the numbering
of the simulation topology (checked: same atom count and element order).

Candidate = bond between two heavy atoms that each have another heavy neighbour.
Each candidate is written to one file:
  <out>.txt         rotatable: single, non-ring, not rigid -> use for clustering
  <out>_rigid.txt   rigid: non-ring double bond, amide-like C(=O/S)-N, ester C(=O)-O,
                    bond next to a triple bond (barely rotates; include if wanted,
                    e.g. the amide omega in the MD-tools tutorial)
  <out>_ring.txt    bonds in non-aromatic rings (not handled further)
Aromatic ring bonds are ignored. Bonds to a 3-fold rotor whose three
substituents are identical atoms with no further heavy neighbours
(CF3, CCl3, C(CH3)3) are dropped: rotating them changes nothing.

With --polar-h also:
  <out>_polarH.txt  rotation of a polar hydrogen: single non-ring bond X-Y, Y = O/N/S
                    with exactly one H and no other heavy neighbour (phenol OH,
                    carboxylic acid OH, SH); torsion = heavy-X-Y-H.

Line format: i j k l   # atom names | class   (cluster_prune.py reads the four indices)

Example (solute-only trajectory -> ligand indices):
  python select_torsions.py --sdf build/FOL.sdf --top build/built.pdb \
      --select "resname FOL" --ligand-indices --polar-h --out REST2-run1/torsions
"""
import argparse

import mdtraj as md
from rdkit import Chem

RIGID_PATTERNS = {
    "amide-like C-N": "[CX3](=[OX1,SX1])-[#7]",     # amide, carbamate, urea, thioamide
    "ester C(=O)-O": "[CX3](=[OX1])-[OX2]",
}


def heavy_nbrs(mol, a, exclude):
    return sorted(n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()
                  if n.GetAtomicNum() > 1 and n.GetIdx() != exclude)


def trivial_rotor_end(mol, center, other):
    """center carries three identical terminal groups (CF3, CCl3, C(CH3)3)."""
    subs = [n for n in mol.GetAtomWithIdx(center).GetNeighbors() if n.GetIdx() != other]
    if len(subs) != 3:
        return False
    if len({n.GetAtomicNum() for n in subs}) != 1 or subs[0].GetAtomicNum() == 1:
        return False
    # identical and terminal (only H besides the center), e.g. F, Cl, CH3
    return all(not heavy_nbrs(mol, n.GetIdx(), center) and n.GetTotalNumHs() == subs[0].GetTotalNumHs()
               for n in subs)


def next_to_triple(mol, a):
    return any(b.GetBondType() == Chem.BondType.TRIPLE for b in mol.GetAtomWithIdx(a).GetBonds())


def classify(mol, ester_rotatable=False):
    patterns = {k: v for k, v in RIGID_PATTERNS.items() if not (ester_rotatable and k.startswith("ester"))}
    rigid_bonds = {}
    for name, sma in patterns.items():
        for m in mol.GetSubstructMatches(Chem.MolFromSmarts(sma)):
            rigid_bonds.setdefault(frozenset((m[0], m[2])), name)

    rows = []
    for bond in mol.GetBonds():
        j, k = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if mol.GetAtomWithIdx(j).GetAtomicNum() == 1 or mol.GetAtomWithIdx(k).GetAtomicNum() == 1:
            continue
        if bond.GetIsAromatic():
            continue
        nj, nk = heavy_nbrs(mol, j, k), heavy_nbrs(mol, k, j)
        if not nj or not nk:
            continue
        ijkl = (nj[0], j, k, nk[0])
        bt = bond.GetBondType()
        if bond.IsInRing():
            cls, note = "ring", ""
        elif bt == Chem.BondType.DOUBLE:
            cls, note = "rigid", "double bond"
        elif bt == Chem.BondType.TRIPLE or next_to_triple(mol, j) or next_to_triple(mol, k):
            cls, note = "rigid", "next to triple bond"
        elif frozenset((j, k)) in rigid_bonds:
            cls, note = "rigid", rigid_bonds[frozenset((j, k))]
        elif trivial_rotor_end(mol, j, k) or trivial_rotor_end(mol, k, j):
            cls, note = "trivial", "3 identical substituents"
        else:
            cls, note = "rotatable", ""
        rows.append(dict(ijkl=ijkl, cls=cls, note=note))
    return rows


def polar_h_rows(mol):
    rows = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        for x, y in ((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
                     (bond.GetEndAtomIdx(), bond.GetBeginAtomIdx())):
            ax, ay = mol.GetAtomWithIdx(x), mol.GetAtomWithIdx(y)
            if ax.GetAtomicNum() == 1 or ay.GetSymbol() not in ("O", "N", "S"):
                continue
            hs = [n.GetIdx() for n in ay.GetNeighbors() if n.GetAtomicNum() == 1]
            other_heavy = [n.GetIdx() for n in ay.GetNeighbors() if n.GetAtomicNum() > 1 and n.GetIdx() != x]
            nx_ = heavy_nbrs(mol, x, y)
            if len(hs) == 1 and not other_heavy and nx_:
                rows.append(dict(ijkl=(nx_[0], x, y, hs[0]), cls="polar-H", note=f"{ay.GetSymbol()}-H"))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sdf", required=True, help="ligand SDF with bond orders and explicit H")
    ap.add_argument("--top", required=True, help="simulation topology, e.g. build/built.pdb")
    ap.add_argument("--select", required=True, help='ligand atoms in --top, e.g. "resname FOL"')
    ap.add_argument("--ligand-indices", action="store_true",
                    help="write ligand indices 0..n-1 (for solute-only trajectories)")
    ap.add_argument("--polar-h", action="store_true", help="also write polar-hydrogen torsions")
    ap.add_argument("--ester-rotatable", action="store_true", help="treat ester C(=O)-O as rotatable")
    ap.add_argument("--out", default="torsions", help="output prefix")
    args = ap.parse_args()

    mol = Chem.MolFromMolFile(args.sdf, removeHs=False)
    if mol is None:
        raise SystemExit(f"RDKit could not read {args.sdf}")

    top = md.load_topology(args.top)
    sel = top.select(args.select)
    sdf_el = [a.GetSymbol().upper() for a in mol.GetAtoms()]
    top_el = [top.atom(int(x)).element.symbol.upper() for x in sel]
    if len(sdf_el) != len(top_el):
        raise SystemExit(f"Atom count differs: SDF {len(sdf_el)} vs topology selection {len(top_el)}")
    if sdf_el != top_el:
        bad = [i for i, (a, b) in enumerate(zip(sdf_el, top_el)) if a != b]
        raise SystemExit(f"Element order differs at ligand positions {bad[:10]}: SDF and topology do not match.")
    print(f"Atom order check passed: {len(sdf_el)} atoms.")

    names = [str(top.atom(int(x))) for x in sel]
    to_out = (lambda x: x) if args.ligand_indices else (lambda x: int(sel[x]))

    rows = classify(mol, args.ester_rotatable)
    files = {"rotatable": f"{args.out}.txt", "rigid": f"{args.out}_rigid.txt", "ring": f"{args.out}_ring.txt"}
    if args.polar_h:
        rows += polar_h_rows(mol)
        files["polar-H"] = f"{args.out}_polarH.txt"

    handles = {c: open(p, "w") for c, p in files.items()}
    for h in handles.values():
        h.write("# i j k l   # atom names | class\n")
    print(f"\n{'class':10s} torsion (atom names)                          note")
    for r in rows:
        label = " ".join(names[x] for x in r["ijkl"])
        print(f"{r['cls']:10s} {label:45s} {r['note']}")
        if r["cls"] in handles:
            i, j, k, l = (to_out(x) for x in r["ijkl"])
            handles[r["cls"]].write(f"{i} {j} {k} {l}   # {label} | {r['cls']} {r['note']}\n")
    for h in handles.values():
        h.close()

    print()
    for c in ("rotatable", "rigid", "ring", "polar-H", "trivial"):
        n = sum(r["cls"] == c for r in rows)
        if c in files:
            print(f"{n:3d} {c:10s} -> {files[c]}")
        elif n:
            print(f"{n:3d} {c:10s} (dropped)")


if __name__ == "__main__":
    main()
