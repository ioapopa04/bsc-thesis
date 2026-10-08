#!/usr/bin/env python
"""
Automatic torsion selection + classification + symmetry numbers for one ligand.

Uses the ligand SDF (has bond orders, e.g. build/<RESNAME>.sdf or build/built.sdf
from MD-tools) and the simulation topology (e.g. build/built.pdb) so that the
written atom indices match the trajectory.

Every candidate bond (single-path bond between two heavy atoms that both have
another heavy neighbour, plus non-aromatic ring bonds) is classified as:

  rotatable  -> torsions.txt          (use for analysis / PCA / clustering)
  rigid      -> torsions_rigid.txt    (non-ring double bond, amide-like C-N,
                                       ester C(=O)-O, bond next to a triple bond;
                                       compute only as a sanity check)
  trivial    -> dropped               (one end is a 3-fold symmetric rotor,
                                       e.g. CF3, C(CH3)3: rotation changes nothing)
  ring       -> torsions_ring.txt     (bonds in non-aromatic rings; small/medium
                                       rings and macrocycles need special treatment,
                                       NOT included in torsions.txt)
Aromatic ring bonds are ignored.

SYMMETRY: the molecule's symmetry operations (graph automorphisms) are written to
<out>_symops.txt, each as the relabelled torsion list. They act on ALL torsions
jointly (e.g. a ring flip shifts both torsions of R1-ring-R2 at once), which
cluster.py / dpca.py use. The per-bond number below is kept only as information.

Per-bond number k_local of a torsion j-k: for each end, if ALL other neighbours of
that atom (hydrogens included) are topologically equivalent, the end has
symmetry = their number (phenyl ipso -> 2, CF3 -> 3), else 1. k = lcm of the
two ends. Equivalence = RDKit canonical ranks with ties kept, computed on the
bare graph (bond orders, aromaticity and charges ignored) so that resonance
partners (carboxylate O's, nitro O's, amidinium N's) count as equivalent.
Only topological symmetry is detected (same limitation as TABS).

Output file format (torsions.txt etc.):
  i j k l symmetry   # atom names | class
Indices are in the numbering of --top (whole system) unless --ligand-indices.
torsion_analysis.py reads columns 1-4. cluster.py and dpca.py read torsions.txt
together with torsions_symops.txt.

Example:
  python select_torsions.py --sdf ../build/PAR.sdf --top ../build/built.pdb \
      --select "resname PAR" --out torsions
"""
import argparse
import os
from math import gcd

import numpy as np
import mdtraj as md
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from rdkit import Chem

RIGID_PATTERNS = {
    "amide-like C-N": "[CX3](=[OX1,SX1])-[#7]",     # amide, carbamate, urea, thioamide
    "ester C(=O)-O": "[CX3](=[OX1])-[OX2]",
}


def lcm(a, b):
    return a * b // gcd(a, b)


def graph_ranks(mol):
    """Canonical ranks with ties kept, on the bare graph (all bonds single,
    no aromaticity, no charges): symmetry-equivalent atoms share a rank."""
    m = Chem.RWMol(mol)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    for a in m.GetAtoms():
        a.SetIsAromatic(False)
        a.SetFormalCharge(0)
        a.SetNoImplicit(True)
    m.UpdatePropertyCache(strict=False)
    return list(Chem.CanonicalRankAtoms(m, breakTies=False))


def end_symmetry(mol, ranks, center, other):
    nbrs = [n.GetIdx() for n in mol.GetAtomWithIdx(center).GetNeighbors() if n.GetIdx() != other]
    if len(nbrs) >= 2 and len({ranks[n] for n in nbrs}) == 1:
        return len(nbrs)
    return 1


def heavy_nbrs(mol, a, exclude):
    return sorted(n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors()
                  if n.GetAtomicNum() > 1 and n.GetIdx() != exclude)


def next_to_triple(mol, a):
    return any(b.GetBondType() == Chem.BondType.TRIPLE for b in mol.GetAtomWithIdx(a).GetBonds())


def smallest_ring(mol, bond_idx):
    sizes = [len(r) for r in mol.GetRingInfo().BondRings() if bond_idx in r]
    return min(sizes) if sizes else 0


def symmetry_ops(mol, torsions, max_iter=50000):
    """Graph automorphisms of the heavy-atom graph (atoms matched by element and
    number of attached H; bond orders ignored so resonance partners are
    equivalent), restricted to those preserving R/S labels. Returned as the
    relabelled torsion quadruples, duplicates removed, identity first."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    G = nx.Graph()
    for a in heavy:
        at = mol.GetAtomWithIdx(a)
        G.add_node(a, key=(at.GetAtomicNum(), at.GetTotalNumHs(includeNeighbors=True)))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in G and j in G:
            G.add_edge(i, j)

    if mol.GetNumConformers() and mol.GetConformer().Is3D():
        Chem.AssignStereochemistryFrom3D(mol)
    cip = {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}

    identity = tuple(tuple(t) for t in torsions)
    ops, seen, n = [identity], {identity}, 0
    gm = GraphMatcher(G, G, node_match=lambda x, y: x["key"] == y["key"])
    for m in gm.isomorphisms_iter():
        n += 1
        if n > max_iter:
            print(f"WARNING: stopped after {max_iter} automorphisms - symmetry list may be incomplete.")
            break
        if any(cip.get(m[a]) != c for a, c in cip.items()):
            continue                                   # would map R onto S: not a real symmetry
        img = tuple(tuple(m[x] for x in t) for t in torsions)
        if img not in seen:
            seen.add(img)
            ops.append(img)
    return ops


def classify(mol):
    ranks = graph_ranks(mol)
    rigid_bonds = {}
    for name, sma in RIGID_PATTERNS.items():
        patt = Chem.MolFromSmarts(sma)
        for match in mol.GetSubstructMatches(patt):
            # the single bond of the pattern is between match[0] (C) and match[2]
            rigid_bonds.setdefault(frozenset((match[0], match[2])), name)

    rows = []
    for bond in mol.GetBonds():
        j, k = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        aj, ak = mol.GetAtomWithIdx(j), mol.GetAtomWithIdx(k)
        if aj.GetAtomicNum() == 1 or ak.GetAtomicNum() == 1:
            continue
        if bond.GetIsAromatic():
            continue
        nj, nk = heavy_nbrs(mol, j, k), heavy_nbrs(mol, k, j)
        if not nj or not nk:
            continue                                         # terminal: no heavy-atom torsion
        i, l = nj[0], nk[0]
        sym = lcm(end_symmetry(mol, ranks, j, k), end_symmetry(mol, ranks, k, j))

        bt = bond.GetBondType()
        if bond.IsInRing():
            cls, note = "ring", f"ring size {smallest_ring(mol, bond.GetIdx())}"
        elif bt == Chem.BondType.DOUBLE:
            cls, note = "rigid", "double bond"
        elif bt == Chem.BondType.TRIPLE or next_to_triple(mol, j) or next_to_triple(mol, k):
            cls, note = "rigid", "linear (triple bond)"
        elif frozenset((j, k)) in rigid_bonds:
            cls, note = "rigid", rigid_bonds[frozenset((j, k))]
        elif sym % 3 == 0:
            cls, note = "trivial", f"{sym}-fold symmetric rotor"
        else:
            cls, note = "rotatable", ""
        rows.append(dict(ijkl=(i, j, k, l), sym=sym, cls=cls, note=note))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sdf", required=True, help="ligand SDF with bond orders and explicit H")
    ap.add_argument("--top", required=True, help="simulation topology, e.g. ../build/built.pdb")
    ap.add_argument("--select", required=True, help='ligand atoms in --top, e.g. "resname PAR"')
    ap.add_argument("--ligand-indices", action="store_true",
                    help="write indices 0..n-1 of the ligand (for solute-only trajectories)")
    ap.add_argument("--ester-rotatable", action="store_true",
                    help="treat ester C(=O)-O as rotatable instead of rigid")
    ap.add_argument("--out", default="torsions", help="output prefix")
    args = ap.parse_args()

    mol = Chem.MolFromMolFile(args.sdf, removeHs=False)
    if mol is None:
        raise SystemExit(f"RDKit could not read {args.sdf}")
    if args.ester_rotatable:
        RIGID_PATTERNS.pop("ester C(=O)-O")

    top = md.load_topology(args.top)
    sel = top.select(args.select)
    sdf_el = [a.GetSymbol() for a in mol.GetAtoms()]
    top_el = [top.atom(int(x)).element.symbol for x in sel]
    if len(sdf_el) != len(top_el):
        raise SystemExit(f"Atom count differs: SDF {len(sdf_el)} vs topology selection {len(top_el)}")
    if [e.upper() for e in sdf_el] != [e.upper() for e in top_el]:
        bad = [i for i, (a, b) in enumerate(zip(sdf_el, top_el)) if a.upper() != b.upper()]
        raise SystemExit(f"Element order differs at ligand positions {bad[:10]} - "
                         "SDF and topology atom orders do not match.")
    print(f"Atom order check passed: {len(sdf_el)} atoms, elements identical.")

    names = [str(top.atom(int(x))) for x in sel]
    to_out = (lambda x: x) if args.ligand_indices else (lambda x: int(sel[x]))

    rows = classify(mol)
    files = {"rotatable": f"{args.out}.txt", "rigid": f"{args.out}_rigid.txt", "ring": f"{args.out}_ring.txt"}
    handles = {c: open(p, "w") for c, p in files.items()}
    for h in handles.values():
        h.write("# i j k l k_local   # atom names | class   (k_local: per-bond symmetry, info only;"
                " symmetry is handled by the _symops.txt file)\n")

    print(f"\n{'class':10s} {'k_loc':>5s}  torsion (atom names)                         note")
    for r in rows:
        ijkl = [to_out(x) for x in r["ijkl"]]
        label = " ".join(names[x] for x in r["ijkl"])
        print(f"{r['cls']:10s} {r['sym']:>5d}  {label:45s} {r['note']}")
        if r["cls"] in handles:
            handles[r["cls"]].write(f"{ijkl[0]} {ijkl[1]} {ijkl[2]} {ijkl[3]} {r['sym']}   "
                                    f"# {label} | {r['cls']} {r['note']}\n")
    for h in handles.values():
        h.close()

    n = {c: sum(r["cls"] == c for r in rows) for c in ("rotatable", "rigid", "trivial", "ring")}
    print(f"\n{n['rotatable']} rotatable -> {files['rotatable']}   "
          f"{n['rigid']} rigid -> {files['rigid']}   {n['ring']} ring -> {files['ring']}   "
          f"{n['trivial']} trivial dropped")
    if n["ring"]:
        print("NOTE: non-aromatic ring bonds found - ring conformations are NOT in torsions.txt.")

    # symmetry operations acting on the rotatable torsions (whole molecule, all torsions jointly)
    rot = [r["ijkl"] for r in rows if r["cls"] == "rotatable"]
    symfile = f"{args.out}_symops.txt"
    if rot:
        ops = symmetry_ops(mol, rot)
        with open(symfile, "w") as f:
            f.write(f"# {len(ops)} symmetry operations; each line = the {len(rot)} torsion quadruples "
                    "after relabelling; first line = identity\n")
            for op in ops:
                f.write(" ".join(str(to_out(x)) for q in op for x in q) + "\n")
        print(f"{len(ops)} symmetry operation(s) on the rotatable torsions -> {symfile}")
        for g, op in enumerate(ops[1:], 1):
            moved = [f"t{t}->" + " ".join(names[x] for x in q) for t, q in enumerate(op) if q != tuple(rot[t])]
            print(f"  op {g}: " + "; ".join(moved))


if __name__ == "__main__":
    main()