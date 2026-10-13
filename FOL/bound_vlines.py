#!/usr/bin/env python
"""Overlay the two protein-bound folic-acid conformers (4LRH, 3QL3) on the full
100 ns state-0 ensemble torsion distributions, as plt.vlines.

The crystal ligands are mapped onto the built topology by enumerating every graph
isomorphism (RDKit AssignBondOrdersFromTemplate + GetSubstructMatches), taking the
mapping that puts the whole conformer closest to the MD ensemble -- the same
mapping validated in FOLn-REST2/REST2-run1/bound_check/assign_bound.py.
"""
import numpy as np
import mdtraj as md
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem

ROOT = "/scratch/ipopa/FOLn-REST2-50ns"
HERE = f"{ROOT}/REST2-run1"
SDF = f"{ROOT}/build/FOL.sdf"
TOP = f"{ROOT}/build/built.solute.pdb"
TRAJ = f"{HERE}/solute_state0_prod1.nc"
TORSFILE = f"{HERE}/ioana_torsions.txt"
BOUND = "/scratch/ipopa/FOLn-REST2/REST2-run1/bound_check"
STRIDE = 5
KEEP = [1, 2, 3, 4, 6, 7]          # torsion indices --drop-free kept (for the title)


def dihedral(p0, p1, p2, p3):
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
    y = np.dot(m, n2) / (np.linalg.norm(m) * np.linalg.norm(n2))
    return float(np.degrees(np.arctan2(y, x)) % 360.0)


# ---- reference ligand: built-topology index of every heavy atom ----
full = Chem.MolFromMolFile(SDF, removeHs=False)
heavy_orig = [a.GetIdx() for a in full.GetAtoms() if a.GetAtomicNum() > 1]
template = Chem.RemoveHs(full)
tmpl_elem = [a.GetSymbol() for a in template.GetAtoms()]


def crystal_mol(pdb, chain, resname, altloc="A"):
    block = []
    for ln in open(pdb):
        if ln.startswith(("HETATM", "ATOM")) and ln[17:20].strip() == resname and ln[21] == chain:
            if ln[16] not in (" ", altloc):
                continue
            block.append(ln[:16] + " " + ln[17:])
    if not block:
        raise SystemExit(f"no {resname} in chain {chain} of {pdb}")
    m = Chem.MolFromPDBBlock("".join(block) + "END\n", sanitize=False, removeHs=False)
    return AllChem.AssignBondOrdersFromTemplate(template, m)


# ---- torsions ----
tors = np.loadtxt(TORSFILE, dtype=int, ndmin=2, comments="#")
quads = tors[:, :4]
n_tor = len(quads)

traj = md.load(TRAJ, top=TOP, stride=STRIDE)
phi_md = np.degrees(md.compute_dihedrals(traj, quads)) % 360.0     # (frames, n_tor)
print(f"ensemble: {traj.n_frames} frames, {n_tor} torsions")

md_heavy = traj.xyz[:, heavy_orig, :]


def crystal_phi(xyz, match):
    """Torsion quadruplets for the crystal, in built-topology atom order."""
    out = np.empty(n_tor)
    for t in range(n_tor):
        a, b, c, d = quads[t]
        ia, ib, ic, id_ = (match[heavy_orig.index(x)] for x in (a, b, c, d))
        out[t] = dihedral(xyz[ia], xyz[ib], xyz[ic], xyz[id_])
    return out


def chord(a, b):
    return 2.0 * np.abs(np.sin(np.radians(a - b) / 2.0))


def best_mapping_phi(pdb, chain, altloc):
    """Enumerate all isomorphisms; keep the one whose conformer is nearest the ensemble."""
    m = crystal_mol(pdb, chain, "FOL", altloc)
    xyz = m.GetConformer().GetPositions()
    matches = [list(x) for x in m.GetSubstructMatches(template, uniquify=False, maxMatches=10000)]
    best = None
    for match in matches:
        phi_cr = crystal_phi(xyz, match)
        # max-chord distance from this conformer to every ensemble frame (raw labels)
        d = np.max(chord(phi_cr[None, :], phi_md), axis=1).min()
        if best is None or d < best[0]:
            best = (d, phi_cr, match)
    d, phi_cr, match = best
    elems = [m.GetAtomWithIdx(c).GetSymbol() for c in match]
    print(f"{pdb.split('/')[-1]} chain {chain} altloc {altloc}: {len(matches)} isomorphisms, "
          f"element match {elems == tmpl_elem}, NN distance {d:.3f}")
    return phi_cr


bound = {
    "4LRH (hFR-alpha)": best_mapping_phi(f"{BOUND}/4LRH.pdb", "A", "A"),
    "3QL3 (DHFR)": best_mapping_phi(f"{BOUND}/3QL3.pdb", "A", "A"),
}
for name, phi in bound.items():
    print(f"  {name} crystal torsions (deg): " + ", ".join(f"{v:.0f}" for v in phi))

# ---- figure: 3x3 panels, one per rotatable torsion ----
ncol = 3
nrow = int(np.ceil(n_tor / ncol))
fig, axes = plt.subplots(nrow, ncol, figsize=(5.0 * ncol, 3.3 * nrow),
                         layout="constrained", squeeze=False)
colors = {"4LRH (hFR-alpha)": "#d62728", "3QL3 (DHFR)": "#1f77b4"}
edges = np.arange(0, 361, 5)
for k in range(n_tor):
    ax = axes.flat[k]
    ax.hist(phi_md[:, k], bins=edges, density=True, color="0.75", edgecolor="none")
    for name, phi in bound.items():
        ax.axvline(phi[k], color=colors[name], lw=2.0, ls="--",
                   label=name if k == 0 else None)
    a, b, c, d = quads[k]
    tag = "kept" if k in KEEP else "free (dropped)"
    ax.set_title(f"torsion {k}: {traj.topology.atom(a)}-{traj.topology.atom(b)}-"
                 f"{traj.topology.atom(c)}-{traj.topology.atom(d)}  [{tag}]", fontsize=9)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 90))
    ax.set_yticks([])
axes.flat[0].legend(fontsize=9, loc="upper right")
for extra in axes.flat[n_tor:]:
    extra.axis("off")
fig.suptitle("State-0 folic-acid torsion distributions (100 ns, 10 000 frames); "
             "dashed lines = protein-bound conformers", fontsize=11)
out = f"{ROOT}/fol-clustering-ioana/bound_vlines.png"
fig.savefig(out, dpi=150)
print(f"\nwritten {out}")
