"""Ligand RMSD to the crystal conformation over time: folic acid in water (vs 4LRH and vs 3QL3)
and in its complexes (4LRH vs 4LRH, 3QL3 vs 3QL3). Four lines in one figure.

  python -I scripts/rmsd_vs_crystal.py <out-prefix>

PLOT: one point per saved frame, no averaging and no lines.

CONFORMATION, NOT POSE: every frame's folate heavy atoms are superposed onto the crystal folate
(Kabsch) before the RMSD, so the number measures how different the ligand's SHAPE is from the
crystal's. Water has no protein to align on, and the complexes are treated the same way so the four
lines are comparable.

SYMMETRY: the crystal has no hydrogens, so atoms that look alike on the heavy-atom graph (the two O
of each carboxyl, the two sides of the benzoyl ring) can be assigned either way. Every graph match
between crystal and simulation is tried and, per frame, the lowest RMSD is kept (as RDKit's
GetBestRMS does). Heavy atoms only (32).

All three trajectories use the SAME folate parameter package and atom order (asserted), state 0.
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
from rdkit import Chem

OUT = sys.argv[1]
WATER = "/scratch/ipopa/thesis/FOLn-REST2-50ns"
CPLX = "/scratch/ipopa/thesis/FOL-complex-REST2"
CRYSTAL = "/scratch/ipopa/thesis/bsc-thesis/crystal"

RUNS = {   # label: (trajectory, topology)
    "water": (f"{WATER}/REST2-run1/solute_state0_prod1.nc", f"{WATER}/build/built.solute.pdb"),
    "4LRH complex": (f"{CPLX}/4LRH/analysis-shuyu/FOL_state0_prod1.nc", f"{WATER}/build/built.solute.pdb"),
    "3QL3 complex": (f"{CPLX}/3QL3/analysis-shuyu/FOL_state0_prod1.nc", f"{WATER}/build/built.solute.pdb"),
}
LINES = [  # (run, crystal, colour, style)
    ("water", "4LRH", "#d62728", "-"),
    ("water", "3QL3", "#1f77b4", "-"),
    ("4LRH complex", "4LRH", "#7f0000", "--"),
    ("3QL3 complex", "3QL3", "#08306b", "--"),
]


def skeleton(mol):
    """Heavy-atom graph, every bond single, no aromaticity or charge: element + connectivity."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rw, new = Chem.RWMol(), {}
    for i in heavy:
        atom = Chem.Atom(mol.GetAtomWithIdx(i).GetAtomicNum())
        atom.SetNoImplicit(True)
        new[i] = rw.AddAtom(atom)
    for b in mol.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in new and c in new:
            rw.AddBond(new[a], new[c], Chem.BondType.SINGLE)
    g = rw.GetMol()
    g.UpdatePropertyCache(strict=False)
    return g, heavy


def kabsch_rmsd(X, Y):
    """RMSD after optimal superposition of every frame of X (F, n, 3) onto Y (n, 3)."""
    Xc = X - X.mean(axis=1, keepdims=True)
    Yc = Y - Y.mean(axis=0)
    H = np.einsum("fni,nj->fij", Xc, Yc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    S[:, -1] *= d                                     # no reflections
    e0 = (Xc ** 2).sum(axis=(1, 2)) + (Yc ** 2).sum()
    msd = np.maximum(e0 - 2 * S.sum(axis=1), 0.0) / X.shape[1]
    return np.sqrt(msd)


sim = Chem.MolFromMolFile(f"{WATER}/build/FOL.sdf", removeHs=False)
sim_g, sim_heavy = skeleton(sim)
ref_atoms = [a.name for a in md.load(f"{WATER}/build/built.solute.pdb").topology.atoms]

crystals = {}
for pdb in ("4LRH", "3QL3"):
    ref = Chem.MolFromMolFile(f"{CRYSTAL}/fol_{pdb}_A.sdf", removeHs=False, sanitize=False)
    ref_g, ref_heavy = skeleton(ref)
    xyz = ref.GetConformer().GetPositions() / 10.0                 # A -> nm, as mdtraj
    matches = sim_g.GetSubstructMatches(ref_g, uniquify=False, useChirality=False, maxMatches=1000)
    assert matches, pdb
    # for each match: crystal coordinates re-ordered into simulation heavy-atom order
    orders = []
    for m in matches:                                  # m[ref_graph_idx] = sim_graph_idx
        inv = {s: r for r, s in enumerate(m)}
        orders.append(np.array([xyz[ref_heavy[inv[k]]] for k in range(len(sim_heavy))]))
    crystals[pdb] = orders
    print(f"{pdb}: {len(matches)} symmetry-equivalent atom assignments tried per frame")

rows, traces = [], {}
for run, (traj, top) in RUNS.items():
    t = md.load(traj, top=top)
    assert [a.name for a in t.topology.atoms] == ref_atoms, f"{run}: atom order differs"
    X = t.xyz[:, sim_heavy, :]
    for run_, pdb, colour, style in LINES:
        if run_ != run:
            continue
        r = np.min([kabsch_rmsd(X, Y) for Y in crystals[pdb]], axis=0) * 10.0    # nm -> A
        ns = t.time / 1000.0
        traces[(run, pdb)] = (ns, r)
        rows.append((run, pdb, len(r), ns[-1], r.mean(), np.median(r), r.min(), np.percentile(r, 5)))

fig, ax = plt.subplots(figsize=(11, 4.8), layout="constrained")
for run, pdb, colour, style in LINES:
    ns, r = traces[(run, pdb)]
    label = (f"folic acid in water vs {pdb} crystal" if run == "water"
             else f"folic acid in {pdb} complex vs {pdb} crystal")
    # one point per saved frame (every 2 ps), no averaging and no connecting lines
    ax.scatter(ns, r, s=1.5, color=colour, alpha=0.35, linewidths=0, rasterized=True,
               marker="o" if run == "water" else "s", label=label)
ax.set_xlabel("time (ns)")
ax.set_ylabel("heavy-atom RMSD to crystal conformation (Å)")
ax.set_xlim(0, max(tr[0][-1] for tr in traces.values()))
ax.set_ylim(bottom=0)
ax.grid(alpha=0.3)
leg = ax.legend(loc="upper right", fontsize=9, markerscale=6)
for h in leg.legend_handles:
    h.set_alpha(1)
ax.set_title("Folic acid: RMSD of the ligand conformation to the crystal over time\n"
             "REST2 state 0; folate superposed on folate, best over symmetry-equivalent atoms; "
             "one point per saved frame (every 2 ps)", fontsize=10)
fig.savefig(f"{OUT}.png", dpi=150)

with open(f"{OUT}.csv", "w") as f:
    f.write("run,crystal,time_ns,rmsd_A\n")
    for (run, pdb), (ns, r) in traces.items():
        for a, b in zip(ns, r):
            f.write(f"{run},{pdb},{a:.3f},{b:.4f}\n")
print("\n| run | vs crystal | frames | length (ns) | mean RMSD (Å) | median | min | 5th pct |")
print("|---|---|---|---|---|---|---|---|")
for run, pdb, n, L, mean, med, mn, p5 in rows:
    print(f"| {run} | {pdb} | {n} | {L:.0f} | {mean:.2f} | {med:.2f} | {mn:.2f} | {p5:.2f} |")
print(f"\nwritten {OUT}.png and {OUT}.csv")
