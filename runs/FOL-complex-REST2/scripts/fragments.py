"""Per-fragment drift of folate at state 0: pteridine, pABA, glutamate. RMSD to the crystal after CA superposition,
per 5 ns block, plus key contacts (3QL3: Arg57 guanidinium to the alpha-carboxyl; Asp27 to the pteridine)."""
import sys
import numpy as np, mdtraj as md
run, ref_top = sys.argv[1], sys.argv[2]
ref = md.load(ref_top); t = md.load(f"{run}/solute_state0_prod1.nc", top=ref_top)
top = ref.topology
ca = top.select("protein and name CA")
t.superpose(ref, atom_indices=ca)
frag = {"pteridine": "C11 C12 C13 C14 C15 C16 N3 N4 N5 N6 N7 O4",
        "pABA+linker": "C4 C5 C6 C7 C8 C9 N2 C10 C3 O3",
        "glutamate": "N1 C2 C1 O1 O2 C17 C18 C19 O5 O6"}
nb = 10; L = t.n_frames // nb
print(f"{run}: {t.n_frames} frames; RMSD (A) per 5 ns block")
print("block(ns)  " + "  ".join(f"{k:>11s}" for k in frag))
sel = {k: top.select("resname FOL and name " + v) for k, v in frag.items()}
for b in range(nb):
    s = slice(b * L, (b + 1) * L)
    row = [np.sqrt(((t.xyz[s][:, i] - ref.xyz[0, i]) ** 2).sum(-1).mean(-1)).mean() * 10 for i in sel.values()]
    print(f"{b*5:3d}-{b*5+5:<3d}     " + "  ".join(f"{v:11.2f}" for v in row))
if len(sys.argv) > 3:
    for label, a, b in [("Arg57 NH1/NH2/NE - FOL O1/O2 (alpha-COOH)", "resid 56 and resname ARG and name NH1 NH2 NE", "resname FOL and name O1 O2"),
                        ("Asp27 OD1/OD2 - FOL N5/N7 (pterin N1/NA2... package names)", "resid 26 and resname ASP and name OD1 OD2", "resname FOL and name N5 N6 N7"),
                        ("NAP (any heavy atom) - FOL pteridine", "resname NAP and not element H", "resname FOL and name C11 C12 C13 C14 C15 C16 N3 N4 N5 N6 N7 O4")]:
        ia, ib = top.select(a), top.select(b)
        if len(ia) == 0 or len(ib) == 0:
            print(f"{label}: selection empty ({len(ia)}, {len(ib)})"); continue
        pairs = np.array([(i, j) for i in ia for j in ib])
        d = md.compute_distances(t, pairs).min(1) * 10
        print(f"{label}: min dist per 5 ns block (A): " + " ".join(f"{d[b*L:(b+1)*L].mean():.1f}" for b in range(nb))
              + f"   crystal {md.compute_distances(ref, pairs).min()*10:.1f}")
