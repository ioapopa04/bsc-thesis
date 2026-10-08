#!/usr/bin/env python
"""
Torsion analysis of REST2 state trajectories (per-state / per-lambda .nc files).

Two modes:

1) List candidate rotatable torsions (atom indices taken from YOUR topology,
   so they match the trajectory exactly):
     python torsion_analysis.py list --top <topology> [--select "resname LIG"]
   -> writes torsions.txt (4 atom indices per line, 0-based). Check it by eye
      and delete lines you don't want (e.g. amide C-N).

2) Analyse:
     python torsion_analysis.py analyse --traj solute_state0_prod1.nc [more state0 files] \
         --top <topology> --torsions torsions.txt --out results \
         [--compare solute_state1_prod1.nc solute_state2_prod1.nc ...] \
         [--skip 0] [--stride 1] [--binwidth 10] [--temp 298]

   Outputs (in --out):
     torsions_hist.png      1D distributions + free-energy profiles per torsion
     jsd_halves.csv         JS divergence 1st vs 2nd half (convergence check)
     jsd_vs_compare.csv     JS divergence reference traj vs each --compare traj
     mi_matrix.csv          mutual information between torsion pairs (nats)
     mi_shuffled.csv        same after shuffling -> finite-sample bias baseline
     torsion_pairs.png      scatter of every torsion pair, one point per frame
     top_mi_pair_2d.png     2D histogram of the most correlated torsion pair
     dihedrals.npy          raw angles (n_frames, n_torsions), degrees 0-360
"""
import argparse
import os

import numpy as np
import mdtraj as md
import networkx as nx
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from sklearn.metrics import mutual_info_score

KB_KCAL = 0.0019872041  # kcal/(mol K)


# ---------------------------------------------------------------- list mode
def list_torsions(top_path, select, out):
    top = md.load_topology(top_path)
    atoms = top.select(select) if select else np.arange(top.n_atoms)
    atoms = set(int(a) for a in atoms)
    g = nx.Graph()
    for b in top.bonds:
        i, j = b[0].index, b[1].index
        if i in atoms and j in atoms:
            g.add_edge(i, j)
    if g.number_of_edges() == 0:
        raise SystemExit("No bonds found in topology - use a format that stores "
                         "bonds (prmtop, psf, PDB with CONECT, OpenMM xml...).")

    heavy = {a.index for a in top.atoms
             if a.element is not None and a.element.symbol != "H"}
    non_ring = {frozenset(e) for e in nx.bridges(g)}   # bridge <=> not in a ring

    def heavy_nbrs(a, exclude):
        return sorted(n for n in g.neighbors(a) if n in heavy and n != exclude)

    torsions = []
    for j, k in g.edges():
        if j not in heavy or k not in heavy or frozenset((j, k)) not in non_ring:
            continue
        nj, nk = heavy_nbrs(j, k), heavy_nbrs(k, j)
        if nj and nk:                       # both ends non-terminal
            torsions.append((nj[0], j, k, nk[0]))

    with open(out, "w") as f:
        f.write("# i j k l (0-based atom indices); central bond j-k\n")
        for t in torsions:
            names = " ".join(str(top.atom(x)) for x in t)
            f.write(f"{t[0]} {t[1]} {t[2]} {t[3]}   # {names}\n")
    print(f"{len(torsions)} candidate torsions written to {out}")
    print("NOTE: bond orders are not checked. Remove non-rotatable ones "
          "(e.g. amide C-N, exocyclic C=C / C=N) by hand.")


# ------------------------------------------------------------- analyse mode
def load_angles(trajs, top, idx, skip, stride):
    """Load one or several trajectory files (same state!), drop `skip` frames
    from the start of EACH file, concatenate, return dihedrals in degrees."""
    if isinstance(trajs, str):
        trajs = [trajs]
    out = []
    for f in trajs:
        t = md.load(f, top=top, stride=stride)[skip:]
        out.append(np.degrees(md.compute_dihedrals(t, idx)) % 360.0)
        print(f"  {f}: {t.n_frames} frames used")
    return np.concatenate(out, axis=0)


def check_same_state(files):
    import re
    states = {m.group(0) for f in files for m in [re.search(r"state\d+", os.path.basename(f))] if m}
    if len(states) > 1:
        print(f"WARNING: --traj mixes {sorted(states)}. Only frames from ONE state "
              "(lambda=1) give Boltzmann populations. Use --compare for other states.")


def hist(x, edges):
    c, _ = np.histogram(x, edges)
    return c / c.sum()


def jsd(p, q):
    return jensenshannon(p, q, base=2) ** 2        # divergence in [0, 1]


def mi_matrix(labels, shuffle=False, seed=0):
    rng = np.random.default_rng(seed)
    n = labels.shape[1]
    m = np.zeros((n, n))
    for a in range(n):
        for b in range(n):
            y = rng.permutation(labels[:, b]) if (shuffle and a != b) else labels[:, b]
            m[a, b] = mutual_info_score(labels[:, a], y)
    return m


def plot_pairs(phi, idx, out_png, lo):
    """Scatter plot of every torsion pair: one point per frame (conformation),
    coloured by frame index so you can see how sampling progressed in time."""
    n = phi.shape[1]
    if n < 2:
        return
    ang = (phi - lo) % 360.0 + lo                  # wrap to [lo, lo+360)
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]
    ncol = min(3, len(pairs))
    nrow = int(np.ceil(len(pairs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 4 * nrow), squeeze=False)
    t = np.arange(len(ang))
    for ax, (a, b) in zip(axes.flat, pairs):
        sc = ax.scatter(ang[:, a], ang[:, b], c=t, s=6, cmap="viridis", alpha=0.7)
        ax.set_xlim(lo, lo + 360)
        ax.set_ylim(lo, lo + 360)
        ax.set_xticks(np.arange(lo, lo + 361, 90))
        ax.set_yticks(np.arange(lo, lo + 361, 90))
        ax.set_xlabel(f"torsion {a} {tuple(int(x) for x in idx[a])} (deg)", fontsize=8)
        ax.set_ylabel(f"torsion {b} {tuple(int(x) for x in idx[b])} (deg)", fontsize=8)
        ax.set_aspect("equal")
        ax.grid(alpha=0.3)
    for ax in list(axes.flat)[len(pairs):]:
        ax.axis("off")
    fig.colorbar(sc, ax=axes, shrink=0.8, label="frame index (time)")
    fig.savefig(out_png, dpi=150, bbox_inches="tight")


def analyse(args):
    os.makedirs(args.out, exist_ok=True)
    idx = np.loadtxt(args.torsions, dtype=int, ndmin=2, comments="#", usecols=(0, 1, 2, 3))
    kT = KB_KCAL * args.temp
    edges = np.arange(0, 360 + args.binwidth, args.binwidth)
    centers = 0.5 * (edges[1:] + edges[:-1])

    check_same_state(args.traj)
    print("Reference ensemble:")
    phi = load_angles(args.traj, args.top, idx, args.skip, args.stride)
    n_frames, n_tor = phi.shape
    np.save(os.path.join(args.out, "dihedrals.npy"), phi)
    print(f"Total: {n_frames} frames, {n_tor} torsions")

    comps = {}
    for c in (args.compare or []):
        print(f"Comparison ensemble:")
        comps[c] = load_angles(c, args.top, idx, args.skip, args.stride)

    # 2D scatter of every torsion pair (reference ensemble only)
    plot_pairs(phi, idx, os.path.join(args.out, "torsion_pairs.png"), args.range_start)

    # 1D distributions and free-energy profiles
    fig, axes = plt.subplots(n_tor, 2, figsize=(9, 2.4 * n_tor), squeeze=False)
    for k in range(n_tor):
        p = hist(phi[:, k], edges)
        axes[k, 0].bar(centers, p, width=args.binwidth, alpha=0.6, label="ref")
        for name, arr in comps.items():
            axes[k, 0].step(centers, hist(arr[:, k], edges), where="mid",
                            label=os.path.basename(name))
        with np.errstate(divide="ignore"):
            G = -kT * np.log(p)
        G -= G[np.isfinite(G)].min()
        axes[k, 1].plot(centers, G, "o-", ms=3)
        axes[k, 0].set_ylabel(f"tor {k}\n{tuple(int(x) for x in idx[k])}", fontsize=8)
        axes[k, 1].set_ylabel("G (kcal/mol)")
    axes[0, 0].legend(fontsize=7)
    axes[-1, 0].set_xlabel("dihedral (deg)")
    axes[-1, 1].set_xlabel("dihedral (deg)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "torsions_hist.png"), dpi=150)

    # convergence: first vs second half
    h = n_frames // 2
    rows = [(k, jsd(hist(phi[:h, k], edges), hist(phi[h:, k], edges))) for k in range(n_tor)]
    np.savetxt(os.path.join(args.out, "jsd_halves.csv"), rows,
               delimiter=",", header="torsion,JSD_1st_vs_2nd_half", fmt=["%d", "%.5f"])

    # reference vs other trajectories
    if comps:
        with open(os.path.join(args.out, "jsd_vs_compare.csv"), "w") as f:
            f.write("torsion," + ",".join(os.path.basename(c) for c in comps) + "\n")
            for k in range(n_tor):
                vals = [jsd(hist(phi[:, k], edges), hist(a[:, k], edges)) for a in comps.values()]
                f.write(f"{k}," + ",".join(f"{v:.5f}" for v in vals) + "\n")

    # mutual information between torsion pairs
    labels = np.digitize(phi, edges[1:-1])
    mi = mi_matrix(labels)
    mi0 = mi_matrix(labels, shuffle=True)
    np.savetxt(os.path.join(args.out, "mi_matrix.csv"), mi, delimiter=",", fmt="%.5f")
    np.savetxt(os.path.join(args.out, "mi_shuffled.csv"), mi0, delimiter=",", fmt="%.5f")

    if n_tor > 1:
        off = mi - np.diag(np.diag(mi))
        a, b = np.unravel_index(np.argmax(off), off.shape)
        fig, ax = plt.subplots(figsize=(4.5, 4))
        H, _, _ = np.histogram2d(phi[:, a], phi[:, b], bins=[edges, edges])
        ax.imshow(H.T, origin="lower", extent=[0, 360, 0, 360], aspect="auto")
        ax.set_xlabel(f"torsion {a}")
        ax.set_ylabel(f"torsion {b}")
        ax.set_title(f"MI = {mi[a, b]:.3f} nats (shuffled {mi0[a, b]:.3f})")
        fig.tight_layout()
        fig.savefig(os.path.join(args.out, "top_mi_pair_2d.png"), dpi=150)
    print(f"Results written to {args.out}/")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    l = sub.add_parser("list")
    l.add_argument("--top", required=True)
    l.add_argument("--select", default=None, help='MDTraj selection, e.g. "resname LIG"')
    l.add_argument("--out", default="torsions.txt")

    a = sub.add_parser("analyse")
    a.add_argument("--traj", required=True, nargs="+",
                   help="one or more trajectory files of the SAME state (e.g. lambda=1)")
    a.add_argument("--top", required=True)
    a.add_argument("--torsions", required=True)
    a.add_argument("--out", default="results")
    a.add_argument("--compare", nargs="*")
    a.add_argument("--skip", type=int, default=0, help="frames to drop (after stride) as equilibration")
    a.add_argument("--stride", type=int, default=1)
    a.add_argument("--binwidth", type=float, default=10.0)
    a.add_argument("--temp", type=float, default=298.0)
    a.add_argument("--range-start", type=float, default=-180.0,
                   help="scatter-plot axes run from this value to +360 (e.g. -180 or 0)")

    args = ap.parse_args()
    if args.mode == "list":
        list_torsions(args.top, args.select, args.out)
    else:
        analyse(args)


if __name__ == "__main__":
    main()
