#!/usr/bin/env python
"""
Cluster ligand conformers on torsions, respecting the molecule's symmetry.

Features: each torsion phi -> cos(phi), sin(phi), scaled by 1/sqrt(n_torsions).
Symmetry: two conformers are compared under every symmetry operation of the
molecule (from torsions_symops.txt, written by select_torsions.py) and the
SMALLEST distance is used:
      d(A, B) = min_g || x(A) - x(g B) ||
Operations act on all torsions jointly (e.g. a ring flip shifts both torsions of
R1-ring-R2 by 180 deg together), so syn/anti-type conformations stay distinct
while relabelled copies of the same conformation merge.
Distance scale (one torsion): d = 0.5 <-> ~29 deg, d = 1.0 <-> 60 deg.

Methods (on the precomputed best-match distances):
  hdbscan (default)  density-based: finds the dense core of each basin (number of
                     clusters not preset); then lower-density frames are assigned to the basin they are
                     connected to by steps < --assign-radius. Only frames cut off by
                     an empty gap (e.g. isolated barrier crossings) stay noise (-1).
  daura              Daura/GROMOS algorithm with a distance --cutoff.

Example:
  python cluster.py --traj whole_state0_prod1.nc --top ../build/built.pdb \
      --torsions torsions.txt --symops torsions_symops.txt --select "resname PAR" --out clusters

Outputs (in --out): clusters.csv, labels.npy, aligned_dihedrals.npy, clusters_pc.png, clusters_torsions.png,
cluster_<c>.pdb (representative = cluster medoid, ligand only if --select).
Torsion means in clusters.csv are circular means after aligning every member to
the cluster representative (so symmetry-related copies don't average out).
"""
import argparse
import os
import warnings

import numpy as np
import mdtraj as md
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import symtools as st

KB_KCAL = 0.0019872041  # kcal/(mol K)


def run_hdbscan(D, mcs, ms):
    try:
        from sklearn.cluster import HDBSCAN
    except ImportError:
        raise SystemExit("HDBSCAN needs scikit-learn >= 1.3 (pip install -U scikit-learn)")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return HDBSCAN(min_cluster_size=mcs, min_samples=ms, metric="precomputed",
                       allow_single_cluster=True).fit_predict(D)


def merge_close_clusters(D, labels, gap):
    ids = [c for c in np.unique(labels) if c >= 0]
    if len(ids) < 2 or gap <= 0:
        return labels
    parent = {c: c for c in ids}

    def find(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    members = {c: np.where(labels == c)[0] for c in ids}
    for a_i, a in enumerate(ids):
        for b in ids[a_i + 1:]:
            if D[np.ix_(members[a], members[b])].min() < gap:
                parent[find(a)] = find(b)
    out = labels.copy()
    for c in ids:
        out[labels == c] = find(c)
    return out


def assign_to_cores(D, labels, radius):
    """Give noise frames the label of their nearest assigned frame if it is closer
    than `radius`, repeatedly, so labels spread outwards from each dense core along
    continuous paths of frames (closest first). Frames separated from every cluster
    by a gap > radius stay noise. HDBSCAN on its own labels only the dense core of a
    basin; its lower-density edges are still part of that basin."""
    labels = labels.copy()
    assigned = labels >= 0
    if not assigned.any() or radius <= 0:
        return labels, 0
    Da = D[:, assigned]
    best_d = Da.min(1)
    best_l = labels[assigned][Da.argmin(1)]
    n_new = 0
    while True:
        cand = np.where(~assigned & (best_d < radius))[0]
        if len(cand) == 0:
            break
        # assign only the closest layer this round, so labels grow outward in order
        layer = cand[best_d[cand] <= best_d[cand].min() + 0.25 * radius]
        labels[layer] = best_l[layer]
        assigned[layer] = True
        n_new += len(layer)
        dn = D[:, layer]
        j = dn.argmin(1)
        closer = dn[np.arange(len(D)), j] < best_d
        best_d = np.where(closer, dn[np.arange(len(D)), j], best_d)
        best_l = np.where(closer, labels[layer][j], best_l)
    return labels, n_new


def run_daura(D, cutoff, max_clusters, coverage):
    n = len(D)
    nb = D < cutoff
    labels = np.full(n, -1)
    alive = np.ones(n, dtype=bool)
    for c in range(max_clusters):
        if (~alive).sum() >= coverage * n or not alive.any():
            break
        counts = (nb[:, alive]).sum(1).astype(float)
        counts[~alive] = -1
        centre = int(np.argmax(counts))
        mem = np.where(nb[centre] & alive)[0]
        labels[np.union1d(mem, [centre])] = c
        alive[labels == c] = False
    return labels


def relabel_by_size(labels):
    ids = sorted([c for c in np.unique(labels) if c >= 0], key=lambda c: -(labels == c).sum())
    new = np.full_like(labels, -1)
    for i, c in enumerate(ids):
        new[labels == c] = i
    return new


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True, nargs="+", help="trajectory file(s) of ONE state (lambda=1)")
    ap.add_argument("--top", required=True)
    ap.add_argument("--torsions", required=True, help="torsions.txt from select_torsions.py")
    ap.add_argument("--symops", default=None,
                    help="torsions_symops.txt from select_torsions.py (omit = no symmetry)")
    ap.add_argument("--select", default=None, help='atoms for representative PDBs, e.g. "resname PAR"')
    ap.add_argument("--skip", type=int, default=0, help="frames to drop from start of each file")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--method", choices=["hdbscan", "daura"], default="hdbscan")
    ap.add_argument("--min-cluster-size", type=float, default=0.02,
                    help="hdbscan: <1 = fraction of frames, >=1 = frames (default 0.02, min 10)")
    ap.add_argument("--min-samples", type=float, default=0.01,
                    help="hdbscan: density smoothing; <1 = fraction of frames, >=1 = frames "
                         "(default 0.01, min 5). Too small -> one basin splits into fragments")
    ap.add_argument("--merge-gap", type=float, default=0.0,
                    help="hdbscan: merge clusters whose closest frames are nearer than this "
                         "(default 0 = off; can wrongly merge basins joined by barrier crossings)")
    ap.add_argument("--assign-radius", type=float, default=0.2,
                    help="hdbscan: noise frames connected to a cluster by steps shorter than this are "
                         "assigned to it (0 = off, keep HDBSCAN's noise)")
    ap.add_argument("--cutoff", type=float, default=0.5, help="daura: distance cutoff (default 0.5)")
    ap.add_argument("--max-clusters", type=int, default=10, help="daura (default 10)")
    ap.add_argument("--coverage", type=float, default=0.95, help="daura (default 0.95)")
    ap.add_argument("--temp", type=float, default=298.0)
    ap.add_argument("--out", default="clusters")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    tors = np.loadtxt(args.torsions, dtype=int, ndmin=2, comments="#")
    ops = st.load_ops(args.symops, tors)
    if args.symops is None:
        print("NOTE: no --symops given - symmetry-equivalent conformations will NOT be merged.")

    parts = []
    for f in args.traj:
        t = md.load(f, top=args.top, stride=args.stride)[args.skip:]
        print(f"  {f}: {t.n_frames} frames")
        parts.append(t)
    traj = md.join(parts) if len(parts) > 1 else parts[0]

    phiops = st.phi_all_ops(traj, ops)                 # (n_ops, frames, n_tor)
    Xops = st.features(phiops)
    n_ops, n_frames, n_tor = phiops.shape
    print(f"{n_frames} frames, {n_tor} torsions, {n_ops} symmetry operation(s), method {args.method}")

    D = st.best_match_distances(Xops)

    if args.method == "hdbscan":
        mcs = int(round(args.min_cluster_size * n_frames)) if args.min_cluster_size < 1 else int(args.min_cluster_size)
        mcs = min(max(mcs, 10), n_frames)
        ms = int(round(args.min_samples * n_frames)) if args.min_samples < 1 else int(args.min_samples)
        ms = min(max(ms, 5), mcs)
        labels = run_hdbscan(D, mcs, ms)
        n_raw = len([c for c in np.unique(labels) if c >= 0])
        labels = merge_close_clusters(D, labels, args.merge_gap)
        n_new = len([c for c in np.unique(labels) if c >= 0])
        n_core = int((labels >= 0).sum())
        labels, n_ass = assign_to_cores(D, labels, args.assign_radius)
        msg = f"HDBSCAN min_cluster_size = {mcs}, min_samples = {ms}; {n_raw} cluster(s)"
        if args.merge_gap > 0:
            msg += f" -> {n_new} after merging clusters closer than {args.merge_gap}"
        print(msg)
        print(f"  {n_core} frames in dense cores; {n_ass} lower-density frames assigned to their basin "
              f"(path of frames closer than {args.assign_radius})")
    else:
        labels = run_daura(D, args.cutoff, args.max_clusters, args.coverage)
        print(f"Daura cutoff d = {args.cutoff}")
    labels = relabel_by_size(labels)
    np.save(os.path.join(args.out, "labels.npy"), labels)

    ids = [c for c in np.unique(labels) if c >= 0]
    if not ids:
        raise SystemExit("No clusters found (all noise) - lower --min-cluster-size or raise --cutoff.")
    noise = (labels == -1).mean()
    kT = KB_KCAL * args.temp
    counts = {c: int((labels == c).sum()) for c in ids}
    pmax = max(counts.values())

    rows = []
    for c in ids:
        mem = np.where(labels == c)[0]
        rep = int(mem[D[np.ix_(mem, mem)].sum(1).argmin()])        # medoid
        _, phi_al, _ = st.align(Xops[:, mem], phiops[:, mem], Xops[0, rep])
        rows.append(dict(c=c, n=counts[c], p=counts[c] / n_frames,
                         dG=-kT * np.log(counts[c] / pmax) + 0.0, rep=rep,
                         mean=st.circ_mean(phi_al), rep_phi=phiops[0, rep]))

    with open(os.path.join(args.out, "clusters.csv"), "w") as f:
        f.write("cluster,frames,population,dG_kcal_mol,rep_frame,"
                + ",".join(f"mean_phi{j}" for j in range(n_tor)) + ","
                + ",".join(f"rep_phi{j}" for j in range(n_tor)) + "\n")
        for r in rows:
            f.write(f"{r['c']},{r['n']},{r['p']:.4f},{r['dG']:.3f},{r['rep']},"
                    + ",".join(f"{v:.1f}" for v in r["mean"]) + ","
                    + ",".join(f"{v:.1f}" for v in r["rep_phi"]) + "\n")
        f.write(f"noise,{int((labels == -1).sum())},{noise:.4f},,,\n")

    print(f"\n{'cluster':>7s} {'frames':>7s} {'pop':>7s} {'dG(kcal/mol)':>13s} {'rep':>6s}  mean torsions (deg)")
    for r in rows:
        print(f"{r['c']:>7d} {r['n']:>7d} {r['p']:>7.1%} {r['dG']:>13.2f} {r['rep']:>6d}  "
              f"{[round(float(v), 1) for v in r['mean']]}")
    print(f"{'noise':>7s} {int((labels == -1).sum()):>7d} {noise:>7.1%}")
    if noise > 0.2:
        print("WARNING: >20% noise frames - sampling may be too sparse or min_cluster_size too large.")

    sub = traj.atom_slice(traj.topology.select(args.select)) if args.select else traj
    for r in rows:
        sub[r["rep"]].save_pdb(os.path.join(args.out, f"cluster_{r['c']}.pdb"))
    print("Representative PDBs written (cluster_<c>.pdb).")

    # ---- plots: each cluster aligned to its own medoid, medoids aligned to the largest cluster's
    Xal, phial = st.align_by_clusters(Xops, phiops, labels, {r["c"]: r["rep"] for r in rows}, rows[0]["rep"])
    np.save(os.path.join(args.out, "aligned_dihedrals.npy"), phial)
    Y = PCA(n_components=min(2, Xal.shape[1])).fit_transform(Xal)
    cmap = plt.get_cmap("tab10")
    colors = [cmap(c % 10) if c >= 0 else (0.7, 0.7, 0.7, 0.5) for c in labels]
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(Y[:, 0], Y[:, 1], c=colors, s=8)
    for c in ids:
        ax.scatter([], [], color=cmap(c % 10), label=f"cluster {c} ({counts[c] / n_frames:.0%})")
    if noise > 0:
        ax.scatter([], [], color=(0.7, 0.7, 0.7), label=f"noise ({noise:.0%})")
    ax.set_xlabel("PC1 (symmetry-aligned)")
    ax.set_ylabel("PC2 (symmetry-aligned)")
    ax.set_title(f"{args.method} clusters")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "clusters_pc.png"), dpi=150)

    edges = np.arange(0, 370, 10)
    fig, axes = plt.subplots(n_tor, 1, figsize=(6, 2.3 * n_tor), squeeze=False)
    for j in range(n_tor):
        data = [phial[labels == c, j] for c in ids] + [phial[labels == -1, j]]
        axes[j, 0].hist(data, bins=edges, stacked=True, color=[cmap(c % 10) for c in ids] + [(0.7, 0.7, 0.7)])
        axes[j, 0].set_ylabel(f"torsion {j}\ncount")
    axes[-1, 0].set_xlabel("dihedral (deg), symmetry-aligned")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "clusters_torsions.png"), dpi=150)
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()