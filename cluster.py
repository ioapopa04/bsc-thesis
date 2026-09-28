#!/usr/bin/env python
"""
Cluster ligand conformers on symmetry-folded torsion features.

Features: each torsion phi_j with symmetry number k_j -> cos(k_j phi_j), sin(k_j phi_j),
all scaled by 1/sqrt(n_torsions). The Euclidean distance between two conformers is then

    d_AB = sqrt( (1/n) * sum_j 2 * (1 - cos(k_j * dphi_j)) )

i.e. exactly the symmetry-corrected dihedral distance of Zivanovic et al. 2020 (their eq 1).
Rough scale for one torsion: d = 0.5 <-> ~29 deg of k*phi, d = 1.0 <-> 60 deg.

Methods:
  hdbscan (default)  density-based; finds the number of clusters itself, any shape;
                     sparse frames (e.g. barrier crossings) get label -1 = noise.
                     Afterwards clusters closer than --merge-gap are merged, because
                     HDBSCAN can fragment one smooth basin into pieces.
  daura              Daura et al. (GROMOS) algorithm as used in Zivanovic 2020: frame with
                     most neighbours within --cutoff = centre, remove it + neighbours,
                     repeat until --coverage of frames is assigned or --max-clusters.

Input (either):
  --traj + --top : dihedrals computed here (needed for representative PDBs)
  --dihedrals    : dihedrals.npy from torsion_analysis.py (no PDBs written)
Symmetry numbers: --torsions torsions.txt (5th column, from select_torsions.py) or --symmetry.

Example:
  python cluster.py --traj whole_state0_prod1.nc --top ../build/built.pdb \
      --torsions torsions.txt --select "resname PAR" --out clusters

Outputs (in --out):
  clusters.csv        per cluster: frames, population, dG (kcal/mol, vs largest cluster),
                      circular mean of each torsion, representative frame
  labels.npy          cluster label per frame (-1 = noise / unassigned)
  clusters_pc.png     PC1 vs PC2 coloured by cluster (noise grey)
  clusters_torsions.png  torsion distributions per cluster
  cluster_<c>.pdb     representative conformer of each cluster (ligand only, if --traj)
"""
import argparse
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

KB_KCAL = 0.0019872041  # kcal/(mol K)


# ------------------------------------------------------------------ helpers
def features(phi_deg, k):
    r = np.radians(phi_deg) * k[None, :]
    X = np.empty((phi_deg.shape[0], 2 * phi_deg.shape[1]))
    X[:, 0::2] = np.cos(r)
    X[:, 1::2] = np.sin(r)
    return X / np.sqrt(phi_deg.shape[1])


def folded_circ_mean(phi_deg, k):
    """circular mean of k*phi, mapped back to phi in [0, 360/k)"""
    r = np.radians(phi_deg) * k[None, :]
    m = np.degrees(np.arctan2(np.sin(r).mean(0), np.cos(r).mean(0)))
    return (m % 360.0) / k


def run_hdbscan(X, min_cluster_size, min_samples):
    try:
        from sklearn.cluster import HDBSCAN
    except ImportError:
        raise SystemExit("HDBSCAN needs scikit-learn >= 1.3 (pip install -U scikit-learn)")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                       allow_single_cluster=True).fit_predict(X)


def merge_close_clusters(X, labels, gap):
    """Merge clusters whose closest members are nearer than `gap` (single linkage
    between clusters). HDBSCAN can split one smooth, continuous basin into pieces;
    pieces not separated by an empty gap are the same conformation."""
    ids = [c for c in np.unique(labels) if c >= 0]
    if len(ids) < 2 or gap <= 0:
        return labels
    parent = {c: c for c in ids}

    def find(c):
        while parent[c] != c:
            parent[c] = parent[parent[c]]
            c = parent[c]
        return c

    nn = {c: NearestNeighbors(n_neighbors=1).fit(X[labels == c]) for c in ids}
    for a_i, a in enumerate(ids):
        for b in ids[a_i + 1:]:
            d, _ = nn[b].kneighbors(X[labels == a])
            if d.min() < gap:
                parent[find(a)] = find(b)
    out = labels.copy()
    for c in ids:
        out[labels == c] = find(c)
    return out


def run_daura(X, cutoff, max_clusters, coverage):
    n = len(X)
    G = NearestNeighbors(radius=cutoff).fit(X).radius_neighbors_graph(X, mode="connectivity").tocsr()
    labels = np.full(n, -1)
    alive = np.ones(n, dtype=bool)
    for c in range(max_clusters):
        if (~alive).sum() >= coverage * n or not alive.any():
            break
        counts = G.dot(alive.astype(float))            # neighbours among remaining frames
        counts[~alive] = -1
        centre = int(np.argmax(counts))
        members = G[centre].indices
        members = members[alive[members]]
        members = np.union1d(members, [centre])
        labels[members] = c
        alive[members] = False
    return labels


def relabel_by_size(labels):
    """cluster 0 = largest; noise stays -1"""
    ids = [c for c in np.unique(labels) if c >= 0]
    order = sorted(ids, key=lambda c: -(labels == c).sum())
    new = np.full_like(labels, -1)
    for i, c in enumerate(order):
        new[labels == c] = i
    return new


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--traj", nargs="+", help="trajectory file(s) of ONE state (lambda=1)")
    src.add_argument("--dihedrals", help="dihedrals.npy (degrees) instead of --traj")
    ap.add_argument("--top", help="topology (needed with --traj)")
    ap.add_argument("--torsions", help="torsions.txt (indices; 5th column = symmetry)")
    ap.add_argument("--symmetry", type=int, nargs="+", help="symmetry numbers if not in --torsions")
    ap.add_argument("--select", default=None, help='atoms for representative PDBs, e.g. "resname PAR"')
    ap.add_argument("--skip", type=int, default=0, help="frames to drop from start of each file")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--method", choices=["hdbscan", "daura"], default="hdbscan")
    ap.add_argument("--min-cluster-size", type=float, default=0.02,
                    help="hdbscan: <1 = fraction of frames, >=1 = number of frames "
                         "(default 0.02, but never below 10 frames)")
    ap.add_argument("--min-samples", type=int, default=5,
                    help="hdbscan: neighbours defining a dense point (default 5); larger = more noise")
    ap.add_argument("--merge-gap", type=float, default=0.2,
                    help="hdbscan: merge clusters whose closest frames are nearer than this d (0 = off)")
    ap.add_argument("--cutoff", type=float, default=0.5, help="daura: distance cutoff d (default 0.5)")
    ap.add_argument("--max-clusters", type=int, default=10, help="daura (default 10)")
    ap.add_argument("--coverage", type=float, default=0.95, help="daura: stop at this fraction (0.95)")
    ap.add_argument("--temp", type=float, default=298.0)
    ap.add_argument("--out", default="clusters")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    # ---- torsions + symmetry
    tors = np.loadtxt(args.torsions, dtype=int, ndmin=2, comments="#") if args.torsions else None
    if args.symmetry:
        k = np.array(args.symmetry, dtype=int)
    elif tors is not None and tors.shape[1] >= 5:
        k = tors[:, 4]
    else:
        raise SystemExit("symmetry numbers needed: --torsions with 5th column, or --symmetry")

    # ---- dihedrals
    traj = None
    if args.traj:
        import mdtraj as md
        if args.top is None or tors is None:
            raise SystemExit("--traj needs --top and --torsions")
        parts, phis = [], []
        for f in args.traj:
            t = md.load(f, top=args.top, stride=args.stride)[args.skip:]
            phis.append(np.degrees(md.compute_dihedrals(t, tors[:, :4])) % 360.0)
            parts.append(t)
            print(f"  {f}: {t.n_frames} frames")
        traj = md.join(parts) if len(parts) > 1 else parts[0]
        phi = np.concatenate(phis)
    else:
        phi = np.load(args.dihedrals)[args.skip:]
    n_frames, n_tor = phi.shape
    if len(k) != n_tor:
        raise SystemExit(f"{len(k)} symmetry numbers for {n_tor} torsions")
    print(f"{n_frames} frames, {n_tor} torsions, symmetry {k.tolist()}, method {args.method}")

    # ---- cluster
    X = features(phi, k)
    if args.method == "hdbscan":
        mcs = int(round(args.min_cluster_size * n_frames)) if args.min_cluster_size < 1 else int(args.min_cluster_size)
        mcs = min(max(mcs, 10), n_frames)
        ms = min(args.min_samples, mcs)
        labels = run_hdbscan(X, mcs, ms)
        n_raw = len([c for c in np.unique(labels) if c >= 0])
        labels = merge_close_clusters(X, labels, args.merge_gap)
        n_new = len([c for c in np.unique(labels) if c >= 0])
        print(f"HDBSCAN min_cluster_size = {mcs}, min_samples = {ms}; {n_raw} raw clusters -> {n_new} after merging "
              f"clusters closer than gap {args.merge_gap}")
    else:
        labels = run_daura(X, args.cutoff, args.max_clusters, args.coverage)
        print(f"Daura cutoff d = {args.cutoff}")
    labels = relabel_by_size(labels)
    np.save(os.path.join(args.out, "labels.npy"), labels)

    ids = [c for c in np.unique(labels) if c >= 0]
    noise = (labels == -1).mean()
    kT = KB_KCAL * args.temp
    if not ids:
        raise SystemExit("No clusters found (all frames noise) - lower --min-cluster-size or raise --cutoff.")

    # ---- per-cluster statistics + representatives
    counts = {c: int((labels == c).sum()) for c in ids}
    pmax = max(counts.values())
    rows = []
    for c in ids:
        m = labels == c
        centre = X[m].mean(0)
        members = np.where(m)[0]
        rep = int(members[np.argmin(((X[m] - centre) ** 2).sum(1))])
        mean_t = folded_circ_mean(phi[m], k)
        rows.append((c, counts[c], counts[c] / n_frames, -kT * np.log(counts[c] / pmax), rep, mean_t, phi[rep]))

    tor_cols = ",".join(f"mean_phi{j}(k={k[j]})" for j in range(n_tor))
    rep_cols = ",".join(f"rep_phi{j}" for j in range(n_tor))
    with open(os.path.join(args.out, "clusters.csv"), "w") as f:
        f.write(f"cluster,frames,population,dG_kcal_mol,rep_frame,{tor_cols},{rep_cols}\n")
        for c, n, p, g, rep, mt, rp in rows:
            f.write(f"{c},{n},{p:.4f},{g:.3f},{rep}," + ",".join(f"{v:.1f}" for v in mt) + ","
                    + ",".join(f"{v:.1f}" for v in rp) + "\n")
        f.write(f"noise,{int((labels == -1).sum())},{noise:.4f},,,\n")

    print(f"\n{'cluster':>7s} {'frames':>7s} {'pop':>7s} {'dG(kcal/mol)':>13s} {'rep':>6s}  mean torsions (deg)")
    for c, n, p, g, rep, mt, rp in rows:
        print(f"{c:>7d} {n:>7d} {p:>7.1%} {g:>13.2f} {rep:>6d}  {np.round(mt, 1).tolist()}")
    print(f"{'noise':>7s} {int((labels == -1).sum()):>7d} {noise:>7.1%}")
    if noise > 0.2:
        print("WARNING: >20% noise frames - sampling may be too sparse or min_cluster_size too large.")

    if traj is not None:
        sub = traj.atom_slice(traj.topology.select(args.select)) if args.select else traj
        for c, n, p, g, rep, mt, rp in rows:
            sub[rep].save_pdb(os.path.join(args.out, f"cluster_{c}.pdb"))
        print(f"Representative PDBs written (cluster_<c>.pdb).")

    # ---- plots
    cmap = plt.get_cmap("tab10")
    colors = np.array([cmap(c % 10) if c >= 0 else (0.7, 0.7, 0.7, 0.5) for c in labels])
    Y = PCA(n_components=min(2, X.shape[1])).fit_transform(X)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(Y[:, 0], Y[:, 1], c=colors, s=8)
    for c in ids:
        ax.scatter([], [], color=cmap(c % 10), label=f"cluster {c} ({counts[c] / n_frames:.0%})")
    if noise > 0:
        ax.scatter([], [], color=(0.7, 0.7, 0.7), label=f"noise ({noise:.0%})")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"{args.method} clusters on torsion features")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "clusters_pc.png"), dpi=150)

    edges = np.arange(0, 370, 10)
    fig, axes = plt.subplots(n_tor, 1, figsize=(6, 2.3 * n_tor), squeeze=False)
    for j in range(n_tor):
        data = [phi[labels == c, j] for c in ids] + [phi[labels == -1, j]]
        cols = [cmap(c % 10) for c in ids] + [(0.7, 0.7, 0.7)]
        axes[j, 0].hist(data, bins=edges, stacked=True, color=cols)
        axes[j, 0].set_ylabel(f"torsion {j}\ncount")
    axes[-1, 0].set_xlabel("dihedral (deg)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "clusters_torsions.png"), dpi=150)
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()