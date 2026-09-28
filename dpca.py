#!/usr/bin/env python
"""
Dihedral PCA (dPCA) with symmetry alignment.

Each torsion phi -> cos(phi), sin(phi). Before PCA, every frame is put in the
symmetry relabelling (from torsions_symops.txt) that brings it closest to a
common reference conformer - the torsion analogue of superposing structures
before an RMSD. Symmetry-equivalent copies of a conformation thus land on the
same point, while genuinely different conformations stay apart.

Example:
  python dpca.py --traj whole_state0_prod1.nc --top ../build/built.pdb \
      --torsions torsions.txt --symops torsions_symops.txt --out pca

Note: without --labels, a conformation whose symmetry copies are equally far from
the reference (e.g. anti vs a syn reference) may be drawn as two blobs. Clustering
is not affected (cluster.py uses best-match distances); give --labels to fix the plot.

Outputs (in --out): pc1_pc2.png (coloured by time), pc1_pc2_free_energy.png,
explained_variance.png/.csv, loadings.csv (feature weights per PC),
pcs.npy (PC coordinates per frame), aligned_dihedrals.npy.
"""
import argparse
import os

import numpy as np
import mdtraj as md
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

import symtools as st

KB_KCAL = 0.0019872041


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True, nargs="+", help="trajectory file(s) of ONE state (lambda=1)")
    ap.add_argument("--top", required=True)
    ap.add_argument("--torsions", required=True)
    ap.add_argument("--symops", default=None, help="torsions_symops.txt (omit = no symmetry)")
    ap.add_argument("--labels", default=None,
                    help="labels.npy from cluster.py (same --traj/--skip/--stride): align each cluster "
                         "to its own medoid - recommended, see note below")
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--n-components", type=int, default=5)
    ap.add_argument("--temp", type=float, default=298.0)
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--out", default="pca")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    tors = np.loadtxt(args.torsions, dtype=int, ndmin=2, comments="#")
    ops = st.load_ops(args.symops, tors)
    parts = [md.load(f, top=args.top, stride=args.stride)[args.skip:] for f in args.traj]
    traj = md.join(parts) if len(parts) > 1 else parts[0]

    phiops = st.phi_all_ops(traj, ops)
    Xops = st.features(phiops)
    if args.labels:
        labels = np.load(args.labels)
        if len(labels) != Xops.shape[1]:
            raise SystemExit("labels.npy length != number of frames (use the same --traj/--skip/--stride)")
        ids = [c for c in np.unique(labels) if c >= 0]
        X0 = Xops[0]
        med = {c: int(np.where(labels == c)[0][((X0[labels == c] - X0[labels == c].mean(0)) ** 2).sum(1).argmin()])
               for c in ids}
        big = max(ids, key=lambda c: (labels == c).sum())
        X, phi_al = st.align_by_clusters(Xops, phiops, labels, med, med[big])
    else:
        X, phi_al, _ = st.iterative_align(Xops, phiops)
    X = X * np.sqrt(phiops.shape[2])                   # undo 1/sqrt(n): PCA on plain cos/sin
    n_frames, n_tor = phi_al.shape
    np.save(os.path.join(args.out, "aligned_dihedrals.npy"), phi_al)

    ncomp = min(args.n_components, X.shape[1], n_frames)
    pca = PCA(n_components=ncomp)
    Y = pca.fit_transform(X)
    evr = pca.explained_variance_ratio_
    print(f"{n_frames} frames, {n_tor} torsions, {len(ops)} symmetry operation(s) -> {X.shape[1]} features")
    for i, (v, c) in enumerate(zip(evr, np.cumsum(evr))):
        print(f"  PC{i+1}: {v:6.1%} of variance (cumulative {c:6.1%})")

    np.save(os.path.join(args.out, "pcs.npy"), Y)
    np.savetxt(os.path.join(args.out, "explained_variance.csv"),
               np.column_stack([np.arange(1, ncomp + 1), evr, np.cumsum(evr)]),
               delimiter=",", header="PC,fraction,cumulative", fmt=["%d", "%.5f", "%.5f"])
    names = [f"{fn}(phi{j})" for j in range(n_tor) for fn in ("cos", "sin")]
    with open(os.path.join(args.out, "loadings.csv"), "w") as f:
        f.write("feature," + ",".join(f"PC{i+1}" for i in range(ncomp)) + "\n")
        for fi, nm in enumerate(names):
            f.write(nm + "," + ",".join(f"{pca.components_[i, fi]:.4f}" for i in range(ncomp)) + "\n")

    fig, ax = plt.subplots(figsize=(5, 3.5))
    x = np.arange(1, ncomp + 1)
    ax.bar(x, evr, label="per PC")
    ax.plot(x, np.cumsum(evr), "ko-", label="cumulative")
    ax.set_xticks(x)
    ax.set_xlabel("principal component")
    ax.set_ylabel("fraction of variance")
    ax.set_ylim(0, 1.05)
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "explained_variance.png"), dpi=150)

    if ncomp < 2:
        return
    l1, l2 = f"PC1 ({evr[0]:.0%})", f"PC2 ({evr[1]:.0%})"
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sc = ax.scatter(Y[:, 0], Y[:, 1], c=np.arange(n_frames), s=8, cmap="viridis", alpha=0.7)
    fig.colorbar(sc, label="frame index (time)")
    ax.set_xlabel(l1)
    ax.set_ylabel(l2)
    ax.set_title("dihedral PCA (symmetry-aligned)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "pc1_pc2.png"), dpi=150)

    H, xe, ye = np.histogram2d(Y[:, 0], Y[:, 1], bins=args.bins)
    with np.errstate(divide="ignore"):
        G = -KB_KCAL * args.temp * np.log(H / H.sum())
    G -= G[np.isfinite(G)].min()
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(np.ma.masked_invalid(G).T, origin="lower", aspect="auto",
                   extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap="viridis_r")
    fig.colorbar(im, label="G (kcal/mol)")
    ax.set_xlabel(l1)
    ax.set_ylabel(l2)
    ax.set_title("free energy on PC1/PC2 (white = never visited)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "pc1_pc2_free_energy.png"), dpi=150)
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()