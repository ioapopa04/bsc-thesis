#!/usr/bin/env python
"""
Dihedral PCA (dPCA) of torsion angles, with user-provided symmetry numbers.

Input: dihedrals.npy written by torsion_analysis.py (shape n_frames x n_torsions,
degrees). Column order = line order in torsions.txt.

Each torsion phi_j with symmetry number k_j becomes two features
    cos(k_j * phi_j), sin(k_j * phi_j)
so that phi and phi + 360/k are the same point (e.g. k=2 for a symmetric phenyl).
PCA is then run on these 2n features (centred, NOT scaled).

Example (paracetamol: torsion 0 = amide, k=1; torsion 1 = aryl C-N, k=2):
  python dpca.py --dihedrals results/dihedrals.npy --symmetry 1 2 --out pca
  python dpca.py --dihedrals results/dihedrals.npy --torsions torsions.txt --out pca

Outputs (in --out):
  pc1_pc2.png            scatter of PC1 vs PC2, one point per frame, coloured by time
  pc1_pc2_free_energy.png  -kT ln P(PC1, PC2) in kcal/mol (2D free-energy map)
  explained_variance.png   variance explained per PC and cumulative
  explained_variance.csv
  loadings.csv           weight of each feature (cos/sin of each torsion) in each PC
  pcs.npy                PC coordinates of every frame (n_frames x n_components)
"""
import argparse
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

KB_KCAL = 0.0019872041  # kcal/(mol K)


def build_features(phi_deg, k):
    """Interleaved features: cos(k0 phi0), sin(k0 phi0), cos(k1 phi1), sin(k1 phi1), ..."""
    r = np.radians(phi_deg) * k[None, :]
    X = np.empty((phi_deg.shape[0], 2 * phi_deg.shape[1]))
    X[:, 0::2] = np.cos(r)
    X[:, 1::2] = np.sin(r)
    names = []
    for j, kj in enumerate(k):
        s = f"{kj}*phi{j}" if kj != 1 else f"phi{j}"
        names += [f"cos({s})", f"sin({s})"]
    return X, names


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dihedrals", required=True, help="dihedrals.npy (degrees)")
    ap.add_argument("--symmetry", type=int, nargs="+", default=None,
                    help="symmetry number per torsion, same order as torsions.txt (1 = none)")
    ap.add_argument("--torsions", default=None,
                    help="torsions.txt from select_torsions.py: symmetry read from its 5th column")
    ap.add_argument("--n-components", type=int, default=5)
    ap.add_argument("--skip", type=int, default=0, help="drop first N frames")
    ap.add_argument("--temp", type=float, default=298.0)
    ap.add_argument("--bins", type=int, default=40, help="bins per axis for the free-energy map")
    ap.add_argument("--out", default="pca")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    phi = np.load(args.dihedrals)[args.skip:]
    if phi.ndim != 2:
        raise SystemExit("dihedrals.npy must be 2D (n_frames, n_torsions)")
    n_frames, n_tor = phi.shape
    if args.torsions:
        tors = np.loadtxt(args.torsions, dtype=int, ndmin=2, comments="#")
        if tors.shape[1] < 5:
            raise SystemExit(f"{args.torsions} has no symmetry column - use --symmetry instead")
        k = tors[:, 4]
    elif args.symmetry:
        k = np.array(args.symmetry, dtype=int)
    else:
        raise SystemExit("give either --torsions (with symmetry column) or --symmetry")
    if len(k) != n_tor:
        raise SystemExit(f"--symmetry needs {n_tor} values (one per torsion), got {len(k)}")
    if np.any(k < 1):
        raise SystemExit("symmetry numbers must be >= 1")

    X, names = build_features(phi, k)
    ncomp = min(args.n_components, X.shape[1], n_frames)
    pca = PCA(n_components=ncomp)
    Y = pca.fit_transform(X)                     # centred automatically, no scaling
    evr = pca.explained_variance_ratio_

    print(f"{n_frames} frames, {n_tor} torsions -> {X.shape[1]} features, symmetry {k.tolist()}")
    for i, (v, c) in enumerate(zip(evr, np.cumsum(evr))):
        print(f"  PC{i+1}: {v:6.1%} of variance (cumulative {c:6.1%})")

    np.save(os.path.join(args.out, "pcs.npy"), Y)
    np.savetxt(os.path.join(args.out, "explained_variance.csv"),
               np.column_stack([np.arange(1, ncomp + 1), evr, np.cumsum(evr)]),
               delimiter=",", header="PC,fraction,cumulative", fmt=["%d", "%.5f", "%.5f"])
    with open(os.path.join(args.out, "loadings.csv"), "w") as f:
        f.write("feature," + ",".join(f"PC{i+1}" for i in range(ncomp)) + "\n")
        for fi, name in enumerate(names):
            f.write(name + "," + ",".join(f"{pca.components_[i, fi]:.4f}" for i in range(ncomp)) + "\n")

    # explained variance
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
        print("Only one PC available - skipping 2D plots.")
        return

    lab1 = f"PC1 ({evr[0]:.0%})"
    lab2 = f"PC2 ({evr[1]:.0%})"

    # scatter
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    sc = ax.scatter(Y[:, 0], Y[:, 1], c=np.arange(n_frames), s=8, cmap="viridis", alpha=0.7)
    fig.colorbar(sc, label="frame index (time)")
    ax.set_xlabel(lab1)
    ax.set_ylabel(lab2)
    ax.set_title("dihedral PCA, one point per conformer")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "pc1_pc2.png"), dpi=150)

    # 2D free-energy map
    H, xe, ye = np.histogram2d(Y[:, 0], Y[:, 1], bins=args.bins)
    P = H / H.sum()
    with np.errstate(divide="ignore"):
        G = -KB_KCAL * args.temp * np.log(P)
    G -= G[np.isfinite(G)].min()
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    im = ax.imshow(np.ma.masked_invalid(G).T, origin="lower", aspect="auto",
                   extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap="viridis_r")
    fig.colorbar(im, label="G (kcal/mol)")
    ax.set_xlabel(lab1)
    ax.set_ylabel(lab2)
    ax.set_title("free energy on PC1/PC2 (white = never visited)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "pc1_pc2_free_energy.png"), dpi=150)
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()