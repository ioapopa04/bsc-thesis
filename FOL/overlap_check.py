#!/usr/bin/env python
"""Quantify cluster overlap for Ioana's vs Shu-Yu's clustering on the 100 ns run.

Overlap measure: for each frame, the fraction of its k nearest neighbours (in the
distance the clustering itself used) that carry the SAME cluster label. 1.0 = a
compact, isolated cluster; 0.33 with 3 clusters = indistinguishable from random.
Also a silhouette score on a subsample of labelled frames.
"""
import sys
import numpy as np
import mdtraj as md

sys.path.insert(0, "/scratch/ipopa/bsc-thesis/ioana")
import symtools as st

ROOT = "/scratch/ipopa/FOLn-REST2-50ns"
HERE = f"{ROOT}/REST2-run1"
TOP = f"{ROOT}/build/built.solute.pdb"
TRAJ = f"{HERE}/solute_state0_prod1.nc"
STRIDE = 5
K = 15


def knn_purity(D, labels, k=K):
    idx = np.argpartition(D, kth=k, axis=1)[:, 1:k + 1]        # k nearest, excluding self
    return (labels[idx] == labels[:, None]).mean(1)


traj = md.load(TRAJ, top=TOP, stride=STRIDE)
n = traj.n_frames
print(f"{n} frames\n")

# ---------------- Ioana ----------------
tors = np.loadtxt(f"{HERE}/ioana_torsions.txt", dtype=int, ndmin=2, comments="#")
ops = st.load_ops(f"{HERE}/ioana_torsions_symops.txt", tors)
phiops = st.phi_all_ops(traj, ops)
free = st.free_rotor_mask(phiops[0], 2.0)
keep = [j for j in range(len(tors)) if not free[j]]
D = st.pairwise_distances(phiops[:, :, keep], "max")
lab_i = np.load(f"{ROOT}/fol-clustering-ioana/clusters/labels.npy")

pur_i = knn_purity(D, lab_i)
print("== Ioana (max metric, 6 torsions, 4 clusters + noise) ==")
for c in sorted(set(lab_i.tolist())):
    m = lab_i == c
    ids = f"cluster {c}" if c >= 0 else "noise"
    print(f"  {ids:>9}: n={m.sum():5d}  kNN purity {pur_i[m].mean():.3f}")
print(f"  overall (labelled frames): {pur_i[lab_i >= 0].mean():.3f}")

# ---------------- Shu-Yu ----------------
from md_tools.analysis import t_hdbscan
names, quads = [], []
for path in ["torsions.txt", "torsions_rigid.txt", "torsions_polarH.txt"]:
    label = path.split(".")[0]
    for line in open(f"{HERE}/{path}"):
        if line.strip().startswith("#") or not line.strip():
            continue
        quads.append(tuple(int(x) for x in line.split("#")[0].split()[:4]))
        names.append(f"{label}_{len(names)}")
th = md.compute_dihedrals(traj, quads)
X = np.concatenate([np.cos(th), np.sin(th)], axis=1)          # MD-tools embedding (uniform)
from sklearn.neighbors import NearestNeighbors
nbr_s = NearestNeighbors(n_neighbors=K + 1).fit(X).kneighbors(X, return_distance=False)[:, 1:K + 1]
lab_s = np.asarray(t_hdbscan(th, names=names, allow_single_cluster=False).fit().labels_)

pur_s = (lab_s[nbr_s] == lab_s[:, None]).mean(1)
print("\n== Shu-Yu (sum metric, 12 torsions, 2 clusters + noise) ==")
for c in sorted(set(lab_s.tolist())):
    m = lab_s == c
    ids = f"cluster {c}" if c >= 0 else f"noise({c})"
    print(f"  {ids:>9}: n={m.sum():5d}  kNN purity {pur_s[m].mean():.3f}")
print(f"  overall (labelled frames): {pur_s[lab_s >= 0].mean():.3f}")


def confusion(labels, neighbour_idx):
    """kNN edge confusion: C[i,j] = fraction of kNN edges leaving label i that land in label j."""
    src = labels[:, None] * np.ones_like(neighbour_idx)
    dst = labels[neighbour_idx]
    keys = sorted(set(labels.tolist()))
    pos = {k: i for i, k in enumerate(keys)}
    C = np.zeros((len(keys), len(keys)))
    for s, d in zip(src.ravel(), dst.ravel()):
        C[pos[int(s)], pos[int(d)]] += 1
    C /= C.sum(1, keepdims=True)
    print("      " + "".join(f"{('c'+str(k)) if k>=0 else 'noise':>9}" for k in keys))
    for k in keys:
        print(f"  {('c'+str(k)) if k>=0 else 'noise':>5} "
              + "".join(f"{C[pos[k], pos[j]]:>9.2f}" for j in keys))


print("\nIoana kNN confusion (rows=source label):")
confusion(lab_i, np.argpartition(D, kth=K, axis=1)[:, 1:K + 1])
print("\nShu-Yu kNN confusion (rows=source label):")
confusion(lab_s, nbr_s)

# ---------------- silhouette on labelled frames ----------------
from sklearn.metrics import silhouette_score
rng = np.random.default_rng(0)
for name, L, neighbours in (("Ioana", lab_i, None), ("Shu-Yu", lab_s, nbr_s)):
    m = np.flatnonzero(L >= 0)
    if len(m) > 2000:
        m = rng.choice(m, 2000, replace=False)
    Dsub = (D[np.ix_(m, m)] if name == "Ioana"
            else np.sqrt(((X[m][:, None, :] - X[m][None, :, :]) ** 2).sum(-1)))
    np.fill_diagonal(Dsub, 0.0)
    sil = silhouette_score(Dsub, L[m], metric="precomputed")
    print(f"{name}: silhouette on {len(m)} labelled frames = {sil:.3f}")
