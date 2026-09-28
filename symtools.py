"""
Shared symmetry tools for cluster.py and dpca.py.

A symmetry operation = a relabelling of atoms that leaves the molecule unchanged
(graph automorphism). Applied to a conformer, it gives the torsion values you
would measure if the atoms had been labelled differently. Two conformers are the
same conformation if SOME operation makes them match.

The operations are stored by select_torsions.py in <prefix>_symops.txt:
one line per operation = the torsion atom quadruples after relabelling
(4 * n_torsions integers). The first line is the identity (= torsions.txt).
"""
import numpy as np


def load_ops(path, torsions):
    """(n_ops, n_tor, 4) array. Without a file: identity only."""
    tors = np.asarray(torsions)[:, :4]
    if path is None:
        return tors[None, :, :]
    raw = np.loadtxt(path, dtype=int, ndmin=2, comments="#")
    n_tor = len(tors)
    if raw.shape[1] != 4 * n_tor:
        raise SystemExit(f"{path}: expected {4 * n_tor} columns (4 x {n_tor} torsions), got {raw.shape[1]}")
    ops = raw.reshape(-1, n_tor, 4)
    if not np.array_equal(ops[0], tors):
        raise SystemExit(f"{path}: first line must be the identity (= torsions.txt). "
                         "Were both files written by the same select_torsions.py run?")
    return ops


def phi_all_ops(traj, ops):
    """Torsions (degrees, 0-360) of every frame under every operation: (n_ops, n_frames, n_tor)."""
    import mdtraj as md
    n_ops, n_tor, _ = ops.shape
    ang = md.compute_dihedrals(traj, ops.reshape(-1, 4))             # (frames, n_ops*n_tor)
    return (np.degrees(ang.astype(np.float64)).reshape(traj.n_frames, n_ops, n_tor).transpose(1, 0, 2)) % 360.0


def features(phi_deg):
    """cos/sin per torsion, scaled by 1/sqrt(n_tor). Works on (..., n_tor) arrays.
    Euclidean distance = Zivanovic 2020 dihedral distance (without S_i)."""
    r = np.radians(phi_deg)
    n = phi_deg.shape[-1]
    X = np.empty(phi_deg.shape[:-1] + (2 * n,))
    X[..., 0::2] = np.cos(r)
    X[..., 1::2] = np.sin(r)
    return X / np.sqrt(n)


def best_match_distances(Xops, max_frames=25000):
    """D[a, b] = min over operations g of || X(a) - X_g(b) ||   (n_frames x n_frames)."""
    from scipy.spatial.distance import cdist
    n = Xops.shape[1]
    if n > max_frames:
        raise SystemExit(f"{n} frames -> distance matrix too large; use --stride to subsample.")
    D = cdist(Xops[0], Xops[0])
    for g in range(1, len(Xops)):
        np.minimum(D, cdist(Xops[0], Xops[g]), out=D)
    np.minimum(D, D.T, out=D)                 # exact symmetry up to small geometric noise
    np.fill_diagonal(D, 0.0)
    return D


def align(Xops, phiops, ref_x):
    """Put every frame in the relabelling that brings it closest to ref_x.
    Returns aligned features (n_frames, 2n), aligned torsions (n_frames, n_tor), chosen op."""
    d = ((Xops - ref_x[None, None, :]) ** 2).sum(-1)                 # (n_ops, n_frames)
    best = d.argmin(0)
    idx = np.arange(Xops.shape[1])
    return Xops[best, idx], phiops[best, idx], best


def align_by_clusters(Xops, phiops, labels, medoids, ref):
    """Consistent alignment for plotting: the reference frame `ref` is used as is;
    each cluster's medoid is put in the relabelling closest to it, and every member
    of that cluster is aligned to its (relabelled) medoid. Noise frames are aligned
    to `ref`. This avoids one conformation being split between equivalent copies."""
    Xa = np.empty(Xops.shape[1:])
    pa = np.empty(phiops.shape[1:])
    ref_x = Xops[0, ref]
    for c, m in medoids.items():
        g = ((Xops[:, m] - ref_x) ** 2).sum(-1).argmin()
        mem = np.where(labels == c)[0]
        Xa[mem], pa[mem], _ = align(Xops[:, mem], phiops[:, mem], Xops[g, m])
    noise = np.where(labels < 0)[0]
    if len(noise):
        Xa[noise], pa[noise], _ = align(Xops[:, noise], phiops[:, noise], ref_x)
    return Xa, pa


def iterative_align(Xops, phiops, n_iter=3):
    """Alignment without a distance matrix: start from frame 0, align, move the
    reference to the frame closest to the aligned mean, repeat."""
    ref = Xops[0, 0]
    for _ in range(n_iter):
        Xa, pa, best = align(Xops, phiops, ref)
        ref = Xa[((Xa - Xa.mean(0)) ** 2).sum(1).argmin()]
    return align(Xops, phiops, ref)


def circ_mean(phi_deg):
    r = np.radians(phi_deg)
    return np.degrees(np.arctan2(np.sin(r).mean(0), np.cos(r).mean(0))) % 360.0