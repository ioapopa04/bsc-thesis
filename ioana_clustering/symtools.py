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

KB_KCAL = 0.0019872041


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


def torsion_pmf(phi_deg, nbins=36, smooth=1, temp=298.0):
    """1D circular free-energy profile of one torsion. Returns (G, centres), both
    length nbins, G in kcal/mol with its minimum subtracted. `smooth` is a
    half-width in bins of a circular boxcar applied to the counts."""
    edges = np.linspace(0.0, 360.0, nbins + 1)
    c, _ = np.histogram(phi_deg, edges)
    c = c + 1.0                                    # keep the log finite
    k = np.ones(2 * smooth + 1) / (2 * smooth + 1)
    cs = np.convolve(np.r_[c[-smooth:], c, c[:smooth]], k, "same")[smooth:-smooth]
    G = -KB_KCAL * temp * np.log(cs)
    return G - G.min(), (edges[:-1] + edges[1:]) / 2


def free_rotor_mask(phi_deg, min_barrier=2.0, nbins=36, smooth=1, temp=298.0):
    """Boolean mask, True for a *near-free rotor*: a torsion whose free-energy
    profile never rises `min_barrier` kcal/mol above its most stable value, so
    every angle is thermally accessible and the torsion defines no metastable
    state. Such torsions only add dimensions (and thermal smear) to the
    distance without separating conformations - see cluster.py --drop-free.

    phi_deg is (n_frames, n_tor); the test is per column."""
    phi = np.asarray(phi_deg, dtype=float)
    if phi.ndim == 1:
        phi = phi[:, None]
    out = np.zeros(phi.shape[1], dtype=bool)
    for j in range(phi.shape[1]):
        G, _ = torsion_pmf(phi[:, j], nbins, smooth, temp)
        out[j] = (G.max() - G.min()) < min_barrier
    return out


def _chord(a, b):
    """Chord length between two sets of angles (deg): 2|sin(dphi/2)|, in [0, 2].
    Equals the Euclidean distance between the (cos, sin) points."""
    from scipy.spatial.distance import cdist
    ca = np.stack([np.cos(np.radians(a)), np.sin(np.radians(a))], axis=-1)
    cb = np.stack([np.cos(np.radians(b)), np.sin(np.radians(b))], axis=-1)
    return cdist(ca, cb)


def pairwise_distances(phiops, metric="max", max_frames=25000):
    """Symmetry-aware distance between conformers, combining the torsions with
    `metric` AFTER choosing the best symmetry operation (the minimum is taken
    over operations on the combined distance, never per torsion - otherwise the
    two frames could be compared under different relabellings).

        d(a, b) = min_g  combine_j  chord(phi_j(a), phi_j(g b))

    metric "max": combine = max over torsions. A single flipped torsion is the
                  whole distance, so a real conformational change is not diluted
                  by thermal jitter in the other torsions. Resilience to noise
                  in one torsion is lower - good for sharp rotamer states.
    metric "rms": combine = sqrt(sum of squares). Each torsion contributes 1/n
                  of the squared distance, so one 180 deg flip is only worth
                  ~uniform small jitter everywhere (this is the old behaviour,
                  kept for comparison).

    Scale (both metrics, one torsion changing by dphi): d = 2|sin(dphi/2)|, so
    d = 0.5 <-> 29 deg, d = 1.0 <-> 60 deg, d = 2.0 <-> 180 deg. Distances are
    in [0, 2] regardless of how many torsions are used.

    phiops is (n_ops, n_frames, n_tor) from phi_all_ops.
    """
    n_ops, n, n_tor = phiops.shape
    if n > max_frames:
        raise SystemExit(f"{n} frames -> distance matrix too large; use --stride to subsample.")
    D = None
    for g in range(n_ops):
        acc = None
        for j in range(n_tor):
            c = _chord(phiops[0, :, j], phiops[g, :, j])
            if metric == "max":
                acc = c if acc is None else np.maximum(acc, c)
            elif metric == "rms":
                acc = c * c if acc is None else acc + c * c
            else:
                raise SystemExit(f"unknown metric {metric!r} (use 'max' or 'rms')")
        if metric == "rms":
            np.sqrt(acc, out=acc)
        D = acc if D is None else np.minimum(D, acc)
    np.minimum(D, D.T, out=D)                       # exact symmetry up to small geometric noise
    np.fill_diagonal(D, 0.0)
    return D


def best_match_distances(Xops, max_frames=25000):
    """Deprecated: kept so older callers still import. This is the old
    Euclidean-on-features (= `rms`) metric; prefer pairwise_distances(phiops, ...),
    which also offers the `max` metric and takes the symmetry minimum on the
    combined distance."""
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