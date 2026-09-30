#!/usr/bin/env python
"""
Cluster ligand conformers on torsions, respecting the molecule's symmetry.

Each selected torsion phi -> cos(phi), sin(phi). Symmetry: two conformers are
compared under every symmetry operation of the molecule (from torsions_symops.txt,
written by select_torsions.py) and the SMALLEST distance is used, taken on the
combined (all-torsion) distance so both frames are always in the same relabelling:
      d(A, B) = min_g combine_j chord(phi_j(A), phi_j(g B))
Operations act on all torsions jointly (e.g. a ring flip shifts both torsions of
R1-ring-R2 by 180 deg together), so syn/anti-type conformations stay distinct
while relabelled copies of the same conformation merge.
Distance scale (one torsion changing by dphi): d = 0.5 <-> 29 deg,
d = 1.0 <-> 60 deg, d = 2.0 <-> 180 deg; d is in [0, 2] whatever the torsion count.

  --metric max (default)  d = the most-changed torsion. A real single-torsion
                          flip is not diluted by thermal jitter in the rest, so
                          distinct rotamer states stay separated. Recommended.
  --metric rms            d = sqrt(sum of squares) - each torsion contributes
                          1/n of the squared distance. Kept for comparison: it
                          makes one 180 deg flip worth about as much as uniform
                          small jitter everywhere, which is why density methods
                          then see a single featureless blob.

Torsion selection: near-free rotors (--drop-free, default on) are removed. A
torsion whose 1D free-energy profile never rises --min-barrier kcal/mol above its
minimum defines no metastable state - it only adds a dimension and thermal smear,
and its contribution to the distance is pure noise for clustering. --torsions-keep
/ --torsions-drop give explicit control.

Methods (on the precomputed best-match distances):
  hdbscan (default)  density-based: finds the dense core of each basin (number of
                     clusters not preset); then lower-density frames are assigned to the basin they are
                     connected to by steps < --assign-radius. Only frames cut off by
                     an empty gap (e.g. isolated barrier crossings) stay noise (-1).
                     NOTE: allow_single_cluster is OFF by default - with it on,
                     a weak/no density contrast is silently reported as one cluster.
  daura              Daura/GROMOS algorithm with a distance --cutoff.

Example:
  python cluster.py --traj whole_state0_prod1.nc --top ../build/built.pdb \
      --torsions torsions.txt --symops torsions_symops.txt --select "resname PAR" --out clusters

Outputs (in --out): clusters.csv, labels.npy, aligned_dihedrals.npy, clusters_pc.png, clusters_torsions.png,
cluster_<c>.pdb (representative = cluster medoid, ligand only if --select).
clusters.csv lists the torsions actually used (torsions column) and gives each
population with a block-bootstrap error (population_err), which is what to read
for convergence - a fine-bin JSD can look small while probability moves between
whole wells. Torsion means are circular means after aligning every member to the
cluster representative (so symmetry-related copies don't average out).
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


def run_hdbscan(D, mcs, ms, allow_single=False):
    try:
        from sklearn.cluster import HDBSCAN
    except ImportError:
        raise SystemExit("HDBSCAN needs scikit-learn >= 1.3 (pip install -U scikit-learn)")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        return HDBSCAN(min_cluster_size=mcs, min_samples=ms, metric="precomputed",
                       allow_single_cluster=allow_single).fit_predict(D)


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


def population_errors(labels, n_blocks, n_boot=200, seed=0):
    """Block-bootstrap standard error of each cluster's population. Resampling
    contiguous blocks (not single frames) keeps the trajectory's autocorrelation
    from making the error look optimistically small. Returns {cluster: (pop, err)};
    err is nan when there are too few blocks to bootstrap."""
    ids = [c for c in np.unique(labels) if c >= 0]
    n = len(labels)
    if n_blocks < 2 or not ids:
        return {c: (float((labels == c).mean()), float("nan")) for c in ids}
    edges = np.linspace(0, n, n_blocks + 1).astype(int)
    rng = np.random.default_rng(seed)
    draws = {c: [] for c in ids}
    for _ in range(n_boot):
        picks = rng.integers(0, n_blocks, n_blocks)
        idx = np.concatenate([np.arange(edges[b], edges[b + 1]) for b in picks])
        lab = labels[idx]
        for c in ids:
            draws[c].append((lab == c).mean())
    return {c: (float((labels == c).mean()), float(np.std(draws[c]))) for c in ids}


def relabel_by_size(labels):
    ids = sorted([c for c in np.unique(labels) if c >= 0], key=lambda c: -(labels == c).sum())
    new = np.full_like(labels, -1)
    for i, c in enumerate(ids):
        new[labels == c] = i
    return new


def cluster_pipeline(D, args, n_frames):
    """The clustering itself, factored out so --stability can re-run it on subsets.
    Returns (labels, note) where note is a human-readable line about the method."""
    if args.method == "hdbscan":
        mcs = int(round(args.min_cluster_size * n_frames)) if args.min_cluster_size < 1 else int(args.min_cluster_size)
        mcs = min(max(mcs, 10), n_frames)
        ms = int(round(args.min_samples * n_frames)) if args.min_samples < 1 else int(args.min_samples)
        ms = min(max(ms, 5), mcs)
        labels = run_hdbscan(D, mcs, ms, args.allow_single_cluster)
        n_raw = len([c for c in np.unique(labels) if c >= 0])
        labels = merge_close_clusters(D, labels, args.merge_gap)
        n_new = len([c for c in np.unique(labels) if c >= 0])
        n_core = int((labels >= 0).sum())
        labels, n_ass = assign_to_cores(D, labels, args.assign_radius)
        note = f"HDBSCAN min_cluster_size = {mcs}, min_samples = {ms}; {n_raw} cluster(s)"
        if args.merge_gap > 0:
            note += f" -> {n_new} after merging clusters closer than {args.merge_gap}"
        note += (f"\n  {n_core} frames in dense cores; {n_ass} lower-density frames assigned to their "
                 f"basin (path of frames closer than {args.assign_radius})")
    else:
        labels = run_daura(D, args.cutoff, args.max_clusters, args.coverage)
        note = f"Daura cutoff d = {args.cutoff}"
    return relabel_by_size(labels), note


def cluster_stability(phiops, base_labels, args, n_resamples=5, frac=0.5, seed=0):
    """Re-cluster independent windows of the trajectory with the SAME settings and
    ask, for each cluster, which fraction of its frames in each window still land
    in one cluster. A value near 1 means the cluster reappears from a subset of the
    data; a low value means its membership (not just its population) is not
    reproducible. Returns {cluster: fraction}."""
    ids = [c for c in np.unique(base_labels) if c >= 0]
    n = phiops.shape[1]
    w = int(round(frac * n))
    if len(ids) == 0 or w < 50 or n_resamples < 1:
        return {c: float("nan") for c in ids}
    rng = np.random.default_rng(seed)
    starts = [0] if w >= n else sorted(rng.choice(n - w, size=min(n_resamples, n - w), replace=False).tolist())
    hits = {c: [] for c in ids}
    for s in starts:
        idx = np.arange(s, s + w)
        sub = st.pairwise_distances(phiops[:, idx], args.metric)
        sub_labels, _ = cluster_pipeline(sub, args, w)
        for c in ids:
            here = base_labels[idx] == c
            if not here.any():
                continue
            # best overlap of this cluster's frames in the window with any sub-cluster
            best = 0.0
            for c2 in np.unique(sub_labels):
                if c2 < 0:
                    continue
                best = max(best, float((sub_labels[here] == c2).mean()))
            hits[c].append(best)
    return {c: (float(np.mean(hits[c])) if hits[c] else float("nan")) for c in ids}


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
    ap.add_argument("--metric", choices=["max", "rms"], default="max",
                    help="how the torsions are combined into one distance: 'max' (default) = the most-changed "
                         "torsion, keeps single-torsion flips visible; 'rms' = sqrt(sum of squares), the old "
                         "behaviour that dilutes one flip into uniform jitter")
    ap.add_argument("--drop-free", action=argparse.BooleanOptionalAction, default=True,
                    help="drop near-free rotors before clustering (default on); see --min-barrier")
    ap.add_argument("--min-barrier", type=float, default=2.0,
                    help="--drop-free: a torsion is near-free if its 1D free-energy profile never rises this "
                         "many kcal/mol above its minimum (default 2.0)")
    ap.add_argument("--torsions-keep", default=None,
                    help="comma-separated torsion indices to use, e.g. '2,3,4' (overrides --drop-free/--torsions-drop)")
    ap.add_argument("--torsions-drop", default=None,
                    help="comma-separated torsion indices to additionally drop, e.g. '0,8'")
    ap.add_argument("--allow-single-cluster", action="store_true",
                    help="hdbscan: allow a single cluster (default off - with it on, absent density contrast "
                         "is reported as one cluster instead of as noise)")
    ap.add_argument("--blocks", type=int, default=5,
                    help="contiguous blocks for the block-bootstrap population error (0 = off, default 5)")
    ap.add_argument("--stability", type=int, default=5,
                    help="re-cluster this many windows of --stability-frac of the trajectory and report how "
                         "often each cluster reappears (0 = off, default 5)")
    ap.add_argument("--stability-frac", type=float, default=0.5,
                    help="window size as a fraction of the frames for --stability (default 0.5)")
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
    ap.add_argument("--cutoff", type=float, default=0.5,
                    help="daura: distance cutoff in chord units (default 0.5; 0.5 <-> 29 deg in the "
                         "most-changed torsion under --metric max, 2.0 = full 180 deg)")
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
    n_ops, n_frames, n_tor_all = phiops.shape
    print(f"{n_frames} frames, {n_tor_all} torsions in file, {n_ops} symmetry operation(s)")

    # ---- torsion selection: drop near-free rotors, honour explicit keep/drop ----
    gmax = np.array([st.torsion_pmf(phiops[0, :, j], temp=args.temp)[0].max() for j in range(n_tor_all)])
    free = st.free_rotor_mask(phiops[0], args.min_barrier, temp=args.temp)
    dropped = np.zeros(n_tor_all, dtype=bool)
    if args.torsions_drop:
        dropped[[int(x) for x in args.torsions_drop.split(",")]] = True
    if args.drop_free:
        dropped |= free
    if args.torsions_keep is not None:
        keep = np.array([int(x) for x in args.torsions_keep.split(",")])
        print(f"torsion indices kept explicitly: {list(keep)}")
        if dropped[keep].any():
            print(f"  NOTE: --torsions-keep overrides --drop-free for indices {list(keep[dropped[keep]])}")
    else:
        keep = np.array([j for j in range(n_tor_all) if not dropped[j]])
    if len(keep) == 0:
        raise SystemExit("No torsions left after selection - lower --min-barrier or drop fewer.")

    print(f"{'tor':>4} {'Gmax':>6} {'kept':>5}   (Gmax = highest point of the 1D PMF, kcal/mol;"
          f" free rotor if < {args.min_barrier})")
    for j in range(n_tor_all):
        why = "free rotor" if free[j] else ""
        print(f"{j:>4} {gmax[j]:>6.2f} {'yes' if j in keep else 'no':>5}   {why}")
    if len(keep) < 2:
        print("WARNING: fewer than 2 torsions - clusters cannot separate conformations.")

    phiops = phiops[:, :, keep]
    tors = tors[keep]
    Xops = st.features(phiops)
    n_tor = len(keep)
    print(f"using {n_tor} torsion(s) {list(keep)}, metric '{args.metric}', method {args.method}")

    D = st.pairwise_distances(phiops, args.metric)

    labels, note = cluster_pipeline(D, args, n_frames)
    print(note)
    if len([c for c in np.unique(labels) if c >= 0]) == 1:
        print("  WARNING: one dense cluster only. Either the ensemble really is a single basin, or the "
              "torsion set/metric does not separate conformations (check the Gmax table above; try "
              "--metric max, --drop-free, or --torsions-keep on the barrier-separated torsions).")
    np.save(os.path.join(args.out, "labels.npy"), labels)

    ids = [c for c in np.unique(labels) if c >= 0]
    if not ids:
        raise SystemExit("No clusters found (all noise) - lower --min-cluster-size or raise --cutoff.")
    noise = (labels == -1).mean()
    kT = KB_KCAL * args.temp
    counts = {c: int((labels == c).sum()) for c in ids}
    pmax = max(counts.values())
    poperr = population_errors(labels, args.blocks)
    if args.stability > 0:
        print(f"checking cluster stability ({args.stability} windows of {args.stability_frac:.0%} of the frames)...")
        stab = cluster_stability(phiops, labels, args, args.stability, args.stability_frac)
    else:
        stab = {c: float("nan") for c in ids}

    rows = []
    for c in ids:
        mem = np.where(labels == c)[0]
        rep = int(mem[D[np.ix_(mem, mem)].sum(1).argmin()])        # medoid
        _, phi_al, _ = st.align(Xops[:, mem], phiops[:, mem], Xops[0, rep])
        rows.append(dict(c=c, n=counts[c], p=counts[c] / n_frames,
                         dG=-kT * np.log(counts[c] / pmax) + 0.0, rep=rep,
                         mean=st.circ_mean(phi_al), rep_phi=phiops[0, rep]))

    # column names carry the real torsion indices, so a later run with a different
    # selection cannot be confused with this one
    cols = [f"mean_phi{j}" for j in keep] + [f"rep_phi{j}" for j in keep]
    with open(os.path.join(args.out, "clusters.csv"), "w") as f:
        f.write(f"cluster,frames,population,population_err,stability,dG_kcal_mol,rep_frame,torsions,"
                + ",".join(cols) + "\n")
        for r in rows:
            err = poperr[r["c"]][1]
            errs = f"{err:.4f}" if err == err else ""        # nan -> empty
            sv = stab[r["c"]]
            svs = f"{sv:.3f}" if sv == sv else ""
            f.write(f"{r['c']},{r['n']},{r['p']:.4f},{errs},{svs},{r['dG']:.3f},{r['rep']},"
                    + " ".join(str(int(j)) for j in keep) + ","
                    + ",".join(f"{v:.1f}" for v in r["mean"]) + ","
                    + ",".join(f"{v:.1f}" for v in r["rep_phi"]) + "\n")
        f.write(f"noise,{int((labels == -1).sum())},{noise:.4f},,,,,\n")

    print(f"\n{'cluster':>7s} {'frames':>7s} {'pop':>7s} {'+-err':>6s} {'stab':>5s} {'dG(kcal/mol)':>13s} "
          f"{'rep':>6s}  mean torsions {list(int(j) for j in keep)} (deg)")
    for r in rows:
        err = poperr[r["c"]][1]
        errs = f"{err:.1%}" if err == err else "-"
        sv = stab[r["c"]]
        svs = f"{sv:.0%}" if sv == sv else "-"
        print(f"{r['c']:>7d} {r['n']:>7d} {r['p']:>7.1%} {errs:>6s} {svs:>5s} {r['dG']:>13.2f} {r['rep']:>6d}  "
              f"{[round(float(v), 1) for v in r['mean']]}")
    print(f"{'noise':>7s} {int((labels == -1).sum()):>7d} {noise:>7.1%}")
    if args.blocks >= 2 and rows:
        worst = max(poperr[r["c"]][1] for r in rows)
        if worst > 0.05:
            print(f"NOTE: largest block-bootstrap population error is {worst:.1%}; populations above a few "
                  "percent are not resolved by this trajectory length.")
    if args.stability > 0 and rows:
        bad = [r["c"] for r in rows if stab[r["c"]] == stab[r["c"]] and stab[r["c"]] < 0.7]
        if bad:
            print(f"NOTE: cluster(s) {bad} reappear in fewer than 70% of the windows - their membership is "
                  "not reproducible from part of the trajectory. Treat them as suggestive, not as states.")
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
        axes[j, 0].set_ylabel(f"torsion {int(keep[j])}\ncount")
    axes[-1, 0].set_xlabel("dihedral (deg), symmetry-aligned")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "clusters_torsions.png"), dpi=150)
    print(f"Results written to {args.out}/")


if __name__ == "__main__":
    main()