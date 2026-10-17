#!/usr/bin/env python
"""
Do REST2 replicates started from DIFFERENT conformers (pick_starts.py) give the same ensemble?

All replicates of one molecule are pooled and clustered ONCE with Shu-Yu's symmetry-first
t_hdbscan (same torsions and settings as shuyu_clustering/cluster_prune.py), so every frame of
every replicate carries a label from one common partition. Clustering each replicate separately
would give labels that cannot be matched across runs.

  --rep rep0:<traj>[:<top>]   (repeatable) state-0 solute trajectory of each replicate; its topology
                              must have the --top atom order (checked: elements and names)
  pooled frames > --max-frames  ->  the partition is fitted on a weighted resample of that many
                              frames (MD-tools resampling=N); EVERY frame is then classified with the
                              same distance and vote. Each replicate contributes equally (same stride,
                              same length).

Outputs (--out):
  populations.csv / .png   cluster mass per replicate (+ unassigned), error = block standard error
                           (--blocks equal blocks per replicate), and the pooled mass
  timeline.png             the cluster label of every frame against time, one row per replicate
  convergence.png          cumulative cluster mass against time, one panel per cluster
  torsions.png             per-torsion distributions per replicate, in the common (aligned) labelling
  agreement.json           Jensen-Shannon divergence between replicates (label distributions,
                           unassigned included, and each aligned torsion), max |mass difference|
                           per cluster vs its block error, start cluster and time spent there first
"""
import argparse
import csv
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
from rdkit import Chem

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shuyu_clustering"))
from cluster_prune import load_torsions  # noqa: E402
from md_tools.analysis import t_hdbscan, TorsionDefinitions, verify_atom_mapping, group_name  # noqa: E402

# categorical slots 1-8 (dataviz reference palette, light mode), assigned in fixed order:
# replicate k always gets slot k+1
REP_COLOURS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
UNASSIGNED_GREY = "#a3a29c"
INK, INK2 = "#0b0b0b", "#52514e"
plt.rcParams.update({"axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": "#e6e5e0", "grid.linewidth": 0.6,
                     "axes.axisbelow": True, "font.size": 9})


def lab(c):
    c = int(c)
    return "unassigned" if c < 0 else group_name(c)


def js(p, q):
    p, q = np.asarray(p, float) + 1e-12, np.asarray(q, float) + 1e-12
    p, q = p / p.sum(), q / q.sum()
    m = 0.5 * (p + q)
    return float(0.5 * np.sum(p * np.log2(p / m)) + 0.5 * np.sum(q * np.log2(q / m)))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rep", action="append", required=True, help="label:traj[:top]")
    ap.add_argument("--top", required=True, help="reference topology with bonds (rep0 built.solute.pdb)")
    ap.add_argument("--sdf", required=True, help="ligand SDF in that atom order (bond orders)")
    ap.add_argument("--torsion", action="append", default=None)
    ap.add_argument("--torsions-file", action="append", default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--no-single-cluster", action="store_true")
    ap.add_argument("--max-frames", type=int, default=10000,
                    help="fit on a resample of this many pooled frames if there are more")
    ap.add_argument("--blocks", type=int, default=5, help="blocks per replicate for the error bars")
    ap.add_argument("--frame-ps", type=float, default=None,
                    help="time between saved frames before striding (default: from the trajectory)")
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    P = lambda name: os.path.join(out, name)
    title = f"{args.name}: " if args.name else ""

    torsions = load_torsions(args.torsions_file, args.torsion)
    names, quads = list(torsions), list(torsions.values())
    reftop = md.load(args.top)
    ref_atoms = [(a.element.symbol, a.name) for a in reftop.topology.atoms]

    # ------------------------------------------------------------- load and pool
    reps, trajs = [], []
    for spec in args.rep:
        parts = spec.split(":")
        label, path = parts[0], parts[1]
        top = parts[2] if len(parts) > 2 else args.top
        if [(a.element.symbol, a.name) for a in md.load(top).topology.atoms] != ref_atoms:
            sys.exit(f"{label}: {top} does not have the reference atom order")
        t = md.load(path, top=top, stride=args.stride)
        t.topology = reftop.topology                      # one topology object, so md.join works
        reps.append(label)
        trajs.append(t)
        print(f"{label}: {t.n_frames} frames from {path}")
    n_per = [t.n_frames for t in trajs]
    if len(set(n_per)) != 1:
        print(f"WARNING: replicates have different frame counts {n_per}; pooled masses weight them unequally")
    pooled = md.join(trajs, check_topology=False)
    tors = md.compute_dihedrals(pooled, quads)
    rep_of = np.repeat(np.arange(len(reps)), n_per)
    # mdtraj's timestep is the spacing of the LOADED (strided) frames, from the file's time stamps
    dt_ns = args.frame_ps * args.stride / 1000.0 if args.frame_ps else float(trajs[0].timestep) / 1000.0
    print(f"pooled: {tors.shape[0]} frames, {tors.shape[1]} torsions, {dt_ns:g} ns between used frames")

    # ------------------------------------------------------------- one common partition
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    definitions = TorsionDefinitions(quads=tuple(quads), names=tuple(names), indexing="zero-based", units="radian")
    bonds = [(b[0].index, b[1].index) for b in pooled.topology.bonds]
    mapping, why = verify_atom_mapping(mol, bonds, [a.element.symbol for a in pooled.topology.atoms])
    resample = int(args.max_frames) if tors.shape[0] > args.max_frames else False
    fit = t_hdbscan(tors, molecule=mol, torsion_definitions=definitions, coordinates=pooled,
                    atom_mapping=mapping, allow_single_cluster=not args.no_single_cluster,
                    resampling=resample, max_precomputed_frames=max(args.max_frames, 2)).fit()
    labels = np.asarray(fit.labels_)
    K = fit.n_clusters_
    print(f"pooled fit: {K} clusters, resampling={resample}, "
          f"masses { {group_name(c): round(float(fit.cluster_mass_[c]), 4) for c in range(K)} }, "
          f"unassigned {fit.unassigned_mass_:.4f}")
    cats = list(range(K)) + [-1]                         # -1 = unassigned (noise + ambiguous)
    lab_u = np.where(labels < 0, -1, labels)

    # ------------------------------------------------------------- per-replicate masses
    def masses(sel):
        n = max(1, int(sel.sum()))
        return np.array([float((lab_u[sel] == c).sum()) / n for c in cats])

    M = np.array([masses(rep_of == r) for r in range(len(reps))])          # reps x cats
    SE = np.zeros_like(M)
    for r in range(len(reps)):
        idx = np.flatnonzero(rep_of == r)
        blocks = np.array_split(idx, args.blocks)
        B = np.array([[float((lab_u[b] == c).mean()) for c in cats] for b in blocks])
        SE[r] = B.std(axis=0, ddof=1) / np.sqrt(len(blocks))
    pooled_m = masses(np.ones_like(rep_of, bool))
    with open(P("populations.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cluster"] + [f"{r}_mass" for r in reps] + [f"{r}_block_se" for r in reps] + ["pooled_mass"])
        for j, c in enumerate(cats):
            w.writerow([lab(c)] + [f"{x:.4f}" for x in M[:, j]] + [f"{x:.4f}" for x in SE[:, j]]
                       + [f"{pooled_m[j]:.4f}"])
    for j, c in enumerate(cats):
        print(f"  {lab(c):10s} " + "  ".join(f"{r} {100 * M[i, j]:5.1f}±{100 * SE[i, j]:4.1f}%"
                                             for i, r in enumerate(reps)))

    # ------------------------------------------------------------- aligned torsions (figures only)
    al = fit.align_all(per_cluster=True)
    aligned = np.asarray(al["selected"])
    edges = np.linspace(-np.pi, np.pi, 37)
    H = np.array([[np.histogram(aligned[rep_of == r, j], bins=edges)[0] for j in range(len(names))]
                  for r in range(len(reps))])                              # reps x torsions x bins

    # ------------------------------------------------------------- agreement numbers
    pair = {}
    for a in range(len(reps)):
        for b in range(a + 1, len(reps)):
            pair[f"{reps[a]}-{reps[b]}"] = {
                "js_labels_bits": js(M[a], M[b]),
                "js_torsions_bits": {n: js(H[a, j], H[b, j]) for j, n in enumerate(names)}}
    per_cluster = {}
    for j, c in enumerate(cats):
        spread = float(M[:, j].max() - M[:, j].min())
        combined_se = float(np.sqrt((SE[:, j] ** 2).max() * 2))
        per_cluster[lab(c)] = {"masses": dict(zip(reps, map(float, M[:, j]))),
                               "block_se": dict(zip(reps, map(float, SE[:, j]))),
                               "max_minus_min": spread,
                               "spread_over_2se": spread / (2 * combined_se) if combined_se > 0 else None}
    start = {}
    for r, label in enumerate(reps):
        l = lab_u[rep_of == r]
        first = int(l[0])
        stay = int(np.argmax(l != first)) if np.any(l != first) else len(l)
        start[label] = {"first_frame_cluster": lab(first), "frames_before_first_change": stay,
                        "ns_before_first_change": stay * dt_ns,
                        "fraction_in_start_cluster_first_10pct":
                            float((l[: max(1, len(l) // 10)] == first).mean())}
    with open(P("agreement.json"), "w") as f:
        json.dump({"name": args.name, "replicates": reps, "frames_per_replicate": n_per,
                   "stride": args.stride, "ns_per_used_frame": dt_ns,
                   "fit": {"n_clusters": int(K), "resampling": resample,
                           "allow_single_cluster": not args.no_single_cluster,
                           "pooled_masses": dict(zip(map(lab, cats), map(float, pooled_m)))},
                   "per_cluster": per_cluster, "pairs": pair, "start": start,
                   "torsions": {n: list(q) for n, q in torsions.items()}}, f, indent=2)
    for k, v in pair.items():
        print(f"  JS(labels) {k}: {v['js_labels_bits']:.4f} bits")

    cols = [REP_COLOURS[r % len(REP_COLOURS)] for r in range(len(reps))]
    cat_names = [lab(c) for c in cats]

    # populations: grouped bars, one bar per replicate, pooled as a black tick
    fig, axis = plt.subplots(figsize=(max(5.5, 1.1 * len(cats) * len(reps) / 2 + 2), 3.6), layout="constrained")
    width = 0.8 / len(reps)
    x = np.arange(len(cats))
    for r, label in enumerate(reps):
        axis.bar(x - 0.4 + width * (r + 0.5), 100 * M[r], width=width * 0.92, color=cols[r],
                 yerr=100 * SE[r], error_kw={"elinewidth": 1, "capsize": 2, "ecolor": INK2}, label=label)
    axis.scatter(x, 100 * pooled_m, marker="_", s=600, color=INK, zorder=3, label="pooled")
    axis.set_xticks(x, cat_names)
    axis.set_ylabel("population (%)")
    axis.set_xlabel("cluster (one partition fitted on all replicates)")
    axis.set_ylim(0, 100 * float(max((M + SE).max(), pooled_m.max())) * 1.2)   # headroom for the legend
    axis.legend(ncol=len(reps) + 1, fontsize=8, frameon=False, loc="upper center")
    axis.set_title(f"{title}cluster populations per replicate (error bars: block SE, {args.blocks} blocks)", fontsize=10)
    fig.savefig(P("populations.png"), dpi=170)
    plt.close(fig)

    # timeline: label of every frame vs time
    fig, axes = plt.subplots(len(reps), 1, figsize=(9, 0.9 * len(reps) + 1.0), sharex=True,
                             layout="constrained", squeeze=False)
    for r, (axis, label) in enumerate(zip(axes[:, 0], reps)):
        l = lab_u[rep_of == r]
        tt = np.arange(len(l)) * dt_ns
        y = np.where(l < 0, K, l)
        axis.scatter(tt, y, s=1.5, color=np.where(l < 0, UNASSIGNED_GREY, cols[r]), linewidths=0)
        axis.set_yticks(range(K + 1), [group_name(c) for c in range(K)] + ["un."], fontsize=7)
        axis.set_ylim(K + 0.6, -0.6)
        axis.set_ylabel(label, rotation=0, ha="right", va="center")
        axis.grid(axis="x")
    axes[-1, 0].set_xlabel("time (ns)")
    fig.suptitle(f"{title}cluster of every frame (state 0); grey = unassigned", fontsize=10)
    fig.savefig(P("timeline.png"), dpi=170)
    plt.close(fig)

    # convergence: cumulative mass vs time, one panel per cluster (+ unassigned)
    ncol = min(4, len(cats))
    nrow = int(np.ceil(len(cats) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.5 * nrow), layout="constrained", squeeze=False)
    for extra in axes.flat[len(cats):]:
        extra.axis("off")
    for j, (axis, c) in enumerate(zip(axes.flat, cats)):
        for r, label in enumerate(reps):
            l = lab_u[rep_of == r]
            cum = np.cumsum(l == c) / np.arange(1, len(l) + 1)
            axis.plot(np.arange(len(l)) * dt_ns, 100 * cum, lw=2, color=cols[r], label=label)
        axis.axhline(100 * pooled_m[j], color=INK2, lw=1, ls=":")
        axis.set_title(lab(c), fontsize=9)
        axis.set_xlabel("time (ns)")
        if j % ncol == 0:
            axis.set_ylabel("cumulative population (%)")
        if j == 0:
            axis.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{title}running cluster populations per replicate (dotted: pooled)", fontsize=10)
    fig.savefig(P("convergence.png"), dpi=170)
    plt.close(fig)

    # torsions: per-torsion distributions per replicate (aligned labelling)
    ncol = min(4, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.5 * nrow), layout="constrained", squeeze=False)
    for extra in axes.flat[len(names):]:
        extra.axis("off")
    centres = np.rad2deg(0.5 * (edges[1:] + edges[:-1]))
    for j, (axis, n) in enumerate(zip(axes.flat, names)):
        for r, label in enumerate(reps):
            d = H[r, j] / max(1, H[r, j].sum()) / np.diff(np.rad2deg(edges))
            axis.step(centres, d, where="mid", lw=1.6, color=cols[r], label=label)
        axis.set_xlim(-180, 180)
        axis.set_xticks([-180, -90, 0, 90, 180])
        axis.set_title(n, fontsize=9)
        if j == 0:
            axis.legend(fontsize=7, frameon=False)
    fig.suptitle(f"{title}torsion distributions per replicate (symmetry-aligned, deg)", fontsize=10)
    fig.savefig(P("torsions.png"), dpi=170)
    plt.close(fig)
    print("written to", out)


if __name__ == "__main__":
    main()
