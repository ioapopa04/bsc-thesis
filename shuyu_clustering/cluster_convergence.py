#!/usr/bin/env python
"""
Convergence of the Shu-Yu clustering with simulated time.

Re-runs exactly the clustering of cluster_prune.py (Shu-Yu's symmetry-first t_hdbscan: symmetry
built into the distance, no merge step) on CUMULATIVE windows 0-5, 0-10, 0-15 ... ns of the
state-0 trajectory and plots the number of clusters against the window end. The run is called
converged once --plateau consecutive windows (default 3) give the same number of clusters.

Windows are cut on the trajectory's own time stamps. Several --traj segments (e.g. REST2-run1
then REST2-run1-ext1) are concatenated in the order given; a segment whose clock restarts is
shifted to follow the previous one.

Example (4ibi, 60 ns = run1 + ext1 already joined):
  python cluster_convergence.py \
      --traj analysis-0-60ns/solute_state0_0-60ns.nc --top build/built.solute.pdb \
      --sdf build/L01.sdf --torsions-file analysis-0-60ns/torsions.txt \
      --torsions-file analysis-0-60ns/torsions_rigid.txt \
      --torsions-file analysis-0-60ns/torsions_polarH.txt \
      --stride 5 --no-single-cluster --step-ns 5 --name 4ibi --rest2-config REST2.config \
      --out analysis-0-60ns/4ibi-convergence

Outputs in --out/: convergence.csv, convergence.png, torsion_cv.json
"""
import argparse
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_prune import build_fit, load_torsions, read_rest2_config, simulated_ns  # noqa: E402


def load_segments(paths, top, stride):
    """Concatenate trajectory segments on one continuous clock (ps)."""
    parts, offset = [], 0.0
    for path in paths:
        t = md.load(path, top=top, stride=stride)
        time = np.asarray(t.time, dtype=float)
        if parts and time[0] <= parts[-1].time[-1]:
            time = time - time[0] + parts[-1].time[-1] + (time[1] - time[0] if len(time) > 1 else 0)
        t.time = time
        parts.append(t)
        print(f"  {path}: {t.n_frames} frames, {time[0] / 1000:.3f}-{time[-1] / 1000:.3f} ns")
    return md.join(parts, check_topology=True) if len(parts) > 1 else parts[0]


def first_plateau(counts, n):
    """Index of the window that completes the first run of n equal counts, or None."""
    for i in range(len(counts) - n + 1):
        if len(set(counts[i:i + n])) == 1:
            return i + n - 1
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", action="append", required=True,
                    help="state-0 solute trajectory; repeat to concatenate segments in order")
    ap.add_argument("--top", required=True)
    ap.add_argument("--sdf", required=True)
    ap.add_argument("--torsion", action="append", default=None)
    ap.add_argument("--torsions-file", action="append", default=None)
    ap.add_argument("--use-torsions", default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--no-single-cluster", action="store_true")
    ap.add_argument("--max-frames", type=int, default=None,
                    help="raise t_hdbscan's max_precomputed_frames (default 10000)")
    ap.add_argument("--step-ns", type=float, default=5.0, help="window increment (ns)")
    ap.add_argument("--plateau", type=int, default=3,
                    help="consecutive windows with the same count that define convergence")
    ap.add_argument("--name", default=None)
    ap.add_argument("--rest2-config", default=None)
    ap.add_argument("--out", default="convergence")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    P = lambda name: os.path.join(out, name)

    torsions = load_torsions(args.torsions_file, args.torsion, args.use_torsions)
    names, quads = list(torsions), list(torsions.values())

    print("loading trajectory")
    cmd = load_segments(args.traj, args.top, args.stride)
    tors = md.compute_dihedrals(cmd, quads)
    time_ns = cmd.time / 1000.0
    total_ns = simulated_ns(cmd, args.stride)
    ends = np.arange(args.step_ns, total_ns + 1e-6, args.step_ns)
    print(f"{cmd.n_frames} frames, {len(names)} torsions, {total_ns:g} ns -> {len(ends)} windows")

    asc = not args.no_single_cluster
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    with open(P("torsion_cv.json"), "w") as f:
        json.dump({"schema_version": 1, "indexing": "zero-based", "units": "radian",
                   "torsions": [{"name": n, "atoms": list(q)} for n, q in torsions.items()]},
                  f, indent=2)

    rows = []
    for end in ends:
        sel = time_ns <= end + 1e-6 if end < ends[-1] else np.ones(len(time_ns), bool)
        fit, _, why = build_fit(tors[sel], cmd[np.flatnonzero(sel)], mol, torsions,
                                allow_single_cluster=asc, max_frames=args.max_frames)
        row = {"t_ns": float(end), "n_frames": int(sel.sum()),
               "n_groups": int(fit.n_clusters_), "noise": float(fit.unassigned_mass_),
               "density_noise": float(fit.noise_mass_), "ambiguous": float(fit.ambiguous_mass_)}
        rows.append(row)
        print(f"  0-{end:g} ns: {row['n_frames']:6d} frames  {row['n_groups']:3d} clusters"
              f"  unassigned {100 * row['noise']:.1f}% (noise {100 * row['density_noise']:.1f}%,"
              f" ambiguous {100 * row['ambiguous']:.1f}%)", flush=True)

    df = pd.DataFrame(rows)
    counts = df["n_groups"].tolist()
    k = first_plateau(counts, args.plateau)
    tail_ok = len(counts) >= args.plateau and len(set(counts[-args.plateau:])) == 1
    df["converged"] = False
    if k is not None:
        df.loc[k:, "converged"] = True
    df.to_csv(P("convergence.csv"), index=False)
    print(df.to_string(index=False))
    if k is None:
        verdict = f"not converged: no {args.plateau} consecutive windows agree"
    else:
        verdict = (f"converged at {df.t_ns[k]:g} ns ({counts[k]} groups over "
                   f"{df.t_ns[k - args.plateau + 1]:g}-{df.t_ns[k]:g} ns)")
        if not tail_ok:
            verdict += f"; but the last {args.plateau} windows disagree ({counts[-args.plateau:]})"
    print(verdict)

    n_rep, tau_max = read_rest2_config(args.rest2_config)
    info = [f"{total_ns:g} ns", f"{n_rep} replicas" if n_rep else None,
            f"τ_max {tau_max}" if tau_max is not None else None, f"stride {args.stride}"]
    fig, ax = plt.subplots(figsize=(7.5, 4.5), layout="constrained")
    if k is not None:
        ax.axvspan(df.t_ns[k - args.plateau + 1] - args.step_ns / 2,
                   df.t_ns[k] + args.step_ns / 2, color="#2ca02c", alpha=0.15,
                   label=f"first {args.plateau} equal windows")
    ax.plot(df.t_ns, df.n_groups, "o-", color="#1f77b4", lw=2, label="symmetry-first clusters")
    ax.set_xlabel("simulated time, cumulative window 0 - t (ns)")
    ax.set_ylabel("number of clusters")
    ax.set_xticks(df.t_ns)
    ax.set_xlim(0, total_ns + args.step_ns / 2)
    ax.yaxis.get_major_locator().set_params(integer=True)
    ax.set_ylim(bottom=0)
    ax.grid(alpha=0.3)
    noise_ax = ax.twinx()
    noise_ax.plot(df.t_ns, 100 * df.noise, "s:", color="#d62728", ms=3, alpha=0.7)
    noise_ax.set_ylabel("unassigned (%)", color="#d62728")
    noise_ax.tick_params(axis="y", colors="#d62728")
    noise_ax.set_ylim(0, max(50, 100 * df.noise.max() * 1.15))
    ax.legend(fontsize=8, loc="upper left")
    ax.set_title(f"{args.name or 'ligand'} — REST2 " + ", ".join(x for x in info if x)
                 + f"\n{verdict}", fontsize=10)
    fig.savefig(P("convergence.png"), dpi=150)
    print(f"\nAll outputs in {out}/")


if __name__ == "__main__":
    main()
