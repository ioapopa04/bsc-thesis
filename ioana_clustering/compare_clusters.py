#!/usr/bin/env python
"""
Cross-tabulate Ioana's clustering against Shu-Yu's, frame by frame.

Shu-Yu's cluster_prune.py does not save per-frame labels, so they are recomputed here
deterministically with the same MD-tools call cluster_prune.py makes:

    t_hdbscan(md.compute_dihedrals(traj, quads), names=..., allow_single_cluster=False)

on the torsion files it used (torsions.txt + torsions_rigid.txt + torsions_polarH.txt) at
the same --stride. Ioana's labels come from clusters/labels.npy. The script prints and writes
the contingency table (counts and row/column fractions), so you can see where Ioana's clusters
fall relative to Shu-Yu's groups -- e.g. whether Ioana's several clusters are all copies of the
same two Shu-Yu groups, or whether they split the frames Shu-Yu left as noise.

Example:
  python compare_clusters.py --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb \
      --stride 5 --labels clusters/labels.npy \
      --torsions-file REST2-run1/torsions.txt --torsions-file REST2-run1/torsions_rigid.txt \
      --torsions-file REST2-run1/torsions_polarH.txt --out crosstab.csv
"""
import argparse
import os
import sys

import numpy as np
import mdtraj as md
import pandas as pd

try:
    from md_tools.analysis import t_hdbscan
except ImportError as e:
    sys.exit(f"md_tools.analysis not available ({e}). Run in the analysis-env environment.")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--labels", required=True, help="Ioana's labels.npy (cluster.py)")
    ap.add_argument("--torsions-file", action="append", required=True,
                    help="repeatable; the torsion files cluster_prune.py was given")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--out", default="crosstab.csv")
    args = ap.parse_args()

    ioana = np.load(args.labels)
    traj = md.load(args.traj, top=args.top, stride=args.stride)
    if len(ioana) != traj.n_frames:
        raise SystemExit(f"Ioana labels {len(ioana)} != {traj.n_frames} frames at stride "
                         f"{args.stride} - pass the same --stride as cluster.py")

    torsions, quads, names = {}, [], []
    for path in args.torsions_file:
        label = os.path.splitext(os.path.basename(path))[0]
        for line in open(path):
            if line.strip().startswith("#") or not line.strip():
                continue
            q = tuple(int(x) for x in line.split("#")[0].split()[:4])
            names.append(f"{label}_{len(torsions)}")
            quads.append(q)
            torsions[names[-1]] = q
    if not quads:
        sys.exit("no torsions read")
    tors = md.compute_dihedrals(traj, quads)
    print(f"{traj.n_frames} frames, {len(quads)} torsions (Shu-Yu set): {names}")

    fit = t_hdbscan(tors, names=names, allow_single_cluster=False).fit()
    shuyu = np.asarray(fit.labels_)
    print(f"Shu-Yu recomputed: n_clusters={fit.summary()['n_clusters']}, "
          f"masses={ {k: round(v, 3) for k, v in fit.summary()['cluster_mass'].items()} }, "
          f"noise={fit.summary()['noise_mass']:.4f}")

    di = pd.Series(ioana, name="ioana").map(lambda x: f"ioana_{x}" if x >= 0 else "ioana_noise")
    ds = pd.Series(shuyu, name="shuyu").map(lambda x: f"shuyu_{x}" if x >= 0 else "shuyu_noise")
    counts = pd.crosstab(di, ds)
    row_frac = counts.div(counts.sum(1), axis=0)          # of each Ioana cluster, where it lands
    col_frac = counts.div(counts.sum(0), axis=1)          # of each Shu-Yu group, where it came from

    print("\ncounts (rows = Ioana clusters, cols = Shu-Yu groups):")
    print(counts.to_string())
    print("\nrow fraction (fraction of each Ioana cluster per Shu-Yu group):")
    print(row_frac.round(3).to_string())
    print("\ncolumn fraction (fraction of each Shu-Yu group per Ioana cluster):")
    print(col_frac.round(3).to_string())

    with open(args.out, "w") as f:
        f.write("# counts\n")
        counts.to_csv(f)
        f.write("\n# row_fraction\n")
        row_frac.round(5).to_csv(f)
        f.write("\n# column_fraction\n")
        col_frac.round(5).to_csv(f)
    print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
