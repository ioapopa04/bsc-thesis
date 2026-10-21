"""Density of ligand-conformation RMSD to the crystal, from the per-frame values that
rmsd_vs_crystal.py writes (<prefix>.csv). Four step histograms, same bins, each normalised to unit
area so runs of different length (water 100 ns, complexes 50 ns) compare directly.

  python -I scripts/rmsd_histogram.py <prefix>        # reads <prefix>.csv, writes <prefix>_hist.png
"""
import csv
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PREFIX = sys.argv[1]
data = defaultdict(list)
with open(f"{PREFIX}.csv") as f:
    for row in csv.DictReader(f):
        data[(row["run"], row["crystal"])].append(float(row["rmsd_A"]))

LINES = [  # same colours as rmsd_vs_crystal.py
    ("water", "4LRH", "#d62728", "-", "folic acid in water vs 4LRH crystal"),
    ("water", "3QL3", "#1f77b4", "-", "folic acid in water vs 3QL3 crystal"),
    ("4LRH complex", "4LRH", "#7f0000", "--", "folic acid in 4LRH complex vs 4LRH crystal"),
    ("3QL3 complex", "3QL3", "#08306b", "--", "folic acid in 3QL3 complex vs 3QL3 crystal"),
]
top = max(max(v) for v in data.values())
edges = np.arange(0.0, np.ceil(top / 0.05) * 0.05 + 0.05, 0.05)

fig, ax = plt.subplots(figsize=(9, 4.8), layout="constrained")
for run, pdb, colour, style, label in LINES:
    r = np.asarray(data[(run, pdb)])
    ax.hist(r, bins=edges, density=True, histtype="step", lw=2.0, color=colour, ls=style,
            label=f"{label} (n = {len(r)}, median {np.median(r):.2f} Å)")
ax.set_xlabel("heavy-atom RMSD to crystal conformation (Å)")
ax.set_ylabel("probability density (1/Å)")
ax.set_xlim(0, edges[-1])
ax.set_ylim(bottom=0)
ax.grid(alpha=0.3)
ax.legend(fontsize=9, loc="upper right")
ax.set_title("Folic acid (MACE-OMOL bespoke torsions, 3H tautomer): distribution of the ligand-conformation RMSD to the crystal\n"
             "REST2 state 0, every saved frame (2 ps); bins 0.05 Å; each curve normalised to area 1",
             fontsize=10)
fig.savefig(f"{PREFIX}_hist.png", dpi=150)
for run, pdb, *_ in LINES:
    r = np.asarray(data[(run, pdb)])
    print(f"{run:13s} vs {pdb}: n={len(r)}  mode bin {edges[np.argmax(np.histogram(r, edges)[0])]:.2f}-"
          f"{edges[np.argmax(np.histogram(r, edges)[0]) + 1]:.2f} A  median {np.median(r):.2f}  "
          f"frac < 1 A {np.mean(r < 1):.1%}")
print(f"written {PREFIX}_hist.png")
