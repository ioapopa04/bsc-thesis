#!/usr/bin/env python
"""
Folic acid in water vs bound to 4LRH (FRalpha) vs bound to 3QL3 (DHFR): one figure, one panel per
torsion, the three state-0 distributions overlaid (all frames, no clustering), plus each complex's
crystal torsion as a vertical line in that complex's colour.

  water (grey, filled)    FOLn-REST2-50ns            REST2 in water, 8 replicas
  4LRH complex (red)      FOL-complex-REST2/4LRH     REST2 bound to FRalpha
  3QL3 complex (blue)     FOL-complex-REST2/3QL3     REST2 bound to DHFR + NADP+

Reading a panel:
  - complex peak on a water peak or tail  -> the protein SELECTS a conformer water already samples
  - complex peak where water has ~nothing -> the protein INDUCES it
  - complex peak away from its own line   -> the complex does not hold the crystal value

Same inputs and conventions as cluster_prune.py (imported from here): the water run's torsion files
(ligand-local indices 0..50, numbered t0..t11 in file order: torsions, torsions_rigid,
torsions_polarH), stride 5, crystal lines matched on heavy-atom element + connectivity (a carboxyl
gets two lines, polar-H torsions none). The complexes use the folate-only trajectories cut out by
FOL-complex-REST2/scripts/shuyu_analysis.sbatch (same 51 atoms, same order as the water run).
Distributions are normalised densities, so the different run lengths (water 100 ns, complexes
50 ns) do not bias the comparison.

Also writes a CSV per torsion x complex: the complex's peak, the water population within
+-WINDOW of that peak (selection vs induction), and the complex population within +-WINDOW of its
own crystal value(s).

Run in analysis-env (needs md_tools.analysis via cluster_prune), through Slurm:
  srun -p defq -x sbrinz,sbrinz1 -c 4 --mem=16G python compare_water_complex.py
"""
import argparse
import io
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import mdtraj as md
import pandas as pd
from rdkit import Chem

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_prune import load_torsions, reference_torsions, molecule_png, simulated_ns

WATER = "/scratch/ipopa/thesis/FOLn-REST2-50ns"
COMPLEX = "/scratch/ipopa/thesis/FOL-complex-REST2"
CRYSTAL = "/scratch/ipopa/thesis/bsc-thesis/crystal"
RESULTS = "/scratch/ipopa/thesis/bsc-thesis/results/folic-acid"


def systems(water, cplx):
    """label -> (trajectory, topology, crystal sdf or None, colour). A complex uses its own
    folate-only PDB when the cut wrote one, else the water run's (same atoms, same order)."""
    out = {"water": (f"{water}/REST2-run1/solute_state0_prod1.nc", f"{water}/build/built.solute.pdb",
                     None, "#9a9a9a")}
    for pdb, colour in (("4LRH", "#d62728"), ("3QL3", "#1f77b4")):
        d = f"{cplx}/{pdb}/analysis-shuyu"
        top = f"{d}/FOL.solute.pdb" if os.path.exists(f"{d}/FOL.solute.pdb") else f"{water}/build/built.solute.pdb"
        out[f"{pdb} complex"] = (f"{d}/FOL_state0_prod1.nc", top, f"{CRYSTAL}/fol_{pdb}_A.sdf", colour)
    return out


def wrap(x):
    return np.angle(np.exp(1j * np.asarray(x)))


def mass_within(values, centre, window):
    """Fraction of frames within +-window (rad) of centre, on the circle."""
    return float(np.mean(np.abs(wrap(values - centre)) <= window))


def peak(values, edges):
    h, _ = np.histogram(values, bins=edges)
    return 0.5 * (edges[h.argmax()] + edges[h.argmax() + 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--water", default=WATER, help="water REST2 run directory")
    ap.add_argument("--complex", default=COMPLEX, help="directory holding 4LRH/ and 3QL3/")
    ap.add_argument("--sdf", default=None, help="ligand SDF (default: <water>/build/FOL.sdf)")
    ap.add_argument("--stride", type=int, default=5)
    ap.add_argument("--bins", type=int, default=72, help="histogram bins over 360 deg")
    ap.add_argument("--window", type=float, default=30.0,
                    help="+- deg around a peak / crystal value for the CSV populations")
    ap.add_argument("--label", default="", help="added to the figure title, e.g. 'v4 (COO-, amide scaled)'")
    ap.add_argument("--out", default=f"{RESULTS}/folic-acid_water_vs_complexes",
                    help="output prefix (.png and .csv are appended)")
    args = ap.parse_args()

    SYSTEMS = systems(args.water, args.complex)
    args.sdf = args.sdf or f"{args.water}/build/FOL.sdf"
    torsions = load_torsions([f"{args.water}/REST2-run1/{f}.txt"
                              for f in ("torsions", "torsions_rigid", "torsions_polarH")])
    names, quads = list(torsions), list(torsions.values())
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    if mol is None:
        sys.exit(f"RDKit could not read {args.sdf}")

    tors, refs, ns = {}, {}, {}
    atom_names = None
    for label, (traj, top, ref, _) in SYSTEMS.items():
        t = md.load(traj, top=top, stride=args.stride)
        if t.n_atoms != mol.GetNumAtoms():
            sys.exit(f"{label}: {t.n_atoms} atoms, ligand has {mol.GetNumAtoms()}")
        here = [[t.topology.atom(i).name for i in q] for q in quads]
        if atom_names is None:
            atom_names = here
        elif here != atom_names:
            sys.exit(f"{label}: torsion atom names differ from the water run - atom order changed")
        tors[label] = md.compute_dihedrals(t, quads)
        ns[label] = simulated_ns(t, args.stride)
        print(f"{label:13s} {t.n_frames:6d} frames (stride {args.stride}), {ns[label]:.1f} ns")
        if ref:
            refs[label] = reference_torsions(ref, mol, quads)

    edges = np.linspace(-np.pi, np.pi, args.bins + 1)
    win = np.deg2rad(args.window)
    colours = {label: s[3] for label, s in SYSTEMS.items()}

    # ---------------------------------------------------------------- table
    rows = []
    for k, name in enumerate(names):
        for label, values in refs.items():
            c = tors[label][:, k]
            p = peak(c, edges)
            row = {"torsion": f"t{k}", "name": name, "atoms": "-".join(atom_names[k]),
                   "complex": label, "complex_peak_deg": round(np.rad2deg(p), 1),
                   "complex_mass_at_peak": round(mass_within(c, p, win), 3),
                   "water_mass_at_complex_peak": round(mass_within(tors["water"][:, k], p, win), 3),
                   "crystal_deg": ", ".join(f"{np.rad2deg(v):.1f}" for v in values[k]) or "n/a"}
            if values[k]:
                row["complex_mass_at_crystal"] = round(max(mass_within(c, v, win) for v in values[k]), 3)
                row["water_mass_at_crystal"] = round(
                    max(mass_within(tors["water"][:, k], v, win) for v in values[k]), 3)
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(args.out + ".csv", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(f"\npopulations within +-{args.window:g} deg (complex peak / crystal value):")
        print(df.drop(columns=["name"]).to_string(index=False))

    # ---------------------------------------------------------------- figure
    ntor = len(names)
    ncol = 4
    nrow = -(-ntor // ncol)
    panel_w, panel_h, mol_w = 2.7, 2.1, 5.0
    fig = plt.figure(figsize=(mol_w + ncol * panel_w, nrow * panel_h + 1.2), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=[mol_w, ncol * panel_w])
    mol_ax = fig.add_subplot(outer[0])
    png = molecule_png(mol, quads, ["#333333"] * ntor, size=(1100, 1100))
    mol_ax.imshow(mpimg.imread(io.BytesIO(png), format="png"))
    mol_ax.axis("off")
    mol_ax.set_title("t<k> marks the central bond of torsion k", fontsize=10)
    grid = outer[1].subgridspec(nrow, ncol)

    xt = np.deg2rad([-180, -90, 0, 90, 180])
    for k, name in enumerate(names):
        ax = fig.add_subplot(grid[divmod(k, ncol)])
        for label, values in tors.items():
            if label == "water":
                ax.hist(values[:, k], bins=edges, density=True, histtype="stepfilled",
                        color=colours[label], alpha=0.55, lw=0, label=label)
            else:
                ax.hist(values[:, k], bins=edges, density=True, histtype="step",
                        color=colours[label], lw=1.6, label=label)
        for label, values in refs.items():
            for n, v in enumerate(values[k]):
                ax.axvline(v, color=colours[label], ls="--", lw=1.3,
                           label=(f"crystal {label.split()[0]}" if n == 0 else None))
        ax.set_xlim(-np.pi, np.pi)
        ax.set_xticks(xt, ["-180", "-90", "0", "90", "180"])
        ax.set_yticks([])
        ax.tick_params(labelsize=8)
        kind = name.rsplit("_", 1)[0].replace("torsions_", "").replace("torsions", "rotatable")
        ax.set_title(f"t{k}  {'-'.join(atom_names[k])}  ({kind})", fontsize=8.5)
        if k == 0:
            ax.legend(fontsize=7, loc="upper left", frameon=False)
    lengths = ", ".join(f"{label} {ns[label]:.0f} ns" for label in tors)
    fig.suptitle(f"Folic acid{' ' + args.label if args.label else ''} torsions: water vs bound (state 0, all frames; "
                 f"{lengths}; dashed = crystal)", fontsize=12)
    fig.savefig(args.out + ".png", dpi=150)
    plt.close(fig)
    print(f"\nwrote {args.out}.png and {args.out}.csv")


if __name__ == "__main__":
    main()
