#!/usr/bin/env python
"""
Shu-Yu-format report for Ioana's clustering.

Reads Ioana's cluster labels (clusters/labels.npy from cluster.py) and writes the SAME
file set that shuyu/cluster_prune.py writes, so the two protocols can be compared
file-for-file:

  torsion_cv.json        torsion definitions (name + atom indices)
  merged_groups.json     symmetry-pruned groups, masses, operations used
  groups.csv             per-group torsion mean/std, n_samples, mass
  representative_<A/B/…>.pdb   one representative ligand frame per group
  group_<i>.png          symmetry-group drawing (RDKit)
  group_torsions.png     per-group symmetry-aligned torsion distributions
  representatives.png    superposed representative structures
  torsions_<i>_<j>_2d.png  joint distribution of two torsions

The clustering itself is IOANA'S (cluster.py: --metric max, --drop-free, assign-to-cores);
only the symmetry pruning and the reporting come from MD-tools, exactly as cluster_prune.py
does it. Labels are all MD-tools' TorsionSymmetry needs -- it never re-fits anything.

Representatives here are the per-group MEDOID (min sum of torsion distances), whereas
cluster_prune.py uses representative_by_vote(); a cosmetic difference, the frame is still an
actual frame of the group.

Frame and torsion selection MUST match cluster.py (--stride/--skip and the
--drop-free/--min-barrier/--torsions-keep/--torsions-drop options), or the labels will not
line up. The defaults match.

Example:
  python group_report.py --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb \
      --sdf build/FOL.sdf --torsions REST2-run1/torsions.txt \
      --labels clusters/labels.npy --stride 5 --out groups
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

try:
    from md_tools.analysis import (load_torsion_json, load_cluster_torsion, TorsionSymmetry,
                                   verify_atom_mapping, group_name, strip_nonpolar_hydrogens,
                                   draw_representative_structures)
except ImportError as e:
    sys.exit(f"md_tools.analysis not available ({e}).\n"
             "Run this in the analysis-env micromamba environment.")

from rdkit import Chem
import pandas as pd

import symtools as st


def circular_mean(values):
    return float(np.arctan2(np.sin(values).mean(), np.cos(values).mean()))


def circular_std(values):
    resultant = np.hypot(np.sin(values).mean(), np.cos(values).mean())
    return float(np.sqrt(-2.0 * np.log(max(resultant, 1e-12))))


def medoid(phi_radians):
    """Index of the frame with the smallest total circular chord distance to the others."""
    if phi_radians.shape[0] == 1:
        return 0
    c = np.stack([np.cos(phi_radians), np.sin(phi_radians)], axis=-1)   # (n, k, 2)
    d = np.sqrt(((c[:, None, :, :] - c[None, :, :, :]) ** 2).sum(-1)).sum(-1)   # (n, n)
    return int(d.sum(1).argmin())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True)
    ap.add_argument("--top", required=True)
    ap.add_argument("--sdf", required=True)
    ap.add_argument("--torsions", required=True, help="torsions file used by cluster.py")
    ap.add_argument("--labels", required=True, help="labels.npy written by cluster.py")
    ap.add_argument("--skip", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--drop-free", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--min-barrier", type=float, default=2.0)
    ap.add_argument("--torsions-keep", default=None)
    ap.add_argument("--torsions-drop", default=None)
    ap.add_argument("--temp", type=float, default=298.0)
    ap.add_argument("--pair", default="3,4",
                    help="0-based positions in the kept torsion list for the 2-D plot (default 3,4)")
    ap.add_argument("--out", default="groups")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    P = lambda name: os.path.join(out, name)

    lab = np.load(args.labels)
    traj = md.load(args.traj, top=args.top, stride=args.stride)[args.skip:]
    if len(lab) != traj.n_frames:
        raise SystemExit(f"labels length {len(lab)} != {traj.n_frames} frames "
                         "(use the same --traj/--skip/--stride as cluster.py)")

    # torsion selection identical to cluster.py / dpca.py: all torsions from file, then
    # honour --torsions-drop / --drop-free / --torsions-keep.
    tors_all = np.loadtxt(args.torsions, dtype=int, ndmin=2, comments="#")
    quads_all = tors_all[:, :4]
    phi_all = np.degrees(md.compute_dihedrals(traj, quads_all)) % 360.0
    n_tor_all = quads_all.shape[0]
    dropped = np.zeros(n_tor_all, dtype=bool)
    if args.torsions_drop:
        dropped[[int(x) for x in args.torsions_drop.split(",")]] = True
    if args.drop_free:
        dropped |= st.free_rotor_mask(phi_all, args.min_barrier, temp=args.temp)
    keep = (np.array([int(x) for x in args.torsions_keep.split(",")]) if args.torsions_keep is not None
            else np.array([j for j in range(n_tor_all) if not dropped[j]]))
    if len(keep) == 0:
        raise SystemExit("no torsions left after selection - match cluster.py's selection options")

    prefix = os.path.splitext(os.path.basename(args.torsions))[0]
    names = [f"{prefix}_{int(j)}" for j in keep]        # name by ORIGINAL file position
    quads = [tuple(int(x) for x in quads_all[j]) for j in keep]
    print(f"trajectory: {traj.n_frames} frames; using {len(keep)} of {n_tor_all} torsions {list(map(int, keep))}")

    # ------------------------------------------------------------ reference molecule + mapping
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    if mol is None:
        sys.exit(f"RDKit could not read {args.sdf}")
    bonds = [(b[0].index, b[1].index) for b in traj.topology.bonds]
    if not bonds:
        sys.exit(f"{args.top} has no bonds - use a topology with CONNECT records / bonds")
    elements = [a.element.symbol for a in traj.topology.atoms]
    mapping, why = verify_atom_mapping(mol, bonds, elements)
    print("atom mapping:", "identity" if mapping == list(range(mol.GetNumAtoms())) else mapping)
    print(" ", why)

    tjson = P("torsion_cv.json")
    with open(tjson, "w") as f:
        json.dump({"schema_version": 1, "indexing": "zero-based", "units": "radian",
                   "torsions": [{"name": n, "atoms": list(q)} for n, q in zip(names, quads)]},
                  f, indent=2)
    torsion_idx = load_torsion_json(tjson, molecule=mol)

    # ------------------------------------------------------ Ioana labels -> MD-tools symmetry
    tors = md.compute_dihedrals(traj, quads)            # radians, (frames, n_tor)
    clusters = load_cluster_torsion(tors, lab, noise_label=-1, angle_unit="radian", names=names)
    print("cluster labels:", clusters.cluster_labels, "| masses:", clusters.mass_table())

    detector = TorsionSymmetry(mol, torsion_idx)
    detector.fit(clusters, coordinates=traj.xyz, atom_mapping=mapping)
    print("operations distinct on these torsions:", len(detector.operations))
    print("groups:", detector.cluster_groups)
    for c in detector.comparisons:
        if c.evidence.get("role") == "transitivity-guard":
            continue
        print(f"  {c.a} vs {c.b}: {c.verdict:9s} op={c.operation} "
              f"E={c.energy if c.energy is None else round(c.energy, 4)}")

    merged = detector.deduplicate(clusters)
    desc = merged.describe()
    print(json.dumps(desc, indent=2))
    with open(P("merged_groups.json"), "w") as f:
        json.dump(desc, f, indent=2)
    total = sum(merged.masses.values()) + merged.noise_mass
    print(f"populations conserved: merged masses + noise = {total:.6f}")

    # ---------------------------------------------------------------- representatives (medoid)
    groups = {group_name(i): g for i, g in enumerate(detector.cluster_groups)}
    chosen = {}
    for i, name in enumerate(groups):
        frames = np.flatnonzero(merged.labels == i)
        chosen[name] = int(frames[medoid(merged.aligned[i])]) if len(frames) else int(merged.entries[i].frames[0])
    aligned = {name: merged.aligned[i] for i, name in enumerate(groups)}
    for name, members in groups.items():
        print(f"group {name}: clusters {members}  mass {merged.masses[list(groups).index(name)]:.4f}"
              f"  representative frame {chosen[name]}")

    # ---------------------------------------------------------------- group drawings (RDKit)
    for g in range(len(detector.cluster_groups)):
        if len(detector.cluster_groups[g]) > 1:          # only merged groups have an operation
            try:
                detector.draw_operation_effect(group=g, output=P(f"operation_group_{g}.png"))
            except Exception as e:
                print(f"  (could not draw operation for group {g}: {e})")
    for i in range(len(detector.cluster_groups)):
        try:
            detector.draw_sym_group(i, output=P(f"group_{i}.png"), cluster_torsions=clusters)
        except Exception as e:
            print(f"  (could not draw group {i}: {e})")

    # ------------------------------------------------- 1-D torsion distributions per group
    ncol = min(3, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.7 * ncol, 3.6 * nrow),
                             layout="constrained", squeeze=False)
    for extra in axes.flat[len(names):]:
        extra.axis("off")
    edges = np.linspace(-np.pi, np.pi, 49)
    palette = plt.get_cmap("tab10")
    for column, torsion in enumerate(names):
        axis = axes.flat[column]
        for k, (name, values) in enumerate(aligned.items()):
            axis.hist(values[:, column], bins=edges, density=True, histtype="step", lw=2.0,
                      color=palette(k), label=f"group {name} ({100 * float(merged.masses[k]):.1f}%)")
        axis.set_xlim(-np.pi, np.pi)
        axis.set_xlabel(f"{torsion} (rad)")
        axis.set_title(torsion)
        if column == 0:
            axis.set_ylabel("density")
            axis.legend(fontsize=9)
    fig.suptitle("Symmetry-aligned torsion distributions, one colour per merged group", fontsize=10)
    fig.savefig(P("group_torsions.png"), dpi=150)

    # ------------------------------------------------- 2-D joint distribution of two torsions
    ix, iy = (int(x) for x in args.pair.split(","))
    if len(names) > max(ix, iy):
        panels = 1 + len(groups)
        ncol2 = min(3, panels)
        nrow2 = int(np.ceil(panels / ncol2))
        fig2, axes2 = plt.subplots(nrow2, ncol2, figsize=(4.5 * ncol2, 4.1 * nrow2),
                                   layout="constrained", squeeze=False)
        for extra in axes2.flat[panels:]:
            extra.axis("off")
        lim = (-np.pi, np.pi, -np.pi, np.pi)

        def _hex2d(axis, xs, ys, title):
            hb = axis.hexbin(xs, ys, gridsize=55, mincnt=1, cmap="viridis", bins="log", extent=lim)
            fig2.colorbar(hb, ax=axis, label="log10(count)")
            axis.set_title(title)
            axis.set_xlabel(f"{names[ix]} (rad)")
            axis.set_ylabel(f"{names[iy]} (rad)")

        _hex2d(axes2.flat[0], tors[:, ix], tors[:, iy], "all frames")
        for k, (name, values) in enumerate(aligned.items()):
            _hex2d(axes2.flat[k + 1], values[:, ix], values[:, iy],
                   f"group {name} ({100 * float(merged.masses[k]):.1f}%)")
        fig2.suptitle(f"Joint distribution: {names[ix]} vs {names[iy]}", fontsize=11)
        fig2.savefig(P(f"torsions_{ix}_{iy}_2d.png"), dpi=150)

    # ---------------------------------------------------------------- groups.csv
    rows = []
    for i, (name, members) in enumerate(groups.items()):
        values = aligned[name]
        row = {"cluster_id": name, "original_clusters": ",".join(str(m) for m in members)}
        for column, torsion in enumerate(names):
            row[f"{torsion}_avg"] = np.rad2deg(circular_mean(values[:, column]))
            row[f"{torsion}_std"] = np.rad2deg(circular_std(values[:, column]))
        row["n_samples"] = int(values.shape[0])
        row["mass"] = float(merged.masses[i])
        rows.append(row)
    df = pd.DataFrame(rows).round(3)
    print(df.to_string(index=False))
    df.to_csv(P("groups.csv"), index=False)

    # ---------------------------------------------------------------- representatives figure + PDBs
    keep_atoms = strip_nonpolar_hydrogens(traj)
    try:
        info = draw_representative_structures(traj, chosen, output=P("representatives.png"),
                                              atom_indices=keep_atoms)
        print(f"{info['n_atoms_drawn']} atoms drawn, superposed: {info['aligned']}")
    except Exception as e:
        print(f"  (could not draw representatives: {e})")
    for name, frame in chosen.items():
        traj[frame].save_pdb(P(f"representative_{name}.pdb"))
    print(f"\nAll outputs in {out}/")


if __name__ == "__main__":
    main()
