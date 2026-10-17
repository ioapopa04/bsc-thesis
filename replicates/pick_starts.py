#!/usr/bin/env python
"""
Pick DIFFERENT starting conformers for REST2 replicates (Shu-Yu: "does a simulation started from a
different RDKit conformer give the same result?").

So far every water run started from ONE structure: the parameter package's SMILES was embedded
(ETKDGv3, 10 conformers, one fixed seed), MMFF94s-minimised, and the lowest-energy conformer was
kept as the package SDF; build-top used those coordinates as given. That run is replicate 0 here.

This script
  1. clusters the EXISTING run (replicate 0) with Shu-Yu's symmetry-first t_hdbscan, exactly as
     shuyu_clustering/cluster_prune.py does (same torsions, stride, allow_single_cluster);
  2. embeds many new conformers of the PACKAGE molecule (same atoms, order, names, hydrogens,
     charges and stereo -- so the build reuses the package and no charge is recomputed), with
     ETKDGv3 + MMFF94s as MD-tools does, and drops any whose stereochemistry changed;
  3. labels every candidate with the replicate-0 clusters (fit.predict: the same symmetry-aware
     distance and 18-of-20 vote) and measures the symmetry-aware distance between candidates;
  4. picks --n-pick starts greedily: first one per cluster NOT yet covered by the starts chosen so
     far (replicate 0's own cluster counts as covered), each time the candidate farthest from all
     chosen starts; once every reachable cluster is covered, farthest-point (maximin) over all.
     Only candidates within --energy-window kcal/mol of the lowest MMFF energy are eligible.

Outputs (--out):
  candidates.csv        every candidate: MMFF energy, converged, stereo ok, cluster, distance to
                        replicate 0 and to the nearest chosen start, torsions (deg)
  starts/rep<k>.sdf     the chosen starts (rep0 = the original start, copied), package atom order
  selection.json        what was chosen and why, cluster masses, coverage, settings
  starts_torsions.png   replicate-0 torsion distributions with the starts as coloured lines
  starts_distance.png   symmetry-aware distance matrix of the starts, with their clusters

Build each start afterwards with build-top -i starts/rep<k>.sdf and `solute.parameters: <package>`
(reuse is required: the input must be the package's exact chemical state).
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
from rdkit.Chem import AllChem, rdDistGeom

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "shuyu_clustering"))
from cluster_prune import load_torsions, build_fit, TORSION_PALETTE  # noqa: E402
from md_tools.analysis import group_name  # noqa: E402

REP_COLOURS = ["black"] + TORSION_PALETTE[:9]


def label_name(c):
    c = int(c)
    return {-2: "noise", -1: "ambiguous"}[c] if c < 0 else group_name(c)


def isomeric_key(mol):
    """Canonical isomeric SMILES with the stereo read from THIS mol's 3D coordinates."""
    m = Chem.Mol(mol)
    Chem.AssignStereochemistryFrom3D(m)
    return Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True)


def embed_candidates(mol, n, seed, prune, max_iters):
    m = Chem.Mol(mol)
    m.RemoveAllConformers()
    params = rdDistGeom.ETKDGv3()
    params.randomSeed = int(seed)
    params.pruneRmsThresh = float(prune)
    params.numThreads = 0
    cids = list(rdDistGeom.EmbedMultipleConfs(m, int(n), params))
    if not cids:
        sys.exit("ETKDGv3 produced no conformers")
    res = AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0, maxIters=int(max_iters),
                                            mmffVariant="MMFF94s")
    return m, cids, [(not_conv == 0, e) for not_conv, e in res]


def closed_embedding(fit, tors, traj):
    """Rows in the fit's symmetry descriptor space (what predict() and the distance use)."""
    desc = fit.descriptor_
    closed = desc.closed_angles(tors, coordinates=traj, tolerance=fit.stored_tolerance,
                                **getattr(fit, "_coordinate_kwargs", {}))
    return desc.embed(closed)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--package-sdf", required=True, help="molecule.sdf of the parameter package")
    ap.add_argument("--start-sdf", required=True, help="the coordinates replicate 0 was built from")
    ap.add_argument("--traj", required=True, help="replicate-0 state-0 solute trajectory")
    ap.add_argument("--top", required=True, help="its topology with bonds (built.solute.pdb)")
    ap.add_argument("--torsion", action="append", default=None)
    ap.add_argument("--torsions-file", action="append", default=None)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--no-single-cluster", action="store_true")
    ap.add_argument("--max-frames", type=int, default=None)
    ap.add_argument("--n-embed", type=int, default=300, help="ETKDGv3 embeddings (before pruning)")
    ap.add_argument("--seed", type=int, default=20261017)
    ap.add_argument("--prune-rms", type=float, default=0.5, help="ETKDG pruneRmsThresh (MD-tools: 0.5)")
    ap.add_argument("--max-iters", type=int, default=2000, help="MMFF94s iterations")
    ap.add_argument("--energy-window", type=float, default=25.0,
                    help="eligible: MMFF energy within this many kcal/mol of the lowest candidate")
    ap.add_argument("--n-pick", type=int, default=3, help="new starts besides replicate 0")
    ap.add_argument("--min-distance", type=float, default=0.1,
                    help="a start must be at least this far (symmetry-aware) from every chosen start; "
                         "fewer than --n-pick are returned if the molecule has no more distinct minima")
    ap.add_argument("--name", default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(os.path.join(out, "starts"), exist_ok=True)
    P = lambda *p: os.path.join(out, *p)
    title = f"{args.name}: " if args.name else ""

    # -------------------------------------------------- molecule, atom order, replicate-0 start
    mol = Chem.SDMolSupplier(args.package_sdf, removeHs=False)[0]
    start = Chem.SDMolSupplier(args.start_sdf, removeHs=False)[0]
    if mol is None or start is None:
        sys.exit("RDKit could not read the package or start SDF")
    top = md.load(args.top).topology
    elems = [a.GetSymbol() for a in mol.GetAtoms()]
    if [a.element.symbol for a in top.atoms] != elems or \
            [a.GetSymbol() for a in start.GetAtoms()] != elems:
        sys.exit("package SDF, start SDF and topology do not share one atom order")
    if mol.HasProp("MDT_ATOM_NAMES") and mol.GetProp("MDT_ATOM_NAMES").split() != [a.name for a in top.atoms]:
        sys.exit("package atom names differ from the topology's")
    ref_key = isomeric_key(start)
    print("molecule:", Chem.MolToSmiles(Chem.RemoveHs(mol)), f"({mol.GetNumAtoms()} atoms)")

    # ------------------------------------------------------------ 1. cluster replicate 0
    torsions = load_torsions(args.torsions_file, args.torsion)
    names, quads = list(torsions), list(torsions.values())
    traj = md.load(args.traj, top=args.top, stride=args.stride)
    tors = md.compute_dihedrals(traj, quads)
    print(f"replicate 0: {tors.shape[0]} frames, {tors.shape[1]} torsions; fitting ...")
    fit, _, _ = build_fit(tors, traj, mol, torsions, allow_single_cluster=not args.no_single_cluster,
                          max_frames=args.max_frames)
    masses = {c: float(fit.cluster_mass_[c]) for c in range(fit.n_clusters_)}
    print(f"n_clusters = {fit.n_clusters_}:",
          {group_name(c): round(m, 4) for c, m in masses.items()},
          f"unassigned {fit.unassigned_mass_:.4f}")

    # ------------------------------------------------------------ 2. candidates
    cand, cids, results = embed_candidates(mol, args.n_embed, args.seed, args.prune_rms, args.max_iters)
    energies = np.array([e for _, e in results])
    converged = np.array([c for c, _ in results])
    stereo_ok = np.array([isomeric_key(Chem.Mol(cand, confId=cid)) == ref_key for cid in cids])
    print(f"{len(cids)} candidates after pruning; {int(converged.sum())} converged; "
          f"{int(stereo_ok.sum())} with unchanged stereo")

    # all structures in one trajectory: row 0 = replicate-0 start, rows 1.. = candidates
    xyz = [start.GetConformer().GetPositions()] + [cand.GetConformer(cid).GetPositions() for cid in cids]
    alltraj = md.Trajectory(np.asarray(xyz) / 10.0, top)
    alltors = md.compute_dihedrals(alltraj, quads)

    # ------------------------------------------------------------ 3. labels and distances
    labels, _ = fit.predict(alltors, coordinates=alltraj)
    X = closed_embedding(fit, alltors, alltraj)
    D = fit.descriptor_.distance_matrix(X)
    print(f"replicate-0 start: cluster {label_name(labels[0])}")
    cand_labels = labels[1:]
    reach = {label_name(c): int((cand_labels == c).sum()) for c in sorted(set(cand_labels.tolist()))}
    print("candidates per cluster:", reach)

    # ------------------------------------------------------------ 4. greedy selection
    emin = energies[converged & stereo_ok].min()
    eligible = converged & stereo_ok & (energies <= emin + args.energy_window)
    print(f"eligible (converged, stereo, E <= {emin:.1f} + {args.energy_window}): {int(eligible.sum())}")
    chosen, why = [0], ["original start (replicate 0)"]           # rows of D / alltraj
    covered = {int(labels[0])}
    for _ in range(args.n_pick):
        dmin = D[:, chosen].min(axis=1)
        pool = [r for r in range(1, len(xyz))
                if eligible[r - 1] and r not in chosen and dmin[r] >= args.min_distance]
        new = [r for r in pool if labels[r] >= 0 and int(labels[r]) not in covered]
        if new:
            r = max(new, key=lambda r: dmin[r])
            why.append(f"first start in cluster {label_name(labels[r])}, farthest from the chosen")
        elif pool:
            r = max(pool, key=lambda r: dmin[r])
            why.append(f"farthest point (every reachable cluster covered); cluster {label_name(labels[r])}")
        else:
            print(f"no eligible candidate left at >= {args.min_distance} from every chosen start")
            break
        chosen.append(r)
        covered.add(int(labels[r]))

    # ------------------------------------------------------------ outputs
    dmin_all = D[:, chosen].min(axis=1)
    with open(P("candidates.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["row", "conf_id", "mmff_kcal_mol", "converged", "stereo_ok", "eligible", "cluster",
                    "d_to_rep0", "d_to_nearest_start", "chosen_as"] + [f"{n}_deg" for n in names])
        for r in range(len(xyz)):
            rep = f"rep{chosen.index(r)}" if r in chosen else ""
            if r == 0:
                w.writerow([0, "start", "", "", True, "", label_name(labels[0]), 0.0, 0.0, rep]
                           + [f"{v:.1f}" for v in np.rad2deg(alltors[0])])
                continue
            k = r - 1
            w.writerow([r, cids[k], f"{energies[k]:.3f}", bool(converged[k]), bool(stereo_ok[k]),
                        bool(eligible[k]), label_name(labels[r]), f"{D[r, 0]:.3f}", f"{dmin_all[r]:.3f}",
                        rep] + [f"{v:.1f}" for v in np.rad2deg(alltors[r])])

    starts = []
    for k, r in enumerate(chosen):
        m = Chem.Mol(start) if r == 0 else Chem.Mol(cand, confId=cids[r - 1])
        m.SetProp("_Name", mol.GetProp("_Name") if mol.HasProp("_Name") else "LIG")
        m.SetProp("REPLICATE_START", f"rep{k}: {why[k]}")
        Chem.MolToMolFile(m, P("starts", f"rep{k}.sdf"))
        starts.append({"replicate": k, "row": int(r), "conf_id": None if r == 0 else int(cids[r - 1]),
                       "mmff_kcal_mol": None if r == 0 else float(energies[r - 1]),
                       "mmff_above_min": None if r == 0 else float(energies[r - 1] - emin),
                       "cluster": label_name(labels[r]),
                       "cluster_mass_rep0": masses.get(int(labels[r])),
                       "d_to_rep0": float(D[r, 0]), "why": why[k],
                       "torsions_deg": {n: round(float(v), 1) for n, v in zip(names, np.rad2deg(alltors[r]))},
                       "sdf": f"starts/rep{k}.sdf"})
        print(f"rep{k}: row {r}  cluster {label_name(labels[r])}  d(rep0) {D[r, 0]:.2f}  -- {why[k]}")
    with open(P("selection.json"), "w") as f:
        json.dump({"name": args.name, "molecule": Chem.MolToSmiles(Chem.RemoveHs(mol)),
                   "replicate0": {"traj": os.path.abspath(args.traj), "stride": args.stride,
                                  "frames": int(tors.shape[0]),
                                  "allow_single_cluster": not args.no_single_cluster,
                                  "n_clusters": int(fit.n_clusters_),
                                  "cluster_masses": {group_name(c): m for c, m in masses.items()},
                                  "unassigned_mass": float(fit.unassigned_mass_)},
                   "torsions": {n: list(q) for n, q in torsions.items()},
                   "candidates": {"n_embedded": len(cids), "n_converged": int(converged.sum()),
                                  "n_stereo_ok": int(stereo_ok.sum()), "n_eligible": int(eligible.sum()),
                                  "mmff_min_kcal_mol": float(emin), "per_cluster": reach},
                   "settings": {k: v for k, v in vars(args).items()},
                   "starts": starts}, f, indent=2)

    # torsion distributions of replicate 0 with the starts
    ncol = min(4, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.8 * ncol, 2.9 * nrow), layout="constrained", squeeze=False)
    edges = np.linspace(-180, 180, 61)
    for extra in axes.flat[len(names):]:
        extra.axis("off")
    for j, (axis, n) in enumerate(zip(axes.flat, names)):
        axis.hist(np.rad2deg(tors[:, j]), bins=edges, density=True, color="0.8")
        for k, r in enumerate(chosen):
            axis.axvline(np.rad2deg(alltors[r, j]), color=REP_COLOURS[k % len(REP_COLOURS)], lw=1.8,
                         ls="-" if k == 0 else "--", label=f"rep{k} ({label_name(labels[r])})")
        axis.set_xlim(-180, 180)
        axis.set_title(n, fontsize=9)
        if j == 0:
            axis.legend(fontsize=7)
    fig.suptitle(f"{title}replicate-0 torsion distributions (grey) and the chosen starting conformers "
                 "(raw labelling)", fontsize=10)
    fig.savefig(P("starts_torsions.png"), dpi=150)
    plt.close(fig)

    sub = D[np.ix_(chosen, chosen)]
    labs = [f"rep{k}\n{label_name(labels[r])}" for k, r in enumerate(chosen)]
    fig, axis = plt.subplots(figsize=(1.2 * len(chosen) + 2.5, 1.1 * len(chosen) + 1.5), layout="constrained")
    im = axis.imshow(sub, cmap="viridis")
    for a in range(len(chosen)):
        for b in range(len(chosen)):
            axis.text(b, a, f"{sub[a, b]:.2f}", ha="center", va="center",
                      color="white" if sub[a, b] < sub.max() / 2 else "black", fontsize=9)
    axis.set_xticks(range(len(chosen)), labs, fontsize=8)
    axis.set_yticks(range(len(chosen)), labs, fontsize=8)
    fig.colorbar(im, label="symmetry-aware torsion distance")
    axis.set_title(f"{title}distance between starts", fontsize=10)
    fig.savefig(P("starts_distance.png"), dpi=150)
    plt.close(fig)
    print("written to", out)


if __name__ == "__main__":
    main()
