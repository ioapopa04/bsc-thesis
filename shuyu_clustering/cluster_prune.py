#!/usr/bin/env python
"""
Terminal version of MD-tools' docs/tutorial/paracetamol/analysis/clustering.ipynb,
run on YOUR trajectory instead of the tutorial's runs.

Same three parts as the notebook:
  1. Clustering        t_hdbscan on the torsions (checked with allow_single_cluster True/False)
  2. Symmetry pruning  TorsionSymmetry: which clusters are the same conformer under a
                       permutation of chemically identical atoms -> merged groups
  3. Summary           per-group masses, circular torsion statistics, representative frames

Needs an MD-tools version that contains md_tools.analysis, plus scikit-learn, rdkit, pandas.

Example (paracetamol, state-0 solute trajectory):
  python para_clustering.py \
      --traj ~/PARA-REST2-long/REST2-run1/solute_state0_prod1.nc \
      --top  ~/PARA-REST2-long/build/built.solute.pdb \
      --sdf  ~/PARA-REST2-long/build/PAR.sdf \
      --out  ~/PARA-REST2-long/tutorial-clustering

Torsions (pick one way; indices must match the trajectory's atom numbering):
  --torsions-file torsions.txt         from select_torsions.py (any molecule, e.g. folic acid);
                                       add torsions_rigid.txt as well to include e.g. the amide
  --torsion name:i,j,k,l  (repeatable)  explicit, e.g. polar-H torsions such as C-O-H of a COOH
  neither                               the tutorial's paracetamol set (omega, aryl, phenol)
  --use-torsions 1,2,3,4                keep only these (0-based, in the order listed)

Folic acid example (solute trajectory, ligand indices):
  python para_clustering.py --traj REST2-run1/solute_state0_prod1.nc --top build/built.solute.pdb \
      --sdf build/FOL.sdf --torsions-file REST2-run1/torsions.txt --stride 5 --out fol-clustering

Labelling and the bioactive reference (all optional):
  --name 4ibi                     compound name for the figure titles
  --rest2-config REST2.config     number of replicas and tau_max for the overview title
  --ref-sdf crystal/4ibi_1D8_A.sdf  crystal ligand (repeatable, e.g. one per chain); its torsions
                                  are drawn as vertical lines. Matched to --sdf on heavy atoms by
                                  element + connectivity only, so torsions through a hydrogen get
                                  no line, and atoms the crystal cannot tell apart (the two O of a
                                  carboxylate) give one line per possible assignment.
  -> clusters_overview.png        molecule on the left, one row per torsion: the full distribution,
                                  then one panel per merged group
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
    from md_tools.analysis import (t_hdbscan, load_torsion_json, load_cluster_torsion,
                                   TorsionSymmetry, verify_atom_mapping, group_name,
                                   representative_by_vote, strip_nonpolar_hydrogens,
                                   draw_representative_structures)
except ImportError as e:
    sys.exit(f"md_tools.analysis not available ({e}).\n"
             "Your MD-tools version does not contain the analysis module - ask which version to install.")
from rdkit import Chem
import pandas as pd

DEFAULT_TORSIONS = ["omega:0,1,3,4", "aryl:1,3,4,5", "phenol:6,7,8,17"]

# One colour per TORSION, all distinct (tab10 without its grey, then tab20's light shades).
# MD-tools' own palette has 5 entries, so with 7 torsions two pairs shared a colour.
TORSION_PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
                   "#e377c2", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a",
                   "#ff9896", "#c5b0d5", "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5"]
REF_STYLES = ["-", "--", ":", "-."]


def circular_mean(values):
    return float(np.arctan2(np.sin(values).mean(), np.cos(values).mean()))


def circular_std(values):
    resultant = np.hypot(np.sin(values).mean(), np.cos(values).mean())
    return float(np.sqrt(-2.0 * np.log(max(resultant, 1e-12))))


def torsion_colours(n):
    if n > len(TORSION_PALETTE):
        sys.exit(f"{n} torsions but only {len(TORSION_PALETTE)} distinct colours defined")
    return TORSION_PALETTE[:n]


def load_torsions(files=None, specs=None, use=None):
    """name -> (i, j, k, l), from select_torsions.py files and/or name:i,j,k,l specs."""
    torsions = {}
    for path in (files or []):
        for line in open(path):
            if line.strip().startswith("#") or not line.strip():
                continue
            q = tuple(int(x) for x in line.split("#")[0].split()[:4])
            label = os.path.splitext(os.path.basename(path))[0]
            torsions[f"{label}_{len(torsions)}"] = q
    for spec in (specs or ([] if torsions else DEFAULT_TORSIONS)):
        name, idx = spec.split(":")
        torsions[name] = tuple(int(x) for x in idx.split(","))
    if use:
        keep_pos = [int(x) for x in use.split(",")]
        items = list(torsions.items())
        torsions = dict(items[i] for i in keep_pos)
    if not torsions:
        sys.exit("no torsions given")
    return torsions


def read_rest2_config(path):
    """(number_of_replicas, tau_max) from a REST2.config, or (None, None)."""
    if not path:
        return None, None
    import yaml
    with open(path) as f:
        rest2 = (yaml.safe_load(f) or {}).get("rest2", {})
    return rest2.get("number_of_replicas"), rest2.get("tau_max")


def simulated_ns(traj, stride=1):
    """Simulated time covered by a trajectory, from its own time stamps (ps).

    A strided load drops up to stride-1 frames at the end (60 ns at stride 5 ends at
    59.992 ns), so they are added back from the unstrided frame interval."""
    t = np.asarray(traj.time, dtype=float)
    if t.size < 2:
        return None
    dt = (t[1] - t[0]) / stride
    return float(t[-1] + (stride - 1) * dt) / 1000.0


def _heavy_graph(mol):
    """Heavy-atom copy with every bond single and no aromaticity: matches on element +
    connectivity only, so it does not depend on how the crystal file assigns bond orders."""
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    rw = Chem.RWMol()
    new = {}
    for i in heavy:
        atom = Chem.Atom(mol.GetAtomWithIdx(i).GetAtomicNum())
        atom.SetNoImplicit(True)
        new[i] = rw.AddAtom(atom)
    for b in mol.GetBonds():
        a, c = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if a in new and c in new:
            rw.AddBond(new[a], new[c], Chem.BondType.SINGLE)
    g = rw.GetMol()
    g.UpdatePropertyCache(strict=False)
    return g, heavy


def reference_torsions(ref_sdf, sim_mol, quads):
    """Torsions of a crystal ligand, in the simulation's atom numbering.

    Returns one list per torsion: the distinct values (rad) over every graph match, or [] when
    the torsion goes through a hydrogen (crystal structures carry none)."""
    ref = Chem.MolFromMolFile(ref_sdf, removeHs=False, sanitize=False)
    if ref is None:
        sys.exit(f"RDKit could not read {ref_sdf}")
    sim_g, sim_heavy = _heavy_graph(sim_mol)
    ref_g, ref_heavy = _heavy_graph(ref)
    if sim_g.GetNumAtoms() != ref_g.GetNumAtoms():
        sys.exit(f"{ref_sdf}: {ref_g.GetNumAtoms()} heavy atoms, simulation ligand has "
                 f"{sim_g.GetNumAtoms()}")
    matches = sim_g.GetSubstructMatches(ref_g, uniquify=False, useChirality=False,
                                        maxMatches=100000)
    if not matches:
        sys.exit(f"{ref_sdf}: heavy-atom graph does not match the simulation ligand")
    xyz = ref.GetConformer().GetPositions()
    sim_pos = {atom: k for k, atom in enumerate(sim_heavy)}

    # The graph match ignores chirality, so for e.g. two identical aryls on two stereocentres it
    # also returns the mirror-image assignment. Keep only matches under which every stereocentre
    # of the simulated ligand has the same handedness in the crystal.
    sim_xyz = sim_mol.GetConformer().GetPositions()
    probe = Chem.Mol(sim_mol)
    Chem.AssignStereochemistryFrom3D(probe)
    centres = []
    for c, _ in Chem.FindMolChiralCenters(probe, includeUnassigned=False):
        nbrs = [n.GetIdx() for n in probe.GetAtomWithIdx(c).GetNeighbors() if n.GetAtomicNum() > 1]
        if len(nbrs) >= 3:
            centres.append((c, nbrs[:3], _handedness(sim_xyz, c, nbrs[:3])))
    if centres:
        def keeps(m):
            inv = {s: r for r, s in enumerate(m)}
            pick = lambda a: ref_heavy[inv[sim_pos[a]]]
            return all(_handedness(xyz, pick(c), [pick(n) for n in nb]) == sign
                       for c, nb, sign in centres)
        kept = [m for m in matches if keeps(m)]
        if not kept:
            print(f"  WARNING {ref_sdf}: crystal ligand is not the simulated stereoisomer; "
                  f"reference torsions ignore chirality")
        else:
            matches = kept
    out = []
    for q in quads:
        if any(a not in sim_pos for a in q):
            out.append([])
            continue
        values = []
        for m in matches:            # m[ref_graph_idx] = sim_graph_idx
            inv = {s: r for r, s in enumerate(m)}
            p = xyz[[ref_heavy[inv[sim_pos[a]]] for a in q]]
            v = _dihedral(p)
            if all(abs(np.angle(np.exp(1j * (v - w)))) > np.deg2rad(2) for w in values):
                values.append(v)
        out.append(values)
    return out


def _handedness(xyz, centre, nbrs):
    """Sign of the volume spanned by three neighbours around a stereocentre."""
    v = xyz[nbrs] - xyz[centre]
    return bool(np.linalg.det(v) > 0)


def _dihedral(p):
    b0, b1, b2 = p[0] - p[1], p[2] - p[1], p[3] - p[2]
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def draw_reference_lines(axis, refs, column, with_label=False):
    """Vertical lines at the crystal value(s) of torsion `column`, one line style per file."""
    for k, (label, values) in enumerate(refs.items()):
        for n, v in enumerate(values[column]):
            axis.axvline(v, color="black", lw=1.3, ls=REF_STYLES[k % len(REF_STYLES)],
                         label=(f"crystal {label}" if with_label and n == 0 else None))


def molecule_png(mol, quads, colours, size=(900, 900)):
    """The ligand (non-polar H stripped) with each torsion's three bonds highlighted in its own
    colour. Bonds or atoms shared by several torsions are drawn in all their colours, so no
    torsion hides another; the central bond carries the torsion number."""
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    keep_h = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1
              and any(n.GetAtomicNum() not in (1, 6) for n in a.GetNeighbors())}
    keep_h |= {a for q in quads for a in q if mol.GetAtomWithIdx(a).GetAtomicNum() == 1}
    rw = Chem.RWMol(mol)
    for a in sorted((a.GetIdx() for a in mol.GetAtoms()
                     if a.GetAtomicNum() == 1 and a.GetIdx() not in keep_h), reverse=True):
        rw.RemoveAtom(a)
    draw = rw.GetMol()
    old_to_new, n = {}, 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 1 or a.GetIdx() in keep_h:
            old_to_new[a.GetIdx()] = n
            n += 1
    draw.UpdatePropertyCache(strict=False)
    AllChem.Compute2DCoords(draw)

    from matplotlib.colors import to_rgb
    atom_hl, bond_hl = {}, {}
    for k, (q, colour) in enumerate(zip(quads, colours)):
        rgb = to_rgb(colour)
        q = [old_to_new[a] for a in q]
        for a in q:
            atom_hl.setdefault(a, []).append(rgb)
        for x, y in zip(q, q[1:]):
            bond = draw.GetBondBetweenAtoms(x, y)
            if bond is not None:
                bond_hl.setdefault(bond.GetIdx(), []).append(rgb)
        centre = draw.GetBondBetweenAtoms(q[1], q[2])
        if centre is not None:
            note = centre.GetProp("bondNote") + "," if centre.HasProp("bondNote") else ""
            centre.SetProp("bondNote", f"{note}t{k}")
    drawer = rdMolDraw2D.MolDraw2DCairo(*size)
    opts = drawer.drawOptions()
    opts.annotationFontScale = 0.9
    opts.fillHighlights = False
    opts.continuousHighlight = False
    drawer.DrawMoleculeWithHighlights(draw, "", atom_hl, bond_hl, {}, {})
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def draw_overview(path, *, mol, names, quads, tors, aligned, masses, noise, refs, title):
    """One figure for the whole clustering: molecule on the left; per torsion row the full
    distribution, then each merged group, all in that torsion's colour, crystal lines on top."""
    import io
    import matplotlib.image as mpimg

    colours = torsion_colours(len(names))
    groups = list(aligned)
    ncol = 1 + len(groups)
    nrow = len(names)
    fig = plt.figure(figsize=(7.0 + 2.3 * ncol, 1.75 * nrow + 1.6), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=[7.0, 2.3 * ncol])
    mol_ax = fig.add_subplot(outer[0])
    mol_ax.imshow(mpimg.imread(io.BytesIO(molecule_png(mol, quads, colours)), format="png"))
    mol_ax.axis("off")
    mol_ax.set_title("t<k> marks the central bond of torsion k", fontsize=10)
    grid = outer[1].subgridspec(nrow, ncol)

    edges = np.linspace(-np.pi, np.pi, 49)
    for row, (name, colour) in enumerate(zip(names, colours)):
        for col in range(ncol):
            ax = fig.add_subplot(grid[row, col])
            if col == 0:
                values, head = tors[:, row], f"all frames (noise {100 * noise:.1f}%)"
            else:
                g = groups[col - 1]
                values, head = aligned[g][:, row], f"group {g} ({100 * float(masses[col - 1]):.1f}%)"
            ax.hist(values, bins=edges, density=True, color=colour,
                    alpha=0.9 if col == 0 else 0.75,
                    histtype="stepfilled", edgecolor=colour)
            draw_reference_lines(ax, refs, row, with_label=(row == 0 and col == 0))
            ax.set_xlim(-np.pi, np.pi)
            ax.set_xticks([-np.pi, 0, np.pi], ["-π", "0", "π"])
            ax.set_yticks([])
            ax.tick_params(labelsize=7)
            if col == 0:
                ax.set_ylabel(f"t{row}\n{name}", color=colour, fontsize=8, rotation=0,
                              ha="right", va="center")
                ax.set_facecolor("#f4f4f4")
            if row == 0:
                ax.set_title(head, fontsize=8)
            if row == 0 and col == 0 and refs:
                ax.legend(fontsize=6, loc="upper left", frameon=False)
    fig.suptitle(title, fontsize=12)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--traj", required=True, help="solute trajectory of the physical state (state 0)")
    ap.add_argument("--top", required=True, help="matching topology WITH bonds (e.g. built.solute.pdb)")
    ap.add_argument("--sdf", required=True, help="ligand SDF from the build (bond orders)")
    ap.add_argument("--torsion", action="append", default=None,
                    help="name:i,j,k,l (repeat); default = tutorial paracetamol set")
    ap.add_argument("--torsions-file", action="append", default=None,
                    help="torsions.txt from select_torsions.py (repeatable; lines 'i j k l ...  # names')")
    ap.add_argument("--use-torsions", default=None,
                    help="comma-separated positions in the final torsion list to keep")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--no-single-cluster", action="store_true",
                    help="run t_hdbscan with allow_single_cluster=False (MD-tools default is True)")
    ap.add_argument("--max-matches", type=int, default=200_000,
                    help="cap on the symmetry (automorphism) enumeration. MD-tools' default is "
                         "200000; molecules with several equivalent methyl/CF3 groups exceed it "
                         "(the group grows factorially with equivalent hydrogens) and "
                         "TorsionSymmetry refuses a truncated search. Raise it above the true "
                         "count, e.g. 8000000, without changing the result.")
    ap.add_argument("--name", default=None, help="compound name for the figure titles")
    ap.add_argument("--rest2-config", default=None,
                    help="REST2.config of the run: number of replicas and tau_max for the title")
    ap.add_argument("--duration-ns", type=float, default=None,
                    help="simulated time for the title (default: from the trajectory time stamps)")
    ap.add_argument("--ref-sdf", action="append", default=None,
                    help="crystal (bioactive) ligand SDF; repeatable (e.g. one per chain)")
    ap.add_argument("--out", default="tutorial-clustering")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    P = lambda name: os.path.join(out, name)

    torsions = load_torsions(args.torsions_file, args.torsion, args.use_torsions)
    names, quads = list(torsions), list(torsions.values())
    label = f"{args.name}: " if args.name else ""

    # ------------------------------------------------------------ 1. clustering
    print("== 1. Clustering ==")
    cmd = md.load(args.traj, top=args.top, stride=args.stride)
    for name, q in torsions.items():
        print(f"  {name:8s} {q}: " + " - ".join(str(cmd.topology.atom(i)) for i in q))
    tors = md.compute_dihedrals(cmd, quads)
    print(f"{tors.shape[0]} frames, {tors.shape[1]} torsions")

    asc = not args.no_single_cluster
    print(f"allow_single_cluster = {asc}")
    fit = t_hdbscan(tors, names=names, allow_single_cluster=asc).fit()
    summary = fit.summary()
    print(f"n_clusters = {summary['n_clusters']}")
    print(f"masses     = { {k: round(v, 3) for k, v in summary['cluster_mass'].items()} }")
    print(f"noise      = {summary['noise_mass']:.4f}")

    print("\nIs one state a possible answer? (allow_single_cluster True vs False)")
    for flag in (True, False):
        s = t_hdbscan(tors, names=names, allow_single_cluster=flag).fit().summary()
        print(f"  allow_single_cluster={flag!s:5s} -> {s['n_clusters']} clusters, "
              f"excluded_by_construction={s['single_cluster_excluded_by_construction']}")
    print("\nCaveats:")
    for caveat in summary["caveats"]:
        print("  -", caveat)

    # ------------------------------------------------------- 2. symmetry pruning
    print("\n== 2. Pruning by symmetry ==")
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    if mol is None:
        sys.exit(f"RDKit could not read {args.sdf}")

    refs = {}
    for path in (args.ref_sdf or []):
        ref_label = os.path.splitext(os.path.basename(path))[0]
        refs[ref_label] = reference_torsions(path, mol, quads)
        print(f"crystal reference {ref_label}:")
        for name, values in zip(names, refs[ref_label]):
            shown = ", ".join(f"{np.rad2deg(v):.1f}" for v in values) or "n/a (through H)"
            print(f"  {name:20s} {shown} deg")
    bonds = [(b[0].index, b[1].index) for b in cmd.topology.bonds]
    if not bonds:
        sys.exit(f"{args.top} has no bonds - use a topology with CONECT records / bonds")
    elements = [a.element.symbol for a in cmd.topology.atoms]
    mapping, why = verify_atom_mapping(mol, bonds, elements)
    print("atom mapping:", "identity" if mapping == list(range(mol.GetNumAtoms())) else mapping)
    print(" ", why)

    definition = {"schema_version": 1, "indexing": "zero-based", "units": "radian",
                  "torsions": [{"name": n, "atoms": list(q)} for n, q in torsions.items()]}
    tjson = P("torsion_cv.json")
    with open(tjson, "w") as f:
        json.dump(definition, f, indent=2)

    torsion_idx = load_torsion_json(tjson, molecule=mol)
    clusters = load_cluster_torsion(tors, fit.labels_, noise_label=-1, angle_unit="radian",
                                    names=list(torsion_idx.names))
    print("cluster labels:", clusters.cluster_labels, "| masses:", clusters.mass_table())

    detector = TorsionSymmetry(mol, torsion_idx, max_matches=args.max_matches)
    detector.fit(clusters, coordinates=cmd.xyz, atom_mapping=mapping)
    print("operations distinct on these torsions:", len(detector.operations))
    print("groups:", detector.cluster_groups)
    for c in detector.comparisons:
        if c.evidence.get("role") == "transitivity-guard":
            continue
        print(f"  {c.a} vs {c.b}: {c.verdict:9s} op={c.operation} "
              f"E={c.energy if c.energy is None else round(c.energy, 4)}")

    for g in range(len(detector.cluster_groups)):
        if len(detector.cluster_groups[g]) > 1:          # only merged groups have an operation to draw
            try:
                detector.draw_operation_effect(group=g, output=P(f"operation_group_{g}.png"))
            except Exception as e:
                print(f"  (could not draw operation for group {g}: {e})")
    for i in range(len(detector.cluster_groups)):
        try:
            detector.draw_sym_group(i, output=P(f"group_{i}.png"), cluster_torsions=clusters)
        except Exception as e:
            print(f"  (could not draw group {i}: {e})")

    merged = detector.deduplicate(clusters)
    desc = merged.describe()
    print(json.dumps(desc, indent=2))
    with open(P("merged_groups.json"), "w") as f:
        json.dump(desc, f, indent=2)
    total = sum(merged.masses.values()) + merged.noise_mass
    print(f"populations conserved: merged masses + noise = {total:.6f}")

    # ---------------------------------------------------------------- 3. summary
    print("\n== 3. Summary ==")
    reps = representative_by_vote(fit)
    groups = {group_name(i): g for i, g in enumerate(detector.cluster_groups)}
    chosen = {name: reps[members[0]] for name, members in groups.items()}
    aligned = {name: merged.aligned[i] for i, name in enumerate(groups)}
    for name, members in groups.items():
        print(f"group {name}: clusters {members}  mass {merged.masses[list(groups).index(name)]:.4f}"
              f"  representative frame {chosen[name]}")

    ncol = min(3, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.7 * ncol, 3.6 * nrow), layout="constrained", squeeze=False)
    for extra in axes.flat[len(names):]:
        extra.axis("off")
    edges = np.linspace(-np.pi, np.pi, 49)
    palette = plt.get_cmap("tab10")
    for column, (axis, torsion) in enumerate(zip(axes.flat, names)):
        for k, (name, values) in enumerate(aligned.items()):
            axis.hist(values[:, column], bins=edges, density=True, histtype="step",
                      lw=2.0, color=palette(k),
                      label=f"group {name} ({100 * float(merged.masses[k]):.1f}%)")
        draw_reference_lines(axis, refs, column, with_label=(column == 0))
        axis.set_xlim(-np.pi, np.pi)
        axis.set_xlabel(f"{torsion} (rad)")
        axis.set_title(torsion)
        if column == 0:
            axis.set_ylabel("density")
            axis.legend(fontsize=9)
    fig.suptitle(f"{label}symmetry-aligned torsion distributions, one colour per merged group"
                 + ("; black lines = crystal (bioactive) conformation" if refs else ""),
                 fontsize=11)
    fig.savefig(P("group_torsions.png"), dpi=150)

    # ------------------------------------------- 3a. everything in one figure
    n_rep, tau_max = read_rest2_config(args.rest2_config)
    duration = args.duration_ns if args.duration_ns is not None else simulated_ns(cmd, args.stride)
    run_info = [f"{duration:g} ns" if duration else None,
                f"{n_rep} replicas" if n_rep else None,
                f"τ_max {tau_max}" if tau_max is not None else None,
                f"stride {args.stride}, {tors.shape[0]} frames"]
    title = (f"{args.name or 'ligand'} — REST2 " + ", ".join(x for x in run_info if x)
             + f"\n{len(groups)} merged groups, noise {100 * merged.noise_mass:.1f}%"
             + ("; black lines = crystal (bioactive) conformation" if refs else ""))
    draw_overview(P("clusters_overview.png"), mol=mol, names=names, quads=quads, tors=tors,
                  aligned=aligned, masses=merged.masses, noise=merged.noise_mass,
                  refs=refs, title=title)

    # ------------------------------------------------- 3b. 2-D torsion marginals
    # The joint distribution of two chosen torsions, which a pair of 1-D histograms
    # cannot show: correlated torsions draw a ridge, independent ones a cross. Panel 0
    # is the whole ensemble (raw torsions); the rest are per merged group (symmetry-
    # aligned, same as the 1-D figure above), each labelled with its population.
    ix, iy = 3, 4                      # "torsion 3" and "torsion 4" (0-based positions)
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
            hb = axis.hexbin(xs, ys, gridsize=55, mincnt=1, cmap="viridis",
                             bins="log", extent=lim)
            fig2.colorbar(hb, ax=axis, label="log10(count)")
            axis.set_title(title)
            axis.set_xlabel(f"{names[ix]} (rad)")
            axis.set_ylabel(f"{names[iy]} (rad)")

        _hex2d(axes2.flat[0], tors[:, ix], tors[:, iy], "all frames")
        for k, (name, values) in enumerate(aligned.items()):
            _hex2d(axes2.flat[k + 1], values[:, ix], values[:, iy],
                   f"group {name} ({100 * float(merged.masses[k]):.1f}%)")
        fig2.suptitle(f"{label}joint distribution: {names[ix]} vs {names[iy]}", fontsize=11)
        fig2.savefig(P(f"torsions_{ix}_{iy}_2d.png"), dpi=150)

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

    keep = strip_nonpolar_hydrogens(cmd)
    try:
        info = draw_representative_structures(cmd, chosen, output=P("representatives.png"), atom_indices=keep)
        print(f"{info['n_atoms_drawn']} atoms drawn, superposed: {info['aligned']}")
    except Exception as e:
        print(f"  (could not draw representatives: {e})")
    for name, frame in chosen.items():
        cmd[frame].save_pdb(P(f"representative_{name}.pdb"))
    print(f"\nAll outputs in {out}/")


if __name__ == "__main__":
    main()
