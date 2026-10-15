#!/usr/bin/env python
"""
Shu-Yu's SYMMETRY-FIRST torsion clustering (MD-tools main, 69c262b; tutorial
docs/tutorial/paracetamol/analysis/symmetry-alignment.ipynb), run on YOUR trajectory, with this
thesis' figures (clusters_overview, group_torsions, 2-D torsion plots, representatives).

Three parts, all MD-tools calls as in the tutorial:
  1. Clustering   t_hdbscan(tors, molecule, torsion_definitions, coordinates, atom_mapping):
                  the molecular symmetry is enumerated FIRST and built into the distance
                  (min over relabellings), then HDBSCAN on that dense matrix, then every frame
                  classified by an 18-of-20 neighbour vote. No merge step: a cluster is final.
                  Labels: -2 density noise, -1 ambiguous vote.
  2. Alignment    fit.align_all(per_cluster=True): every frame written in one labelling for the
                  figures (each cluster aligned to its own density peak, each peak to the density
                  peak of cluster A). Changes no label and no population.
  3. Summary      fit.representatives() (highest vote, then symmetry-aware medoid), circular
                  statistics on the aligned torsions, figures.

The previous cluster-then-merge version (t_hdbscan + TorsionSymmetry) is in git history.
Needs MD-tools main with the symmetry-first analysis (env analysis-env-main), scikit-learn, rdkit, pandas.

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
  -> clusters_overview_unaligned.png  the same layout with each cluster's torsions AS STORED (not
                                  relabelled): one conformer in two labellings shows both images
  -> clusters_overview.png        molecule on the left; to its right one column per torsion: top row
                                  the full distribution, then one row per cluster; many
                                  torsions wrap into 2-3 bands, many groups into 2-3 side-by-side
                                  blocks, to keep the figure ~16:9
"""
import argparse
import importlib.metadata
import json
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mdtraj as md

try:
    import md_tools
    from md_tools.analysis import (t_hdbscan, TorsionDefinitions, verify_atom_mapping, group_name,
                                   strip_nonpolar_hydrogens, draw_representative_structures)
    from md_tools.analysis import THDBSCAN
    if not hasattr(THDBSCAN, "align_all"):
        raise ImportError("this md_tools has no symmetry-first THDBSCAN.align_all")
except ImportError as e:
    sys.exit(f"md_tools.analysis not available ({e}).\n"
             "Needs Shu-Yu's MD-tools main (symmetry-first t_hdbscan): env analysis-env-main.")
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


def draw_overview(path, *, mol, names, quads, tors, aligned, masses, noise, refs, title,
                  row_labels=None):
    """One figure for the whole clustering: molecule on the left; to its right one column per
    torsion and one row per population -- the full distribution first, then each merged group --
    all in the torsion's colour, crystal lines on top. The grid is wrapped to keep the whole
    figure near 16:9: many torsions with few groups -> the torsions wrap into 1-3 bands stacked
    vertically (each band repeats the population rows); many groups with few torsions -> the
    populations wrap into 1-3 blocks side by side (each block repeats the torsion columns)."""
    import io
    import matplotlib.image as mpimg

    colours = torsion_colours(len(names))
    groups = list(aligned)
    nrow = 1 + len(groups)
    ntor = len(names)
    gap, target = 0.5, 16 / 9

    def size(nband, nblock, panel_h, panel_w=2.4):
        ncol = -(-ntor // nband)            # torsion columns per band
        rows = -(-nrow // nblock)           # population rows per block
        grid_w = nblock * (1.3 + panel_w * ncol) + (nblock - 1) * gap
        grid_h = nband * panel_h * rows + (nband - 1) * gap
        mol_w = min(max(4.5, 0.2 * grid_w), max(3.0, 0.9 * grid_h))   # larger molecules get more room
        return ncol, rows, grid_w, grid_h, mol_w

    def fit(nband, nblock):
        """Panel height (1.3-2.6 in) and, if the tallest panels are still too wide a figure, a
        narrower panel width (2.4 down to 1.6 in), that bring this split closest to 16:9."""
        panel_h, panel_w = 1.55, 2.4
        for _ in range(6):   # mol_w depends on the grid size, so settle in a few passes
            ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
            want_h = (mol_w + grid_w) / target - 1.0
            panel_h = float(np.clip((want_h - (nband - 1) * gap) / (nband * rows), 1.3, 2.6))
            ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
            if (mol_w + grid_w) / (grid_h + 1.0) > target:   # still too wide: narrow the panels
                want_grid_w = target * (grid_h + 1.0) - mol_w
                per_block = (want_grid_w - (nblock - 1) * gap) / nblock - 1.3
                panel_w = float(np.clip(per_block / ncol, 1.6, 2.4))
        ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
        return panel_h, panel_w, (mol_w + grid_w) / (grid_h + 1.0)

    # one split direction at a time; ties go to the simpler layout (fewer splits come first)
    splits = [(1, 1)] + [(b, 1) for b in (2, 3) if b <= ntor] + [(1, k) for k in (2, 3) if k <= nrow]
    nband, nblock = min(splits, key=lambda bk: round(abs(np.log(fit(*bk)[2] / target)), 3))
    panel_h, panel_w, _ = fit(nband, nblock)
    ncol, rows, grid_w, grid_h, mol_w = size(nband, nblock, panel_h, panel_w)
    # text scales with the figure: unchanged up to ~20 in wide, larger on big grids so it stays
    # readable when the figure is shrunk to page width
    fs = max(1.0, (mol_w + grid_w) / 20.0)
    fig = plt.figure(figsize=(mol_w + grid_w, grid_h + 1.0), layout="constrained")
    outer = fig.add_gridspec(1, 2, width_ratios=[mol_w, grid_w])
    mol_ax = fig.add_subplot(outer[0])
    png = molecule_png(mol, quads, colours, size=(int(220 * mol_w), int(220 * mol_w)))
    mol_ax.imshow(mpimg.imread(io.BytesIO(png), format="png"))
    mol_ax.axis("off")
    mol_ax.set_title("t<k> marks the central bond of torsion k", fontsize=10 * fs)
    cells = outer[1].subgridspec(nband, nblock, hspace=gap / (panel_h * rows),
                                 wspace=gap / (1.3 + panel_w * ncol))
    grids = [[cells[b, k].subgridspec(rows, ncol) for k in range(nblock)] for b in range(nband)]
    labels = [f"all frames\n(noise {100 * noise:.1f}%)"] + (list(row_labels) if row_labels else [
        f"group {g}\n({100 * float(masses[i]):.1f}%)" for i, g in enumerate(groups)])

    edges = np.linspace(-np.pi, np.pi, 49)
    for t, (name, colour) in enumerate(zip(names, colours)):
        band, col = divmod(t, ncol)
        for pop in range(nrow):
            block, row = divmod(pop, rows)
            ax = fig.add_subplot(grids[band][block][row, col])
            values = tors[:, t] if pop == 0 else aligned[groups[pop - 1]][:, t]
            ax.hist(values, bins=edges, density=True, color=colour,
                    alpha=0.9 if pop == 0 else 0.75,
                    histtype="stepfilled", edgecolor=colour)
            draw_reference_lines(ax, refs, t, with_label=(pop == 0 and t == 0))
            ax.set_xlim(-np.pi, np.pi)
            ax.set_xticks([-np.pi, 0, np.pi], ["-π", "0", "π"])
            ax.set_yticks([])
            ax.tick_params(labelsize=7 * fs)
            if row == 0:
                ax.set_title(f"t{t}\n{name}", color=colour, fontsize=8 * fs)
            if pop == 0:
                ax.set_facecolor("#f4f4f4")
            if col == 0:
                ax.set_ylabel(labels[pop], fontsize=8 * fs, rotation=0, ha="right", va="center")
            if pop == 0 and t == 0 and refs:
                ax.legend(fontsize=6 * fs, loc="upper left", frameon=False)
    fig.suptitle(title, fontsize=12 * fs)
    fig.savefig(path, dpi=130)
    plt.close(fig)


def build_fit(tors, cmd, mol, torsions, *, allow_single_cluster=True, max_frames=None):
    """Shu-Yu's symmetry-first t_hdbscan (MD-tools main, docs/tutorial/paracetamol/analysis/
    symmetry-alignment.ipynb): TorsionDefinitions on the build SDF, the atom mapping checked
    against the trajectory topology, and the trajectory itself so the symmetry images of the
    selected torsions can be recomputed. Returns (fit, mapping, why)."""
    definitions = TorsionDefinitions(quads=tuple(torsions.values()), names=tuple(torsions),
                                     indexing="zero-based", units="radian")
    bonds = [(b[0].index, b[1].index) for b in cmd.topology.bonds]
    if not bonds:
        sys.exit("the topology has no bonds - use one with CONECT records / bonds")
    mapping, why = verify_atom_mapping(mol, bonds, [a.element.symbol for a in cmd.topology.atoms])
    extra = {} if max_frames is None else {"max_precomputed_frames": int(max_frames)}
    fit = t_hdbscan(tors, molecule=mol, torsion_definitions=definitions, coordinates=cmd,
                    atom_mapping=mapping, allow_single_cluster=allow_single_cluster, **extra).fit()
    return fit, mapping, why


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
                    help="allow_single_cluster=False (Shu-Yu's default, and this script's, is True)")
    ap.add_argument("--max-frames", type=int, default=None,
                    help="raise t_hdbscan's max_precomputed_frames (default 10000): the symmetry-aware "
                         "distance is a dense N x N matrix, 8 N^2 bytes per copy, about three copies")
    ap.add_argument("--name", default=None, help="compound name for the figure titles")
    ap.add_argument("--rest2-config", default=None,
                    help="REST2.config of the run: number of replicas and tau_max for the title")
    ap.add_argument("--duration-ns", type=float, default=None,
                    help="simulated time for the title (default: from the trajectory time stamps)")
    ap.add_argument("--ref-sdf", action="append", default=None,
                    help="crystal (bioactive) ligand SDF; repeatable (e.g. one per chain)")
    ap.add_argument("--pair-2d", action="append", default=None, metavar="X,Y",
                    help="draw the joint distribution of torsion Y (vertical axis) against torsion X "
                         "(horizontal), 0-based positions; repeatable. Default: 3,4")
    ap.add_argument("--out", default="tutorial-clustering")
    args = ap.parse_args()
    out = os.path.abspath(args.out)
    os.makedirs(out, exist_ok=True)
    P = lambda name: os.path.join(out, name)

    torsions = load_torsions(args.torsions_file, args.torsion, args.use_torsions)
    names, quads = list(torsions), list(torsions.values())
    label = f"{args.name}: " if args.name else ""
    print("md-tools", importlib.metadata.version("md-tools"), "from",
          os.path.dirname(md_tools.__file__))

    # ------------------------------------------------- 1. symmetry-first clustering
    print("== 1. Symmetry-first clustering (Shu-Yu, MD-tools main) ==")
    cmd = md.load(args.traj, top=args.top, stride=args.stride)
    for name, q in torsions.items():
        print(f"  {name:8s} {q}: " + " - ".join(str(cmd.topology.atom(i)) for i in q))
    tors = md.compute_dihedrals(cmd, quads)
    print(f"{tors.shape[0]} frames, {tors.shape[1]} torsions")
    mol = Chem.SDMolSupplier(args.sdf, removeHs=False)[0]
    if mol is None:
        sys.exit(f"RDKit could not read {args.sdf}")
    with open(P("torsion_cv.json"), "w") as f:
        json.dump({"schema_version": 1, "indexing": "zero-based", "units": "radian",
                   "torsions": [{"name": n, "atoms": list(q)} for n, q in torsions.items()]},
                  f, indent=2)

    asc = not args.no_single_cluster
    fit, mapping, why = build_fit(tors, cmd, mol, torsions, allow_single_cluster=asc,
                                  max_frames=args.max_frames)
    print("atom mapping:", "identity" if mapping == list(range(mol.GetNumAtoms())) else mapping)
    print(" ", why)
    summary = fit.summary()
    desc = fit.descriptor_
    sym = summary["symmetry"]
    print(f"allow_single_cluster = {asc}; {summary['n_fit_frames']} frames fitted, "
          f"min_cluster_size {summary['settings']['min_cluster_size']}, "
          f"min_samples {summary['settings']['min_samples']}")
    print(f"symmetry operations acting on the torsions: {len(desc.operations)}")
    print("descriptor columns (selected + symmetry images):")
    for n, q, w, o in zip(desc.names, desc.quads, desc.column_weights, desc.origin):
        origin = "selected" if o["kind"] == "selected" else f"image of {o['image_of']}"
        print(f"  {n:28s} {'-'.join(map(str, q)):14s} weight {w:.3f}  {origin}")
    masses = {i: float(fit.cluster_mass_[i]) for i in range(fit.n_clusters_)}
    noise = float(fit.unassigned_mass_)
    print(f"n_clusters = {fit.n_clusters_}")
    print(f"masses     = { {group_name(k): round(v, 4) for k, v in masses.items()} }")
    print(f"density noise (-2) = {fit.noise_mass_:.4f}, ambiguous vote (-1) = "
          f"{fit.ambiguous_mass_:.4f}, unassigned total = {noise:.4f}")
    print("Caveats:")
    for caveat in summary["caveats"]:
        print("  -", caveat)
    with open(P("summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=str)

    refs = {}
    for path in (args.ref_sdf or []):
        ref_label = os.path.splitext(os.path.basename(path))[0]
        refs[ref_label] = reference_torsions(path, mol, quads)
        print(f"crystal reference {ref_label}:")
        for name, values in zip(names, refs[ref_label]):
            shown = ", ".join(f"{np.rad2deg(v):.1f}" for v in values) or "n/a (through H)"
            print(f"  {name:20s} {shown} deg")

    # ------------------------- 2. one labelling for the figures (two-level Dirichlet alignment)
    print("\n== 2. Alignment: align_all(per_cluster=True) ==")
    al = fit.align_all(per_cluster=True)
    groups = {group_name(c): c for c in range(fit.n_clusters_)}
    aligned = {g: al["selected"][fit.labels_ == c] for g, c in groups.items()}
    raw = {g: tors[fit.labels_ == c] for g, c in groups.items()}
    print(f"reference: frame {al['reference']} (density peak of cluster A)")
    for g, c in groups.items():
        m = fit.labels_ == c
        print(f"  cluster {g}: density peak frame {int(al['cluster_references'][c])}, "
              f"relabelled {100 * float((al['operations'][m] != 0).mean()):.1f}%, "
              f"margin < 0.1 for {100 * float((al['margin'][m] < 0.1).mean()):.1f}%")
    reps = fit.representatives()
    chosen = {g: int(reps[c]) for g, c in groups.items()}

    # ---------------------------------------------------------------- 3. summary + figures
    print("\n== 3. Summary ==")
    for i, (g, c) in enumerate(groups.items()):
        print(f"group {g}: cluster {c}  mass {masses[i]:.4f}  representative frame {chosen[g]}")
    group_list = list(groups)
    with open(P("merged_groups.json"), "w") as f:     # same file name as before, for downstream
        json.dump({"route": summary["route"], "groups": {g: [c] for g, c in groups.items()},
                   "masses": {g: masses[i] for i, g in enumerate(group_list)},
                   "noise_mass": float(fit.noise_mass_), "ambiguous_mass": float(fit.ambiguous_mass_),
                   "representatives": chosen, "alignment_reference": int(al["reference"]),
                   "cluster_density_peaks": {g: int(al["cluster_references"][c])
                                             for g, c in groups.items()}}, f, indent=2)

    ncol = min(3, len(names))
    nrow = int(np.ceil(len(names) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.7 * ncol, 3.6 * nrow), layout="constrained", squeeze=False)
    for extra in axes.flat[len(names):]:
        extra.axis("off")
    edges = np.linspace(-np.pi, np.pi, 49)
    palette = plt.get_cmap("tab10")
    for column, (axis, torsion) in enumerate(zip(axes.flat, names)):
        for k, (g, values) in enumerate(aligned.items()):
            axis.hist(values[:, column], bins=edges, density=True, histtype="step",
                      lw=2.0, color=palette(k % 10), label=f"group {g} ({100 * masses[k]:.1f}%)")
        draw_reference_lines(axis, refs, column, with_label=(column == 0))
        axis.set_xlim(-np.pi, np.pi)
        axis.set_xlabel(f"{torsion} (rad)")
        axis.set_title(torsion)
        if column == 0:
            axis.set_ylabel("density")
            axis.legend(fontsize=9)
    fig.suptitle(f"{label}symmetry-aligned torsion distributions (two-level alignment), one colour "
                 "per cluster" + ("; black lines = crystal (bioactive) conformation" if refs else ""),
                 fontsize=11)
    fig.savefig(P("group_torsions.png"), dpi=150)
    plt.close(fig)

    n_rep, tau_max = read_rest2_config(args.rest2_config)
    duration = args.duration_ns if args.duration_ns is not None else simulated_ns(cmd, args.stride)
    run_info = [f"{duration:g} ns" if duration else None,
                f"{n_rep} replicas" if n_rep else None,
                f"τ_max {tau_max}" if tau_max is not None else None,
                f"stride {args.stride}, {tors.shape[0]} frames"]
    head = f"{args.name or 'ligand'} — REST2 " + ", ".join(x for x in run_info if x)
    crystal = "; black lines = crystal (bioactive) conformation" if refs else ""
    unassigned = (f"unassigned {100 * noise:.1f}% (density noise {100 * fit.noise_mass_:.1f}%, "
                  f"ambiguous {100 * fit.ambiguous_mass_:.1f}%)")
    draw_overview(P("clusters_overview.png"), mol=mol, names=names, quads=quads, tors=tors,
                  aligned=aligned, masses=masses, noise=noise, refs=refs,
                  title=f"{head}\n{len(groups)} symmetry-first clusters, {unassigned}{crystal}")

    # The same overview with each cluster's torsions AS STORED (no relabelling): a cluster that
    # holds one conformer in two labellings shows both images here (e.g. peaks shifted by pi).
    draw_overview(P("clusters_overview_unaligned.png"), mol=mol, names=names, quads=quads,
                  tors=tors, aligned=raw, masses=masses, noise=noise, refs=refs,
                  title=f"{head}\n{len(groups)} symmetry-first clusters, torsions AS STORED "
                        f"(not relabelled), {unassigned}{crystal}")

    # ------------------------------------------------- 2-D torsion marginals
    pairs = [tuple(int(v) for v in p.split(",")) for p in (args.pair_2d or ["3,4"])]
    for ix, iy in pairs:               # 0-based torsion positions: x axis, y axis
        if len(names) <= max(ix, iy):
            print(f"  (skipping 2-D plot {ix},{iy}: only {len(names)} torsions)")
            continue
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
        for k, (g, values) in enumerate(aligned.items()):
            _hex2d(axes2.flat[k + 1], values[:, ix], values[:, iy],
                   f"group {g} ({100 * masses[k]:.1f}%)")
        fig2.suptitle(f"{label}joint distribution: {names[iy]} as a function of {names[ix]}", fontsize=11)
        fig2.savefig(P(f"torsions_{ix}_{iy}_2d.png"), dpi=150)
        plt.close(fig2)

    rows = []
    for i, (g, c) in enumerate(groups.items()):
        values = aligned[g]
        row = {"cluster_id": g, "original_clusters": str(c)}
        for column, torsion in enumerate(names):
            row[f"{torsion}_avg"] = np.rad2deg(circular_mean(values[:, column]))
            row[f"{torsion}_std"] = np.rad2deg(circular_std(values[:, column]))
        row["n_samples"] = int(values.shape[0])
        row["mass"] = masses[i]
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
    for g, frame in chosen.items():
        cmd[frame].save_pdb(P(f"representative_{g}.pdb"))
    print(f"\nAll outputs in {out}/")


if __name__ == "__main__":
    main()
