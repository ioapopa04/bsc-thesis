"""Health check of one complex REST2 run: exchange acceptance and whether folate stays bound.

  python -I scripts/check_run.py <run-dir> [--states 0,7]

  acceptance   per neighbouring pair, from remd_records/restart_prod1.json (completion_report);
               while the run is still going, from the live REST2.runstate.json if it has counts
  binding      for each listed state: folate heavy-atom RMSD to the crystal pose after
               superposing the protein CA atoms (no ligand fit, so drift out of the pocket shows),
               and the folate centre-of-mass displacement. The reference is build/built.solute.pdb,
               whose ligand heavy atoms are the deposited coordinates.
"""
import argparse
import json
from pathlib import Path

import mdtraj as md
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("run")
ap.add_argument("--states", default="0,7")
ap.add_argument("--ligand", default="FOL")
args = ap.parse_args()
run = Path(args.run).resolve()

rep = run / "remd_records/restart_prod1.json"
if rep.exists() and "completion_report" in json.loads(rep.read_text()):
    cr = json.loads(rep.read_text())["completion_report"]
    print("acceptance (completed run)")
    for p in cr["by_neighbouring_pair"]:
        flag = "" if p["acceptance"] >= 0.2 else "   <-- below 0.2"
        print(f"  {p['state_pair']}  tau {p['tau_pair'][0]:.3f}-{p['tau_pair'][1]:.3f}  "
              f"{p['accepted']:5d}/{p['proposed']:5d}  {p['acceptance']:.3f}{flag}")
    ov = cr.get("overall")
    print(f"  overall {ov.get('acceptance', ov) if isinstance(ov, dict) else ov}")
else:
    print("acceptance: no completion_report yet")

top = run.parent / "build/built.solute.pdb"
ref = md.load(str(top))
ca = ref.topology.select("protein and name CA")
# by name, not the selection language: it cannot parse residue names that start with a digit (1D8)
lig = np.array([a.index for a in ref.topology.atoms
                if a.residue.name == args.ligand and a.element.symbol != "H"])
print(f"\nbinding: {len(lig)} {args.ligand} heavy atoms, superposed on {len(ca)} CA")
for s in map(int, args.states.split(",")):
    nc = run / f"solute_state{s}_prod1.nc"
    if not nc.exists():
        print(f"  state {s}: no {nc.name}")
        continue
    t = md.load(str(nc), top=str(top))
    t.superpose(ref, atom_indices=ca)
    d = t.xyz[:, lig] - ref.xyz[0, lig]
    rmsd = np.sqrt((d ** 2).sum(-1).mean(-1)) * 10
    com = np.linalg.norm(t.xyz[:, lig].mean(1) - ref.xyz[0, lig].mean(0), axis=-1) * 10
    ca_rmsd = md.rmsd(t, ref, atom_indices=ca) * 10
    q = lambda a: f"mean {a.mean():.2f}  last {a[-1]:.2f}  max {a.max():.2f}"
    print(f"  state {s} ({t.n_frames} frames, {t.time[-1] - t.time[0]:.0f} ps)")
    print(f"    {args.ligand} RMSD to crystal (A)  {q(rmsd)}")
    print(f"    {args.ligand} COM shift (A)        {q(com)}")
    print(f"    protein CA RMSD (A)        {q(ca_rmsd)}")
