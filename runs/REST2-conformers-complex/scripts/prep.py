"""Complex inputs for the five REST2-conformers ligands, from the deposited crystal structures.

  python -I scripts/prep.py <root>

For each entry: chain A protein (altloc ' '/'A', altloc column blanked) + the chain-A copy of the
ligand, written to <root>/<name>/prep/complex.pdb. Crystallisation additives (DMS, GOL, SO4, NO3),
the chain-B copies and crystal waters are dropped. The ligand stays heavy-atom only: md-tools puts
the WATER-RUN package's hydrogens on it (the package is reused unchanged, see packages/).
Every chain-A ligand was checked to sit in chain A's own pocket (4IBI/4IBK: one chain-B residue
within 4 A, a crystal contact).
"""
import sys
from pathlib import Path

ROOT = Path(sys.argv[1]).resolve()
ENTRIES = {"4ibi": ("4IBI", "1D8"), "4ibk": ("4IBK", "1DE"), "4j7d": ("4J7D", "I31"),
           "5mwh": ("5MWH", "UWX"), "5ov8": ("5OV8", "AXN")}

for name, (pdb, lig) in ENTRIES.items():
    prot, het = [], []
    for line in (ROOT / "crystal" / f"{pdb}.pdb").read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")) or line[21] != "A" or line[16] not in " A":
            continue
        line = line[:16] + " " + line[17:]
        if line.startswith("ATOM"):
            prot.append(line)
        elif line[17:20] == lig:
            het.append(line)
    out = ROOT / name / "prep" / "complex.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join([f"REMARK   1 {pdb} chain A + {lig}; additives, chain B, HOH dropped"]
                             + prot + ["TER"] + het + ["END"]) + "\n")
    print(f"{name}: {out}  {len(prot)} protein atoms, {len(het)} {lig} atoms")
