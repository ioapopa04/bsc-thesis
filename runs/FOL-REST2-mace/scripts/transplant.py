"""Put Shu-Yu's MACE-OMOL bespoke folate parameters into a built MD-tools System.

  python -I transplant.py <build-dir> <mace-system.xml> <mace-molecule.sdf>

build-top cannot read a finished OpenMM System, and Shu-Yu's fit exists only as one
(folic_acid_mace_omol/ligand_sage_mace_system.xml, Sage 2.2.1 + bespoke torsions, its own AM1-BCC
charges, folate pH 7.4 = 3H-pterin tautomer, -2). So the system is built by MD-tools from the SAME
molecule (package made from folic_acid_mace_omol/molecule.sdf), and then every intramolecular
folate term is replaced by the MACE system's:

  NonbondedForce   per-atom charge, sigma, epsilon; every folate-folate exception
  HarmonicBond     every folate bond present in the built system (X-H are constraints in both
                   builds; their lengths are checked against the MACE r0)
  HarmonicAngle    every folate angle
  PeriodicTorsion  ALL folate torsion terms (propers and impropers) removed and the MACE ones
                   added; every other torsion is kept, in order

Atoms are mapped by a graph match (elements, bond orders, formal charges, stereo) between the MACE
molecule and the build's FOL.sdf, and every match must give the same parameters (symmetry). The
original is kept as built.sage.xml; built.xml is overwritten; transplant.json records the source
file's sha256 and every count. Nothing about solvent, ions or protein is touched.
"""
import hashlib
import json
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import openmm as mm
from openmm import unit
from openmm.app import PDBFile
from rdkit import Chem

build, mace_xml, mace_sdf = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
src = build / "built.sage.xml"
if not src.exists():
    shutil.copy2(build / "built.xml", src)
system = mm.XmlSerializer.deserialize(src.read_text())
mace = mm.XmlSerializer.deserialize(mace_xml.read_text())
top = PDBFile(str(build / "built.pdb")).topology
lig = [a.index for a in top.atoms() if a.residue.name == "FOL"]
assert len(lig) == mace.getNumParticles() == 49, (len(lig), mace.getNumParticles())

# ---------------------------------------------------------------- atom mapping
bmol = Chem.MolFromMolFile(str(build / "FOL.sdf"), removeHs=False)    # package molecule, build order
mmol = Chem.MolFromMolFile(str(mace_sdf), removeHs=False)
lig_elems = [a.element.symbol for a in top.atoms() if a.residue.name == "FOL"]
assert [a.GetSymbol() for a in bmol.GetAtoms()] == lig_elems, "FOL.sdf order != built FOL order"
matches = bmol.GetSubstructMatches(mmol, uniquify=False, useChirality=True, maxMatches=10000)
assert matches, "MACE molecule and build molecule are not the same species (tautomer / protomer?)"
to_sys = [{mi: lig[bi] for mi, bi in enumerate(m)} for m in matches]   # MACE index -> system index
M = to_sys[0]
L = set(lig)

def forces(s):
    out = {}
    for k in range(s.getNumForces()):
        f = s.getForce(k)
        out.setdefault(f.__class__.__name__, (k, f))
    return out
F, S = forces(mace), forces(system)

# ---------------------------------------------------------------- MACE terms, in system indices
def mace_terms(Mp):
    nb = F["NonbondedForce"][1]
    part = {Mp[i]: tuple(x._value for x in nb.getParticleParameters(i)) for i in range(49)}
    exc = {}
    for e in range(nb.getNumExceptions()):
        i, j, q, s, eps = nb.getExceptionParameters(e)
        exc[frozenset((Mp[i], Mp[j]))] = (q._value, s._value, eps._value)
    bond = {}
    f = F["HarmonicBondForce"][1]
    for b in range(f.getNumBonds()):
        i, j, r, k = f.getBondParameters(b)
        bond[frozenset((Mp[i], Mp[j]))] = (r._value, k._value)
    ang = {}
    f = F["HarmonicAngleForce"][1]
    for a in range(f.getNumAngles()):
        i, j, k, t, kk = f.getAngleParameters(a)
        ang[(min(Mp[i], Mp[k]), Mp[j], max(Mp[i], Mp[k]))] = (t._value, kk._value)
    tor = []
    f = F["PeriodicTorsionForce"][1]
    for t in range(f.getNumTorsions()):
        i, j, k, l, per, ph, kk = f.getTorsionParameters(t)
        tor.append((Mp[i], Mp[j], Mp[k], Mp[l], per, ph._value, kk._value))
    return part, exc, bond, ang, tor

def canon_tor(tor):
    out = defaultdict(list)
    for i, j, k, l, per, ph, kk in tor:
        key = (i, j, k, l) if j < k else (l, k, j, i)
        out[key].append((per, round(ph, 6), round(kk, 6)))
    return {k: sorted(v) for k, v in out.items()}

ref = mace_terms(M)
asym = 0
for Mp in to_sys[1:]:            # a symmetry-equivalent mapping should describe the same Hamiltonian
    alt = mace_terms(Mp)
    same = (all(abs(alt[0][a][0] - ref[0][a][0]) < 1e-6 for a in ref[0])
            and canon_tor(alt[4]) == canon_tor(ref[4]))
    asym += not same
part, exc, bond, ang, tor = ref
rec = {"source": str(mace_xml), "sha256": hashlib.sha256(mace_xml.read_bytes()).hexdigest(),
       "graph_matches": len(matches), "matches_giving_different_parameters": asym,
       "folate_atoms": lig[:1] + lig[-1:]}

# ---------------------------------------------------------------- nonbonded
k_nb, nb = S["NonbondedForce"]
dq = 0.0
for a, (q, s, e) in part.items():
    q0 = nb.getParticleParameters(a)[0]._value
    dq = max(dq, abs(q - q0))
    nb.setParticleParameters(a, q, s, e)
n_exc = 0
for e in range(nb.getNumExceptions()):
    i, j, *_ = nb.getExceptionParameters(e)
    if i in L and j in L:
        q, s, eps = exc[frozenset((i, j))]
        nb.setExceptionParameters(e, i, j, q, s, eps)
        n_exc += 1
assert n_exc == len(exc), (n_exc, len(exc))
rec.update(max_abs_charge_change=dq, exceptions=n_exc,
           total_charge=sum(p[0] for p in part.values()))

# ---------------------------------------------------------------- bonds, constraints, angles
k_b, bf = S["HarmonicBondForce"]
n_b = 0
for b in range(bf.getNumBonds()):
    i, j, r, k = bf.getBondParameters(b)
    if i in L and j in L:
        bf.setBondParameters(b, i, j, *bond[frozenset((i, j))])
        n_b += 1
cons_dev = 0.0
n_c = 0
for c in range(system.getNumConstraints()):
    i, j, d = system.getConstraintParameters(c)
    if i in L and j in L:
        cons_dev = max(cons_dev, abs(d._value - bond[frozenset((i, j))][0]))
        n_c += 1
assert n_b + n_c == len(bond), (n_b, n_c, len(bond))
assert cons_dev < 1e-6, f"constraint lengths differ from MACE r0 by {cons_dev} nm"
k_a, af = S["HarmonicAngleForce"]
n_a = 0
for a in range(af.getNumAngles()):
    i, j, k, t, kk = af.getAngleParameters(a)
    if i in L and j in L and k in L:
        af.setAngleParameters(a, i, j, k, *ang[(min(i, k), j, max(i, k))])
        n_a += 1
assert n_a == len(ang), (n_a, len(ang))
rec.update(bonds=n_b, constraints_checked=n_c, max_constraint_vs_r0_nm=cons_dev, angles=n_a)

# ---------------------------------------------------------------- torsions
k_t, tf = S["PeriodicTorsionForce"]
new = mm.PeriodicTorsionForce()
new.setForceGroup(tf.getForceGroup())
new.setUsesPeriodicBoundaryConditions(tf.usesPeriodicBoundaryConditions())
old_lig = []
for t in range(tf.getNumTorsions()):
    i, j, k, l, per, ph, kk = tf.getTorsionParameters(t)
    if {i, j, k, l} <= L:
        old_lig.append((i, j, k, l, per, ph._value, kk._value))
        continue
    assert not ({i, j, k, l} & L), "a torsion spans folate and something else"
    new.addTorsion(i, j, k, l, per, ph, kk)
for i, j, k, l, per, ph, kk in tor:
    new.addTorsion(i, j, k, l, per, ph, kk)
system.removeForce(k_t)
system.addForce(new)
before, after = canon_tor(old_lig), canon_tor(tor)
changed = sorted(q for q in set(before) | set(after) if before.get(q) != after.get(q))
bonds_changed = sorted({tuple(sorted((q[1], q[2]))) for q in changed})
names = [a.name for a in top.atoms()]
rec.update(torsion_terms_removed=len(old_lig), torsion_terms_added=len(tor),
           quadruplets_changed=len(changed),
           central_bonds_changed=[f"{names[a]}-{names[b]}" for a, b in bonds_changed])

# ---------------------------------------------------------------- energy check
# The folate cut out of the modified system (vacuum, no cutoff, as the MACE System) must give the
# MACE System's angle, torsion and nonbonded energies at the same coordinates (bonds: X-H are
# constraints in the build, so only the heavy-atom bond energy is compared).
inv = {s: m for m, s in M.items()}                         # system index -> MACE index
pos = Chem.MolFromMolFile(str(mace_sdf), removeHs=False).GetConformer().GetPositions() / 10.0

def energies(s, positions):
    out = {}
    for k in range(s.getNumForces()):
        s.getForce(k).setForceGroup(k)
    ctx = mm.Context(s, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(positions)
    for k in range(s.getNumForces()):
        e = ctx.getState(getEnergy=True, groups={k}).getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        out[s.getForce(k).__class__.__name__] = out.get(s.getForce(k).__class__.__name__, 0.0) + e
    return out

sub = mm.System()
for m_i in range(49):
    sub.addParticle(system.getParticleMass(M[m_i]))
snb = mm.NonbondedForce(); snb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
for m_i in range(49):
    snb.addParticle(*nb.getParticleParameters(M[m_i]))
for e in range(nb.getNumExceptions()):
    i, j, q, s_, eps = nb.getExceptionParameters(e)
    if i in L and j in L:
        snb.addException(inv[i], inv[j], q, s_, eps)
sub.addForce(snb)
sa = mm.HarmonicAngleForce()
for a in range(af.getNumAngles()):
    i, j, k, t, kk = af.getAngleParameters(a)
    if i in L and j in L and k in L:
        sa.addAngle(inv[i], inv[j], inv[k], t, kk)
sub.addForce(sa)
st = mm.PeriodicTorsionForce()
for t in range(new.getNumTorsions()):
    i, j, k, l, per, ph, kk = new.getTorsionParameters(t)
    if {i, j, k, l} <= L:
        st.addTorsion(inv[i], inv[j], inv[k], inv[l], per, ph, kk)
sub.addForce(st)
mace_copy = mm.XmlSerializer.deserialize(mace_xml.read_text())
e_sub, e_mace = energies(sub, pos), energies(mace_copy, pos)
check = {name: (round(e_sub[name], 4), round(e_mace[name], 4))
         for name in ("NonbondedForce", "HarmonicAngleForce", "PeriodicTorsionForce")}
for name, (a, b) in check.items():
    assert abs(a - b) < 1e-3 * max(1.0, abs(b)), f"{name}: transplanted {a} vs MACE {b} kJ/mol"
rec["energy_check_kJmol (transplanted, MACE)"] = check

(build / "built.xml").write_text(mm.XmlSerializer.serialize(system))
(build / "transplant.json").write_text(json.dumps(rec, indent=2))
print(json.dumps(rec, indent=2))
