"""Shu-Yu's MACE-OMOL bespoke folate System (folic_acid_mace_omol/ligand_sage_mace_system.xml)
vs the v4 folate in a built system (FOL-REST2-v4/water/build/built.xml): atom mapping by graph
match (elements, bond orders, charges), then every parameter class compared on mapped atoms."""
import sys
from collections import defaultdict
import numpy as np
import openmm as mm
from openmm import unit
from openmm.app import PDBFile
from rdkit import Chem

T = "/scratch/ipopa/thesis"
mace_sys = mm.XmlSerializer.deserialize(open(f"{T}/folic_acid_mace_omol/ligand_sage_mace_system.xml").read())
mace_mol = Chem.MolFromMolFile(f"{T}/folic_acid_mace_omol/molecule.sdf", removeHs=False)
v4_sys = mm.XmlSerializer.deserialize(open(f"{T}/FOL-REST2-v4/water/build/built.xml").read())
v4_top = PDBFile(f"{T}/FOL-REST2-v4/water/build/built.pdb").topology
v4_mol = Chem.MolFromMolFile(f"{T}/FOL-REST2-v4/water/build/FOL.sdf", removeHs=False)
lig = [a.index for a in v4_top.atoms() if a.residue.name == "FOL"]
assert lig == list(range(len(lig))), "FOL is expected first"
# sdf order == built order for the ligand (checked by element)
assert [a.GetSymbol() for a in v4_mol.GetAtoms()] == [a.element.symbol for a in v4_top.atoms()][:len(lig)]
matches = v4_mol.GetSubstructMatches(mace_mol, uniquify=False, useChirality=True, maxMatches=1000)
print(f"{mace_sys.getNumParticles()} MACE atoms, {len(lig)} v4 atoms, {len(matches)} graph matches")
m = matches[0]                                   # m[mace_idx] = v4_idx
inv = {v: k for k, v in enumerate(m)}

def forces(s):
    return {f.__class__.__name__: f for f in s.getForces()}
F, G = forces(mace_sys), forces(v4_sys)
print("MACE forces:", sorted(F)); print("v4 forces  :", sorted(G))

nbm, nbv = F["NonbondedForce"], G["NonbondedForce"]
dq, dsig, deps = [], [], []
for i in range(mace_sys.getNumParticles()):
    q1, s1, e1 = nbm.getParticleParameters(i); q2, s2, e2 = nbv.getParticleParameters(m[i])
    dq.append(q1.value_in_unit(unit.elementary_charge) - q2.value_in_unit(unit.elementary_charge))
    dsig.append(s1._value - s2._value); deps.append(e1._value - e2._value)
print(f"charges: max |dq| {np.max(np.abs(dq)):.4f} e, rms {np.sqrt(np.mean(np.square(dq))):.4f}; total MACE "
      f"{sum(nbm.getParticleParameters(i)[0]._value for i in range(49)):+.3f}")
print(f"LJ: max |dsigma| {np.max(np.abs(dsig)):.2e} nm, max |deps| {np.max(np.abs(deps)):.2e} kJ/mol")

def bonds(f, n, remap=None):
    out = {}
    for b in range(f.getNumBonds()):
        i, j, r, k = f.getBondParameters(b)
        if remap: i, j = remap[i], remap[j]
        if i < n and j < n: out[tuple(sorted((i, j)))] = (r._value, k._value)
    return out
bm = bonds(F["HarmonicBondForce"], 10**9, m); bv = bonds(G["HarmonicBondForce"], len(lig))
cons = {tuple(sorted(v4_sys.getConstraintParameters(c)[:2])) for c in range(v4_sys.getNumConstraints())}
common = set(bm) & set(bv)
print(f"bonds: MACE {len(bm)}, v4 {len(bv)} (+{len([c for c in cons if max(c) < len(lig)])} ligand constraints); "
      f"max |dr| {max(abs(bm[b][0]-bv[b][0]) for b in common):.2e}, max |dk| {max(abs(bm[b][1]-bv[b][1]) for b in common):.3g}")
def angles(f, n, remap=None):
    out = {}
    for a in range(f.getNumAngles()):
        i, j, k, t, kk = f.getAngleParameters(a)
        if remap: i, j, k = remap[i], remap[j], remap[k]
        if max(i, j, k) < n: out[(min(i, k), j, max(i, k))] = (t._value, kk._value)
    return out
am = angles(F["HarmonicAngleForce"], 10**9, m); av = angles(G["HarmonicAngleForce"], len(lig))
ca = set(am) & set(av)
print(f"angles: MACE {len(am)}, v4 {len(av)}, common {len(ca)}; max |dtheta| {max(abs(am[a][0]-av[a][0]) for a in ca):.2e}, "
      f"max |dk| {max(abs(am[a][1]-av[a][1]) for a in ca):.3g}")

def torsions(f, n, remap=None):
    out = defaultdict(list)
    for t in range(f.getNumTorsions()):
        i, j, k, l, per, ph, kk = f.getTorsionParameters(t)
        if remap: i, j, k, l = (remap[x] for x in (i, j, k, l))
        if max(i, j, k, l) >= n: continue
        key = (i, j, k, l) if j < k else (l, k, j, i)
        out[key].append((per, round(ph._value, 6), round(kk._value, 6)))
    return {k: sorted(v) for k, v in out.items()}
tm = torsions(F["PeriodicTorsionForce"], 10**9, m); tv = torsions(G["PeriodicTorsionForce"], len(lig))
diff = sorted(k for k in set(tm) | set(tv) if tm.get(k) != tv.get(k))
print(f"torsions: MACE {len(tm)} quadruplets, v4 {len(tv)}; {len(diff)} differ")
by_bond = defaultdict(list)
for k in diff: by_bond[tuple(sorted((k[1], k[2])))].append(k)
names = [a.name for a in v4_top.atoms()]
for b, ks in sorted(by_bond.items()):
    print(f"  central bond {names[b[0]]}-{names[b[1]]} ({b[0]}-{b[1]}): {len(ks)} quadruplet(s) differ, e.g. "
          f"{'-'.join(names[x] for x in ks[0])}: v4 {tv.get(ks[0])} -> MACE {tm.get(ks[0])}")
