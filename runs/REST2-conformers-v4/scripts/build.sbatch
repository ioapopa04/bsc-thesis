#!/bin/bash
#SBATCH --job-name=rc4-build
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=/scratch/ipopa/thesis/REST2-conformers-v4/logs/build-%A_%a.out
# REST2 v4 (Shu-Yu's protonation branch, env openmm-env-protonation) for the three small molecules
# whose protomer Open Babel changes at pH 7.4 (REST2-conformers-complex/protonation-check):
#   0 4ibi (L01, -2)   1 4ibk (L02, -2)   2 4j7d (L03, +1)
# One new package per compound, shared by the water and the complex build. Settings are those of
# the original runs: water = REST2-conformers/<n> (6 replicas, tau 0.4, now 60 ns in one run),
# complex = REST2-conformers-complex/<n> (6 replicas, tau 0.4, ligand-only hot region, 40 ns).
set -uo pipefail
ROOT=/scratch/ipopa/thesis/REST2-conformers-v4
OLDW=/scratch/ipopa/thesis/REST2-conformers
OLDC=/scratch/ipopa/thesis/REST2-conformers-complex
NAMES=(4ibi 4ibk 4j7d); RESS=(L01 L02 L03); LIGS=(1D8 1DE I31)
n=${NAMES[$SLURM_ARRAY_TASK_ID]}; res=${RESS[$SLURM_ARRAY_TASK_ID]}; lig=${LIGS[$SLURM_ARRAY_TASK_ID]}
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scratch/ipopa/md-data
echo "== $n ($res) md-tools $(python -c 'import importlib.metadata as m; print(m.version("md-tools"))')" \
     "$(git -C /scratch/ipopa/software/MD-tools-protonation log --oneline -1)  on $(hostname)  $(date)"
fail() { echo "FAILED: $1"; exit 1; }

echo "== 1. parameterise at pH 7.4"
cd "$ROOT/packages"
if [ ! -f "$n/parameter.config" ]; then
  md-openmm build-top --parameterize -i "$n.input.sdf" --resname "$res" --config "$n.parameterize.config" \
      -op "$n/$res.pdb" -os "$n/$res.xml" -log "$n.parameterize.log" > "$n.parameterize.out" 2>&1 \
    || { grep -v '^$' "$n.parameterize.out" | tail -20; fail parameterise; }
fi
grep -E 'net_formal_charge|canonical_smiles|parameter_id' "$n/parameter.config"
grep -iE 'protonation:' "$n.parameterize.out" | head -3
# the complex atom_map names heavy atoms of the OLD package: they must be the same in the new one
python -I - "$OLDC/packages/$n-water/$res.pdb" "$n/$res.pdb" <<'EOF' || fail "heavy-atom names differ from the old package"
import sys
from openmm.app import PDBFile
heavy = lambda p: [a.name for a in PDBFile(p).topology.atoms() if a.element.symbol != "H"]
old, new = heavy(sys.argv[1]), heavy(sys.argv[2])
print(f"   heavy atoms: old {len(old)}, new {len(new)}, names identical: {old == new}")
sys.exit(0 if old == new else 1)
EOF

echo "== 2. water"
b=$ROOT/$n-water/build
if [ -x "$ROOT/$n-water/REST2-run1/run.sh" ]; then
  echo "   already built (REST2-run1/run.sh exists)"
else
( cd "$b" && md-openmm build-top -i "../../packages/$n/$res.sdf" --config build-top.config \
    -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) || { tail -20 "$b/built.out"; fail "water build-top"; }
grep -E 'reused|net charge|atoms  |protonation|status' "$b/built.out" | head -8
cp -p "$ROOT/packages/$n/molecule.sdf" "$b/$res.sdf"
cp -p "$OLDW/$n/build/scaler.config" "$b/scaler.config"
( cd "$b" && md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1 ) \
  || { tail -20 "$b/scaler.out"; fail "water scaler"; }
grep -E 'selection|unscaled torsions|amide' "$b/scaler.out" | head -5
( cd "$ROOT/$n-water" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
  || { tail -20 "$ROOT/$n-water/build-md.out"; fail "water build-md"; }
fi
grep -hE 'number_of_replicas|number_of_exchanges' "$ROOT/$n-water/input/REST2.in"

echo "== 3. complex"
b=$ROOT/$n-complex/build
[ -s "$b/built.xml" ] && grep -q 'status: completed' "$b/built.out" && echo "   build-top already completed" || \
( cd "$b" && md-openmm build-top -i ../prep/complex.pdb --config build-top.config \
    -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) || { grep -v '^$' "$b/built.out" | tail -25; fail "complex build-top"; }
grep -E 'missing atoms|ligand |reused|net charge|atoms  |WARNING|status' "$b/built.out" | head -12
read -r IDX RNAME < <(python -I - "$b/built.pdb" "$res" "$lig" <<'EOF'
import sys
from openmm.app import PDBFile
hits = [(i, r.name) for i, r in enumerate(PDBFile(sys.argv[1]).topology.residues(), 1)
        if r.name in (sys.argv[2], sys.argv[3])]
assert len(hits) == 1, hits
print(*hits[0])
EOF
) || fail "ligand residue not found"
[ -n "${IDX:-}" ] || fail "ligand residue not found"
echo "   ligand residue :$IDX ($RNAME)"
cp -p "$ROOT/packages/$n/molecule.sdf" "$b/$RNAME.sdf"
sed -E "s/mask: \":[0-9]+\"/mask: \":$IDX\"/; s/residue [0-9]+ =/residue $IDX =/" "$OLDC/$n/build/scaler.config" > "$b/scaler.config"
( cd "$b" && md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1 ) \
  || { tail -20 "$b/scaler.out"; fail "complex scaler"; }
grep -E 'selection|unscaled torsions|amide' "$b/scaler.out" | head -5
sed -E "s/mask: \":[0-9]+\"/mask: \":$IDX\"/" "$OLDC/$n/REST2.config" > "$ROOT/$n-complex/REST2.config"
( cd "$ROOT/$n-complex" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
  || { tail -20 "$ROOT/$n-complex/build-md.out"; fail "complex build-md"; }
grep -hE 'number_of_replicas|number_of_exchanges' "$ROOT/$n-complex/input/REST2.in"
echo "== $n built  $(date)"
