#!/bin/bash
#SBATCH --job-name=rep-build
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ipopa/thesis/REST2-replicates/logs/build-%A_%a.out
# Build one REST2 replicate from a start chosen by pick_starts.py (env openmm-env-protonation, as v4).
# Everything is the v4 water run's except the starting coordinates and the velocity seed:
#   same parameter package (reused, so no charges are recomputed), same scaler, same REST2.config.
#   sbatch --array=0-3 replicates/build.sbatch
#   0 paracetamol rep1     1-3 folate water rep1..rep3
# Output: /scratch/ipopa/thesis/REST2-replicates/<MOL>/rep<k>/{build,REST2.config,REST2-run1}
set -uo pipefail
TH=/scratch/ipopa/thesis
PICK=$TH/bsc-thesis/replicates
case ${SLURM_ARRAY_TASK_ID} in
0) MOL=PARA; K=1; V4=$TH/PARA-REST2-v4; RES=TYL; PKG=$V4/build/parameter
   STARTS=$PICK/paracetamol-v4/starts; V4BUILD=$V4/build; V4CFG=$V4/REST2.config ;;
1|2|3) MOL=FOL; K=${SLURM_ARRAY_TASK_ID}; V4=$TH/FOL-REST2-v4; RES=FOL; PKG=$V4/packages/FOL
   STARTS=$PICK/folic-acid-v4-water/starts; V4BUILD=$V4/water/build; V4CFG=$V4/water/REST2.config ;;
esac
SEED=$((1000 + K))           # rep0 ran with dynamics.seed 1
D=$TH/REST2-replicates/$MOL/rep$K
B=$D/build
fail() { echo "FAILED: $1"; exit 1; }

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scratch/ipopa/md-data
echo "== $MOL rep$K (seed $SEED) on $(hostname)  $(date)"
echo "   md-tools $(python -c 'import importlib.metadata as m; print(m.version("md-tools"))')," \
     "branch $(git -C /scratch/ipopa/software/MD-tools-protonation log --oneline -1)"
[ -e "$D/REST2-run1/run.sh" ] && fail "$D already built (remove it to rebuild)"

echo "== 0. set up $D"
mkdir -p "$B"
cp -p "$STARTS/rep$K.sdf" "$B/start.sdf"
cat > "$B/build-top.config" <<EOF
# REST2 replicate rep$K of the v4 water run: build-top.config of $V4BUILD with the package by path.
# Input = start.sdf (pick_starts.py rep$K); the package is REUSED, its coordinates are not used.
solute:
  kind: ligand
  residue_name: $RES
  parameters: $PKG
protonation:
  ligand_method: openbabel
  ph: 7.4
EOF
cp -p "$V4BUILD/scaler.config" "$B/scaler.config"
{ echo "# Replicate rep$K: $V4CFG with only dynamics.seed changed (rep0: 1)."
  cat "$V4CFG"
  printf "\ndynamics:\n  seed: %s\n" "$SEED"; } > "$D/REST2.config"
grep -q '^dynamics:' "$V4CFG" && fail "$V4CFG already has a dynamics block; merge the seed by hand"

echo "== 1. build-top"
( cd "$B" && md-openmm build-top -i start.sdf --config build-top.config \
    -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) \
  || { grep -v '^$' "$B/built.out" | tail -20; fail build-top; }
grep -E 'reused|net charge|atoms  |solvation|protonation|status' "$B/built.out" | head -10
grep -q 'reused' "$B/built.out" || fail "the package was not reused"
[ -s "$B/$RES.sdf" ] || cp -p "$PKG/molecule.sdf" "$B/$RES.sdf"   # bond orders for the scaler

echo "== 1b. the built solute is the chosen start, and the parameters are rep0's"
python -I - "$B/start.sdf" "$B/built.solute.pdb" "$B/built.xml" "$V4BUILD/built.xml" <<'EOF' || fail "start/parameter check"
import sys
import numpy as np
import mdtraj as md
from rdkit import Chem
import openmm as mm
start = Chem.SDMolSupplier(sys.argv[1], removeHs=False)[0].GetConformer().GetPositions()
built = md.load(sys.argv[2]).xyz[0] * 10.0
n = len(start)
a, b = start - start.mean(0), built[:n] - built[:n].mean(0)
u, _, vt = np.linalg.svd(b.T @ a)
r = u @ np.diag([1, 1, np.sign(np.linalg.det(u @ vt))]) @ vt
rmsd = np.sqrt((((b @ r) - a) ** 2).sum(1).mean())
print(f"   RMSD built solute vs start.sdf: {rmsd:.4f} A ({n} atoms)")
def charges(path):
    s = mm.XmlSerializer.deserialize(open(path).read())
    f = [f for f in s.getForces() if isinstance(f, mm.NonbondedForce)][0]
    return np.array([f.getParticleParameters(i)[0]._value for i in range(n)])
dq = np.abs(charges(sys.argv[3]) - charges(sys.argv[4])).max()
print(f"   max |charge difference| vs rep0 build: {dq:.2e} e")
sys.exit(0 if rmsd < 0.05 and dq < 1e-6 else 1)
EOF

echo "== 2. scaler"
( cd "$B" && md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config \
    > scaler.out 2>&1 ) || { tail -20 "$B/scaler.out"; fail scaler; }
grep -E 'selection|unscaled torsions|ligand amides|amide' "$B/scaler.out" | head -6

echo "== 3. build-md"
( cd "$D" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
  || { tail -20 "$D/build-md.out"; fail build-md; }
grep -E 'status' "$D/build-md.out"
grep -hE 'number_of_replicas|number_of_exchanges' "$D/input/REST2.in"
grep -m1 -A9 '^#   dynamics:' "$D/REST2-run1/build-md.log" | grep seed
echo "== done  $(date)"
