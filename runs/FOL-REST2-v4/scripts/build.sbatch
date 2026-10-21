#!/bin/bash
#SBATCH --job-name=fol4-build
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ipopa/thesis/FOL-REST2-v4/logs/build-%j.out
# Folate REST2 with Shu-Yu's md-tools branch (ligand pH model + REST2 v4 amide scaling), env
# openmm-env-protonation. One new folate package (COO-, pH 7.4) shared by water, 4LRH and 3QL3.
set -uo pipefail
ROOT=/scratch/ipopa/thesis/FOL-REST2-v4
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scratch/ipopa/md-data
echo "== md-tools $(python -c 'import importlib.metadata as m; print(m.version("md-tools"))') from $(python -c 'import md_tools,os;print(os.path.dirname(md_tools.__file__))')"
echo "   branch $(git -C /scratch/ipopa/software/MD-tools-protonation log --oneline -1)  on $(hostname)  $(date)"
fail() { echo "FAILED: $1"; exit 1; }

echo "== 1. parameterise folate"
cd "$ROOT/packages"
if [ ! -f FOL/parameter.config ]; then
  md-openmm build-top --parameterize -i FOL.input.sdf --resname FOL --config FOL.parameterize.config \
      -op FOL/FOL.pdb -os FOL/FOL.xml -log FOL.parameterize.log > FOL.parameterize.out 2>&1 \
    || { grep -v '^$' FOL.parameterize.out | tail -20; fail parameterise; }
fi
grep -E 'net_formal_charge|canonical_smiles|parameter_id' FOL/parameter.config
grep -iE 'protonat|openbabel|charge' FOL.parameterize.out | head -8

build_one() {   # build_one <dir> <input> <resid or ""> <n_states>
  local d=$1 inp=$2 idx=$3 n=$4 b=$ROOT/$1/build
  echo "== 2. $d: build-top"
  ( cd "$b" && md-openmm build-top -i "$inp" --config build-top.config \
      -os built.xml -op built.pdb -log built.log > built.out 2>&1 ) \
    || { grep -v '^$' "$b/built.out" | tail -20; fail "$d build-top"; }
  grep -E 'reused|net charge|atoms  |solvation|protonation|status' "$b/built.out" | head -10
  cp -p "$ROOT/packages/FOL/molecule.sdf" "$b/FOL.sdf"     # bond orders for the scaler
  if [ -n "$idx" ]; then
    python -I - "$b/built.pdb" "$idx" <<'EOF' || fail "$d: residue check"
import sys
from openmm.app import PDBFile
r = list(PDBFile(sys.argv[1]).topology.residues())[int(sys.argv[2]) - 1]
print(f"   residue {sys.argv[2]} = {r.name}")
sys.exit(0 if r.name == "FOL" else 1)
EOF
    cat > "$b/scaler.config" <<EOF
# Selective REST2 v4: folate only (residue $idx), $n states, tau 0..0.4 linear. Protein/NADP+ cold.
method: REST2
schedule: {kind: linear, n_states: $n, tau_min: 0.0, tau_max: 0.4}
ligand_scaling_dict:
  FOL:
    mask: ":$idx"
    torsion_exclusions: auto
EOF
  fi
  echo "== 3. $d: scaler"
  ( cd "$b" && md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config \
      > scaler.out 2>&1 ) || { tail -20 "$b/scaler.out"; fail "$d scaler"; }
  grep -E 'selection|unscaled torsions|proline-like|amide' "$b/scaler.out" | head -6
  echo "== 4. $d: build-md"
  ( cd "$ROOT/$d" && md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1 ) \
    || { tail -20 "$ROOT/$d/build-md.out"; fail "$d build-md"; }
  grep -E 'claimed region|status' "$ROOT/$d/build-md.out"
  grep -hE 'number_of_replicas|number_of_exchanges' "$ROOT/$d/input/REST2.in"
}

build_one water ../../packages/FOL/FOL.sdf "" 8
build_one 4LRH ../prep/complex.pdb 208 10
build_one 3QL3 ../prep/complex.pdb 160 10
echo "== all built  $(date)"
