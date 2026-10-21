#!/bin/bash
#SBATCH --job-name=para4-build
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ipopa/thesis/PARA-REST2-v4/logs/build-%j.out
# Paracetamol with Shu-Yu's branch (env openmm-env-protonation): parameterise, solvate, scale, build-md.
set -euo pipefail
ROOT=/scratch/ipopa/thesis/PARA-REST2-v4
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation
export MD_DATA=/scratch/ipopa/md-data
echo "== md-tools from $(python -c 'import md_tools,os;print(os.path.dirname(md_tools.__file__))')  on $(hostname)  $(date)"
cd "$ROOT/build"
md-openmm build-top --parameterize -i paracetamol.sdf --resname TYL --config parameterize.config \
    -op parameter/TYL.pdb -os parameter/TYL.xml -log parameterize.log > parameterize.out 2>&1
grep -E 'net_formal_charge|parameter_id|canonical_smiles' parameter/parameter.config
grep -iE 'protonation' parameterize.out | head -3
md-openmm build-top -i parameter/TYL.sdf --config build-top.config -os built.xml -op built.pdb -log built.log > built.out 2>&1
grep -E 'reused|ligand protonation|solvation|atoms  |status' built.out
md-openmm build-top --rest2-scaler -s built.xml -p built.pdb --config scaler.config > scaler.out 2>&1
grep -E 'selection|unscaled torsions|ligand amides|state[0-9]+\.xml' scaler.out
cd "$ROOT"
md-openmm build-md -odir ./REST2-run1 --config REST2.config > build-md.out 2>&1
grep -E 'status' build-md.out; grep -hE 'number_of_replicas|number_of_exchanges' input/REST2.in
echo "== done $(date)"
