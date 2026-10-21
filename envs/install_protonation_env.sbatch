#!/bin/bash
#SBATCH --job-name=env-protonation
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/ipopa/software/logs/%x-%j.out
# Shu-Yu's feature/ligand-protonation-and-ligand-amides (ligand pH model + REST2 v4 amide scaling),
# installed side by side: env openmm-env-protonation, editable from the worktree MD-tools-protonation.
# The original openmm-env (md-tools 0.6.2) and the MD-tools checkout on main are left untouched.
set -euo pipefail
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
SRC=/scratch/ipopa/software/MD-tools-protonation
echo "== $(git -C $SRC log --oneline -1)  on $(hostname)  $(date)"
micromamba create -y -n openmm-env-protonation -f "$SRC/environment.yml"
micromamba activate openmm-env-protonation
pip install -e "$SRC"
echo "== check"
python -c "import md_tools, importlib.metadata as m, os; print('md-tools', m.version('md-tools'), 'from', os.path.dirname(md_tools.__file__))"
python -c "from md_tools.openmm.ligand_protonation import openbabel_version; print('openbabel', openbabel_version())"
python -c "import openmm; print('openmm', openmm.__version__)"
which md-openmm
echo "== done $(date)"
