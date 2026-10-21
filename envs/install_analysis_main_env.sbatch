#!/bin/bash
#SBATCH --job-name=env-analysis-main
#SBATCH --partition=defq
#SBATCH --exclude=sbrinz,sbrinz1
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/ipopa/software/logs/%x-%j.out
# Shu-Yu's latest main (symmetry-first t_hdbscan, density_peak, align_all) for ANALYSIS only:
# env analysis-env-main = clone of analysis-env + md-tools from the MD-tools-main worktree (pinned,
# non-editable). analysis-env (md-tools 0.6.4) is left untouched.
set -euo pipefail
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
SRC=/scratch/ipopa/software/MD-tools-main
echo "== $(git -C $SRC log --oneline -1)  on $(hostname)  $(date)"
micromamba create -y -n analysis-env-main --clone analysis-env
micromamba activate analysis-env-main
pip install "$SRC[analysis]"
python - <<'PY'
import md_tools, importlib.metadata as m, os
from md_tools.analysis import t_hdbscan, THDBSCAN, TorsionDefinitions
print("md-tools", m.version("md-tools"), "from", os.path.dirname(md_tools.__file__))
print("align_all:", hasattr(THDBSCAN, "align_all"), " density_peak:", hasattr(THDBSCAN, "density_peak"))
PY
echo "== done $(date)"
