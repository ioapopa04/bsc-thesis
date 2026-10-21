# Shared env for the FOL REST2 workflow
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)" 2>/dev/null
micromamba activate openmm-env
export SCRIPTS=/scratch/ipopa/thesis/bsc-thesis
export ROOT=/scratch/ipopa/thesis/FOLn-REST2
export RES=FOL
export MD_DATA=/scratch/ipopa/md-data
