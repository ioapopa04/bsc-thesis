#!/usr/bin/env bash
# Parameterise one package: ./parameterize.sh FOL|NAP  (OpenFF + AM1-BCC/sqm, md-tools defaults)
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env
r=$1
md-openmm build-top --parameterize -i "$r.sdf" --resname "$r" --config "$r.parameterize.config" \
  -op "$r/$r.pdb" -os "$r/$r.xml" -log "$r.parameterize.log"
