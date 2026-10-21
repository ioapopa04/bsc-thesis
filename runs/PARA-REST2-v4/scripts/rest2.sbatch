#!/bin/bash
#SBATCH --job-name=para4-rest2
#SBATCH --partition=defq
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:2
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/ipopa/thesis/PARA-REST2-v4/logs/rest2-%x-%j.out
# Paracetamol v4 REST2 run (min, eq_1..3, production: 4 replicas on 2 GPUs, 2 per GPU under MPS).
#   sbatch --exclude=sbrinz1 scripts/rest2.sbatch REST2-run1
#
# 10 replicas on 5 GPUs = 2 ranks per GPU, which md-tools allows only under NVIDIA MPS. The MPS
# daemon is private to this job (pipe dir keyed by job id) and is started on exactly the GPUs Slurm
# allocated; clients then see them as 0..n-1. It is stopped when the job ends.
# run.sh resumes an interrupted stage from its checkpoint, so resubmitting after a crash is safe.
set -uo pipefail
ROOT=/scratch/ipopa/thesis/PARA-REST2-v4
RUN="$ROOT/$1"
[ -x "$RUN/run.sh" ] || { echo "no $RUN/run.sh"; exit 2; }

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env-protonation   # Shu-Yu branch b3603e0 (ligand pH model, REST2 v4)
export MD_DATA=/scratch/ipopa/md-data

ALLOC="${CUDA_VISIBLE_DEVICES:?Slurm set no CUDA_VISIBLE_DEVICES}"
NG=$(awk -F, '{print NF}' <<<"$ALLOC")
echo "host $(hostname)  job $SLURM_JOB_ID  GPUs (Slurm) $ALLOC  run $1  $(date)"
nvidia-smi -L

export CUDA_MPS_PIPE_DIRECTORY="/tmp/mps-$SLURM_JOB_ID/pipe"
export CUDA_MPS_LOG_DIRECTORY="/tmp/mps-$SLURM_JOB_ID/log"
mkdir -p "$CUDA_MPS_PIPE_DIRECTORY" "$CUDA_MPS_LOG_DIRECTORY"
stop_mps() { echo quit | nvidia-cuda-mps-control 2>/dev/null; rm -rf "/tmp/mps-$SLURM_JOB_ID"; }
trap stop_mps EXIT
export CUDA_DEVICE_ORDER=PCI_BUS_ID
nvidia-cuda-mps-control -d || { echo "could not start MPS"; exit 1; }
echo "MPS daemon on $ALLOC ($CUDA_MPS_PIPE_DIRECTORY)"

export CUDA_VISIBLE_DEVICES=$(seq -s, 0 $((NG - 1)))   # the daemon's devices, as clients see them
# md-tools takes the MPI world size from OMPI_COMM_WORLD_SIZE, else PMI_SIZE, else SLURM_NTASKS. The
# serial stages (min, eq_1..3) are not under mpirun, so SLURM_NTASKS=10 made them refuse ("launcher
# world size 10, MPI communicator size 1"; jobs 8024/8025). mpirun still takes its 10 slots from
# SLURM_TASKS_PER_NODE, and its ranks get OMPI_COMM_WORLD_SIZE, which is read first.
unset SLURM_NTASKS SLURM_NPROCS
cd "$RUN"
./run.sh
rc=$?
echo "run.sh rc=$rc  $(date)"
grep -E 'throughput|run_status' remd_records/REST2_prod1.out 2>/dev/null
exit $rc
