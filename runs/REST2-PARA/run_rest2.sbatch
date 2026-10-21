#!/bin/bash
#SBATCH --job-name=para-rest2
#SBATCH --partition=defq
#SBATCH -N 1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:2
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/ipopa/thesis/REST2-PARA/logs/rest2-%x-%j.out
# Paracetamol REST2 (4 replicas) inside a Slurm GPU allocation; resumes an interrupted run.
#   sbatch run_rest2.sbatch REST2-run1
#   sbatch --export=ALL,RUNSCRIPT=run_resume.sh run_rest2.sbatch REST2-run1     (continue an interrupted run)
#   An interrupted REMD stage needs --resume (md-run refuses to write into existing output without it:
#   job 9486), but the serial min/eq stages REJECT --resume (job 11466). run_resume.sh is run.sh with
#   --resume added to the REST2 line only.
#
# 4 replicas on 2 GPUs = 2 ranks per GPU, which md-tools allows only under NVIDIA MPS. The MPS
# daemon is private to this job (pipe dir keyed by job id) and is started on exactly the GPUs Slurm
# allocated; clients then see them as 0..n-1. It is stopped when the job ends.
# run.sh resumes an interrupted stage from its checkpoint, so resubmitting after a crash is safe.
set -uo pipefail
ROOT=/scratch/ipopa/thesis/REST2-PARA
RUN="$ROOT/$1"
[ -x "$RUN/run.sh" ] || { echo "no $RUN/run.sh"; exit 2; }

export MAMBA_ROOT_PREFIX=/scratch/ipopa/software/md-stack
eval "$(/scratch/ipopa/software/md-stack/bin/micromamba shell hook -s bash)"
micromamba activate openmm-env
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
# serial stages (min, eq_1..3) are not under mpirun, so SLURM_NTASKS=4 made them refuse ("launcher
# world size N, MPI communicator size 1"; FOL-complex jobs 8024/8025). mpirun still takes its 4 slots from
# SLURM_TASKS_PER_NODE, and its ranks get OMPI_COMM_WORLD_SIZE, which is read first.
unset SLURM_NTASKS SLURM_NPROCS
cd "$RUN"
./${RUNSCRIPT:-run.sh} "${@:2}"   # RUNSCRIPT=run_resume.sh: --resume on the REST2 stage only
rc=$?
echo "run.sh rc=$rc  $(date)"
grep -E 'throughput|run_status' remd_records/REST2_prod1.out 2>/dev/null
exit $rc
