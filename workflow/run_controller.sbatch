#!/usr/bin/env bash
#SBATCH --partition=htc_daemon
#SBATCH --account=juno
#SBATCH --licenses=sps
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=/sps/juno/yqiyu/GiBUUGEN/minerva_useCdA_2025p5/snakemake_controller_%j.slurm.log

set -euo pipefail

stage=${1:-production}
case "$stage" in
  smoke|production) ;;
  *)
    printf 'Expected stage smoke or production, got: %s\n' "$stage" >&2
    exit 2
    ;;
esac

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=${SNAKEMAKE_REPO_ROOT:-${SLURM_SUBMIT_DIR:-}}
if [[ -z "$repo_root" ]]; then
  repo_root=$(cd -- "$script_dir/.." && pwd)
fi
controller_log=${SNAKEMAKE_CONTROLLER_LOG:-/sps/juno/yqiyu/GiBUUGEN/minerva_useCdA_2025p5/snakemake_${stage}_controller.log}
controller_env=${SNAKEMAKE_ENV_SCRIPT:-/sps/juno/yqiyu/GiBUUGEN/env.sh}
controller_python=${SNAKEMAKE_PYTHON:-/cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/bin/python3}
snakemake_python=${SNAKEMAKE_VENV_PYTHON:-$repo_root/.venv/bin/python}

if [[ ! -f "$controller_env" ]]; then
  printf 'Controller environment script is unavailable: %s\n' "$controller_env" >&2
  exit 1
fi

set +u
source "$controller_env"
set -u

if ! env -u PYTHONPATH "$controller_python" --version; then
  printf 'Controller Python is unavailable after loading %s: %s\n' \
    "$controller_env" "$controller_python" >&2
  exit 1
fi
if ! env -u PYTHONPATH "$snakemake_python" -m snakemake --version; then
  printf 'Snakemake venv is unavailable: %s\n' "$snakemake_python" >&2
  exit 1
fi

cd "$repo_root"

{
  printf 'controller_job_id=%s\n' "${SLURM_JOB_ID:-unknown}"
  printf 'controller_host=%s\n' "$(hostname)"
  printf 'controller_stage=%s\n' "$stage"
  printf 'repo_root=%s\n' "$repo_root"
  printf 'controller_env=%s\n' "$controller_env"
  printf 'controller_python=%s\n' "$controller_python"
  printf 'snakemake_python=%s\n' "$snakemake_python"
  printf 'controller_started_at=%s\n' "$(date --iso-8601=seconds)"
} >> "$controller_log"

if [[ ${SNAKEMAKE_UNLOCK_FIRST:-0} == 1 ]]; then
  env -u PYTHONPATH "$snakemake_python" -m snakemake \
    --snakefile workflow/Snakefile \
    --configfile workflow/config/minerva_usecda.yaml \
    --unlock >> "$controller_log" 2>&1
fi

exec env -u PYTHONPATH "$snakemake_python" -m snakemake \
  --snakefile workflow/Snakefile \
  --configfile workflow/config/minerva_usecda.yaml \
  --profile workflow/profiles/slurm \
  "$stage" >> "$controller_log" 2>&1
