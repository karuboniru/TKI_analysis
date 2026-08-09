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
repo_root=$(cd -- "$script_dir/.." && pwd)
controller_log=${SNAKEMAKE_CONTROLLER_LOG:-/sps/juno/yqiyu/GiBUUGEN/minerva_useCdA_2025p5/snakemake_${stage}_controller.log}

cd "$repo_root"

{
  printf 'controller_job_id=%s\n' "${SLURM_JOB_ID:-unknown}"
  printf 'controller_host=%s\n' "$(hostname)"
  printf 'controller_stage=%s\n' "$stage"
  printf 'controller_started_at=%s\n' "$(date --iso-8601=seconds)"
} >> "$controller_log"

if [[ ${SNAKEMAKE_UNLOCK_FIRST:-0} == 1 ]]; then
  env -u PYTHONPATH uv run snakemake \
    --snakefile workflow/Snakefile \
    --configfile workflow/config/minerva_usecda.yaml \
    --unlock >> "$controller_log" 2>&1
fi

exec env -u PYTHONPATH uv run snakemake \
  --snakefile workflow/Snakefile \
  --configfile workflow/config/minerva_usecda.yaml \
  --profile workflow/profiles/slurm \
  "$stage" >> "$controller_log" 2>&1
