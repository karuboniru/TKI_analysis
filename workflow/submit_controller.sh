#!/usr/bin/env bash
set -euo pipefail

stage=${1:-production}
wait_for_job_name=${2:-}
case "$stage" in
  smoke|production) ;;
  *)
    printf 'Usage: %s [smoke|production] [orphan-worker-job-name]\n' "$0" >&2
    exit 2
    ;;
esac

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
controller_name="smk-minerva-usecda-${stage}"
if [[ -n $(squeue -h -u "$USER" -n "$controller_name" -o '%A') ]]; then
  printf 'A controller named %s is already queued or running.\n' "$controller_name" >&2
  exit 1
fi

sbatch_args=(
  --parsable
  --job-name="$controller_name"
)
unlock_first=0

if [[ -n "$wait_for_job_name" ]]; then
  mapfile -t worker_ids < <(
    squeue -h -u "$USER" -n "$wait_for_job_name" -o '%A' | sort -u
  )
  if (( ${#worker_ids[@]} > 0 )); then
    dependency=$(IFS=:; printf '%s' "${worker_ids[*]}")
    sbatch_args+=(
      --dependency="afterany:${dependency}"
    )
    unlock_first=1
    printf 'Waiting for %d orphan worker jobs named %s.\n' \
      "${#worker_ids[@]}" "$wait_for_job_name"
  else
    printf 'No active orphan workers named %s; resuming immediately.\n' \
      "$wait_for_job_name"
    unlock_first=1
  fi
fi

sbatch_args+=(--export="ALL,SNAKEMAKE_UNLOCK_FIRST=${unlock_first}")
job_id=$(sbatch "${sbatch_args[@]}" "$script_dir/run_controller.sbatch" "$stage")
printf 'Submitted %s controller as Slurm job %s on htc_daemon.\n' \
  "$stage" "$job_id"
