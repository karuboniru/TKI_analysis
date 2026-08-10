# GiBUU MINERvA nuclear-model workflow

This workflow compares the default GiBUU nuclear momentum distribution with
the Ciofi degli Atti--Simula parametrization selected by
`&initNucleus_in_PS useCdA = T`. The default configuration changes only that
switch, uses common random seeds for each F/T pair, and produces MINERvA 0-pion
and neutral-pion TKI analyses plus direct T/F comparisons.

## Remote bootstrap

All repository edits are made locally and synchronized through Git. On the
cluster, use an ordinary SSH shell; do not run an AI agent there.

```bash
cd /sps/juno/yqiyu/TKI_analysis
git fetch origin
git switch dev/minerva-usecda-snakemake
git pull --ff-only

bash workflow/bootstrap_gibuu.sh \
  /sps/juno/yqiyu/GiBUU2025_p5 \
  /sps/juno/yqiyu/GiBUUGEN/env.sh

source /sps/juno/yqiyu/GiBUUGEN/env.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4 --target \
  build_jobcard gibuu2root minerva_all_gibuu_plot \
  minerva_nuclear_model_compare

UV_CACHE_DIR=/sps/juno/yqiyu/uv/cache \
  uv sync --frozen --python "$(command -v python3)"
```

The bootstrap is idempotent. It downloads the two official 2025 archives,
records their SHA256 sums, builds with four compiler jobs, and refuses to
continue unless `version.txt` identifies Release 2025 patch 5.

## Validate and run

The smoke test runs one 100-ensemble job for each nuclear model. It still uses
10 same-energy runs and all 300 FSI time steps, so it exercises the complete
physics and analysis path.

```bash
source /sps/juno/yqiyu/GiBUUGEN/env.sh
cd /sps/juno/yqiyu/TKI_analysis

env -u PYTHONPATH uv run snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config/minerva_usecda.yaml --lint

env -u PYTHONPATH uv run snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config/minerva_usecda.yaml \
  --profile workflow/profiles/slurm --dry-run smoke

workflow/submit_controller.sh smoke
workflow/submit_controller.sh production
```

Production submits 256 GiBUU jobs for each of `useCdA=F` and `useCdA=T`.
Snakemake retries failed jobs twice, retains failed raw event files, and removes
successful `FinalEvents.dat` files only after validating the converted ROOT
tree with uproot. The submission wrapper runs the Snakemake controller itself
as an `sbatch` job on the infinite-time `htc_daemon` partition, so losing the
SSH connection does not stop DAG scheduling. The controller loads the shared
GiBUUGEN/LCG environment so its Python 3.11 interpreter also works on the
daemon node, clears LCG's `PYTHONPATH` before invoking Snakemake, and uses the
frozen lockfile. Prepare that shared environment from a login node with:

```bash
set +u
source /sps/juno/yqiyu/GiBUUGEN/env.sh
set -u
env -u PYTHONPATH uv sync --frozen \
  --python /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/bin/python3
```

If a controller was accidentally run in an interactive SSH session and that
session was lost, first get the orphan workers' common Slurm job name from
`squeue`, then submit a dependency-safe recovery controller:

```bash
workflow/submit_controller.sh production ORPHAN_WORKER_JOB_NAME
```

The recovery controller stays pending on `htc_daemon` until every captured
orphan worker is terminal. It then removes the stale Snakemake lock and resumes
from the outputs already present, avoiding duplicate GiBUU jobs.

Useful monitoring and recovery commands are:

```bash
squeue -u "$USER"
env -u PYTHONPATH uv run snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config/minerva_usecda.yaml \
  --profile workflow/profiles/slurm --summary
env -u PYTHONPATH uv run snakemake --snakefile workflow/Snakefile \
  --configfile workflow/config/minerva_usecda.yaml \
  --profile workflow/profiles/slurm --rerun-incomplete production
```

## Outputs and matrix expansion

The default output root is
`/sps/juno/yqiyu/GiBUUGEN/minerva_useCdA_2025p5/minerva_usecda_2025p5`.
Each stage contains paired `use_cda_F` and `use_cda_T` sample directories with
jobcards, compressed logs, validated ROOT files, analysis products, and a
provenance manifest. The scenario-level `comparison` directory contains four
two-panel TKI comparison plots in PDF/SVG/EPS, `comparison.root`, and
`comparison.json`.

`workflow/config/minerva_full_factorial.yaml` exposes all legacy GiBUUGEN axes:
T, T2p2h, Oset width, in-medium NN mode, 2-pion background, A dependence, FSI,
and useCdA. It expands to 384 samples and 98,304 production GiBUU jobs at the
default statistics. Submission is intentionally blocked until
`matrix.allow_large_matrix` is explicitly changed to `true`; inspect a dry-run
before doing so.
