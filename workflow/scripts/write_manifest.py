import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


jobcards = [Path(path) for path in snakemake.input.jobcards]
records = [Path(path) for path in snakemake.input.records]
validations = [Path(path) for path in snakemake.input.validations]
preflight = json.loads(Path(snakemake.input.preflight).read_text())

manifest = {
    "study": snakemake.params.study,
    "stage": snakemake.params.stage,
    "scenario": snakemake.params.scenario_config,
    "use_cda": snakemake.params.use_cda,
    "replicates": len(records),
    "num_ensembles": snakemake.params.num_ensembles,
    "num_runs_same_energy": snakemake.params.num_runs_same_energy,
    "seeds": snakemake.params.seeds,
    "preflight": preflight,
    "jobcards": [
        {"path": str(path), "sha256": sha256(path)} for path in jobcards
    ],
    "records": [
        {"path": str(path), "bytes": path.stat().st_size} for path in records
    ],
    "validations": [json.loads(path.read_text()) for path in validations],
    "analysis": {
        "root": str(snakemake.input.analysis_root),
        "root_sha256": sha256(Path(snakemake.input.analysis_root)),
        "chi2": str(snakemake.input.analysis_json),
        "chi2_sha256": sha256(Path(snakemake.input.analysis_json)),
    },
}

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(manifest, indent=2) + "\n")
