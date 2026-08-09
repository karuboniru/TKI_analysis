import json
from pathlib import Path

import uproot


record_path = Path(snakemake.input.record)
events_path = Path(snakemake.input.events)
if not record_path.is_file() or record_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty ROOT output: {record_path}")
if not events_path.is_file() or events_path.stat().st_size == 0:
    raise RuntimeError(f"Missing or empty GiBUU event output: {events_path}")

required_branches = {
    "StdHepN",
    "StdHepPdg",
    "StdHepStatus",
    "StdHepP4",
    "StdHepX4",
    "weight",
    "channel",
}
with uproot.open(record_path) as root_file:
    if "out_tree" not in root_file:
        raise RuntimeError(f"out_tree is missing from {record_path}")
    tree = root_file["out_tree"]
    entries = tree.num_entries
    if entries <= 0:
        raise RuntimeError(f"out_tree has no entries in {record_path}")
    branches = set(tree.keys())
    missing = sorted(required_branches - branches)
    if missing:
        raise RuntimeError(f"Missing branches in {record_path}: {missing}")

result = {
    "stage": snakemake.params.stage,
    "scenario": snakemake.params.scenario,
    "use_cda": snakemake.params.use_cda,
    "replicate": int(snakemake.params.replicate),
    "seed": int(snakemake.params.seed),
    "record": str(record_path),
    "record_bytes": record_path.stat().st_size,
    "entries": entries,
    "branches": sorted(branches),
}
output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=2) + "\n")
