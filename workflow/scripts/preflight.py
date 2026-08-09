import hashlib
import json
import os
from pathlib import Path
import subprocess


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


paths = {name: Path(value) for name, value in snakemake.params.paths.items()}
required_files = [
    "env_script",
    "gibuu_executable",
    "gibuu_version",
    "release_archive",
    "buuinput_archive",
    "build_jobcard",
    "gibuu2root",
    "analyzer",
    "comparator",
]
for name in required_files:
    path = paths[name]
    if not path.is_file():
        raise RuntimeError(f"Required file does not exist: {name}={path}")

for name in ["gibuu_executable", "build_jobcard", "gibuu2root", "analyzer", "comparator"]:
    if not os.access(paths[name], os.X_OK):
        raise RuntimeError(f"Required executable is not executable: {name}={paths[name]}")

if not paths["buuinput"].is_dir():
    raise RuntimeError(f"BUUInput directory does not exist: {paths['buuinput']}")

version = paths["gibuu_version"].read_text().strip()
if "Release 2025, patch 5" not in version:
    raise RuntimeError(f"Expected GiBUU Release 2025 patch 5, got: {version}")

repo_root = paths["repo_root"]
git_commit = subprocess.run(
    ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
    check=True,
    capture_output=True,
    text=True,
).stdout.strip()
git_status = subprocess.run(
    ["git", "-C", str(repo_root), "status", "--porcelain"],
    check=True,
    capture_output=True,
    text=True,
).stdout.strip()
if git_status:
    raise RuntimeError("Remote analysis repository is not clean:\n" + git_status)

manifest = {
    "study": snakemake.params.study,
    "git_commit": git_commit,
    "gibuu_version": version,
    "paths": {name: str(path) for name, path in paths.items()},
    "sha256": {
        "release_archive": sha256(paths["release_archive"]),
        "buuinput_archive": sha256(paths["buuinput_archive"]),
        "gibuu_executable": sha256(paths["gibuu_executable"]),
        "build_jobcard": sha256(paths["build_jobcard"]),
        "gibuu2root": sha256(paths["gibuu2root"]),
        "analyzer": sha256(paths["analyzer"]),
        "comparator": sha256(paths["comparator"]),
    },
}

output = Path(snakemake.output[0])
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(manifest, indent=2) + "\n")
