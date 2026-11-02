"""Create a checksum lock before final labels are opened, or verify an existing lock.

Keys are stored relative to the lock file's directory so a checkout at a different
absolute path still verifies (absolute keys elsewhere, e.g. inside run.json, remain
path-bound and force a recut on a clone -- see benchmarks/README.md).
"""
import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from sirna_consensus.fasta import sha256_file


def key_for(lock: Path, path: Path) -> str:
    return os.path.relpath(Path(path).resolve(), lock.parent.resolve())


def resolve_key(lock: Path, key: str) -> Path:
    candidate = Path(key)
    return candidate if candidate.is_absolute() else lock.parent / candidate


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--lock", type=Path, required=True)
    p.add_argument("--verify", action="store_true")
    p.add_argument("files", nargs="*", type=Path)
    a = p.parse_args(argv)
    if a.verify:
        try:
            d = json.loads(a.lock.read_text())
        except json.JSONDecodeError as exc:
            raise ValueError("freeze not parseable: " + str(a.lock)) from exc
        for name, digest in d["artifacts"].items():
            target = resolve_key(a.lock, name)
            if not target.is_file():
                raise ValueError("freeze mismatch (missing): " + name)
            if sha256_file(target) != digest:
                raise ValueError("freeze mismatch: " + name)
        print("freeze verified")
    else:
        if a.lock.exists():
            raise ValueError("refusing to overwrite an existing freeze")
        if not a.files:
            raise ValueError("provide code, config, splits, model manifests and label-blind predictions")
        a.lock.write_text(json.dumps({
            "artifacts": {key_for(a.lock, f): sha256_file(f) for f in a.files}}, indent=2) + "\n")


if __name__ == "__main__":
    main()
