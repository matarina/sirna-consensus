"""Verify freeze locks whose pinned artifacts are present in this checkout.

The locks pin gitignored, generated artifacts (dev results, holdout runs), so on a
bare clone most entries are legitimately absent. This entry point reports, per lock,
how many entries are missing vs. mismatched, and fails only on a MISMATCH of a file
that is present. Keys are relative to the lock directory (see freeze.py).
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from sirna_consensus.fasta import sha256_file  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]   # benchmarks/
LOCK_DIR = ROOT


def check(lock: Path) -> tuple[int, list[str], int]:
    try:
        data = json.loads(lock.read_text())["artifacts"]
    except (json.JSONDecodeError, KeyError, TypeError, AttributeError) as exc:
        raise ValueError(f"not parseable ({exc})") from exc
    if not isinstance(data, dict):
        raise ValueError(f"not parseable ('artifacts' is {type(data).__name__})")
    present = missing = 0
    mismatched = []
    for key, digest in data.items():
        candidate = Path(key)
        target = candidate if candidate.is_absolute() else lock.parent / candidate
        if not target.is_file():
            missing += 1
        elif sha256_file(target) != digest:
            mismatched.append(key)
        else:
            present += 1
    return present, mismatched, missing


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("locks", nargs="*", type=Path, default=[
        LOCK_DIR / "development_freeze.json", LOCK_DIR / "final_freeze.json"])
    a = p.parse_args(argv)
    failed = False
    # Never abort the whole run on one bad lock: report it and check the rest.
    for lock in a.locks:
        if not lock.is_file():
            print(f"{lock.name}: absent (skip)")
            continue
        try:
            present, mismatched, missing = check(lock)
        except ValueError as exc:
            print(f"{lock.name}: NOT PARSEABLE — {exc}")
            failed = True
            continue
        total = present + len(mismatched)
        status = "MISMATCH" if mismatched else "ok"
        print(f"{lock.name}: {status} — {present}/{total} present entries verify, {missing} absent")
        for key in mismatched:
            print(f"  mismatch: {key}")
        failed = failed or bool(mismatched)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
