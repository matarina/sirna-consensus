"""Evaluate a sirna-consensus run against an experimental benchmark.

This is a comparison harness, not bundled experimental data. It answers two
questions for each run:

* recovery@K: how many known effective sites are matched by the Top-K?
* precision@K: how many Top-K guides match a known effective site?

Benchmark TSV columns required:
    transcript_id, target_start, target_end, effective

`effective` may be 1/0 or a boolean-ish string. Optionally use `knockdown` with
`--effective-column knockdown --threshold 0.5`.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from sirna_consensus.consensus import agree  # noqa: E402
from sirna_consensus.types import Core  # noqa: E402


def _as_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y", "effective", "high"}


def _as_float(value: str) -> float:
    try:
        return float(value)
    except ValueError:
        return 0.0


def _effective(row: dict[str, str], column: str, threshold: float) -> bool:
    if column not in row:
        raise ValueError(f"benchmark missing column {column!r}")
    value = row[column]
    if column == "knockdown":
        return _as_float(value) >= threshold
    return _as_bool(value)


def _load_benchmark(path: Path, column: str, threshold: float) -> list[dict]:
    items = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            transcript = row.get("transcript_id") or row.get("transcript")
            start = row.get("target_start") or row.get("start")
            end = row.get("target_end") or row.get("end")
            if not transcript or not start or not end:
                continue
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                continue
            length = end_i - start_i + 1
            items.append({
                "transcript_id": transcript,
                "start": start_i,
                "end": end_i,
                "length": length,
                "effective": _effective(row, column, threshold),
            })
    return items


def _run_top(path: Path, k: int) -> list[Core]:
    record = json.loads(path.read_text())
    candidates = record.get("candidates", [])
    top_keys = record.get("top", [])
    by_key = {
        f"{c['transcript_id']}:{c['start']}:{c['length']}": c
        for c in candidates
    }
    out = []
    for key in top_keys[:k]:
        c = by_key.get(key)
        if not c:
            continue
        out.append(Core(c["transcript_id"], c["start"], c["length"],
                        c.get("sense_dna", "A" * c["length"])))
    return out


def _matches(core: Core, item: dict, tolerance: int = 2) -> bool:
    bench = Core(item["transcript_id"], item["start"], item["length"],
                 "A" * item["length"])
    return agree(core, bench, tolerance)


def evaluate(run_path: Path, benchmark: list[dict], k: int,
             tolerance: int) -> dict:
    top = _run_top(run_path, k)
    effective = [item for item in benchmark if item["effective"]]
    matched_effective = sum(
        1 for item in effective
        if any(_matches(core, item, tolerance) for core in top))
    matched_top = sum(
        1 for core in top
        if any(_matches(core, item, tolerance) for item in effective))
    return {
        "run": str(run_path),
        "k": k,
        "top_count": len(top),
        "effective_sites": len(effective),
        "recovery_at_k": matched_effective / len(effective) if effective else 0.0,
        "precision_at_k": matched_top / len(top) if top else 0.0,
        "matched_effective": matched_effective,
        "matched_top": matched_top,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True)
    parser.add_argument("--run", required=True, action="append",
                        help="run.json path; repeat to compare multiple runs")
    parser.add_argument("--top-k", default="5,10")
    parser.add_argument("--tolerance", type=int, default=2)
    parser.add_argument("--effective-column", default="effective")
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--output", default="")
    args = parser.parse_args(argv)

    benchmark = _load_benchmark(Path(args.benchmark), args.effective_column,
                                args.threshold)
    ks = [int(k) for k in args.top_k.split(",") if k.strip()]
    results = []
    for run in args.run:
        for k in ks:
            results.append(evaluate(Path(run), benchmark, k, args.tolerance))

    text = json.dumps(results, indent=2)
    if args.output:
        Path(args.output).write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
