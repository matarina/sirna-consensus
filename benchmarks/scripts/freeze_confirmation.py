"""Freeze the completed strict round; original locks are never rewritten."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PKG / "src"))
from sirna_consensus.confirmation import dump, load_protocol  # noqa: E402
from sirna_consensus.fasta import sha256_file  # noqa: E402
from sirna_consensus.strict_evaluation import verify_artifacts  # noqa: E402


def main():
    directory = PKG / "benchmarks/confirmation_20261106"
    destination = directory / "final_freeze.json"
    if destination.exists():
        raise ValueError("final lock already exists; preserve frozen evidence")
    _, _, baseline, registration = load_protocol(directory / "protocol.json")
    summary = json.loads((directory / "summary.json").read_text())
    output = Path(summary["evaluation"]).parent
    candidate_lock = json.loads((output / "candidate_freeze.json").read_text())
    verify_artifacts(candidate_lock)
    decision = json.loads(Path(summary["decision"]).read_text())
    verification = json.loads(Path(summary["verification"]).read_text())
    if not (
        verification["all_passed"]
        and summary["test_evidence"]["all_passed"]
        and not decision["promote"]
        and summary["default_policy_unchanged"]
        and summary["confirmation_candidate_metrics"] is None
    ):
        raise ValueError("not a verified blocked result")
    if verification["evaluation_sha256"] != sha256_file(Path(summary["evaluation"])):
        raise ValueError("verification is stale")
    evidence = directory / "evidence"
    evidence.mkdir(exist_ok=True)
    for name in (
        "audit.json",
        "evaluation.json",
        "verification.json",
        "decision.json",
        "failure_records.json",
        "test_evidence.json",
        "engineering_evidence.json",
        "development_diagnostics.json",
        "environment.json",
        "candidate_freeze.json",
    ):
        target = evidence / name
        if target.exists():
            raise ValueError("preserve existing evidence copy: " + str(target))
        shutil.copyfile(output / name, target)
    shutil.copyfile(
        Path(summary["development_verification"]), evidence / "development_verification.json"
    )
    paths = set()
    for folder in [
        "src",
        "tests",
        "docs",
        "benchmarks/confirmation_20261106",
        "models/huesken_only_20261106",
        "benchmarks/results/confirmation_20261106",
    ]:
        paths.update(
            p
            for p in (PKG / folder).rglob("*")
            if p.is_file() and "__pycache__" not in p.parts and ".pytest_cache" not in p.parts
        )
    paths.update((PKG / "benchmarks/scripts").glob("*.py"))
    for name in [
        "README.md",
        "pyproject.toml",
        ".github/workflows/ci.yml",
        ".gitignore",
        "benchmarks/evaluate.py",
        "benchmarks/three_benchmark_report_zh.html",
        "benchmarks/computational_defensibility_report_zh.html",
        "models/default_policy.json",
    ]:
        paths.add(PKG / name)
    paths.update(Path(p) for p in candidate_lock["artifacts"])
    dump(
        destination,
        dict(
            schema_version=1,
            kind="strict_confirmation_round",
            baseline_commit=baseline["commit"],
            registration_commit=registration,
            implementation_commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=PKG, text=True
            ).strip(),
            status=summary["status"],
            action="retain_baseline",
            independent_validation_passed=False,
            eligible_confirmation_records=0,
            artifacts={os.path.relpath(p, directory): sha256_file(p) for p in sorted(paths)},
        ),
    )
    print(
        json.dumps(
            dict(lock=str(destination), artifacts=len(paths), independent_validation_passed=False)
        )
    )


if __name__ == "__main__":
    main()
