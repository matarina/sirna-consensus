"""Reproduce the frozen 2026-11-06 attempt without acquiring or reusing test labels.

This registration contains zero eligible confirmation IDs. It deliberately cannot
accept new cohorts after fitting: register a separate future round before exposure.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

PKG = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PKG / "src"))
from sirna_consensus.confirmation import decision, dump, load_protocol  # noqa: E402
from sirna_consensus.fasta import sha256_file  # noqa: E402
from sirna_consensus.strict_evaluation import verify_artifacts, verify_baseline  # noqa: E402


def read(path):
    return json.loads(Path(path).read_text())


def create(path, value):
    if path.exists():
        raise ValueError("preserve first artifact; use a new attempt: " + str(path))
    dump(path, value)


def audit(protocol_path, out):
    p, registry, baseline, registration = load_protocol(protocol_path)
    verify_baseline(baseline, protocol_path.parent)
    counts = []
    for cohort in registry["cohorts"]:
        path = Path(cohort["path"])
        if sha256_file(path) != cohort["sha256"]:
            raise ValueError("cohort source changed: " + str(path))
        counts.append({k: cohort.get(k) for k in ("name", "role", "records", "exclusion_reasons")})
    if p["confirmation_ids"] or registry["eligible_confirmation_records"]:
        raise ValueError("this frozen round has no eligible confirmation IDs")
    create(
        out / "audit.json",
        dict(
            status="no_eligible_independent_confirmation",
            registration_commit=registration,
            baseline_commit=baseline["commit"],
            baseline_artifacts_verified=len(baseline["artifacts"]),
            protocol_sha256=sha256_file(protocol_path),
            cohorts=counts,
            original_data_preserved=True,
            downloads=0,
            eligible_confirmation_records=0,
            reason="Prior exposure cannot be undone by changing the training source.",
        ),
    )


def freeze(protocol_path, model_dir, out):
    p, _, baseline, registration = load_protocol(protocol_path)
    result = read(model_dir / "results.json")
    verification = read(model_dir / "independent_verification.json")
    if not verification["all_passed"] or verification["results_sha256"] != sha256_file(
        model_dir / "results.json"
    ):
        raise ValueError("development arithmetic verification failed or stale")
    if result["protocol_sha256"] != sha256_file(protocol_path) or result["training_sources"] != [
        p["training"]["source"]
    ]:
        raise ValueError("candidate trained outside registered source/protocol")
    artifacts = {
        str(path.resolve()): sha256_file(path) for path in model_dir.rglob("*") if path.is_file()
    }
    source_paths = [
        *list((PKG / "src/sirna_consensus").glob("*.py")),
        PKG / "benchmarks/evaluate.py",
        Path(__file__).resolve(),
        PKG / "benchmarks/scripts/verify_guarded.py",
        PKG / "benchmarks/scripts/verify_platform_v3.py",
    ]
    for path in source_paths:
        artifacts[str(path)] = sha256_file(path)
    for path in (
        protocol_path,
        protocol_path.with_name("cohorts.json"),
        protocol_path.with_name("baseline.json"),
        Path(p["training"]["source"]),
        Path(p["training"]["hgnc"]),
    ):
        artifacts[str(path.resolve())] = sha256_file(path)
    manifest = read(model_dir / "scorer_manifest.json")
    for name, sha in manifest["artifacts"].items():
        path = Path(name) if Path(name).is_absolute() else model_dir / name
        if sha256_file(path) != sha:
            raise ValueError("scorer implementation changed since fitting")
    baseline_model = read(PKG / "models/development_20261106/scorer_manifest.json")
    create(
        out / "candidate_freeze.json",
        dict(
            schema_version=1,
            protocol_sha256=sha256_file(protocol_path),
            registration_commit=registration,
            implementation_commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=PKG, text=True
            ).strip(),
            model_dir=str(model_dir),
            checkpoint_sha256=result["checkpoint_sha256"],
            baseline_checkpoint_sha256=baseline_model["checkpoint_sha256"],
            baseline_commit=baseline["commit"],
            artifacts=artifacts,
            development_only=True,
            default_replaced=False,
        ),
    )


def predict(protocol_path, out):
    p, _, _, _ = load_protocol(protocol_path)
    lock = out / "candidate_freeze.json"
    verify_artifacts(read(lock))
    if p["confirmation_ids"]:
        raise ValueError("cannot add confirmation data to the zero-ID registration after fitting")
    create(
        out / "predictions.json",
        dict(
            protocol_sha256=sha256_file(protocol_path),
            candidate_freeze_sha256=sha256_file(lock),
            scores=[],
            generated_runs=[],
            status="abstained_no_eligible_confirmation",
            labels_read_for_prediction=False,
            excluded_cohorts_scored=False,
        ),
    )


def evaluate(protocol_path, out):
    subprocess.run(
        [
            sys.executable,
            str(PKG / "benchmarks/evaluate.py"),
            "--phase",
            "confirmation",
            "--protocol",
            str(protocol_path),
            "--candidate-freeze",
            str(out / "candidate_freeze.json"),
            "--predictions",
            str(out / "predictions.json"),
            "--output",
            str(out / "evaluation.json"),
        ],
        check=True,
        cwd=PKG,
    )


def verify(protocol_path, out):
    subprocess.run(
        [
            sys.executable,
            str(PKG / "benchmarks/scripts/verify_platform_v3.py"),
            "--evaluation",
            str(out / "evaluation.json"),
            "--protocol",
            str(protocol_path),
            "--candidate-freeze",
            str(out / "candidate_freeze.json"),
            "--predictions",
            str(out / "predictions.json"),
            "--output",
            str(out / "verification.json"),
        ],
        check=True,
        cwd=PKG,
    )


def report(protocol_path, model_dir, out):
    p, registry, baseline, _ = load_protocol(protocol_path)
    evaluation, verification = read(out / "evaluation.json"), read(out / "verification.json")
    if not verification["all_passed"] or verification["evaluation_sha256"] != sha256_file(
        out / "evaluation.json"
    ):
        raise ValueError("independent verification absent, failed or stale")
    evidence = dict(evaluation["gate_evidence"])
    evidence["independent_recomputation_passed"] = verification["all_passed"]
    verdict = decision(p, evidence)
    if verdict["promote"]:
        raise ValueError("zero-ID registration must never authorize promotion")
    verdict.update(
        status="rollback_to_frozen_baseline_without_switching",
        default_policy_sha256=sha256_file(Path(baseline["default_policy_path"])),
        original_default_policy_sha256=baseline["default_policy_sha256"],
        evaluation_sha256=sha256_file(out / "evaluation.json"),
        verification_sha256=sha256_file(out / "verification.json"),
    )
    create(out / "decision.json", verdict)
    model_result = read(model_dir / "results.json")
    summary = dict(
        status=evaluation["status"],
        action=verdict["action"],
        promoted=False,
        independent_validation_passed=False,
        confirmation_candidate_metrics=None,
        confirmation_baseline_metrics=None,
        development_diagnostics=dict(
            records=model_result["records"],
            genes=model_result["genes"],
            components=model_result["components"],
            metrics=model_result["metrics"],
            selected_params=model_result["selected_params"],
            independent_efficacy_confirmation=False,
        ),
        cohorts=registry["cohorts"],
        loso=p["loso"],
        failed_gates=verdict["failures"],
        baseline_commit=baseline["commit"],
        default_policy_unchanged=verdict["default_policy_sha256"]
        == baseline["default_policy_sha256"],
        protocol_sha256=sha256_file(protocol_path),
        evaluation=str(out / "evaluation.json"),
        verification=str(out / "verification.json"),
        decision=str(out / "decision.json"),
        development_verification=str(model_dir / "independent_verification.json"),
        model_dir=str(model_dir),
    )
    for artifact in ("test_evidence.json", "engineering_evidence.json"):
        path = out / artifact
        if path.is_file():
            summary[artifact.removesuffix(".json")] = read(path)
    create(out / "summary.json", summary)
    dump(protocol_path.with_name("summary.json"), summary)
    create(
        out / "failure_records.json",
        dict(
            failures=evaluation["failures"],
            failed_gates=verdict["failures"],
            training_failures=[
                dict(
                    reason="unique case-only approved HGNC symbol Ufc1/UFC1 unresolved before fitting",
                    log="benchmarks/confirmation_20261106/failures/initial_gene_case_mapping.log",
                    correction_commit="06b197c",
                    fitted_models=0,
                )
            ],
            evidence_not_available=[
                "independent NDCG@5 and gain",
                "independent paired CI",
                "actual generated Top5 budget gain",
                "multi-study LOSO under Huesken-only constraint",
            ],
            action="retain existing default; no claim of confirmed performance improvement",
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "stage", choices=["audit", "freeze", "predict", "evaluate", "verify", "report", "all"]
    )
    parser.add_argument(
        "--protocol", type=Path, default=PKG / "benchmarks/confirmation_20261106/protocol.json"
    )
    parser.add_argument("--model-dir", type=Path, default=PKG / "models/huesken_only_20261106")
    parser.add_argument(
        "--out", type=Path, default=PKG / "benchmarks/results/confirmation_20261106"
    )
    args = parser.parse_args()
    protocol, model, out = args.protocol.resolve(), args.model_dir.resolve(), args.out.resolve()
    stages = (
        ["audit", "freeze", "predict", "evaluate", "verify", "report"]
        if args.stage == "all"
        else [args.stage]
    )
    for stage in stages:
        if stage in ("freeze", "report"):
            globals()[stage](protocol, model, out)
        else:
            globals()[stage](protocol, out)
        print(json.dumps(dict(stage=stage, completed=True)))


if __name__ == "__main__":
    main()
