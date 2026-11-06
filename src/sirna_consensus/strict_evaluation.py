"""Confirmation stage: frozen inputs, fail-closed eligibility, fixed-budget utility."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

from .confirmation import (
    assay_metrics,
    cross_split_audit,
    decision,
    dump,
    load_protocol,
    macro_gene,
    paired_bootstrap,
    paired_budget_bounds,
    panel_bounds,
    require_confirmation_ids,
    summarize,
)
from .fasta import sha256_file


def verify_artifacts(lock):
    for name, sha in lock["artifacts"].items():
        path = Path(name)
        if not path.is_file() or sha256_file(path) != sha:
            raise ValueError("missing or changed frozen artifact: " + name)


def verify_baseline(baseline, directory):
    mapping = json.loads((directory / "historical_snapshots.json").read_text())["locks"][
        "platform_v3_freeze.json"
    ]
    bench = directory.parent
    for name, sha in baseline["artifacts"].items():
        original = Path(name)
        key = os.path.relpath(original, bench)
        path = bench / mapping[key] if key in mapping else original
        if not path.is_file() or sha256_file(path) != sha:
            raise ValueError("baseline artifact absent or mismatched: " + name)


def canonical_records(registry, protocol):
    rows = []
    for cohort in registry["cohorts"]:
        if cohort.get("role") != "confirmation" or not cohort.get("independent_eligible"):
            continue
        if sha256_file(Path(cohort["path"])) != cohort["sha256"]:
            raise ValueError("confirmation source checksum mismatch")
        with Path(cohort["path"]).open() as h:
            for raw in csv.DictReader(h, delimiter="\t"):
                if raw["record_id"] not in protocol["confirmation_ids"]:
                    continue
                row = dict(raw)
                row["guide"] = row["guide_dna"].upper().replace("U", "T")
                row["value"] = float(row["value"])
                row["active_threshold"] = protocol["active_thresholds"][row["endpoint"]]
                row["core"] = row["guide"][:19].translate(str.maketrans("ACGT", "TGCA"))[::-1]
                rows.append(row)
    return rows


def evaluate(protocol_path, candidate_freeze, prediction_path, output):
    protocol_path, candidate_freeze, prediction_path, output = map(
        Path, (protocol_path, candidate_freeze, prediction_path, output)
    )
    if output.exists():
        raise ValueError("preserve first evaluation; choose a new attempt, never overwrite")
    protocol, registry, baseline, registration = load_protocol(protocol_path)
    verify_baseline(baseline, protocol_path.parent)
    freeze = json.loads(candidate_freeze.read_text())
    verify_artifacts(freeze)
    if freeze["protocol_sha256"] != sha256_file(protocol_path):
        raise ValueError("candidate was not frozen against this protocol")
    predictions = json.loads(prediction_path.read_text())
    if predictions["candidate_freeze_sha256"] != sha256_file(candidate_freeze) or predictions[
        "protocol_sha256"
    ] != sha256_file(protocol_path):
        raise ValueError("predictions do not match frozen model/protocol")
    result = dict(
        schema_version=3,
        phase="confirmation",
        protocol_sha256=sha256_file(protocol_path),
        registration_commit=registration,
        candidate_freeze_sha256=sha256_file(candidate_freeze),
        predictions_sha256=sha256_file(prediction_path),
        baseline_commit=baseline["commit"],
        candidate_metrics=None,
        baseline_metrics=None,
        paired_ndcg=None,
        deployment=[],
        independent_validation_passed=False,
        development_metrics_used_for_gate=False,
    )
    evidence = dict(
        independence_verified=False,
        provenance_verified=True,
        independent_studies=0,
        independent_components=0,
        candidate_ndcg5=None,
        ndcg_gain=None,
        paired_ci95_lower=None,
        no_completion_loss=False,
        no_prediction_coverage_loss=False,
        actual_budget_gain_lower=None,
        independent_recomputation_passed=False,
    )
    if not protocol["confirmation_ids"]:
        if predictions.get("scores") or predictions.get("generated_runs"):
            raise ValueError("predictions supplied for excluded/nonexistent confirmation set")
        result.update(
            status="blocked_no_independent_confirmation",
            failures=[
                "all_existing_cohorts_excluded_or_unknown",
                "LOSO_not_permitted_under_Huesken_only",
                "actual_generated_panel_budget_gain_unconfirmed",
            ],
        )
    else:
        records = canonical_records(registry, protocol)
        require_confirmation_ids(records, registry, protocol)
        candidate_dir = Path(freeze["model_dir"])
        development = json.loads((candidate_dir / "development_records.json").read_text())
        baseline_records = next(
            Path(p)
            for p in baseline["artifacts"]
            if p.endswith("/development_20261106/development_records.json")
        )
        exposure = cross_split_audit(
            json.loads(baseline_records.read_text()) + development, records
        )
        if not exposure["all_clear"]:
            raise ValueError("cross-split leakage detected: " + json.dumps(exposure["findings"]))
        ids = {r["record_id"] for r in records}
        scoremap = {r["record_id"]: r for r in predictions["scores"]}
        if len(scoremap) != len(predictions["scores"]) or set(scoremap) != ids:
            raise ValueError("predictions must cover exact confirmation IDs including abstentions")
        a = assay_metrics(records, [scoremap[r["record_id"]]["candidate"] for r in records])
        b = assay_metrics(records, [scoremap[r["record_id"]]["baseline"] for r in records])
        from .learning import group_records

        groups, _ = group_records(records)
        gene_groups = {r["gene_id"]: int(g) for r, g in zip(records, groups)}
        delta = paired_bootstrap(
            a, b, gene_groups, protocol["bootstrap_replicates"], protocol["seed"]
        )
        result.update(
            status="evaluated_pending_independent_recomputation",
            candidate_metrics=summarize(a),
            baseline_metrics=summarize(b),
            paired_ndcg=delta,
            cross_split_audit=exposure,
        )
        by_unit = {r["unit"]: [x for x in records if x["unit"] == r["unit"]] for r in records}
        runmap = {r["unit"]: r for r in predictions["generated_runs"]}
        completion = []
        if set(runmap) == set(by_unit) and len(runmap) == len(predictions["generated_runs"]):
            for unit, measured in by_unit.items():
                panels = {}
                for arm in ("candidate", "baseline"):
                    entry = runmap[unit][arm]
                    if sha256_file(Path(entry["path"])) != entry["sha256"]:
                        raise ValueError("generated run changed after prediction freeze")
                    run = json.loads(Path(entry["path"]).read_text())
                    first = measured[0]
                    if (
                        run["transcript_id"] != first["transcript_id"]
                        or run["input_sha256"] != first["target_sha256"]
                    ):
                        raise ValueError("generated run target/assay mismatch")
                    expected = (
                        freeze["checkpoint_sha256"]
                        if arm == "candidate"
                        else freeze["baseline_checkpoint_sha256"]
                    )
                    if run["provenance"]["scorer"]["checkpoint_sha256"] != expected:
                        raise ValueError("wrong model used for generated panel")
                    by_key = {r["key"]: r for r in run["candidates"]}
                    panels[arm] = {}
                    for name in protocol["deployment"]["no_completion_loss_panels"]:
                        keys = [k for group in run[name].values() for k in group]
                        if len(keys) != len(set(keys)) or len(keys) > 5:
                            raise ValueError("invalid fixed-budget generated panel")
                        panel = []
                        for k in keys:
                            c = by_key[k]
                            panel.append(
                                dict(
                                    guide=c["guide_dna"],
                                    passenger=c["passenger"],
                                    chemistry=c["chemistry"],
                                    guide_overhang=c["guide_overhang"],
                                    passenger_overhang=c["passenger_overhang"],
                                    mapping_status="exact",
                                )
                            )
                        panels[arm][name] = panel
                name = protocol["deployment"]["panel"]
                bounds = paired_budget_bounds(
                    panels["candidate"][name], panels["baseline"][name], measured
                )
                completion.append(
                    all(
                        len(panels["candidate"][n]) >= len(panels["baseline"][n])
                        for n in panels["candidate"]
                    )
                )
                row = dict(
                    unit=unit,
                    gene_id=measured[0]["gene_id"],
                    budget_lower=bounds[0],
                    budget_upper=bounds[1],
                )
                for arm in panels:
                    stats = panel_bounds(panels[arm][name], measured)
                    row[arm] = {
                        k: v for k, v in stats.items() if k not in ("active", "unknown", "inactive")
                    }
                result["deployment"].append(row)
        evidence.update(
            independence_verified=True,
            independent_studies=len({r["study"] for r in records}),
            independent_components=delta["components"],
            candidate_ndcg5=result["candidate_metrics"]["ndcg5"],
            ndcg_gain=delta["delta"],
            paired_ci95_lower=delta["ci95"][0],
            no_completion_loss=bool(completion)
            and len(completion) == len(by_unit)
            and all(completion),
            no_prediction_coverage_loss=all(
                x["prediction_coverage"] >= y["prediction_coverage"] for x, y in zip(a, b)
            ),
            actual_budget_gain_lower=macro_gene(result["deployment"], "budget_lower"),
        )
    result["gate_evidence"] = evidence
    result["provisional_decision"] = decision(protocol, evidence)
    dump(output, result)
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--phase", choices=["confirmation"], default="confirmation")
    p.add_argument("--protocol", required=True)
    p.add_argument("--candidate-freeze", required=True)
    p.add_argument("--predictions", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args(argv)
    result = evaluate(a.protocol, a.candidate_freeze, a.predictions, a.output)
    print(json.dumps(dict(status=result["status"], independent_validation_passed=False)))
    return 0
