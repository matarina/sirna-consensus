"""Independent verifier for the preregistered Huesken-only attempt."""

from __future__ import annotations

import csv
import itertools
import json
import math
import pickle
from collections import defaultdict
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

from sirna_consensus.fasta import sha256_file


def encode(guides, version):
    rows = []
    for value in guides:
        g = value.upper().replace("U", "T")[:19]
        if len(g) != 19 or set(g) - set("ACGT"):
            raise ValueError("invalid verification guide")
        composition = [g.count(b) / 19 for b in "ACGT"]
        dinucs = ["".join(x) for x in itertools.product("ACGT", repeat=2)]
        di = [sum(g[i : i + 2] == m for i in range(18)) / 18 for m in dinucs]

        def gc(s):
            return sum(s.count(b) for b in "GC") / len(s)

        run = max(len(list(v)) for _, v in itertools.groupby(g)) / 19
        row = (
            [int(g[i] == b) for i in range(19) for b in "ACGT"]
            + composition
            + di
            + [gc(g), gc(g[1:8]), gc(g[:4]) - gc(g[-4:]), run]
        )
        if version == "sequence_v2":
            row += [int(g[i : i + 2] == m) for i in range(18) for m in dinucs]
            row += [
                sum(g[i : i + 3] == "".join(m) for i in range(17)) / 17
                for m in itertools.product("ACGT", repeat=3)
            ]
        rows.append(row)
    return np.asarray(rows, dtype=float)


def metrics(records, scores):
    units = defaultdict(list)
    for r, s in zip(records, scores):
        units[r["unit"]].append((r, s))
    result = []
    for unit, items in sorted(units.items()):
        vals = [r["value"] for r, _ in items]
        relevance = {
            r["record_id"]: (
                sum(v < r["value"] for v in vals) + 0.5 * (sum(v == r["value"] for v in vals) - 1)
            )
            / (len(vals) - 1)
            if len(vals) > 1
            else 0.0
            for r, _ in items
        }
        ordered = sorted(
            [(r, float(s)) for r, s in items if s is not None and math.isfinite(float(s))],
            key=lambda p: (-p[1], p[0]["guide"], p[0].get("molecule_id", ""), p[0]["record_id"]),
        )
        top = [r for r, _ in ordered[:5]]
        ndcg = None
        if len(items) >= 5 and len(set(vals)) > 1:
            dcg = sum(
                math.expm1(math.log(2) * relevance[r["record_id"]]) / math.log2(i + 2)
                for i, r in enumerate(top)
            )
            ideal = sum(
                math.expm1(math.log(2) * v) / math.log2(i + 2)
                for i, v in enumerate(sorted(relevance.values(), reverse=True)[:5])
            )
            ndcg = dcg / ideal
        active = sum(r["value"] >= r["active_threshold"] for r in top)
        rho = None
        if (
            len(ordered) > 1
            and len({r["value"] for r, _ in ordered}) > 1
            and len({s for _, s in ordered}) > 1
        ):
            rho = float(
                spearmanr([r["value"] for r, _ in ordered], [s for _, s in ordered]).statistic
            )
        result.append(
            dict(
                unit=unit,
                gene_id=items[0][0]["gene_id"],
                ndcg5=ndcg,
                precision5=active / len(top) if top else None,
                top5_active_fraction=active / len(top) if top else None,
                top3_at_least_one_active=float(
                    any(r["value"] >= r["active_threshold"] for r in top[:3])
                ),
                fill_rate=len(top) / 5,
                prediction_coverage=len(ordered) / len(items),
                budget_yield=active / 5,
                measurement_coverage=1.0 if top else None,
                spearman=rho,
            )
        )
    summary = {}
    for key in (
        "ndcg5",
        "precision5",
        "top5_active_fraction",
        "top3_at_least_one_active",
        "fill_rate",
        "prediction_coverage",
        "budget_yield",
        "measurement_coverage",
        "spearman",
    ):
        gene = defaultdict(list)
        for r in result:
            if r[key] is not None:
                gene[r["gene_id"]].append(r[key])
        summary[key] = mean(mean(v) for v in gene.values()) if gene else None
    return summary


def verify_development(model_dir):
    p = Path(model_dir)
    data = json.loads((p / "development_records.json").read_text())
    split = json.loads((p / "split_manifest.json").read_text())
    predictions = json.loads((p / "oof_predictions.json").read_text())
    fits = json.loads((p / "fit_audit.json").read_text())
    result = json.loads((p / "results.json").read_text())
    index = {r["record_id"]: i for i, r in enumerate(data)}
    scores = {r["record_id"]: r["score"] for r in predictions}
    checks = dict(
        huesken_only=all(r["study"] == "huesken2005" for r in data),
        exact_ids=len(scores) == len(predictions) == len(data) == len(index)
        and set(scores) == set(index),
        training_ledger_huesken_only=all(
            f["studies"] == ["huesken2005"] and set(f["train_ids"]) <= set(index) for f in fits
        ),
    )
    protocol_path = Path(__file__).resolve().parents[1] / "confirmation_20261106/protocol.json"
    protocol = json.loads(protocol_path.read_text())
    source = Path(protocol["training"]["source"])
    with source.open() as handle:
        original = {r["record_id"]: r for r in csv.DictReader(handle, delimiter="\t")}
    approved, alternatives = {}, defaultdict(set)
    with Path(protocol["training"]["hgnc"]).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["status"] != "Approved":
                continue
            approved[row["symbol"]] = row["hgnc_id"]
            alternatives[row["symbol"].casefold()].add(row["hgnc_id"])
            for field in ("prev_symbol", "alias_symbol"):
                for alias in row[field].split("|"):
                    if alias:
                        alternatives[alias.casefold()].add(row["hgnc_id"])
    checks["pinned_training_sources"] = (
        sha256_file(source) == protocol["training"]["source_sha256"]
        and sha256_file(Path(protocol["training"]["hgnc"])) == protocol["training"]["hgnc_sha256"]
        and result["protocol_sha256"] == sha256_file(protocol_path)
    )
    checks["authoritative_training_ids_values_sequences_genes"] = set(original) == set(
        index
    ) and all(
        r["value"] == float(original[r["record_id"]]["value"])
        and r["guide"] == original[r["record_id"]]["guide_dna"].upper().replace("U", "T")
        and r["study"] == original[r["record_id"]]["study"]
        and r["active_threshold"]
        == protocol["active_thresholds"][original[r["record_id"]]["endpoint"]]
        and (
            {approved[original[r["record_id"]]["gene"]]}
            if original[r["record_id"]]["gene"] in approved
            else alternatives[original[r["record_id"]]["gene"].casefold()]
        )
        == {r["gene_id"]}
        for r in data
    )
    replay_error = weight_error = feature_error = 0.0
    covered = []
    groups = split["groups"]
    bad17 = []
    bad_hamming = []
    words = {}
    signatures = {}
    for r in data:
        s = r["core"]
        rc = s.translate(str.maketrans("ACGT", "TGCA"))[::-1]
        for seq in (s, rc):
            for i in range(3):
                word = seq[i : i + 17]
                if word in words and words[word] != groups[r["record_id"]]:
                    bad17.append(r["record_id"])
                words[word] = groups[r["record_id"]]
            for a, b in itertools.combinations(range(19), 2):
                key = (a, b, seq[:a] + seq[a + 1 : b] + seq[b + 1 :])
                if key in signatures and signatures[key] != groups[r["record_id"]]:
                    bad_hamming.append(r["record_id"])
                signatures[key] = groups[r["record_id"]]
    checks.update(shared17_grouped=not bad17, Hamming2_grouped=not bad_hamming)
    gene_groups = defaultdict(set)
    for r in data:
        gene_groups[r["gene_id"]].add(groups[r["record_id"]])
    checks["same_gene_grouped"] = all(len(x) == 1 for x in gene_groups.values())
    from sirna_consensus.challenger import features

    for version in ("sequence_v1", "sequence_v2"):
        feature_error = max(
            feature_error,
            float(
                np.max(
                    np.abs(
                        encode([r["guide"] for r in data], version)
                        - features([r["guide"] for r in data], version)
                    )
                )
            ),
        )
    with threadpool_limits(limits=4):
        for fold in split["folds"]:
            tr, te = fold["train_ids"], fold["test_ids"]
            covered += te
            checks["outer_" + str(fold["fold"]) + "_disjoint"] = not (set(tr) & set(te)) and not (
                {groups[i] for i in tr} & {groups[i] for i in te}
            )
            checks["outer_" + str(fold["fold"]) + "_ledger"] = all(
                set(f["train_ids"]) <= set(tr)
                for f in fits
                if f["context"] == "outer_" + str(fold["fold"])
            )
            rows = [data[index[i]] for i in tr]
            genes = {r["gene_id"] for r in rows}
            units = defaultdict(list)
            conditions = defaultdict(set)
            for r in rows:
                units[r["unit"]].append(r)
                conditions[r["gene_id"]].add(r["unit"])
            w = np.array(
                [
                    1 / (len(genes) * len(conditions[r["gene_id"]]) * len(units[r["unit"]]))
                    for r in rows
                ]
            )
            w /= w.mean()
            weight_error = max(
                weight_error, float(np.max(np.abs(w - np.asarray(fold["train_weights"]))))
            )
            ck = p / fold["checkpoint"]
            assert sha256_file(ck) == fold["checkpoint_sha256"]
            bundle = pickle.loads(ck.read_bytes())
            matrix = encode([data[index[i]]["guide"] for i in te], bundle["params"]["feature"])
            if bundle["scaler"] is not None:
                matrix = bundle["scaler"].transform(matrix)
            estimator = bundle["model"]
            out = (
                estimator.decision_function(matrix)
                if bundle["params"]["family"] == "pairwise"
                else estimator.predict(matrix)
            )
            replay_error = max(
                replay_error, max(abs(float(v) - scores[i]) for i, v in zip(te, out))
            )
    checks["exact_outer_coverage"] = len(covered) == len(set(covered)) == len(data) and set(
        covered
    ) == set(index)
    independently = metrics(data, [scores[r["record_id"]] for r in data])
    errors = {}
    for k, v in independently.items():
        expected = result["metrics"][k]
        errors[k] = (
            abs(v - expected)
            if v is not None and expected is not None
            else (0.0 if v is expected else float("inf"))
        )
    checks.update(
        independent_features=feature_error < 1e-12,
        equal_gene_assay_weights=weight_error < 1e-12,
        checkpoint_reproduction=replay_error < 1e-12,
        independent_metrics=max(errors.values()) < 1e-12,
    )
    report = dict(
        kind="independent_development_arithmetic_verification",
        all_passed=all(checks.values()),
        checks=checks,
        feature_max_error=feature_error,
        weight_max_error=weight_error,
        checkpoint_max_error=replay_error,
        metric_errors=errors,
        independent_metrics=independently,
        fit_operations=len(fits),
        independent_efficacy_confirmation=False,
        results_sha256=sha256_file(p / "results.json"),
    )
    (p / "independent_verification.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


def verify_confirmation(evaluation_path, protocol_path, freeze_path, predictions_path, output_path):
    """Independent verification of this round's zero-ID, blocked result.

    Arithmetic correctness is a separate fact from efficacy confirmation. A future
    nonempty registration must supply its own independent paired-metric verifier;
    this routine refuses it rather than accepting a production evaluator's claims.
    """
    evaluation_path, protocol_path, freeze_path, predictions_path, output_path = map(
        Path, (evaluation_path, protocol_path, freeze_path, predictions_path, output_path)
    )
    if output_path.exists():
        raise ValueError("preserve first independent verification")

    def read(p):
        return json.loads(p.read_text())

    p, e, f, predictions = map(
        read, (protocol_path, evaluation_path, freeze_path, predictions_path)
    )
    if p["confirmation_ids"]:
        raise ValueError(
            "independent confirmation verifier is specific to the zero-ID preregistration"
        )
    registry = read(protocol_path.with_name("cohorts.json"))
    baseline = read(protocol_path.with_name("baseline.json"))
    checks = dict(
        frozen_protocol_metadata=sha256_file(protocol_path.with_name("cohorts.json"))
        == p["cohort_registry_sha256"]
        and sha256_file(protocol_path.with_name("baseline.json")) == p["baseline_manifest_sha256"],
        candidate_artifacts_unchanged=all(
            Path(path).is_file() and sha256_file(Path(path)) == sha
            for path, sha in f["artifacts"].items()
        ),
        evaluation_links=e["protocol_sha256"] == sha256_file(protocol_path)
        and e["candidate_freeze_sha256"] == sha256_file(freeze_path)
        and e["predictions_sha256"] == sha256_file(predictions_path),
        prediction_links=predictions["candidate_freeze_sha256"] == sha256_file(freeze_path)
        and predictions["protocol_sha256"] == sha256_file(protocol_path),
        authoritative_cohort_sources_unchanged=all(
            sha256_file(Path(c["path"])) == c["sha256"] for c in registry["cohorts"]
        ),
        zero_eligible_records=registry["eligible_confirmation_records"] == 0
        and not any(c.get("independent_eligible") for c in registry["cohorts"])
        and registry["confirmation_ids"] == [],
        no_excluded_cohort_predictions=predictions["scores"] == []
        and predictions["generated_runs"] == [],
        no_manufactured_performance=all(
            e[k] is None for k in ("candidate_metrics", "baseline_metrics", "paired_ndcg")
        )
        and e["deployment"] == [],
        correct_blocked_status=e["status"] == "blocked_no_independent_confirmation",
        no_independent_efficacy_claim=e["independent_validation_passed"] is False,
        development_not_used_for_confirmation=e["development_metrics_used_for_gate"] is False,
        original_default_retained=sha256_file(Path(baseline["default_policy_path"]))
        == baseline["default_policy_sha256"],
        no_false_LOSO=p["loso"]["status"] == "not_permitted_single_development_study",
    )
    expected = dict(
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
    checks["independently_derived_gate_evidence"] = e["gate_evidence"] == expected
    provisional = e["provisional_decision"]
    checks["gate_rejects_without_data"] = (
        provisional["promote"] is False and provisional["action"] == "retain_baseline"
    )
    report = dict(
        kind="independent_zero_confirmation_verification",
        all_passed=all(checks.values()),
        checks=checks,
        independent_efficacy_confirmation=False,
        eligible_confirmation_records=0,
        evaluation_sha256=sha256_file(evaluation_path),
        protocol_sha256=sha256_file(protocol_path),
        limitation="No independent efficacy metric exists. Verification pass confirms a blocked result only.",
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    return report
