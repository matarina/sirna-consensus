"""Adversarial end-to-end checks of blocked confirmation and independent arithmetic."""

import copy
import importlib.util
import json
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")

from sirna_consensus import strict_evaluation  # noqa: E402
from sirna_consensus.confirmation import assay_metrics, dump, summarize  # noqa: E402
from sirna_consensus.fasta import sha256_file  # noqa: E402

SCRIPT = Path(__file__).resolve().parents[1] / "benchmarks/scripts/verify_guarded.py"
spec = importlib.util.spec_from_file_location("independent_guarded", SCRIPT)
independent = importlib.util.module_from_spec(spec)
spec.loader.exec_module(independent)


def test_independent_metrics_agree_with_tied_multigene_assays_and_abstentions():
    rng = np.random.default_rng(49)
    records, scores = [], []
    for gene, assays in [("A", 3), ("B", 1), ("C", 2)]:
        for assay in range(assays):
            for i in range(11):
                records.append(
                    dict(
                        record_id=f"{gene}:{assay}:{i}",
                        study="toy",
                        gene_id=gene,
                        unit=f"{gene}:{assay}",
                        value=float(rng.integers(0, 5)),
                        active_threshold=3.0,
                        guide=f"{i:019d}",
                    )
                )
                scores.append(None if i % 4 == 0 else float(rng.integers(0, 7)))
    production = summarize(assay_metrics(records, scores))
    separately = independent.metrics(records, scores)
    for key, value in separately.items():
        assert production[key] == pytest.approx(value, abs=1e-12)


@pytest.fixture
def blocked(tmp_path, monkeypatch):
    protocol_path = tmp_path / "protocol.json"
    cohort_source = tmp_path / "source.tsv"
    cohort_source.write_text("record_id\n")
    policy = tmp_path / "policy.json"
    policy.write_text('{"default":"existing"}\n')
    baseline = dict(
        commit="frozen",
        artifacts={},
        default_policy_path=str(policy),
        default_policy_sha256=sha256_file(policy),
    )
    registry = dict(
        cohorts=[
            dict(
                name="exposed",
                path=str(cohort_source),
                sha256=sha256_file(cohort_source),
                independent_eligible=False,
            )
        ],
        eligible_confirmation_records=0,
        confirmation_ids=[],
    )
    dump(tmp_path / "baseline.json", baseline)
    dump(tmp_path / "cohorts.json", registry)
    protocol = dict(
        confirmation_ids=[],
        baseline_manifest_sha256=sha256_file(tmp_path / "baseline.json"),
        cohort_registry_sha256=sha256_file(tmp_path / "cohorts.json"),
        loso=dict(status="not_permitted_single_development_study"),
        gate=dict(
            minimum_independent_studies=2,
            minimum_independent_components=5,
            ndcg_strictly_above=0.8,
            gain_strictly_above=0.0,
            paired_ci95_lower_strictly_above=0.0,
            actual_generated_budget_gain_lower_strictly_above=0.0,
        ),
    )
    dump(protocol_path, protocol)
    dump(tmp_path / "historical_snapshots.json", dict(locks={"platform_v3_freeze.json": {}}))
    checkpoint = tmp_path / "checkpoint"
    checkpoint.write_bytes(b"frozen candidate")
    freeze = tmp_path / "freeze.json"
    dump(
        freeze,
        dict(
            artifacts={str(checkpoint): sha256_file(checkpoint)},
            protocol_sha256=sha256_file(protocol_path),
        ),
    )
    predictions = tmp_path / "predictions.json"
    dump(
        predictions,
        dict(
            candidate_freeze_sha256=sha256_file(freeze),
            protocol_sha256=sha256_file(protocol_path),
            scores=[],
            generated_runs=[],
        ),
    )
    monkeypatch.setattr(
        strict_evaluation, "load_protocol", lambda _: (protocol, registry, baseline, "registered")
    )
    return protocol_path, freeze, predictions, tmp_path / "evaluation.json"


def test_empty_confirmation_remains_unknown_and_cannot_promote(blocked):
    protocol, freeze, predictions, output = blocked
    result = strict_evaluation.evaluate(protocol, freeze, predictions, output)
    assert result["candidate_metrics"] is None
    assert result["provisional_decision"]["promote"] is False
    verification = independent.verify_confirmation(
        output, protocol, freeze, predictions, output.with_name("verified.json")
    )
    assert verification["all_passed"]
    assert verification["independent_efficacy_confirmation"] is False
    with pytest.raises(ValueError, match="preserve first evaluation"):
        strict_evaluation.evaluate(protocol, freeze, predictions, output)


@pytest.mark.parametrize(
    "attack", ["manufactured_gain", "fake_independence", "false_promote", "default_switch"]
)
def test_independent_verifier_catches_forged_success(blocked, attack):
    protocol, freeze, predictions, output = blocked
    result = strict_evaluation.evaluate(protocol, freeze, predictions, output)
    if attack == "manufactured_gain":
        result["candidate_metrics"] = dict(ndcg5=0.99)
    elif attack == "fake_independence":
        result["gate_evidence"]["independence_verified"] = True
    elif attack == "false_promote":
        result["provisional_decision"]["promote"] = True
    else:
        baseline = json.loads(protocol.with_name("baseline.json").read_text())
        Path(baseline["default_policy_path"]).write_text("changed")
    dump(output, result)
    verification = independent.verify_confirmation(
        output, protocol, freeze, predictions, output.with_name("verified.json")
    )
    assert not verification["all_passed"]


def test_excluded_predictions_and_model_mutation_fail_before_evaluation(blocked):
    protocol, freeze, predictions, output = blocked
    payload = json.loads(predictions.read_text())
    original = copy.deepcopy(payload)
    payload["scores"] = [dict(record_id="excluded", candidate=0.99, baseline=0.0)]
    dump(predictions, payload)
    with pytest.raises(ValueError, match="excluded/nonexistent"):
        strict_evaluation.evaluate(protocol, freeze, predictions, output)
    assert not output.exists()
    dump(predictions, original)
    artifact = next(iter(json.loads(freeze.read_text())["artifacts"]))
    Path(artifact).write_text("different model")
    with pytest.raises(ValueError, match="frozen artifact"):
        strict_evaluation.evaluate(protocol, freeze, predictions, output)
