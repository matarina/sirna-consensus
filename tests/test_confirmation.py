"""Adversarial protocol checks and hand-computable evaluation cases."""

import copy
import json

import pytest

from sirna_consensus.challenger import features, fit, predict_matrix
from sirna_consensus.confirmation import (
    assay_metrics,
    cross_split_audit,
    decision,
    ensure_training_records,
    macro_gene,
    molecule_identity,
    paired_bootstrap,
    paired_budget_bounds,
    panel_bounds,
    require_confirmation_ids,
)
from sirna_consensus.types import reverse_complement

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")


def records(gene="G1", count=6):
    return [
        dict(
            record_id=f"{gene}_{i}",
            gene_id=gene,
            study="huesken2005",
            unit=gene + "|assay",
            guide=("A" * 18 + "ACGT"[i % 4]),
            core=reverse_complement("A" * 18 + "ACGT"[i % 4]),
            value=float(i),
            active_threshold=4.0,
            target_rank=i / (count - 1),
            condition="assay",
        )
        for i in range(count)
    ]


def protocol():
    return dict(
        training={"studies": ["huesken2005"]},
        seed=20261106,
        gate=dict(
            minimum_independent_studies=2,
            minimum_independent_components=5,
            ndcg_strictly_above=0.8,
            gain_strictly_above=0.0,
            paired_ci95_lower_strictly_above=0.0,
            actual_generated_budget_gain_lower_strictly_above=0.0,
        ),
    )


def molecule(base="A", chemistry="all_RNA", value=None, unit="assay"):
    sense = base * 19
    result = dict(
        guide=reverse_complement(sense) + "AC",
        passenger=sense + "GT",
        guide_overhang="AC",
        passenger_overhang="GT",
        chemistry=chemistry,
        mapping_status="exact",
        unit=unit,
    )
    if value is not None:
        result.update(value=value, active_threshold=70.0)
    return result


def test_ndcg_ideal_abstention_and_fixed_assay_pool():
    rows = records()
    ideal = assay_metrics(rows, list(range(6)))[0]
    assert ideal["ndcg5"] == pytest.approx(1.0)
    assert ideal["precision5"] == pytest.approx(2 / 5)
    assert ideal["top3_at_least_one_active"] == 1.0
    empty = assay_metrics(rows, [None] * 6)[0]
    assert empty["ndcg5"] == 0.0 and empty["fill_rate"] == 0.0 and empty["abstentions"] == 6
    assert empty["precision5"] is None
    assert json.loads(json.dumps(empty))["ndcg_eligible"] is True


def test_gene_macro_does_not_overweight_genes_with_more_assays():
    rows = [
        dict(gene_id="A", score=1.0),
        dict(gene_id="A", score=1.0),
        dict(gene_id="B", score=0.0),
    ]
    assert macro_gene(rows, "score") == 0.5


def test_constant_endpoint_is_not_invented_perfect_ranking():
    rows = records()
    for r in rows:
        r["value"] = 1.0
    assert assay_metrics(rows, list(range(6)))[0]["ndcg5"] is None


def test_duplicate_and_nonfinite_truth_rejected():
    rows = records()
    rows[1]["record_id"] = rows[0]["record_id"]
    with pytest.raises(ValueError, match="duplicate"):
        assay_metrics(rows, list(range(6)))
    rows = records()
    rows[0]["value"] = float("nan")
    with pytest.raises(ValueError, match="nonfinite"):
        assay_metrics(rows, list(range(6)))


def test_ties_do_not_use_experimental_labels():
    rows = records()
    a = assay_metrics(rows, [1.0] * 6)[0]["top_ids"]
    assert a == assay_metrics(list(reversed(rows)), [1.0] * 6)[0]["top_ids"]


def test_unknown_activity_and_chemistry_are_not_negative_labels():
    a = molecule()
    truth = molecule(value=80.0)
    assert panel_bounds([a], [truth])["precision5"] == 1.0
    a["chemistry"] = "unknown"
    stats = panel_bounds([a], [truth])
    assert stats["precision5"] is None and stats["yield_bounds"] == [0.0, 0.2]
    assert stats["top3_at_least_one_active"] is None


def test_full_strand_mismatch_cannot_inherit_core_activity():
    a = molecule()
    truth = molecule(value=80.0)
    a["passenger"] = a["passenger"][:-2] + "CC"
    a["passenger_overhang"] = "CC"
    assert panel_bounds([a], [truth])["yield_bounds"] == [0.0, 0.2]
    assert molecule_identity(molecule(chemistry="unknown")) is None


def test_unknown_shared_molecule_cancels_but_unknown_identity_does_not():
    a = molecule()
    assert paired_budget_bounds([a], [a], []) == [0.0, 0.0]
    a = molecule(chemistry="unknown")
    assert paired_budget_bounds([a], [a], []) == [-0.2, 0.2]
    assert paired_budget_bounds([molecule()], [molecule("C")], []) == [-0.2, 0.2]


def test_observed_positive_can_prove_positive_budget_change():
    truth = [molecule(value=80.0), molecule("C", value=0.0)]
    assert paired_budget_bounds([molecule()], [molecule("C")], truth) == [0.2, 0.2]
    with pytest.raises(ValueError, match="assay"):
        panel_bounds([molecule()], [*truth, molecule("G", value=20.0, unit="other")])


def test_panel_budget_and_duplicate_checks():
    with pytest.raises(ValueError):
        panel_bounds([molecule()] * 6, [])
    with pytest.raises(ValueError):
        panel_bounds([molecule(), molecule()], [])


def test_huesken_allowlist_rejects_external_fitting():
    rows = records()
    rows[0]["study"] = "sciabola2013"
    with pytest.raises(ValueError, match="unapproved"):
        ensure_training_records(rows, protocol())
    with pytest.raises(ValueError, match="confirmation"):
        fit(
            dict(feature="sequence_v1", family="ridge", alpha=1.0),
            {"sequence_v1": np.zeros((6, 100))},
            rows,
            np.arange(6),
            protocol(),
        )


def test_excluded_ids_cannot_be_promoted_to_confirmation():
    with pytest.raises(ValueError, match="eligible"):
        require_confirmation_ids(records(), {"cohorts": []}, {"confirmation_ids": []})


def test_sequence_audit_detects_reverse_complement_shift_and_gene(random_dna):
    s = random_dna(19)
    training = [dict(core=s, gene_id="G")]
    tests = [
        dict(record_id="rc", core=reverse_complement(s), gene_id="H"),
        dict(record_id="shift", core="A" + s[:18], gene_id="J"),
        dict(record_id="gene", core=random_dna(19), gene_id="G"),
    ]
    result = cross_split_audit(training, tests)
    assert not result["all_clear"] and {x["record_id"] for x in result["findings"]} == {
        "rc",
        "shift",
        "gene",
    }


def test_gate_requires_every_condition_and_strictly_above_point_eight():
    p = protocol()
    e = dict(
        independence_verified=True,
        provenance_verified=True,
        independent_studies=2,
        independent_components=5,
        candidate_ndcg5=0.81,
        ndcg_gain=0.03,
        paired_ci95_lower=0.01,
        no_completion_loss=True,
        no_prediction_coverage_loss=True,
        actual_budget_gain_lower=0.01,
        independent_recomputation_passed=True,
    )
    assert decision(p, e)["promote"]
    assert not decision(p, {**e, "candidate_ndcg5": 0.8})["promote"]
    for key in e:
        v = copy.deepcopy(e)
        v.pop(key)
        assert not decision(p, v)["promote"], key
    assert not decision(p, {})["independent_validation_passed"]


def test_bootstrap_pairs_assays_and_keeps_all_genes():
    a = [dict(unit="a", gene_id="A", ndcg5=0.9), dict(unit="b", gene_id="B", ndcg5=0.8)]
    b = [dict(unit="a", gene_id="A", ndcg5=0.7), dict(unit="b", gene_id="B", ndcg5=0.6)]
    result = paired_bootstrap(a, b, {"A": 0, "B": 1}, replicates=100)
    assert result["delta"] == pytest.approx(0.2) and result["ci95"] == pytest.approx([0.2, 0.2])
    with pytest.raises(ValueError, match="identical"):
        paired_bootstrap(a, b[:1], {"A": 0})


def test_both_feature_groups_and_pairwise_model_are_real(random_dna):
    rows = records(count=10)
    for r in rows:
        r["guide"] = random_dna(19)
        r["core"] = reverse_complement(r["guide"])
    assert features([r["guide"] for r in rows], "sequence_v1").shape == (10, 100)
    x = {"sequence_v2": features([r["guide"] for r in rows], "sequence_v2")}
    assert x["sequence_v2"].shape == (10, 452)
    model = fit(
        dict(feature="sequence_v2", family="pairwise", C=0.1), x, rows, np.arange(10), protocol()
    )
    assert np.isfinite(predict_matrix(model, x["sequence_v2"])).all()
