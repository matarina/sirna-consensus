"""Small synthetic checks: training weights and sequence isolation, not fitted accuracy."""

import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("sklearn")
from sirna_consensus.learning import features, folds, group_records, weights  # noqa: E402


def record(study, gene, condition, seq):
    return {
        "study": study,
        "gene_id": gene,
        "condition": condition,
        "unit": "|".join((study, gene, condition)),
        "core": seq,
        "record_id": study + gene + condition + seq,
    }


def test_study_gene_condition_weights_are_balanced():
    rows = [record("A", "g1", "c1", "A" * 19)] * 10 + [record("A", "g2", "c1", "C" * 19)] * 2
    rows += [record("B", "g3", "c1", "G" * 19)] * 3 + [record("B", "g3", "c2", "T" * 19)] * 5
    w = weights(rows)
    assert w[:12].sum() == pytest.approx(w[12:].sum())
    assert w[:10].sum() == pytest.approx(w[10:12].sum())
    assert w[12:15].sum() == pytest.approx(w[15:].sum())


def test_gene_and_shifted_sequence_neighbors_cannot_cross_folds(random_dna):
    a = random_dna(19, seed=12)
    b = random_dna(19, seed=13)
    c = random_dna(19, seed=14)
    rows = [
        record("A", "g1", "c1", a),
        record("B", "g1", "c2", b),
        record("A", "g2", "c1", a[1:] + "A"),
        record("A", "g3", "c1", c),
    ]
    groups, _ = group_records(rows)
    assert groups[0] == groups[1] == groups[2] and groups[3] != groups[0]
    for tr, te in folds(np.arange(len(rows)), groups, 5):
        assert not set(groups[tr]) & set(groups[te])


def test_features_reject_incomplete_or_ambiguous_paired_regions():
    with pytest.raises(ValueError):
        features(["A" * 18])
    with pytest.raises(ValueError):
        features(["A" * 18 + "N"])
    assert features(["A" * 19 + "TT"]).shape == features(["A" * 19 + "CC"]).shape
    assert np.array_equal(features(["A" * 19 + "TT"]), features(["A" * 19 + "CC"]))
