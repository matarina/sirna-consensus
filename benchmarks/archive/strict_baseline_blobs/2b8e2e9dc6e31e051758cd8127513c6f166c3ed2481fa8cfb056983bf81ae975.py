"""Final-run summariser: gene-weighted complete Top-5 aggregation."""
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from run_final import gene_weighted, needs_recut  # noqa: E402


def _result(gene, k, complete, precision):
    return {"gene": gene, "requested_k": k, "complete_measured_panel": complete,
            "precision_at_k": precision}


def test_gene_weighted_means_complete_panels_per_gene(tmp_path):
    path = tmp_path / "metrics.json"
    path.write_text(json.dumps({"results": [
        _result("A", 5, True, 1.0), _result("A", 5, True, 0.0),   # mean 0.5
        _result("B", 5, True, 0.4),
        _result("A", 10, True, 0.9),                              # different k -> ignored
        _result("C", 5, False, 0.7),                              # incomplete -> ignored
    ]}))
    assert gene_weighted(path) == {"A": 0.5, "B": 0.4}


def test_gene_weighted_empty(tmp_path):
    path = tmp_path / "m.json"
    path.write_text(json.dumps({"results": []}))
    assert gene_weighted(path) == {}


def test_needs_recut_on_summary_content_change(tmp_path):
    import hashlib
    import json as _json
    import os
    summary = tmp_path / "REVIEW.md"
    summary.write_text("v1")
    key = os.path.relpath(summary.resolve(), tmp_path.resolve())
    lock = tmp_path / "lock.json"

    # absent lock -> recut
    assert needs_recut(lock, {key}, summary) is True

    # lock covering the file with the CURRENT digest -> no recut
    lock.write_text(_json.dumps({"artifacts": {key: hashlib.sha256(b"v1").hexdigest()}}))
    assert needs_recut(lock, {key}, summary) is False

    # content changed (same key) -> recut again
    summary.write_text("v2")
    assert needs_recut(lock, {key}, summary) is True

    # a new key appears -> recut
    assert needs_recut(lock, {key, "other"}, summary) is True


def test_needs_recut_handles_absent_and_malformed(tmp_path):
    summary = tmp_path / "REVIEW.md"
    summary.write_text("v1")
    lock = tmp_path / "lock.json"
    key = os.path.relpath(summary.resolve(), tmp_path.resolve())
    # absent summary -> recut
    assert needs_recut(lock, {key}, tmp_path / "missing.md") is True
    # malformed lock -> recut, no crash
    lock.write_text("{ not json")
    assert needs_recut(lock, {key}, summary) is True
    # lock missing the 'artifacts' key -> recut
    lock.write_text("{}")
    assert needs_recut(lock, {key}, summary) is True


def test_needs_recut_on_absent_summary_with_matching_keys(tmp_path):
    # The set matches, so the absent-summary guard (not the set diff) must fire.
    summary = tmp_path / "REVIEW.md"
    lock = tmp_path / "lock.json"
    key = os.path.relpath(summary.resolve(), tmp_path.resolve())
    lock.write_text(json.dumps({"artifacts": {key: "deadbeef"}}))
    assert not summary.exists()
    assert needs_recut(lock, {key}, summary) is True
