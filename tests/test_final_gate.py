"""The final gate must reject an out-of-cohort or empty evaluation, not just a bad lock."""
import json
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "benchmarks"))
sys.path.insert(0, str(REPO / "src"))
import evaluate  # noqa: E402

from sirna_consensus.fasta import sha256_file  # noqa: E402

FMT = "guide21-passenger_unknown"
GUIDE = "AGCTTACGATCGTACCTGATA"


def _selection(tmp: Path) -> Path:
    run = tmp / "selection.json"
    run.write_text(json.dumps({
        "schema_version": 2, "kind": "exact_guide_selection", "config": {"top_k": 1},
        "candidates": [{"key": "0", "guide_dna": GUIDE, "duplex_format": FMT, "eligible": True,
                        "group": FMT + "/native", "score": 1.0}],
        "provenance": {"transcript_id": "TX", "target_sha256": "H"},
        "research_top": {FMT + "/native": ["0"]},
        "screened_top": {FMT + "/native": ["0"]},
    }))
    return run


def _measurements(tmp: Path, study: str, guide: str = GUIDE) -> Path:
    path = tmp / "measurements.tsv"
    cols = ["record_id", "study", "gene", "assay", "transcript_id", "target_sha256",
            "guide_dna", "duplex_format", "endpoint", "value"]
    row = ["r1", study, "G", "a", "TX", "H", guide, FMT, "percent_inhibition", "80"]
    path.write_text("\t".join(cols) + "\n" + "\t".join(row) + "\n")
    return path


def _cohort(tmp: Path, studies: list[str]) -> Path:
    path = tmp / "cohort.json"
    path.write_text(json.dumps({"independent_holdout_eligible": True,
                                "checkpoint_training_audit": "verified",
                                "holdout_studies": studies}))
    return path


def _protocol(tmp: Path) -> Path:
    path = tmp / "protocol.json"
    path.write_text(json.dumps({"seed": 2026}))
    return path


def _lock(tmp: Path, files: list[Path]) -> Path:
    lock = tmp / "lock.json"
    base = tmp.resolve()
    lock.write_text(json.dumps({"artifacts": {os.path.relpath(f.resolve(), base): sha256_file(f)
                                              for f in files}}))
    return lock


def _argv(tmp, lock, cohort, protocol, benchmark, run, out) -> list[str]:
    return ["--phase", "final", "--freeze-lock", str(lock), "--cohort-manifest", str(cohort),
            "--protocol", str(protocol), "--benchmark", str(benchmark), "--run", str(run),
            "--output", str(out), "--top-k", "5"]


def test_final_gate_accepts_in_cohort_benchmark(tmp_path):
    run, cohort, protocol = _selection(tmp_path), _cohort(tmp_path, ["S1"]), _protocol(tmp_path)
    bench = _measurements(tmp_path, "S1")
    lock = _lock(tmp_path, [bench, cohort, protocol, run])
    out = tmp_path / "final.json"
    assert evaluate.main(_argv(tmp_path, lock, cohort, protocol, bench, run, out)) == 0
    assert json.loads(out.read_text())["results"]


def test_final_gate_rejects_out_of_cohort_benchmark(tmp_path):
    run, cohort, protocol = _selection(tmp_path), _cohort(tmp_path, ["S1"]), _protocol(tmp_path)
    bench = _measurements(tmp_path, "S2")            # study not in the cohort
    lock = _lock(tmp_path, [bench, cohort, protocol, run])
    with pytest.raises(SystemExit):
        evaluate.main(_argv(tmp_path, lock, cohort, protocol, bench, run, tmp_path / "f.json"))


def test_final_gate_rejects_empty_evaluation(tmp_path):
    run, cohort, protocol = _selection(tmp_path), _cohort(tmp_path, ["S1"]), _protocol(tmp_path)
    empty = tmp_path / "empty.tsv"
    empty.write_text("\t".join(["record_id", "study", "gene", "assay", "transcript_id",
                                "target_sha256", "guide_dna", "duplex_format", "endpoint", "value"]) + "\n")
    lock = _lock(tmp_path, [empty, cohort, protocol, run])
    with pytest.raises(SystemExit):
        evaluate.main(_argv(tmp_path, lock, cohort, protocol, empty, run, tmp_path / "f.json"))


def test_final_gate_rejects_flag_only_cohort(tmp_path):
    run, protocol = _selection(tmp_path), _protocol(tmp_path)
    cohort = tmp_path / "cohort.json"          # two-field cohort, no holdout_studies
    cohort.write_text(json.dumps({"independent_holdout_eligible": True,
                                  "checkpoint_training_audit": "verified"}))
    bench = _measurements(tmp_path, "S1")
    lock = _lock(tmp_path, [bench, cohort, protocol, run])
    with pytest.raises(SystemExit):
        evaluate.main(_argv(tmp_path, lock, cohort, protocol, bench, run, tmp_path / "f.json"))


def test_final_gate_tolerates_absent_extra_lock_entry(tmp_path):
    # A pinned entry that is absent here (e.g. a gitignored run.json on a bare clone)
    # must not crash the gate, unlike a mismatched present file.
    run, cohort, protocol = _selection(tmp_path), _cohort(tmp_path, ["S1"]), _protocol(tmp_path)
    bench = _measurements(tmp_path, "S1")
    lock = _lock(tmp_path, [bench, cohort, protocol, run])
    data = json.loads(lock.read_text())
    data["artifacts"]["ghost/run.json"] = "deadbeef"
    lock.write_text(json.dumps(data))
    out = tmp_path / "f.json"
    assert evaluate.main(_argv(tmp_path, lock, cohort, protocol, bench, run, out)) == 0
    assert json.loads(out.read_text())["results"]
