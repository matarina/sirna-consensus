"""Exposure audit: training vs evaluation exposure must not be conflated."""
import copy
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
from audit_exposure import build_audit, core, load_molecules, normalize, write_outputs  # noqa: E402

TRAIN = "AAATCAATTAACATATTAG"
DEV = "CTAATATGTTAATTGATTT"
EVAL = "GCTTACGATCGTACCTGAT"
FREE = "TTTTGGGGCCCCAAAATTT"
EVAL2 = "ACGACGACGACGACGACGA"


def snap(study, gene, sequence):
    return {"Authors": study, "Gene": gene, "Antisense_21mer": sequence + "at"}


def records():
    return [
        snap("TrainStudy", "G1", TRAIN), snap("TrainStudy", "G1", TRAIN),
        snap("DevStudy", "G2", DEV),
        snap("EvalStudy", "G3", EVAL),
        snap("FreeStudy", "G4", FREE),
        snap("FreeStudy", "G5", EVAL2),
    ]


def audit(**kw):
    kw.setdefault("training", {TRAIN})
    kw.setdefault("development", {DEV})
    kw.setdefault("upstream_eval", {EVAL})
    return build_audit(records(), kw["training"], kw["development"], kw["upstream_eval"],
                       criterion=kw.get("criterion", "weight"),
                       minimum_complete_paired_genes=kw.get("minimum_complete_paired_genes", 2))


def test_core_normalization():
    assert normalize("acgu") == "ACGT"
    assert core("CUAAUAUGUUAAUUGAUUUau") == "CTAATATGTTAATTGATTT"


def test_training_and_development_exposure_exclude_a_study():
    result = audit()
    assert result["per_study"]["TrainStudy"]["classification"] == "training_exposed"
    assert result["per_study"]["DevStudy"]["classification"] == "development_exposed"
    assert "TrainStudy" not in result["holdout_studies"]
    assert "DevStudy" not in result["holdout_studies"]


def test_eval_exposure_kept_under_weight_criterion_but_not_fully():
    weight = audit(criterion="weight")
    fully = audit(criterion="fully")
    assert "EvalStudy" in weight["holdout_studies"]
    assert "EvalStudy" not in fully["holdout_studies"]
    assert "FreeStudy" in fully["holdout_studies"]
    assert weight["per_study"]["EvalStudy"]["classification"] == "upstream_evaluation_exposed"


def test_eligibility_requires_enough_genes():
    ok = audit(minimum_complete_paired_genes=2)          # Free(2 genes)+Eval(1 gene) >= 2
    tight = audit(minimum_complete_paired_genes=99)
    assert ok["independent_holdout_eligible"] is True
    assert tight["independent_holdout_eligible"] is False


def test_absence_of_training_source_degrades_to_unverified():
    result = audit(training=set())
    assert result["checkpoint_training_audit"] == "unverified"
    # dev exposure still excludes DevStudy; the training-only study becomes independent
    assert "TrainStudy" in result["holdout_studies"]
    assert "DevStudy" not in result["holdout_studies"]


def test_stray_training_file_is_unverified_until_a_majority_matches():
    # 1 molecule overlapping a snapshot study out of 10 -> not the training file
    stray = {TRAIN} | {f"{'ACGT'[i % 4]}" * 19 for i in range(9)}
    assert audit(training=stray)["checkpoint_training_audit"] == "unverified"
    # a file that is overwhelmingly a snapshot study verifies
    assert audit(training={TRAIN})["checkpoint_training_audit"] == "verified"


def test_audit_is_deterministic():
    assert audit() == audit()


def test_overlap_is_exact_molecule_not_gene():
    # a study sharing a gene but no molecule must stay independent
    result = audit(training={EVAL2})
    assert result["per_study"]["FreeStudy"]["overlap_checkpoint_training"] == 1
    assert result["per_study"]["EvalStudy"]["overlap_checkpoint_training"] == 0


def test_training_exposure_takes_precedence_over_development():
    result = build_audit(records(), {TRAIN, DEV}, {DEV}, set(), minimum_complete_paired_genes=1)
    assert result["per_study"]["DevStudy"]["classification"] == "training_exposed"


def test_zero_molecule_study_is_not_eligible():
    result = build_audit([{"Authors": "ShortStudy", "Gene": "G1", "Antisense_21mer": "ACG"}],
                         set(), set(), set(), minimum_complete_paired_genes=1)
    assert result["per_study"]["ShortStudy"]["unique_molecules"] == 0
    assert result["holdout_studies"] == []


def test_fully_criterion_can_leave_no_holdout():
    exposed = {EVAL, EVAL2}  # every non-training study is also upstream-eval-exposed
    weights = build_audit(records(), {TRAIN}, {DEV}, exposed, criterion="weight", minimum_complete_paired_genes=1)
    fully = build_audit(records(), {TRAIN}, {DEV}, exposed, criterion="fully", minimum_complete_paired_genes=1)
    assert weights["independent_holdout_eligible"] is True
    assert weights["fully_independent_available"] is False
    assert fully["holdout_studies"] == [] and fully["independent_holdout_eligible"] is False


def test_load_molecules_rejects_missing_column(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text(f"siRNA,label\n{TRAIN},1\n")
    assert load_molecules(path, "siRNA") == {TRAIN}
    with pytest.raises(ValueError):
        load_molecules(path, "guide")


def test_write_outputs_does_not_mutate_cohort_and_creates_dir(tmp_path):
    cohort = audit()
    before = copy.deepcopy(cohort)
    out = tmp_path / "nested" / "dir"
    write_outputs(out, {"snapshot": {"path": "x"}}, cohort, 6)
    assert "per_study" in cohort and cohort == before
    assert (out / "cohort_manifest.json").is_file()
    assert (out / "exposure_audit.json").is_file()
