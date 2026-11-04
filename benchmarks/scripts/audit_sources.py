"""Audit the publisher-linked snapshot and its checkpoint/development exposure.

Snapshot provenance (record counts, publisher MD5s, licenses) is recorded here as
before.  The exposure verdict is no longer asserted: it is computed by
``audit_exposure`` from the local assets and written to ``dataset_manifest.json``,
``exposure_audit.json`` and ``cohort_manifest.json``.
"""
import argparse
import collections
import hashlib
import json
import sys
from pathlib import Path

from audit_exposure import (
    build_audit,
    build_inputs,
    load_development,
    load_molecules,
    parse_snapshot,
    write_outputs,
)

ROOT = Path(__file__).resolve().parents[1]
_raw = ROOT / "data/raw"
p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--snapshot", type=Path, default=_raw / "sirna_data.xlsx")
p.add_argument("--train-source", type=Path,
               default=ROOT.parent.parent / "sirnaforge-platform/data/raw/oligoformer/data/Hu.csv")
p.add_argument("--upstream-eval", type=Path, nargs="*", default=[
    ROOT.parent.parent / "sirnaforge-platform/data/raw/oligoformer/data/Mix.csv",
    ROOT.parent.parent / "sirnaforge-platform/data/raw/oligoformer/data/Taka.csv"])
p.add_argument("--development", type=Path, default=ROOT / "data/huesken/measurements.tsv")
p.add_argument("--criterion", choices=("weight", "fully"), default="weight")
a = p.parse_args()
raw = a.snapshot.parent

meta = json.loads((raw / "figshare_25908841_v1.json").read_text())
files = {"46578241": "target_gene_sequences.fa", "46578250": "sirna_data.xlsx"}
artifacts = {}
for f in meta["files"]:
    p = raw / files[str(f["id"])]
    data = p.read_bytes()
    if hashlib.md5(data).hexdigest() != f["computed_md5"]:
        raise ValueError("publisher MD5 mismatch")
    artifacts[p.name] = {"sha256": hashlib.sha256(data).hexdigest(), "publisher_md5": f["computed_md5"],
                         "source_url": f["download_url"], "bytes": len(data)}

records = parse_snapshot(a.snapshot)
studies = collections.Counter(r["Authors"] for r in records)

if not a.train_source.is_file():
    print(f"warning: training source missing ({a.train_source}); exposure unverified", file=sys.stderr)
training = load_molecules(a.train_source, "siRNA") if a.train_source.is_file() else set()
development = load_development(a.development) if a.development.is_file() else set()
upstream_eval = {m for path in a.upstream_eval if path.is_file()
                 for m in load_molecules(path, "siRNA")}
cohort = build_audit(records, training, development, upstream_eval, criterion=a.criterion)

inputs = build_inputs(a.snapshot, a.train_source, list(a.upstream_eval), a.development,
                      training, development)
write_outputs(ROOT, inputs, cohort, len(records))

report = {
    "figshare_article": 25908841,
    "version": meta["version"], "snapshot_records": len(records), "studies": dict(studies),
    "artifacts": artifacts, "figshare_license": meta["license"],
    "paper_license": "CC BY-NC 2.0 UK; provenance discrepancy retained",
    "redistribution_policy": "raw snapshot local only pending license reconciliation",
    "independent_holdout_eligible": cohort["independent_holdout_eligible"],
    "checkpoint_training_audit": cohort["checkpoint_training_audit"],
    "external_holdout_records": cohort["holdout_records"],
    "holdout_studies": cohort["holdout_studies"], "holdout_genes": cohort["holdout_genes"],
    "fully_independent_available": cohort["fully_independent_available"],
    "criteria": {"criterion": cohort["criterion"],
                 "minimum_complete_paired_genes": cohort["minimum_complete_paired_genes"],
                 "identity": "exact 19-nt antisense core"},
    "reasons": [
        "Exposure is computed from local assets by scripts/audit_exposure.py, not asserted.",
        "Checkpoint training data is upstream data/Hu.csv (Huesken); the released split is an upstream assumption, recorded.",
        "Studies upstream kept as evaluation collections (Mix/Taka) are weight-independent but may have shaped upstream model selection.",
        "Database Huesken metadata (HeLa/luciferase) conflicts with original-source H1299/fluorescent-reporter curation.",
        "Database percent-inhibition field rescales Huesken normalized activity; do not interpret it as percent knockdown.",
    ],
    "details": "exposure_audit.json", "cohort": "cohort_manifest.json",
    "independent_validation_passed": False,
}
(ROOT / "dataset_manifest.json").write_text(json.dumps(report, indent=2) + "\n")
print(json.dumps({"snapshot_records": len(records), "studies": dict(studies),
                  "independent_holdout_eligible": cohort["independent_holdout_eligible"],
                  "checkpoint_training_audit": cohort["checkpoint_training_audit"]}))
