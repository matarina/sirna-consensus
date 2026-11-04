"""Compute checkpoint/development molecular exposure before any independent claim.

The blocker was never measured: audit_sources.py *asserted* that all seven database
studies sit in OligoFormer's development/evaluation collections and set
``independent_holdout_eligible=False``.  This module actually separates the two
kinds of exposure, because only one of them invalidates a frozen-checkpoint
benchmark:

* **training exposure** - the molecule was in the checkpoint's fitting data
  (upstream trains on Huesken / ``data/Hu.csv``).  This changes the weights.
* **upstream evaluation exposure** - the molecule was in an upstream *evaluation*
  collection (``Mix.csv`` / ``Taka.csv``).  This did not fit the weights, but may
  have shaped upstream model selection.  Disclosed, not automatically fatal.

A study is ``weight_independent`` when it shares no molecule with the checkpoint
training set and none with this project's development cohort.  With
``--criterion fully`` a study must additionally avoid every upstream evaluation
collection.  Identity is exact 19-nt antisense-core match; nothing is inferred.

Assets are read-only and hashed.  Missing assets degrade to ``unverified`` rather
than guessing.  The module is importable so the logic is unit-testable without the
gitignored raw downloads.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
_XLSX_NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
CORE_LEN = 19


def normalize(value: str) -> str:
    return re.sub(r"[^ACGTU]", "", (value or "").upper().replace("U", "T"))


def core(value: str) -> str:
    """Canonical 19-nt antisense binding core; empty when too short."""
    dna = normalize(value)
    return dna[:CORE_LEN] if len(dna) >= CORE_LEN else ""


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_snapshot(path: Path) -> list[dict]:
    """Return snapshot rows keyed by header name (shared-string xlsx)."""
    with zipfile.ZipFile(path) as archive:
        strings = ["".join(t.itertext()) for t in
                   ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("x:si", _XLSX_NS)]
        sheet = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        raw = []
        for row in sheet.findall("x:sheetData/x:row", _XLSX_NS):
            cells = {}
            for cell in row:
                value = cell.findtext("x:v", "", _XLSX_NS)
                column = "".join(c for c in cell.attrib["r"] if c.isalpha())
                cells[column] = strings[int(value)] if cell.attrib.get("t") == "s" else value
            raw.append(cells)
    header = raw.pop(0)
    return [{header[k]: v for k, v in row.items()} for row in raw]


def load_molecules(path: Path, column: str) -> set[str]:
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or column not in reader.fieldnames:
            raise ValueError(f"{path}: missing column {column!r}")
        return {core(row[column]) for row in reader if core(row[column])}


def load_development(path: Path) -> set[str]:
    with Path(path).open(newline="") as handle:
        return {core(row["guide_dna"]) for row in csv.DictReader(handle, delimiter="\t")
                if row.get("guide_dna")}


def _study_rows(records: list[dict]) -> dict[str, list[dict]]:
    groups: dict[str, list[dict]] = collections.defaultdict(list)
    for row in records:
        groups[row.get("Authors", "unknown")].append(row)
    return groups


def _gene(row: dict) -> str:
    return (row.get("Gene", "") or "").strip()


def build_audit(snapshot_records: list[dict], training: set[str], development: set[str],
                upstream_eval: set[str], *, criterion: str = "weight",
                minimum_complete_paired_genes: int = 5) -> dict:
    """Pure core: classify every study and decide holdout eligibility."""
    if criterion not in ("weight", "fully"):
        raise ValueError("criterion must be 'weight' or 'fully'")
    per_study = {}
    for study, rows in sorted(_study_rows(snapshot_records).items()):
        molecules = {core(r.get("Antisense_21mer", "")) for r in rows} - {""}
        overlap_train = len(molecules & training)
        overlap_dev = len(molecules & development)
        overlap_eval = len(molecules & upstream_eval)
        if overlap_train:
            label = "training_exposed"
        elif overlap_dev:
            label = "development_exposed"
        elif overlap_eval:
            label = "upstream_evaluation_exposed"
        else:
            label = "independent"
        weight_independent = overlap_train == 0 and overlap_dev == 0 and len(molecules) > 0
        fully_independent = weight_independent and overlap_eval == 0
        per_study[study] = {
            "records": len(rows),
            "unique_genes": len({_gene(r) for r in rows} - {""}),
            "unique_molecules": len(molecules),
            "overlap_checkpoint_training": overlap_train,
            "overlap_project_development": overlap_dev,
            "overlap_upstream_evaluation": overlap_eval,
            "classification": label,
            "weight_independent": weight_independent,
            "fully_independent": fully_independent,
        }
    holdout = {s: v for s, v in per_study.items()
               if v["fully_independent" if criterion == "fully" else "weight_independent"]}
    genes = {_gene(r) for s in holdout
             for r in _study_rows(snapshot_records)[s]} - {""}
    records = sum(v["records"] for v in holdout.values())
    fully = {s: v for s, v in per_study.items() if v["fully_independent"]}
    fully_genes = {_gene(r) for s in fully
                   for r in _study_rows(snapshot_records)[s]} - {""}
    exposed_to_training = [s for s, v in per_study.items() if v["overlap_checkpoint_training"]]
    eval_exposed = [s for s, v in per_study.items() if v["overlap_upstream_evaluation"]]
    # "verified" means the located training file actually IS a snapshot study, not a
    # stray non-empty file: a majority of its molecules must fall in one study. A file
    # that merely shares one molecule with the snapshot stays unverified.
    match_fraction = (max((v["overlap_checkpoint_training"] for v in per_study.values()), default=0)
                      / len(training)) if training else 0.0
    audit = "verified" if match_fraction >= 0.5 else "unverified"
    eligible = bool(holdout) and len(genes) >= minimum_complete_paired_genes
    return {
        "criterion": criterion,
        "minimum_complete_paired_genes": minimum_complete_paired_genes,
        "checkpoint_training_audit": audit,
        "training_match_fraction": round(match_fraction, 4),
        "training_split_basis": "checkpoint training set is upstream data/Hu.csv; verified when >=50% of "
                                "its molecules fall in one snapshot study (released split is an upstream assumption)",
        "independent_holdout_eligible": eligible,
        "holdout_studies": sorted(holdout),
        "excluded_studies": {s: v["classification"] for s, v in per_study.items() if s not in holdout},
        "holdout_records": records,
        "holdout_genes": len(genes),
        "holdout_gene_names": sorted(genes),
        "fully_independent_studies": sorted(fully),
        "fully_independent_records": sum(v["records"] for v in fully.values()),
        "fully_independent_genes": len(fully_genes),
        "fully_independent_available": bool(fully),
        "upstream_evaluation_exposed_studies": sorted(eval_exposed),
        "training_exposed_studies": sorted(exposed_to_training),
        "per_study": per_study,
    }


def build_inputs(snapshot: Path, train_source: Path, upstream_eval: list[Path], development: Path,
                 training: set[str], development_set: set[str]) -> dict:
    """The pinned-input record, identical whichever entry point calls it."""
    return {
        "snapshot": {"path": str(snapshot), "sha256": sha256_file(snapshot)},
        "checkpoint_training_source": ({"path": str(train_source), "sha256": sha256_file(train_source)}
                                       if train_source.is_file() else {"path": str(train_source)}),
        "upstream_evaluation_sources": [{"path": str(p), "sha256": sha256_file(p)}
                                        for p in upstream_eval if p.is_file()],
        "development": ({"path": str(development), "sha256": sha256_file(development)}
                        if development_set else {"path": str(development)}),
    }


def _reason(cohort: dict) -> str:
    if cohort["checkpoint_training_audit"] != "verified":
        return ("checkpoint training data (upstream Hu.csv) unavailable; exposure unverified")
    if not cohort["holdout_studies"]:
        if cohort["criterion"] == "fully":
            return ("no study is additionally free of upstream-evaluation exposure; "
                    "only a weight-independent cohort exists")
        return "no study is free of checkpoint-training and development molecular overlap"
    if not cohort["independent_holdout_eligible"]:
        return (f"holdout has {cohort['holdout_genes']} genes, below the "
                f"{cohort['minimum_complete_paired_genes']}-gene floor")
    reason = (f"{len(cohort['holdout_studies'])} studies / {cohort['holdout_records']} records / "
              f"{cohort['holdout_genes']} genes share no molecule with checkpoint training or "
              "the development cohort")
    if cohort["criterion"] == "weight" and not cohort["fully_independent_available"]:
        reason += "; all are upstream-evaluation-exposed, so this is weight-independent, not fully independent"
    return reason


def _residual(cohort: dict) -> list[str]:
    notes = ["Exact 19-nt antisense-core identity only; no coordinates, isoforms or near-duplicates.",
             "Checkpoint training set is taken as upstream data/Hu.csv (the Huesken file its "
             "trainer consumes); the released split itself is an upstream assumption, recorded not assumed away.",
             "Gene names are source strings; synonyms (e.g. 'Firefly luciferase' vs 'FireflyLuc') are not merged."]
    if cohort["criterion"] == "weight":
        notes.append("Studies upstream kept as evaluation collections (Mix/Taka) are retained at weight level. "
                     "Upstream trains on Hu and selects epochs on Mix, so those molecules did not fit the "
                     "weights but did steer early stopping: exposure is stronger than pure 'published eval'.")
    notes.append("Strategy/policy selection is confined to the Huesken development cohort, so the "
                 "holdout is never used to choose a winner here.")
    if cohort["criterion"] == "weight" and not cohort["fully_independent_available"]:
        notes.append("No study is free of upstream-evaluation exposure either; run --criterion fully to see it.")
    return notes


def write_outputs(base: Path, inputs: dict, cohort: dict, records_note: str) -> None:
    # No created_utc: these files are pinned in freeze locks, so a timestamp would make
    # any regeneration (and therefore freeze --verify) fail for no semantic reason.
    base = Path(base)
    base.mkdir(parents=True, exist_ok=True)
    audit = {
        "schema_version": 1,
        "purpose": "separate checkpoint-training exposure from evaluation exposure before an independent claim",
        "inputs": inputs, "snapshot_records": records_note,
        "per_study": cohort["per_study"],
        "summary": {"training_exposed_studies": cohort["training_exposed_studies"],
                    "upstream_evaluation_exposed_studies": cohort["upstream_evaluation_exposed_studies"],
                    "fully_independent_available": cohort["fully_independent_available"],
                    "checkpoint_training_audit": cohort["checkpoint_training_audit"]},
    }
    (base / "exposure_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    manifest = {
        "schema_version": 1,
        "source": "exposure_audit.json", "criterion": cohort["criterion"],
        "checkpoint_training_audit": cohort["checkpoint_training_audit"],
        "training_match_fraction": cohort["training_match_fraction"],
        "training_split_basis": cohort["training_split_basis"],
        "independent_holdout_eligible": cohort["independent_holdout_eligible"],
        "minimum_complete_paired_genes": cohort["minimum_complete_paired_genes"],
        "holdout_studies": cohort["holdout_studies"],
        "excluded_studies": cohort["excluded_studies"],
        "holdout_records": cohort["holdout_records"],
        "holdout_genes": cohort["holdout_genes"],
        "holdout_gene_names": cohort["holdout_gene_names"],
        "fully_independent_studies": cohort["fully_independent_studies"],
        "fully_independent_records": cohort["fully_independent_records"],
        "fully_independent_available": cohort["fully_independent_available"],
        "reason": _reason(cohort),
        "residual_exposure": _residual(cohort),
        "buildable_holdout": {"dataset_manifest": "data/holdout/dataset_manifest.json",
                              "note": "records whose accession or strands do not place exactly on a "
                                      "reference are excluded there, so the built cohort can be smaller "
                                      "than holdout_records (e.g. an accession absent from the snapshot FASTA)"},
        "inputs": inputs,
    }
    (base / "cohort_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--snapshot", type=Path, default=ROOT / "data/raw/sirna_data.xlsx")
    # The sibling project lives beside this repo, not inside it (benchmarks -> repo -> parent).
    sibling = ROOT.parent.parent / "sirnaforge-platform/data/raw/oligoformer/data"
    p.add_argument("--train-source", type=Path, default=sibling / "Hu.csv")
    p.add_argument("--upstream-eval", type=Path, nargs="*",
                   default=[sibling / "Mix.csv", sibling / "Taka.csv"])
    p.add_argument("--development", type=Path, default=ROOT / "data/huesken/measurements.tsv")
    p.add_argument("--out", type=Path, default=ROOT)
    p.add_argument("--criterion", choices=("weight", "fully"), default="weight")
    p.add_argument("--min-genes", type=int, default=5)
    a = p.parse_args(argv)

    records = parse_snapshot(a.snapshot)
    if not a.train_source.is_file():
        print(f"warning: training source missing ({a.train_source}); exposure cannot be verified", file=sys.stderr)
    if not any(path.is_file() for path in a.upstream_eval):
        print("warning: no upstream evaluation collection found; 'fully' criterion is under-informed", file=sys.stderr)
    if not a.development.is_file():
        print(f"warning: development file missing ({a.development}); overlap not assessed",
              file=sys.stderr)
    training = load_molecules(a.train_source, "siRNA") if a.train_source.is_file() else set()
    development = load_development(a.development) if a.development.is_file() else set()
    upstream_eval = {m for path in a.upstream_eval if path.is_file()
                     for m in load_molecules(path, "siRNA")}
    inputs = build_inputs(a.snapshot, a.train_source, list(a.upstream_eval), a.development,
                          training, development)
    cohort = build_audit(records, training, development, upstream_eval,
                         criterion=a.criterion, minimum_complete_paired_genes=a.min_genes)
    write_outputs(a.out, inputs, cohort, len(records))
    print(json.dumps({"holdout_studies": cohort["holdout_studies"],
                      "holdout_genes": cohort["holdout_genes"],
                      "independent_holdout_eligible": cohort["independent_holdout_eligible"],
                      "checkpoint_training_audit": cohort["checkpoint_training_audit"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

# touch
