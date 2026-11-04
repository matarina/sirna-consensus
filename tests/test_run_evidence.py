"""--run-evidence attaches native votes onto supplied candidates by 19-nt site."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from sirna_consensus import select_cli  # noqa: E402
from sirna_consensus.fasta import sha256, sha256_file  # noqa: E402
from sirna_consensus.types import reverse_complement  # noqa: E402

TX = "ACGT" * 20
SENSE = TX[10:29]
GUIDE = reverse_complement(SENSE) + "AT"
TID = "ACC1.1"


def _prepare(tmp: Path):
    fasta = tmp / "t.fa"
    sequence = TX
    fasta.write_text(f">{TID}\n{sequence}\n")
    candidates = tmp / "candidates.tsv"
    candidates.write_text("\t".join(["transcript_id", "start", "length", "sense", "ordered_guide",
                                     "guide_overhang", "passenger", "passenger_overhang"]) + "\n"
                          + "\t".join([TID, str(sequence.find(SENSE) + 1), "19", SENSE, GUIDE, "AT", "", ""]) + "\n")
    run = tmp / "run.json"
    run.write_text(json.dumps({
        "schema_version": 2, "transcript_id": TID, "input_sha256": sha256(sequence),
        "normalized_hits": [{"tool": "sirnadesigner", "role": "generator", "vote_eligible": True,
                             "score": 0.9,
                             "core": {"transcript_id": TID, "start": sequence.find(SENSE) + 1, "length": 19,
                                      "sense": SENSE, "ordered_guide": GUIDE, "guide_overhang": "AT",
                                      "passenger": "", "passenger_overhang": ""}}]}))
    return fasta, candidates, run


def _args(**kw):
    base = dict(fasta=None, candidates=None, run=None, enumerate=False, run_evidence=None, out=None,
                strategy="native", filter_policy="soft", disable_components="", top_k=5,
                scorer_manifest=None, policy=None, screen_manifest=None, target_gene=None, bowtie="bowtie")
    base.update(kw)
    return argparse.Namespace(**base)


def test_run_evidence_attaches_and_reports(tmp_path):
    fasta, candidates, run = _prepare(tmp_path)
    out = tmp_path / "out"
    select_cli.command(_args(fasta=str(fasta), candidates=str(candidates), run_evidence=str(run), out=str(out)))
    record = json.loads((out / "selection.json").read_text())
    assert record["provenance"]["evidence_attached"] == 1
    assert record["candidates"][0]["tools"] == ["sirnadesigner"]


def test_run_evidence_is_deterministic(tmp_path):
    fasta, candidates, run = _prepare(tmp_path)
    out = tmp_path / "out"
    select_cli.command(_args(fasta=str(fasta), candidates=str(candidates), run_evidence=str(run), out=str(out)))
    first = (out / "selection.json").read_text()
    select_cli.command(_args(fasta=str(fasta), candidates=str(candidates), run_evidence=str(run), out=str(out)))
    assert (out / "selection.json").read_text() == first
    assert sha256_file(out / "selection.json")
