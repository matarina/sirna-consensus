"""Holdout preparation: placement, exclusion and cohort gating."""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "benchmarks/scripts"))
import prepare_holdout  # noqa: E402

from sirna_consensus.types import reverse_complement  # noqa: E402

TX = "ACGT" * 20
SENSE = TX[10:29]
GUIDE = reverse_complement(SENSE) + "AT"


def _row(**kw):
    row = {"ID": "MC1", "Authors": "StudyA", "Accession_number": "ACC1.1", "Gene": "GENE1",
           "Sense_19mer": SENSE, "Antisense_21mer": GUIDE, "%Inhibition": "80", "Cell": "H",
           "Concentration": "50nM", "Hours": "48h", "PMID": "1"}
    row.update(kw)
    return row


def test_read_fasta_sanitizes_ambiguous_bases(tmp_path):
    path = tmp_path / "r.fa"
    path.write_text(">ACC1.1 desc\nACGTNNNacgu\n")
    assert prepare_holdout.read_fasta(path) == {"ACC1.1": "ACGTAAAACGT"}


def test_safe_name():
    assert prepare_holdout.safe("NM_002046.7") == "NM_002046.7"


def test_build_places_and_excludes(tmp_path, monkeypatch):
    rows = [
        _row(),
        _row(ID="MC2", Accession_number="MISSING"),           # no reference sequence
        _row(ID="MC3", Authors="Huesken"),                    # not a holdout study
        _row(ID="MC4", Sense_19mer="A" * 19, Antisense_21mer="T" * 19 + "AT"),  # sense not in transcript
    ]
    monkeypatch.setattr(prepare_holdout, "parse_snapshot", lambda _p: rows)
    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">ACC1.1\n{TX}\n")
    cohort = tmp_path / "cohort.json"
    cohort.write_text(json.dumps({"holdout_studies": ["StudyA"]}))

    records, candidates, transcripts, exclusions, _ = prepare_holdout.build(tmp_path / "x.xlsx", fasta, cohort)
    assert len(records) == 1 and records[0]["endpoint"] == "percent_inhibition"
    assert records[0]["transcript_id"] == "ACC1.1" and records[0]["value"] == 80.0
    assert len(candidates) == 1 and candidates[0]["start"] == TX.find(SENSE) + 1
    assert set(transcripts) == {"ACC1.1"}
    reasons = {e["record_id"]: e["reason"] for e in exclusions}
    assert reasons == {"snapshot:MC2": "accession 'MISSING' not in reference",
                       "snapshot:MC4": "sense not found in accession"}
