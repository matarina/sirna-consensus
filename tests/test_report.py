"""Output writers.

The pipeline test only asserts the files exist; these pin what is *in* them,
which is what a wet-lab reader actually consumes.
"""
import csv
import json
import re
from pathlib import Path

import pytest

from sirna_consensus import tools as tools_module
from sirna_consensus.pipeline import Config, run
from sirna_consensus.report import TSV_COLUMNS
from sirna_consensus.tools.base import Tool
from sirna_consensus.types import RawHit, reverse_complement

COMP = str.maketrans("ACGT", "TGCA")


class FakeTool(Tool):
    def __init__(self, name, sites, env=None, available=True):
        super().__init__(env)
        self.name = name
        self.sites = sites
        self._available = available

    def availability(self):
        return (self._available, "fake")

    def version(self):
        return "fake-1"

    def run(self, fasta, workdir, transcript_id):
        lines = Path(fasta).read_text().splitlines()[1:]
        sequence = "".join(line.strip() for line in lines)
        hits = [RawHit(tool=self.name, sense=sequence[s - 1:s + 20], length=21,
                       start=s, score=float(s)) for s in self.sites]
        return self.wrap("ok", "fake", hits)


@pytest.fixture
def run_dir(tmp_path, fasta_file, monkeypatch):
    monkeypatch.setattr(tools_module, "REGISTRY", {
        "alpha": lambda env=None: FakeTool("alpha", [1, 42, 70], env),
        "beta": lambda env=None: FakeTool("beta", [1, 70], env),
    })
    out = tmp_path / "out"
    run(Config(fasta=str(fasta_file), out_dir=str(out),
               tools=("alpha", "beta"), min_votes=2, include_singletons=True))
    return out


def rows(run_dir: Path):
    return list(csv.DictReader((run_dir / "candidates.tsv").open(), delimiter="\t"))


def test_tsv_has_the_advertised_columns(run_dir):
    assert list(rows(run_dir)[0]) == TSV_COLUMNS


def test_tsv_sense_and_guide_are_complementary(run_dir):
    for row in rows(run_dir):
        assert row["guide_sequence_dna"] == reverse_complement(row["target_sequence_dna"])
        assert row["guide_sequence_rna"] == row["guide_sequence_dna"].replace("T", "U")


def test_tsv_is_sorted_by_descending_score(run_dir):
    scores = [float(r["score"]) for r in rows(run_dir)]
    assert scores == sorted(scores, reverse=True)


def test_votes_dominate_score(run_dir):
    """A 2-vote site must outrank a 1-vote site however good the latter looks."""
    by_votes: dict[int, list[float]] = {}
    for row in rows(run_dir):
        by_votes.setdefault(int(row["votes"]), []).append(float(row["score"]))
    assert max(by_votes[1]) < min(by_votes[2])


def test_fasta_guides_match_the_table(run_dir):
    text = (run_dir / "top_guides.fa").read_text().splitlines()
    headers = [line for line in text if line.startswith(">")]
    sequences = [line for line in text if not line.startswith(">")]
    table = {(r["target_start"], r["target_end"]): r for r in rows(run_dir)}
    assert headers and len(headers) == len(sequences)
    for header, guide in zip(headers, sequences):
        site = re.search(r"site=\S+?:(\d+)-(\d+)", header)
        assert site, header
        key = (site.group(1), site.group(2))
        assert key in table
        assert table[key]["guide_sequence_rna"] == guide
        assert not set(guide) - set("ACGU")


def test_html_is_escaped_and_complete(run_dir):
    html = (run_dir / "report.html").read_text()
    for section in ("How to read this", "External tools", "Limitations"):
        assert section in html
    assert "<svg" not in html          # this package ships no figures
    assert "&lt;" not in html.split("<style>")[0]
    assert not re.search(r">\s*None\s*<", html)


def test_html_separates_vote_threshold_from_filters(run_dir):
    record = json.loads((run_dir / "run.json").read_text())
    counts = record["counts"]
    assert counts["sites"] == counts["candidates"] + counts["below_min_votes"]
    assert counts["candidates"] == counts["passing"] + counts["rejected_by_filters"]
    html = (run_dir / "report.html").read_text()
    assert f"{counts['rejected_by_filters']} dropped by sequence filters" in html
    assert f"{counts['below_min_votes']} below the vote threshold" in html


def test_run_record_is_self_consistent(run_dir):
    record = json.loads((run_dir / "run.json").read_text())
    assert len(record["candidates"]) == record["counts"]["candidates"]
    assert record["transcript_length"] > 0
    assert record["input_sha256"]
    for entry in record["tools"]:
        assert entry["status"] in {"ok", "unavailable", "failed", "empty"}
