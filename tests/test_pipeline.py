import json
from pathlib import Path

import pytest

from sirna_consensus import filters
from sirna_consensus import tools as tools_module
from sirna_consensus.pipeline import Config, run
from sirna_consensus.tools.base import Tool
from sirna_consensus.types import RawHit


class FakeTool(Tool):
    """Emits fixed sites; no external process involved."""

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
        hits = [RawHit(tool=self.name, sense=sequence[s - 1: s - 1 + 21], length=21,
                       start=s, score=0.5) for s in self.sites]
        return self.wrap("ok", "fake", hits)


@pytest.fixture
def registered(monkeypatch):
    def install(mapping):
        monkeypatch.setattr(tools_module, "REGISTRY", mapping)
    return install


def make_config(tmp_path, fasta_file, tools, **kw):
    return Config(fasta=str(fasta_file), out_dir=str(tmp_path / "out"),
                  tools=tuple(tools), **kw)


def test_pipeline_runs_and_writes_outputs(tmp_path, fasta_file, registered):
    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1, 42, 70], env),
        "beta": lambda env=None: FakeTool("beta", [1, 42], env),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "beta"]))
    assert result["counts"]["consensus_unanimous"] == 2
    counts = result["counts"]
    assert counts["sites"] == counts["consensus_candidates"] + counts["below_min_votes"]
    assert counts["candidates"] == (
        counts["consensus_candidates"] + counts["exploratory_singletons"])
    assert counts["candidates"] == counts["passing"] + counts["rejected_by_filters"]
    assert result["counts"]["sites"] >= 2
    out = tmp_path / "out"
    for name in ("candidates.tsv", "top_guides.fa", "all_candidates.fa",
                 "run.json", "report.html"):
        assert (out / name).is_file(), name
    record = json.loads((out / "run.json").read_text())
    assert record["input_sha256"] == result["input_sha256"]


def test_singletons_can_be_excluded(tmp_path, fasta_file, registered):
    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1, 42], env),
        "beta": lambda env=None: FakeTool("beta", [1], env),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "beta"], min_votes=2,
                             include_singletons=False))
    assert {c.core.start for c in result["candidates"]} == {1}


def test_unavailable_tool_is_recorded_not_fatal(tmp_path, fasta_file, registered):
    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1, 42], env),
        "broken": lambda env=None: FakeTool("broken", [], env, available=False),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "broken"]))
    statuses = {t["tool"]: t["status"] for t in result["tools"]}
    assert statuses["broken"] == "unavailable"
    assert result["config"]["acting_tools"] == ["alpha"]


def test_failing_tool_does_not_kill_the_run(tmp_path, fasta_file, registered):
    class Boom(FakeTool):
        def run(self, fasta, workdir, transcript_id):
            raise RuntimeError("kaboom")

    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1, 42], env),
        "boom": lambda env=None: Boom("boom", [], env),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "boom"]))
    statuses = {t["tool"]: t["status"] for t in result["tools"]}
    assert statuses["boom"] == "failed"
    assert result["counts"]["passing"] >= 1


def test_tool_returning_nothing_reports_empty(tmp_path, fasta_file, registered):
    registered({"alpha": lambda env=None: FakeTool("alpha", [], env)})
    with pytest.raises(RuntimeError, match="no candidate could be placed"):
        run(make_config(tmp_path, fasta_file, ["alpha"], min_votes=1))


def test_no_usable_tool_raises(tmp_path, fasta_file, registered):
    registered({"nope": lambda env=None: FakeTool("nope", [], env, available=False)})
    with pytest.raises(RuntimeError, match="no external tool is usable"):
        run(make_config(tmp_path, fasta_file, ["nope"]))


def test_normalisation_stats_are_recorded(tmp_path, fasta_file, registered):
    registered({"alpha": lambda env=None: FakeTool("alpha", [1, 42], env)})
    result = run(make_config(tmp_path, fasta_file, ["alpha"], min_votes=1))
    stats = result["normalization"]["alpha"]
    assert stats["placed"] == 2 and stats["total"] == 2


def test_reusing_output_dir_uses_unique_work_dirs(tmp_path, fasta_file, registered):
    registered({"alpha": lambda env=None: FakeTool("alpha", [1, 42], env)})
    out = tmp_path / "out"
    first = run(Config(fasta=str(fasta_file), out_dir=str(out), tools=("alpha",)))
    second = run(Config(fasta=str(fasta_file), out_dir=str(out), tools=("alpha",)))
    assert first["run_id"] != second["run_id"]
    work = out / "work"
    assert work.is_dir()
    assert len([p for p in work.iterdir() if p.is_dir()]) == 2


def test_reusing_output_dir_with_different_transcripts(tmp_path, fasta_file,
                                                      registered):
    registered({"alpha": lambda env=None: FakeTool("alpha", [1, 42], env)})
    out = tmp_path / "out"
    other = tmp_path / "other.fa"
    other.write_text(">OTHER_0002 synthetic\n" + "ACGT" * 20 + "\n")
    first = run(Config(fasta=str(fasta_file), out_dir=str(out), tools=("alpha",)))
    second = run(Config(fasta=str(other), out_dir=str(out), tools=("alpha",)))
    assert first["transcript_id"] != second["transcript_id"]
    assert first["run_id"] != second["run_id"]
    assert len([p for p in (out / "work").iterdir() if p.is_dir()]) == 2
    record = json.loads((out / "run.json").read_text())
    assert record["transcript_id"] == second["transcript_id"]


def test_scores_survive_coordinate_relocation(tmp_path, fasta_file, registered):
    """A hit whose reported start is wrong must keep its score after placement."""
    class Shifted(FakeTool):
        def run(self, fasta, workdir, transcript_id):
            lines = Path(fasta).read_text().splitlines()[1:]
            sequence = "".join(line.strip() for line in lines)
            hits = [
                RawHit(tool=self.name, sense=sequence[0:21], length=21,
                       start=1, score=1.0),
                RawHit(tool=self.name, sense=sequence[41:62], length=21,
                       start=99, score=5.0),
            ]
            return self.wrap("ok", "shifted", hits)

    registered({"alpha": lambda env=None: Shifted("alpha", [], env)})
    result = run(make_config(tmp_path, fasta_file, ["alpha"], min_votes=1))
    by_start = {c.core.start: c for c in result["candidates"]}
    assert 42 in by_start
    assert by_start[42].tool_scores["alpha"] == 1.0
    assert by_start[1].tool_scores["alpha"] == 0.0


def test_min_votes_excludes_singletons_from_top(tmp_path, fasta_file, registered):
    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1, 42], env),
        "beta": lambda env=None: FakeTool("beta", [1], env),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "beta"],
                             min_votes=2, include_singletons=True))
    candidate_starts = {c.core.start for c in result["candidates"]}
    top_starts = {c.core.start for c in result["top"]}
    assert 42 in candidate_starts       # exploratory singleton
    assert top_starts == {1}
    assert all(c.votes >= 2 for c in result["top"])


def test_consensus_tier_uses_acting_tool_denominator(tmp_path, fasta_file, registered):
    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1], env),
        "beta": lambda env=None: FakeTool("beta", [42], env),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "beta"], min_votes=1))
    candidate = next(c for c in result["candidates"] if c.core.start == 1)
    assert candidate.acting_tools == 2
    assert candidate.consensus_tier == "1of_2"


def test_run_record_contains_effective_provenance(tmp_path, fasta_file, registered):
    registered({"alpha": lambda env=None: FakeTool("alpha", [1, 42], env)})
    result = run(make_config(
        tmp_path, fasta_file, ["alpha"],
        site_mode="bounded", include_singletons=False,
        strict_specificity=True, siren_max_offtargets=1,
        filter_config=filters.FilterConfig(gc_min=0.2, gc_max=0.6)))
    config = result["config"]
    for key in ("site_mode", "include_singletons", "strict_specificity",
                "siren_max_offtargets", "filter_config", "transcriptome",
                "transcriptome_sha256"):
        assert key in config
    assert config["site_mode"] == "bounded"
    assert config["include_singletons"] is False
    assert config["filter_config"]["gc_min"] == 0.2
    record = json.loads((tmp_path / "out" / "run.json").read_text())
    assert record["config"]["site_mode"] == "bounded"
    assert record["tools"][0]["role"] == "generator"
    assert "selection_limits" in record["tools"][0]


def test_evaluator_only_site_is_not_a_design_vote(tmp_path, fasta_file, registered):
    class Evaluator(FakeTool):
        role = "evaluator"

        def run(self, fasta, workdir, transcript_id):
            lines = Path(fasta).read_text().splitlines()[1:]
            sequence = "".join(line.strip() for line in lines)
            hits = [
                RawHit(tool=self.name, sense=sequence[0:21], length=21,
                       start=1, score=0.5, vote_eligible=False,
                       extra={"offtargets": 2, "screen_state": "input_only",
                              "transcriptome_backed": False}),
            ]
            return self.wrap("ok", "evaluator-only", hits)

    registered({"eye": lambda env=None: Evaluator("eye", [], env)})
    result = run(make_config(tmp_path, fasta_file, ["eye"], min_votes=1))
    assert result["config"]["acting_tools"] == []
    assert result["candidates"] == []
    assert result["status"] == "screening_only"
    assert result["counts"]["evaluator_hits"] == 1


def test_strict_specificity_rejects_unscreened_candidates(tmp_path, fasta_file, registered):
    class Evaluator(FakeTool):
        role = "evaluator"

        def run(self, fasta, workdir, transcript_id):
            lines = Path(fasta).read_text().splitlines()[1:]
            sequence = "".join(line.strip() for line in lines)
            hits = [
                RawHit(tool=self.name, sense=sequence[0:21], length=21,
                       start=1, score=0.5, vote_eligible=False,
                       extra={"offtargets": 0, "screen_state": "input_only",
                              "transcriptome_backed": False}),
            ]
            return self.wrap("ok", "evaluator-only", hits)

    registered({
        "alpha": lambda env=None: FakeTool("alpha", [1], env),
        "eye": lambda env=None: Evaluator("eye", [], env),
    })
    result = run(make_config(tmp_path, fasta_file, ["alpha", "eye"],
                             min_votes=1, strict_specificity=True))
    candidate = result["candidates"][0]
    assert candidate.specificity_status == "input_only"
    assert "specificity_unverified" in candidate.filter_flags
    assert result["top"] == []
