import json
from pathlib import Path

import pytest

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
    assert counts["sites"] == counts["candidates"] + counts["below_min_votes"]
    assert counts["candidates"] == counts["passing"] + counts["rejected_by_filters"]
    assert result["counts"]["sites"] >= 2
    out = tmp_path / "out"
    for name in ("candidates.tsv", "top_guides.fa", "run.json", "report.html"):
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
