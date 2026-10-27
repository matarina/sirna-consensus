"""CLI surface. Hermetic: a fake tool registry stands in for the real tools."""
import csv
import json
from pathlib import Path

import pytest

from sirna_consensus import tools as tools_module
from sirna_consensus.cli import main
from sirna_consensus.tools.base import Tool
from sirna_consensus.types import RawHit


class FakeTool(Tool):
    def __init__(self, name, sites, env=None, length=21):
        super().__init__(env)
        self.name = name
        self.sites = sites
        self.length = length

    def availability(self):
        return True, "fake"

    def version(self):
        return "fake-1"

    def run(self, fasta, workdir, transcript_id):
        lines = Path(fasta).read_text().splitlines()[1:]
        sequence = "".join(line.strip() for line in lines)
        n = self.length
        hits = [RawHit(tool=self.name, sense=sequence[s - 1:s - 1 + n], length=n,
                       start=s, score=1.0) for s in self.sites]
        return self.wrap("ok", "fake", hits)


@pytest.fixture
def registry(monkeypatch):
    monkeypatch.setattr(tools_module, "REGISTRY", {
        "alpha": lambda env=None: FakeTool("alpha", [1, 42], env),
        "beta": lambda env=None: FakeTool("beta", [1, 42], env),
    })


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert capsys.readouterr().out.strip()


def test_tools_command_lists_the_registry(registry, capsys):
    assert main(["tools"]) == 0
    out = capsys.readouterr().out
    assert "alpha" in out and "beta" in out


def test_tools_command_fails_when_nothing_is_usable(monkeypatch, capsys):
    monkeypatch.setattr(tools_module, "REGISTRY", {
        "dead": lambda env=None: FakeTool("dead", [])})

    class Dead(FakeTool):
        def availability(self):
            return False, "nope"

    monkeypatch.setattr(tools_module, "REGISTRY",
                        {"dead": lambda env=None: Dead("dead", [])})
    assert main(["tools"]) == 1
    assert "no tool is usable" in capsys.readouterr().err


def test_run_writes_outputs_and_reports_counts(registry, tmp_path, fasta_file, capsys):
    out = tmp_path / "out"
    code = main(["run", "--fasta", str(fasta_file), "--out", str(out),
                 "--tools", "alpha,beta", "--top-k", "2"])
    assert code == 0
    printed = capsys.readouterr().out
    assert "placed" in printed
    assert "below vote threshold" in printed
    for name in ("candidates.tsv", "top_guides.fa", "run.json", "report.html"):
        assert (out / name).is_file()
    counts = json.loads((out / "run.json").read_text())["counts"]
    assert counts["top"] <= 2
    rows = list(csv.DictReader((out / "candidates.tsv").open(), delimiter="\t"))
    assert counts["candidates"] == len(rows)


def test_exclude_singletons_keeps_only_agreed_sites(registry, tmp_path, fasta_file):
    out = tmp_path / "out"
    main(["run", "--fasta", str(fasta_file), "--out", str(out),
          "--tools", "alpha,beta", "--exclude-singletons", "--min-votes", "2"])
    rows = list(csv.DictReader((out / "candidates.tsv").open(), delimiter="\t"))
    assert rows and all(int(r["votes"]) >= 2 for r in rows)


def test_lengths_option_is_honoured(monkeypatch, tmp_path, fasta_file):
    monkeypatch.setattr(tools_module, "REGISTRY", {
        "short": lambda env=None: FakeTool("short", [1, 42], env, length=19)})
    out = tmp_path / "out"
    main(["run", "--fasta", str(fasta_file), "--out", str(out),
          "--tools", "short", "--lengths", "19"])
    rows = list(csv.DictReader((out / "candidates.tsv").open(), delimiter="\t"))
    assert rows and {r["core_length"] for r in rows} == {"19"}


def test_fuzzy_site_mode_is_reachable(registry, tmp_path, fasta_file):
    out = tmp_path / "out"
    main(["run", "--fasta", str(fasta_file), "--out", str(out),
          "--tools", "alpha,beta", "--site-mode", "fuzzy"])
    record = json.loads((out / "run.json").read_text())
    assert record["config"]["offset_tolerance"] >= 0
    assert (out / "candidates.tsv").is_file()


def test_missing_fasta_argument_exits():
    with pytest.raises(SystemExit):
        main(["run", "--out", "x"])


def test_unknown_tool_is_a_runtime_error(registry, tmp_path, fasta_file):
    from sirna_consensus.tools.base import ToolError

    with pytest.raises(ToolError, match="unknown tool"):
        main(["run", "--fasta", str(fasta_file), "--out", str(tmp_path / "o"),
              "--tools", "ghost"])
