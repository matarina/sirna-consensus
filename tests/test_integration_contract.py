import json
from pathlib import Path

from sirna_consensus.normalize import place
from sirna_consensus.pipeline import Config, run
from sirna_consensus.types import RawHit


def test_complete_duplex_without_chemistry_is_rejected():
    sense = "ACGT" * 5 + "A"
    hit = RawHit("tool", sense, 21, start=1,
                 extra={"ordered_guide": "T" * 21, "passenger": sense})
    placement = place(hit, "TX", sense)
    assert not placement.ok and placement.reason == "missing_chemistry"


def test_coordinate_and_strand_contracts_are_rejected():
    sense = "ACGT" * 5 + "A"
    malformed = RawHit("tool", sense, 21, start=1.5)
    assert place(malformed, "TX", sense).reason == "malformed_coordinate"
    inconsistent = RawHit(
        "tool", sense, 21, start=1, chemistry="unmodified_rna_21nt",
        extra={"ordered_guide": "A" * 21, "passenger": sense})
    assert not place(inconsistent, "TX", sense).ok


def test_unmodified_source_chemistry_is_explicitly_supported():
    sense = "ACGT" * 5 + "A"
    hit = RawHit("tool", sense, 21, start=1, chemistry="unmodified_rna_21nt",
                 extra={"ordered_guide": "".join({"A": "T", "C": "G", "G": "C", "T": "A"}[b] for b in sense[::-1]),
                        "passenger": sense, "source_guide_modifications": "none"})
    assert place(hit, "TX", sense).ok


def test_run_records_native_sets_and_validated_union(tmp_path, fasta_file, monkeypatch):
    from sirna_consensus.types import ToolRun

    class Emitting:
        name = "emitting"
        role = "generator"
        acceptance = "fixture"

        def __init__(self, env=None):
            self.env = env

        def availability(self):
            return True, "fixture"

        def version(self):
            return "fixture-1"

        def hash_outputs(self, workdir):
            return {p.name: "fixture-hash" for p in workdir.iterdir() if p.is_file()}

        def run(self, fasta, workdir, transcript_id):
            (workdir / "native.csv").write_text("native\n")
            sequence = "".join(Path(fasta).read_text().splitlines()[1:])
            hit = RawHit("emitting", sequence[:21], 21, start=1, score=0.7)
            return ToolRun("emitting", "ok", [hit], version=self.version(), role=self.role)

    from sirna_consensus import tools as tools_module
    monkeypatch.setattr(tools_module, "REGISTRY", {"emitting": Emitting})
    result = run(Config(str(fasta_file), str(tmp_path / "out"), tools=("emitting",)))
    assert result["native_candidates"]["emitting"][0]["native_score"] == 0.7
    assert result["candidate_union"][0]["source_tools"] == ["emitting"]
    assert result["format_compatibility"]["status"] == "compatible"
    assert result["site_agreement"]["multi_tool_sites"] == 0
    saved = json.loads((tmp_path / "out" / "run.json").read_text())
    assert saved["tools"][0]["version"] == "fixture-1"
    assert "native.csv" in saved["tools"][0]["raw_output_hashes"]
    from sirna_consensus.selection import load_candidates
    rows, _ = load_candidates(tmp_path / "out" / "run.json")
    assert rows[0]["core"].chemistry == "unknown"


def test_legacy_adapter_hit_list_is_wrapped(tmp_path, fasta_file, monkeypatch):
    class Legacy:
        name = "legacy"
        role = "generator"
        acceptance = "legacy fixture"

        def __init__(self, env=None):
            self.env = env

        def availability(self):
            return True, "fixture"

        def version(self):
            return "legacy-1"

        def run(self, fasta, workdir, transcript_id):
            sequence = "".join(Path(fasta).read_text().splitlines()[1:])
            return [RawHit("legacy", sequence[:21], 21, start=1)]

    from sirna_consensus import tools as tools_module
    monkeypatch.setattr(tools_module, "REGISTRY", {"legacy": Legacy})
    result = run(Config(str(fasta_file), str(tmp_path / "out"), tools=("legacy",)))
    assert result["tools"][0]["detail"].startswith("legacy adapter")
