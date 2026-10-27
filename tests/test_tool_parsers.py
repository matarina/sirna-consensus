"""Adapter output decoding, tested against synthetic files.

These are the parts most likely to rot: every tool's output format is an upstream
implementation detail. Pinning them with fixtures means a format change fails
here rather than silently producing zero candidates.
"""
from pathlib import Path

from sirna_consensus.tools import get_tool, survey

SEQ = "ACGTACGTACGTACGTACGTA"


def test_survey_never_raises():
    rows = survey()
    assert rows
    for name, usable, detail in rows:
        assert isinstance(name, str) and isinstance(usable, bool)
        assert isinstance(detail, str)


def test_unknown_tool_is_rejected():
    import pytest
    from sirna_consensus.tools.base import ToolError

    with pytest.raises(ToolError, match="unknown tool"):
        get_tool("nope")


def test_sirnadesigner_parses_its_csv(tmp_path: Path):
    path = tmp_path / "_TX.csv"
    # start,end,24nt,score,guide,passenger
    path.write_text(f"100,122,{SEQ},7,UCUCCUUCATATTTCCTGACG,{SEQ}\n")
    hits = get_tool("sirnadesigner")._parse(path)
    assert len(hits) == 1
    assert hits[0].start == 102          # passenger starts 2 nt into the 24-mer
    assert hits[0].sense == SEQ
    assert hits[0].score == 7.0


def test_sirnadesigner_skips_short_rows(tmp_path: Path):
    path = tmp_path / "bad.csv"
    path.write_text("1,2,ACGT\n")
    assert get_tool("sirnadesigner")._parse(path) == []


def test_sirnaforge_parses_comma_delimited_table(tmp_path: Path):
    path = tmp_path / "out.csv"
    path.write_text("id,transcript_id,position,guide_sequence,passenger_sequence,"
                    "composite_score,off_target_count,passes_filters\n"
                    f"SIRNAF_1,TX,100,{SEQ},{SEQ},82.5,0,True\n")
    hits = get_tool("sirnaforge")._parse(path)
    assert len(hits) == 1
    assert hits[0].start == 100 and hits[0].score == 82.5


def test_sirnaforge_parses_tab_delimited_table(tmp_path: Path):
    path = tmp_path / "out.tsv"
    path.write_text("position\tpassenger_sequence\tcomposite_score\n"
                    f"250\t{SEQ}\t12.0\n")
    hits = get_tool("sirnaforge")._parse(path)
    assert hits[0].start == 250


def test_sirnaforge_ignores_foreign_table(tmp_path: Path):
    path = tmp_path / "other.csv"
    path.write_text("a,b,c\n1,2,3\n")
    assert get_tool("sirnaforge")._parse(path) == []


def test_sirnaforge_handles_missing_file(tmp_path: Path):
    assert get_tool("sirnaforge")._parse(tmp_path / "absent.csv") == []


def test_siren_parses_sirnas_and_inverts_offtarget_burden(tmp_path: Path):
    other = tmp_path / "other_files"
    other.mkdir()
    (other / "sirnas.fa").write_text(
        f">sirna_1-21\n{SEQ}\n>sirna_1-21_r\nTACGTACGTACGTACGTACGT\n"
        f">sirna_5-25\n{SEQ}\n>sirna_5-25_r\nTACGTACGTACGTACGTACGT\n")
    (tmp_path / "off_targets_summary.tsv").write_text(
        "Off target\tsiRNA number\tsiRNA names\n"
        "ENST1\t2\tsirna_1-21, sirna_5-25_r\n")

    tool = get_tool("siren")
    burden = tool._offtarget_burden(tmp_path)
    assert burden == {"sirna_1-21": 1, "sirna_5-25": 1}

    hits = tool._parse_sirnas(tmp_path, burden)
    assert [h.start for h in hits] == [1, 5]      # the _r records are guides
    assert hits[0].score == 0.5                   # 1 / (1 + burden)


def test_siren_returns_nothing_without_sirnas_file(tmp_path: Path):
    assert get_tool("siren")._parse_sirnas(tmp_path, {}) == []


def test_siren_drops_records_whose_name_has_no_coordinates(tmp_path: Path):
    path = tmp_path / "sirnas.fa"
    path.write_text(f">unnamed\n{SEQ}\n")
    assert get_tool("siren")._parse_sirnas(tmp_path, {}) == []
