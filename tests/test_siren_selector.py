"""Selector resolution for SIREN's ``-g``, against a synthetic transcriptome.

Hermetic: a handful of fake cDNA records, no network, no docker, no real 1.4 GB
Ensembl file.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from sirna_consensus.tools.base import ToolEnv, ToolError
from sirna_consensus.tools.siren import ENV_SELECTOR, SirenTool, resolve_selector

# Deliberately distinct sequences so the digest mapping is unambiguous.
SEQ_A = "ACGT" * 40                       # the transcript we feed in as input
SEQ_B = "TTGCAAGCTT" * 12
SEQ_C = ("AAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACC"
         "CTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTATTTGAAGGAGTCTAGCAGCCGCAG"
         "TAAGGCACAATACCTCGTCCGTGTTACCAGACCAAACAAGACGTCCTCTTCAATGTTTAAA"
         "TGACCCTCTCGTCATAAAACCTTTCTACTATGTGTT")


@pytest.fixture
def transcriptome(tmp_path: Path) -> Path:
    path = tmp_path / "cdna.fa"
    path.write_text(
        f">ENST00000000001.1 cdna chromosome:GRCh38:1:1:100:1 "
        f"gene:ENSG00000000001.1 gene_symbol:GENEONE "
        f"description:thing [Source:HGNC Symbol;Acc:HGNC:1]\n{SEQ_A}\n"
        f">ENST00000000002.2 cdna chromosome:GRCh38:1:1:100:1 "
        f"gene:ENSG00000000002.1 gene_symbol:GENETWO\n{SEQ_B}\n"
        f">ENST00000000003.3 cdna chromosome:GRCh38:2:1:100:1 "
        f"gene:ENSG00000000003.1 gene_symbol:GENEONE\n{SEQ_C}\n")
    return path


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("SIRNA_CONSENSUS_CACHE", str(tmp_path / "cache"))
    monkeypatch.setattr("sirna_consensus.tools.siren._INDEX_CACHE", {})


def _input(tmp_path: Path, header: str, sequence: str) -> Path:
    path = tmp_path / "input.fa"
    path.write_text(f">{header}\n{sequence}\n")
    return path


def test_exact_id_wins(tmp_path, transcriptome):
    fasta = _input(tmp_path, "ENST00000000001.1 whatever", SEQ_A)
    res = resolve_selector(fasta, "ENST00000000001.1", transcriptome)
    assert (res.selector, res.source, res.matched) == (
        "ENST00000000001.1", "exact-id", True)


def test_version_suffix_is_dropped_when_unique(tmp_path, transcriptome):
    fasta = _input(tmp_path, "ENST00000000002.7", SEQ_B)
    res = resolve_selector(fasta, "ENST00000000002.7", transcriptome)
    assert (res.selector, res.source) == ("ENST00000000002", "versionless-id")


def test_gene_symbol_from_input_header_resolves_to_the_id(tmp_path, transcriptome):
    # GENETWO names exactly one transcript.
    fasta = _input(tmp_path, "XM_1.1 gene_symbol:GENETWO", SEQ_B)
    res = resolve_selector(fasta, "XM_1.1", transcriptome)
    assert (res.selector, res.source, res.gene) == (
        "ENST00000000002.2", "gene-symbol", "GENETWO")


def test_gene_ambiguity_is_broken_by_the_mapped_transcript(tmp_path, transcriptome):
    # A refseq id with no gene token: GENEONE names two transcripts, so only the
    # transcript our sequence actually matches can decide.
    fasta = _input(tmp_path, "NM_009.4 amphiregulin-like (GENEONE), mRNA", SEQ_C)
    res = resolve_selector(fasta, "NM_009.4", transcriptome)
    assert res.source == "gene-symbol"
    assert res.selector == "ENST00000000003.3"      # the one we actually map to
    assert res.mapped_id == "ENST00000000003.3"


def test_common_token_is_the_last_header_based_resort(tmp_path, transcriptome):
    # A foreign accession with no gene token and no sequence match: the longest
    # token it shares with exactly one transcriptome id is used.
    fasta = _input(tmp_path, "XM_00000000002.1 x", "GGGGCCCC" * 10)
    res = resolve_selector(fasta, "XM_00000000002.1", transcriptome)
    assert (res.selector, res.source, res.matched) == (
        "00000000002", "common-token", True)


def test_a_foreign_id_does_not_resolve(tmp_path, transcriptome):
    fasta = _input(tmp_path, "ZZZ.1", "GGGGCCCC" * 10)
    res = resolve_selector(fasta, "ZZZ.1", transcriptome)
    assert res.matched is False
    assert res.selector == "ZZZ.1"


def test_override_is_used_verbatim(tmp_path, transcriptome):
    fasta = _input(tmp_path, "NM_1.1", SEQ_A)
    res = resolve_selector(fasta, "NM_1.1", transcriptome,
                           override="ENST00000000002")
    assert (res.selector, res.source, res.matched) == (
        "ENST00000000002", "override", True)


def test_override_that_matches_nothing_is_reported_as_unmatched(tmp_path,
                                                                transcriptome):
    fasta = _input(tmp_path, "NM_1.1", SEQ_A)
    res = resolve_selector(fasta, "NM_1.1", transcriptome, override="NOSUCHGENE")
    assert (res.selector, res.matched) == ("NOSUCHGENE", False)


# -- run() wiring: which scope the detail string claims ---------------------

class _Fake(SirenTool):
    """SIREN stand-in: records invocations and writes a fake off-target summary."""

    def __init__(self, env, fail_on_transcriptome=False):
        super().__init__(env)
        self.calls: list[tuple[str, str]] = []
        self.fail_on_transcriptome = fail_on_transcriptome

    def availability(self):
        return True, "fake"

    def _invoke(self, targets, selector, workdir):
        self.calls.append((str(targets), selector))
        if self.fail_on_transcriptome and targets.name == "cdna.fa":
            raise ToolError("fake SIREN: gene not found")
        workdir.mkdir(parents=True, exist_ok=True)
        other = workdir / "other_files"
        other.mkdir(exist_ok=True)
        (other / "sirnas.fa").write_text(f">sirna_1-21\n{SEQ_A[:21]}\n")
        (workdir / "off_targets_summary.tsv").write_text(
            "Off target\tsiRNA number\tsiRNA names\n"
            "ENST9\t1\tsirna_1-21_r\n")


def _env(transcriptome: Path | None, selector: str = "") -> ToolEnv:
    env = ToolEnv()
    env.transcriptome = str(transcriptome) if transcriptome else ""
    env.siren_gene_selector = selector
    return env


def test_detail_says_screened_against_the_transcriptome(tmp_path, transcriptome):
    fasta = _input(tmp_path, "NM_009.4 (GENEONE)", SEQ_C)
    tool = _Fake(_env(transcriptome))
    run = tool.run(fasta, tmp_path / "wd", "NM_009.4")
    assert run.status == "ok"
    assert "screened against the transcriptome" in run.detail
    assert "ENST00000000003.3" in run.detail
    assert tool.calls == [(str(transcriptome), "ENST00000000003.3")]
    assert run.hits[0].extra["offtargets"] == 1


def test_detail_admits_the_last_resort_fallback_when_nothing_resolves(
        tmp_path, transcriptome):
    fasta = _input(tmp_path, "ZZZ.1", "GGGGCCCC" * 10)
    tool = _Fake(_env(transcriptome))
    run = tool.run(fasta, tmp_path / "wd", "ZZZ.1")
    assert run.status == "ok"
    assert "fell back to the input FASTA only" in run.detail
    assert "no transcriptome selector matched" in run.detail
    assert tool.calls == [(str(fasta), "ZZZ.1")]


def test_detail_admits_fallback_when_siren_rejects_a_resolved_selector(
        tmp_path, transcriptome):
    fasta = _input(tmp_path, "NM_009.4 (GENEONE)", SEQ_C)
    tool = _Fake(_env(transcriptome), fail_on_transcriptome=True)
    run = tool.run(fasta, tmp_path / "wd", "NM_009.4")
    assert run.status == "ok"
    assert "fell back to the input FASTA only" in run.detail
    assert "no transcriptome-backed screen happened" in run.detail
    assert [call[0] for call in tool.calls] == [str(transcriptome), str(fasta)]


def test_fallback_parse_never_reads_the_aborted_attempt_summary(
        tmp_path, transcriptome):
    """run3 (2026-10-27) hazard, as a regression test.

    The transcriptome attempt wrote a real (partial) transcriptome-wide
    off_targets_summary.tsv and then died. The input-only retry finished
    without producing its own summary. Parsing must take the *retry's*
    directory only: input-scoped means no off-targets -- never a silent
    adoption of the aborted screen's burden mislabelled as input-scoped.
    """
    class _Leaky(SirenTool):
        def __init__(self, env):
            super().__init__(env)
            self.calls = []

        def availability(self):
            return True, "fake"

        def _invoke(self, targets, selector, workdir):
            self.calls.append((str(targets), selector))
            workdir.mkdir(parents=True, exist_ok=True)
            other = workdir / "other_files"
            other.mkdir(exist_ok=True)
            (other / "sirnas.fa").write_text(f">sirna_1-21\n{SEQ_A[:21]}\n")
            if targets.name == "cdna.fa":
                # aborted transcriptome screen: partial summary, then death
                (workdir / "off_targets_summary.tsv").write_text(
                    "Off target\tsiRNA number\tsiRNA names\n"
                    "ENST9\t3\tsirna_1-21_r\n")
                raise ToolError("fake SIREN: timed out mid screen")
            # the input-only retry completes but screens nothing beyond the
            # input itself: it deliberately writes no summary at all

    fasta = _input(tmp_path, "NM_009.4 (GENEONE)", SEQ_C)
    tool = _Leaky(_env(transcriptome))
    run = tool.run(fasta, tmp_path / "wd", "NM_009.4")
    assert run.status == "ok"
    assert "fell back to the input FASTA only" in run.detail
    # burden comes from the retry's own directory, which has no summary
    assert run.hits[0].extra["offtargets"] == 0
    # the aborted attempt's artifacts are preserved, but isolated
    assert (tmp_path / "wd" / "transcriptome" /
            "off_targets_summary.tsv").is_file()
    assert not (tmp_path / "wd" / "input" /
                "off_targets_summary.tsv").exists()


def test_no_transcriptome_configured_says_so(tmp_path):
    fasta = _input(tmp_path, "NM_009.4", SEQ_A)
    tool = _Fake(_env(None))
    run = tool.run(fasta, tmp_path / "wd", "NM_009.4")
    assert run.status == "ok"
    assert "no transcriptome configured" in run.detail
    assert tool.calls == [(str(fasta), "NM_009.4")]


def test_env_selector_is_read_from_the_environment(monkeypatch):
    monkeypatch.setenv(ENV_SELECTOR, "MYGENE")
    assert ToolEnv().siren_gene_selector == "MYGENE"


def test_coordinates_from_the_mapped_accession_still_place_on_the_input(
        tmp_path, transcriptome):
    """The resolved transcript is a *different* accession; place() must cope."""
    from sirna_consensus.normalize import place
    from sirna_consensus.types import RawHit

    input_id = "NM_009.4"
    fasta = _input(tmp_path, input_id, SEQ_C)
    res = resolve_selector(fasta, input_id, transcriptome)
    assert res.selector != input_id            # SIREN coordinates are on another id

    sense = SEQ_C[100:121]
    # A start SIREN would report on ENST00000000003.3 that disagrees with ours.
    hit = RawHit(tool="siren", sense=sense, length=21, start=7, score=1.0)
    placement = place(hit, input_id, SEQ_C)
    assert placement.ok
    assert placement.core.start == 101         # found by sequence, not by start


def test_siren_threads_default_and_override(monkeypatch):
    from sirna_consensus.tools.base import ToolEnv

    monkeypatch.delenv("SIREN_THREADS", raising=False)
    assert ToolEnv().siren_threads > 0

    monkeypatch.setenv("SIREN_THREADS", "48")
    assert ToolEnv().siren_threads == 48

    monkeypatch.setenv("SIREN_THREADS", "nonsense")
    assert ToolEnv().siren_threads > 0          # falls back, never raises

    monkeypatch.setenv("SIREN_THREADS", "0")
    assert ToolEnv().siren_threads > 0
