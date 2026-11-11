"""Behavioral contracts for tool-driven molecular identity and provenance."""

import pytest

from sirna_consensus.cli import main
from sirna_consensus.pipeline import Config, run
from sirna_consensus.reference import utr3
from sirna_consensus.selection import attach_evidence, enumerate_duplexes
from sirna_consensus.types import Core, reverse_complement


def test_duplex_geometry_edges_and_chemistry(random_dna):
    tx = random_dna(90)
    rows = enumerate_duplexes("T", tx)
    assert len(rows) == len(tx) - 22
    assert rows[0]["core"].start == 3 and rows[-1]["core"].end == 88
    for row in rows:
        c = row["core"]
        s = c.start - 1
        assert c.guide == reverse_complement(tx[s - 2 : s + 19])
        assert c.passenger == tx[s : s + 21]
        assert c.guide[:19] == reverse_complement(c.sense)
        assert c.chemistry == "all_RNA"
    assert enumerate_duplexes("T", "A" * 22) == []
    assert len(enumerate_duplexes("T", "A" * 23)) == 1
    assert not enumerate_duplexes("T", "A" * 11 + "N" + "A" * 11)


def test_chemistry_changes_identity_without_changing_legacy_id():
    old = Core("T", 1, 19, "A" * 19)
    assert old.molecule_id == Core("T", 1, 19, "A" * 19, chemistry="unknown").molecule_id
    assert old.molecule_id != Core("T", 1, 19, "A" * 19, chemistry="all_RNA").molecule_id


def test_longer_tool_guide_supports_suffix_and_never_votes_for_unknown_chemistry(random_dna):
    tx = random_dna(70)
    rows = enumerate_duplexes("T", tx)
    hit = Core("T", 8, 21, tx[7:28])
    assert (
        attach_evidence(rows, [{"core": hit, "tools": ["forge"], "tool_scores": {"forge": 0.9}}])
        == 1
    )
    attached = [r for r in rows if r["site_evidence"]]
    assert attached[0]["core"].start == 10
    assert attached[0]["core"].sense == hit.sense[-19:]
    assert attached[0]["tools"] == []
    assert rows[0]["site_evidence"] == []


def test_repeat_at_another_coordinate_does_not_receive_evidence():
    rows = enumerate_duplexes("T", "ACGT" * 30)
    c = rows[1]["core"]
    attach_evidence(rows, [{"core": c, "tools": ["a"], "tool_scores": {}}])
    assert sum(bool(r["site_evidence"]) for r in rows) == 1
    assert rows[1]["tools"] == ["a"]


def test_run_requires_an_external_generator(tmp_path, fasta_file):
    with pytest.raises(ValueError, match="at least one tool"):
        run(Config(str(fasta_file), str(tmp_path / "out"), tools=()))


def test_cli_default_and_legacy_are_explicit(tmp_path, fasta_file):
    with pytest.raises(SystemExit):
        main(["run", "--profile", "rna19p2", "--fasta", str(fasta_file), "--out", str(tmp_path / "new")])
    with pytest.raises(ValueError, match="at least one tool"):
        run(Config(str(fasta_file), str(tmp_path / "old"), tools=(), profile="legacy-v2"))


def test_removed_model_screen_options_are_rejected(tmp_path, fasta_file):
    config = Config(
        str(fasta_file), str(tmp_path / "out"), tools=("sirnadesigner",), screen_manifest=str(tmp_path / "missing")
    )
    with pytest.raises(ValueError, match="do not accept"):
        run(config)


def test_spliced_utr_coordinates_on_both_strands():
    tx = "ABCDEFGHIJKLMNO"
    plus = {"exons": [(100, 110), (200, 205)], "coding": [(103, 110), (200, 202)], "strand": "+"}
    minus = {"exons": [(100, 105), (200, 210)], "coding": [(103, 105), (200, 208)], "strand": "-"}
    assert utr3(tx, plus) == "MNO"
    assert utr3(tx, minus) == "MNO"
    assert utr3(tx + "X", plus) is None


def test_removed_replay_option_is_rejected(tmp_path, fasta_file):
    with pytest.raises(ValueError, match="do not accept"):
        run(Config(str(fasta_file), str(tmp_path / "out"), tools=("sirnadesigner",), evidence_run="old.json"))


def test_reference_alternate_loci_share_pinned_hgnc_identity(tmp_path):
    from sirna_consensus.reference import canonicalize_genes
    from sirna_consensus.screening import _resolve_gene

    hgnc = tmp_path/'hgnc.tsv'
    hgnc.write_text('symbol\thgnc_id\tstatus\tprev_symbol\talias_symbol\n'
                    'TARGET\tHGNC:1\tApproved\tOLD\tALIAS\n')
    mapping,names,aliases = canonicalize_genes(
        {'TA':'ENSG1.1','TA_ALT':'ENSG2.3','OTHER':'ENSG3.2'},
        {'ENSG1.1':'TARGET','ENSG2.3':'OLD','ENSG3.2':'Unannotated'},hgnc)
    assert mapping['TA']==mapping['TA_ALT']=='HGNC:1'
    assert mapping['OTHER']=='ENSG3.2'
    assert _resolve_gene('ENSG2',mapping,names,aliases)=='HGNC:1'
    assert _resolve_gene('TARGET',mapping,names,aliases)=='HGNC:1'
