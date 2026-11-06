import pytest

from sirna_consensus.fasta import InputError, read_transcript
from sirna_consensus.types import Core, RawHit, reverse_complement


def test_reverse_complement_roundtrip():
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement(reverse_complement("ACGTTGCA")) == "ACGTTGCA"


def test_core_guide_is_reverse_complement():
    core = Core("TX", 1, 21, "ACGTACGTACGTACGTACGTA")
    assert core.guide == reverse_complement(core.sense)
    assert core.end == 21


def test_core_rejects_length_mismatch():
    with pytest.raises(ValueError):
        Core("TX", 1, 19, "ACGT")
    with pytest.raises(ValueError):
        Core("TX", 0, 4, "ACGT")


def test_core_key_is_stable():
    assert Core("TX", 5, 19, "A" * 19).key == "TX:5:19"


def test_raw_hit_normalises_rna_and_length():
    hit = RawHit(tool="t", sense="acguacguacguacguacgua", length=0).normalise_sequence()
    assert hit.sense == "ACGTACGTACGTACGTACGTA"
    assert hit.length == 21


def test_read_transcript_ok(fasta_file, transcript, transcript_id):
    assert read_transcript(fasta_file) == (transcript_id, transcript)


def test_read_transcript_rejects_multiple_records(tmp_path):
    path = tmp_path / "two.fa"
    path.write_text(">a\nACGTACGTACGTACGTACGTA\n>b\nACGTACGTACGTACGTACGTA\n")
    with pytest.raises(InputError, match="exactly one"):
        read_transcript(path)


def test_read_transcript_rejects_ambiguous(tmp_path):
    path = tmp_path / "n.fa"
    path.write_text(">a\nACGTACGTACGTACGTACGTN\n")
    with pytest.raises(InputError, match="unambiguous"):
        read_transcript(path)


def test_read_transcript_rejects_empty(tmp_path):
    path = tmp_path / "e.fa"
    path.write_text(">a\n\n")
    with pytest.raises(InputError):
        read_transcript(path)
