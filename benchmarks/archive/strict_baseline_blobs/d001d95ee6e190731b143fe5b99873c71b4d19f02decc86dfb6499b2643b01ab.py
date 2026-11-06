from sirna_consensus.normalize import normalize, place
from sirna_consensus.types import RawHit

TID = "TX"
SEQ = "ACGTACGTACGTACGTACGTAAAA" + "TTTTGGGGCCCC" + "ACGTACGTACGTACGTACGT"


def test_places_hit_using_reported_start():
    hit = RawHit(tool="t", sense=SEQ[10:31], length=21, start=11)
    placement = place(hit, TID, SEQ)
    assert placement.ok
    assert placement.core.start == 11
    assert placement.core.sense == SEQ[10:31]


def test_places_hit_without_start_by_search():
    hit = RawHit(tool="t", sense=SEQ[10:31], length=21)
    placement = place(hit, TID, SEQ)
    assert placement.ok and placement.core.start == 11


def test_wrong_start_falls_back_to_sequence():
    hit = RawHit(tool="t", sense=SEQ[10:31], length=21, start=99)
    placement = place(hit, TID, SEQ)
    assert placement.ok
    assert placement.core.start == 11


def test_missing_sequence_is_reported():
    hit = RawHit(tool="t", sense="GCGCGCGCGCGCGCGCGCG", length=21)
    placement = place(hit, TID, SEQ)
    assert not placement.ok
    assert placement.reason in {"not_found_in_transcript", "ambiguous_repeat"}


def test_unsupported_length_is_reported():
    hit = RawHit(tool="t", sense=SEQ[0:25], length=25, start=1)
    placement = place(hit, TID, SEQ, lengths=(19, 20, 21))
    assert not placement.ok
    assert placement.reason == "length_25_not_supported"


def test_ambiguous_repeat_needs_coordinates():
    seq = "ACGTACGTACGTACGTACGTA" * 2
    hit = RawHit(tool="t", sense=seq[:21], length=21)
    assert place(hit, TID, seq).reason == "ambiguous_repeat"
    with_start = RawHit(tool="t", sense=seq[:21], length=21, start=22)
    assert place(with_start, TID, seq).ok


def test_normalize_deduplicates_per_tool():
    hit = RawHit(tool="t", sense=SEQ[10:31], length=21, start=11)
    cores, stats = normalize([hit, hit], TID, SEQ)
    assert len(cores) == 1
    assert stats.duplicates == 1
    assert stats.placed == 1


def test_normalize_keeps_two_tools_on_the_same_site():
    a = RawHit(tool="a", sense=SEQ[10:31], length=21, start=11)
    b = RawHit(tool="b", sense=SEQ[10:31], length=21, start=11)
    cores, stats = normalize([a, b], TID, SEQ)
    assert len(cores) == 2 and stats.placed == 2
