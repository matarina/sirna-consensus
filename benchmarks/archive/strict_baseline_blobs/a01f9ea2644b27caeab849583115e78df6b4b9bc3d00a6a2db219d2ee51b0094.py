import pytest

from sirna_consensus import filters
from sirna_consensus.consensus import cluster
from sirna_consensus.rank import build_candidates, diverse_top, diverse_top_balanced
from sirna_consensus.types import Core


def core(sense: str, start: int = 1) -> Core:
    return Core("TX", start, len(sense), sense)


def test_gc_window_rejects_extremes():
    assert "gc_outside_window" in filters.hard_flags(core("A" * 20 + "C"))
    assert "gc_outside_window" in filters.hard_flags(core("GC" * 10 + "C"))


def test_homopolymer_detected():
    assert "homopolymer" in filters.hard_flags(core("ACGTAAAACGTACGTACGTAA"))


def test_low_complexity_detected():
    assert "low_complexity" in filters.hard_flags(core("AT" * 10 + "A"))


def test_guide_5p_rule():
    # guide is the reverse complement of the sense; sense ending in C/G gives a
    # guide starting with G/C, which the rule rejects.
    bad = core("ACGTACGTACGTACGTACGTA"[:-1] + "C")
    assert "guide_5p_not_AU" in filters.hard_flags(bad)
    good = core("ACGTACGTACGTACGTACGTA")
    assert "guide_5p_not_AU" not in filters.hard_flags(good)


def test_internal_palindrome_flag():
    flags = filters.hard_flags(core("ACGTACGTACGTACGTACGTA"))
    assert isinstance(flags, list)


def test_soft_flags_annotate_seed():
    flags = filters.soft_flags(core("GCGCGCGCGCGCGCGCGCGC"))
    assert any(flag.startswith("seed_gc") for flag in flags)


def test_gc_fraction_and_seed_gc_bounds():
    c = core("ACGTACGTACGTACGTACGTA")
    assert 0.0 <= filters.gc_fraction(c.sense) <= 1.0
    assert 0.0 <= filters.seed_gc(c) <= 1.0


def test_votes_dominate_ranking():
    a = Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA")
    b = Core("TX", 300, 21, "ACGTACGTACGTACGTACGTA")
    clusters = cluster([("t1", a, 0.1), ("t2", a, 0.1), ("t1", b, 0.99)])
    candidates = build_candidates(clusters, min_votes=1)
    assert candidates[0].votes == 2
    assert candidates[0].core.start == 100
    assert candidates[1].votes == 1


def test_rank_numbers_are_contiguous():
    a = Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA")
    b = Core("TX", 300, 21, "ACGTACGTACGTACGTACGTA")
    candidates = build_candidates(cluster([("t1", a, None), ("t1", b, None)]),
                                  min_votes=1)
    assert [c.rank for c in candidates] == [1, 2]


def test_min_votes_filters_singletons():
    a = Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA")
    b = Core("TX", 300, 21, "ACGTACGTACGTACGTACGTA")
    clusters = cluster([("t1", a, None), ("t2", a, None), ("t3", b, None)])
    candidates = build_candidates(clusters, min_votes=2, include_singletons=False)
    assert [c.core.start for c in candidates] == [100]


def test_diverse_top_enforces_spacing():
    candidates = build_candidates(
        cluster([("t1", Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA"), None),
                 ("t2", Core("TX", 110, 21, "ACGTACGTACGTACGTACGTA"), None),
                 ("t3", Core("TX", 400, 21, "ACGTACGTACGTACGTACGTA"), None)]),
        min_votes=1)
    kept = diverse_top(candidates, top_k=2, min_spacing=30)
    assert [c.core.start for c in kept] == [100, 400]


def test_diverse_top_enforces_min_votes():
    candidates = build_candidates(
        cluster([("t1", Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA"), None),
                 ("t2", Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA"), None),
                 ("t3", Core("TX", 400, 21, "ACGTACGTACGTACGTACGTA"), None)]),
        min_votes=1, include_singletons=True)
    kept = diverse_top(candidates, top_k=2, min_spacing=30, min_votes=2)
    assert [c.core.start for c in kept] == [100]


def test_consensus_tier_denominator_is_acting_tools():
    candidates = build_candidates(
        cluster([("t1", Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA"), None)]),
        min_votes=1, acting_tools=2)
    assert candidates[0].consensus_tier == "1of_2"


def test_filter_config_is_customisable():
    c = core("ACGTACGTACGTACGTACGTA")
    default = filters.hard_flags(c)
    custom = filters.hard_flags(c, filters.FilterConfig(gc_min=0.0, gc_max=1.0,
                                                       max_homopolymer=2))
    assert "gc_outside_window" not in custom
    assert isinstance(default, list)


def test_filter_config_rejects_invalid_ranges():
    with pytest.raises(ValueError):
        filters.FilterConfig(gc_min=0.7, gc_max=0.5)
    with pytest.raises(ValueError):
        filters.FilterConfig(max_homopolymer=1)


def test_diverse_top_balanced_caps_per_bin():
    candidates = build_candidates(
        cluster([
            ("t1", Core("TX", 100, 21, "ACGTACGTACGTACGTACGTA"), None),
            ("t1", Core("TX", 130, 21, "ACGTACGTACGTACGTACGTA"), None),
            ("t1", Core("TX", 320, 21, "ACGTACGTACGTACGTACGTA"), None),
            ("t1", Core("TX", 520, 21, "ACGTACGTACGTACGTACGTA"), None),
        ]),
        min_votes=1)
    balanced = diverse_top_balanced(candidates, top_k=4, min_spacing=30,
                                    bin_size=200, max_per_bin=1)
    starts = [c.core.start for c in balanced]
    assert starts == [100, 320, 520]
