from sirna_consensus.consensus import agree, bounded_cluster, bucket, by_votes, cluster
from sirna_consensus.types import Core


def core(start: int, length: int = 21, tid: str = "TX") -> Core:
    return Core(tid, start, length, "ACGTACGTACGTACGTACGTA"[:length])


def test_same_site_agrees_across_lengths_and_offsets():
    assert agree(core(100, 21), core(102, 19))
    assert agree(core(102, 19), core(100, 21))


def test_far_apart_cores_do_not_agree():
    assert not agree(core(100), core(160))
    assert not agree(core(100), core(104))


def test_different_transcript_never_agrees():
    assert not agree(core(100), Core("OTHER", 100, 21, "ACGTACGTACGTACGTACGTA"))


def test_cluster_votes_count_distinct_tools():
    items = [("a", core(100), 0.5), ("b", core(101), 0.2), ("a", core(100), 0.9)]
    clusters = cluster(items)
    assert len(clusters) == 1
    assert clusters[0].votes == 2  # tool "a" twice still counts once
    assert clusters[0].tools["a"] == 0.9


def test_cluster_splits_disjoint_sites():
    items = [("a", core(100), None), ("b", core(300), None)]
    clusters = cluster(items)
    assert len(clusters) == 2
    assert all(c.votes == 1 for c in clusters)


def test_representative_is_modal_geometry():
    items = [("a", core(100, 21), None), ("b", core(102, 19), None),
             ("c", core(102, 19), None)]
    cluster_ = cluster(items)[0]
    assert cluster_.representative() == core(102, 19)


def test_bucket_groups_only_identical_sites():
    items = [("a", core(100), None), ("b", core(100), None), ("c", core(101), None)]
    buckets = bucket(items)
    assert sorted(b.votes for b in buckets) == [1, 2]


def test_bucket_is_not_defeated_by_a_saturated_tool():
    # a tool proposing every position must not collapse the transcript to one site
    items = [("a", core(s), None) for s in range(1, 51)]
    items += [("b", core(10), None)]
    buckets = bucket(items)
    assert len(buckets) == 50
    assert max(b.votes for b in buckets) == 2


def test_fuzzy_cluster_chains_the_whole_transcript():
    items = [("a", core(s), None) for s in range(1, 51)]
    assert len(cluster(items)) == 1


def test_by_votes_sorts_descending():
    items = [("a", core(100), None), ("b", core(100), None), ("c", core(300), None)]
    tiers = by_votes(cluster(items))
    assert list(tiers) == [2, 1]


def test_bounded_cluster_joins_shifted_19_and_21_mers():
    items = [("a", core(100, 21), None), ("b", core(102, 19), None)]
    clusters = bounded_cluster(items)
    assert len(clusters) == 1
    assert clusters[0].votes == 2
    assert clusters[0].anchor().start == 100


def test_bounded_cluster_does_not_chain_dense_overlaps():
    items = [("a", core(100, 21), None), ("b", core(102, 19), None),
             ("c", core(104, 21), None)]
    clusters = bounded_cluster(items)
    assert len(clusters) > 1
    for cluster_ in clusters:
        starts = [member.start for member in cluster_.members]
        assert max(starts) - min(starts) <= 2


def test_bounded_cluster_counts_distinct_tools_only_once():
    items = [("a", core(100), None), ("a", core(101), None),
             ("b", core(102), None)]
    clusters = bounded_cluster(items)
    assert len(clusters) == 1
    assert clusters[0].votes == 2     # tool a twice, tool b once
