"""Ranking.

Votes come first. A candidate that several independent tools picked beats one
that a single tool loved — that is the whole point of the exercise. Only inside
an equal-vote group does the quality tie-break apply, and that tie-break is
deliberately rule-based rather than another opaque model.
"""
from __future__ import annotations

from . import filters
from .consensus import Cluster
from .types import Candidate, Core

IDEAL_GC = 0.45


def _tool_score_term(scores: dict[str, float | None]) -> float:
    """Mean of the normalised in-[0,1] tool scores, if the tools provide any."""
    usable = [s for s in scores.values() if s is not None]
    if not usable:
        return 0.0
    usable = [min(max(s, 0.0), 1.0) for s in usable]
    return sum(usable) / len(usable)


def _quality_term(core: Core) -> float:
    gc = filters.gc_fraction(core.sense)
    gc_term = max(0.0, 1.0 - abs(gc - IDEAL_GC) / 0.25)
    seed_term = 1.0 - abs(filters.seed_gc(core) - 0.5) * 2
    end_term = 1.0 if core.guide[0] in "AT" else 0.0
    return 0.6 * gc_term + 0.25 * max(seed_term, 0.0) + 0.15 * end_term


def score_components(cluster: Cluster) -> dict[str, float]:
    core = cluster.representative()
    return {
        "vote_component": 2.0 * cluster.votes,
        "tool_score_component": 0.5 * _tool_score_term(cluster.tools),
        "quality_component": 0.5 * _quality_term(core),
    }


def score_candidate(cluster: Cluster) -> float:
    return sum(score_components(cluster).values())


def build_candidates(clusters: list[Cluster], min_votes: int = 2,
                     include_singletons: bool = True,
                     acting_tools: int = 0,
                     filter_config: filters.FilterConfig = filters.DEFAULT_FILTERS
                     ) -> list[Candidate]:
    """Turn clusters into filtered, ranked candidates.

    ``include_singletons`` controls the *reporting* list: singletons are kept for
    inspection, but they are never promoted into the default Top-K because
    :func:`diverse_top` enforces ``min_votes``.
    """
    candidates: list[Candidate] = []
    for cluster in clusters:
        if cluster.votes < min_votes and not (include_singletons and cluster.votes == 1):
            continue
        core = cluster.representative()
        flags = filters.hard_flags(core, filter_config)
        parts = score_components(cluster)
        candidates.append(Candidate(
            core=core,
            votes=cluster.votes,
            tools=sorted(cluster.tools),
            tool_scores={k: v for k, v in cluster.tools.items() if v is not None},
            acting_tools=acting_tools,
            filter_flags=flags,
            soft_flags=filters.soft_flags(core, filter_config),
            score=sum(parts.values()),
            score_parts=parts,
        ))
    # A vote is worth 2.0 and everything else tops out at 1.0, so sorting by
    # score alone is the same lexicographic order as (votes, tool score,
    # quality) — the weighting makes votes strictly dominant.
    candidates.sort(key=lambda c: (-c.score, c.core.start))
    for position, candidate in enumerate(candidates, start=1):
        candidate.rank = position
    return candidates


def diverse_top(candidates: list[Candidate], top_k: int,
                min_spacing: int = 30, min_votes: int = 0) -> list[Candidate]:
    """Greedy spacer so the Top-K are not all one clump on the transcript.

    ``min_votes`` is a hard floor for the exported Top-K. Singletons remain in
    ``candidates.tsv`` / ``all_candidates.fa`` for inspection but do not enter
    the order-ready set unless the caller sets ``min_votes <= 1``.
    """
    kept: list[Candidate] = []
    for candidate in candidates:
        if candidate.votes < min_votes:
            continue
        if all(abs(candidate.core.start - k.core.start) >= min_spacing for k in kept):
            kept.append(candidate)
        if len(kept) == top_k:
            break
    return kept


def diverse_top_balanced(candidates: list[Candidate], top_k: int,
                         min_spacing: int = 30, min_votes: int = 0,
                         bin_size: int = 200, max_per_bin: int = 2
                         ) -> list[Candidate]:
    """Alternative diversity-aware Top-K.

    Same spacing floor as :func:`diverse_top`, plus a cap on how many selected
    guides may come from one transcript bin. This is a comparison mode, not the
    default, so the exported order-ready set remains the deterministic greedy one.
    """
    if bin_size <= 0 or max_per_bin <= 0:
        raise ValueError("bin_size and max_per_bin must be positive")
    kept: list[Candidate] = []
    bin_counts: dict[int, int] = {}
    for candidate in candidates:
        if candidate.votes < min_votes:
            continue
        if all(abs(candidate.core.start - k.core.start) >= min_spacing for k in kept):
            bin_id = candidate.core.start // bin_size
            if bin_counts.get(bin_id, 0) >= max_per_bin:
                continue
            kept.append(candidate)
            bin_counts[bin_id] = bin_counts.get(bin_id, 0) + 1
        if len(kept) == top_k:
            break
    return kept
