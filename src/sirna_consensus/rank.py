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


def score_candidate(cluster: Cluster) -> float:
    core = cluster.representative()
    return (2.0 * cluster.votes
            + 0.5 * _tool_score_term(cluster.tools)
            + 0.5 * _quality_term(core))


def build_candidates(clusters: list[Cluster], min_votes: int = 2,
                     include_singletons: bool = True) -> list[Candidate]:
    """Turn clusters into filtered, ranked candidates."""
    candidates: list[Candidate] = []
    for cluster in clusters:
        # `min_votes` is the reporting threshold; `include_singletons` is the
        # escape hatch that still surfaces single-tool sites (tagged by tier).
        if cluster.votes < min_votes and not (include_singletons and cluster.votes == 1):
            continue
        core = cluster.representative()
        flags = filters.hard_flags(core)
        candidates.append(Candidate(
            core=core,
            votes=cluster.votes,
            tools=sorted(cluster.tools),
            tool_scores={k: v for k, v in cluster.tools.items() if v is not None},
            filter_flags=flags,
            soft_flags=filters.soft_flags(core),
            score=score_candidate(cluster),
        ))
    # A vote is worth 2.0 and everything else tops out at 1.0, so sorting by
    # score alone is the same lexicographic order as (votes, tool score,
    # quality) — the weighting makes votes strictly dominant.
    candidates.sort(key=lambda c: (-c.score, c.core.start))
    for position, candidate in enumerate(candidates, start=1):
        candidate.rank = position
    return candidates


def diverse_top(candidates: list[Candidate], top_k: int,
                min_spacing: int = 30) -> list[Candidate]:
    """Greedy spacer so the Top-K are not all one clump on the transcript."""
    kept: list[Candidate] = []
    for candidate in candidates:
        if all(abs(candidate.core.start - k.core.start) >= min_spacing for k in kept):
            kept.append(candidate)
        if len(kept) == top_k:
            break
    return kept
