"""Put every tool's output on one coordinate axis.

This is the stage that replaces "enumerate a window for every position": we only
place the candidates we actually received, onto the input transcript.

A hit is rejected (with a reason) when its sequence cannot be placed
unambiguously. Repeats are the only genuinely hard case: if a sequence occurs at
several positions and the tool did not report coordinates, we cannot tell which
copy it meant, so we drop it rather than guess.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from .types import Core, RawHit, to_dna

OFFSET_TOLERANCE = 2


@dataclass
class Placement:
    core: Core | None = None
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.core is not None


@dataclass
class PlacedHit:
    """A placed hit with its score and evidence carried through normalisation."""
    tool: str
    core: Core
    score: float | None = None
    vote_eligible: bool = True
    extra: dict[str, Any] = field(default_factory=dict)


@dataclass
class NormalizeStats:
    total: int = 0
    placed: int = 0
    duplicates: int = 0
    reasons: dict[str, int] = field(default_factory=dict)

    def reject(self, reason: str) -> None:
        self.reasons[reason] = self.reasons.get(reason, 0) + 1


def _positions(sequence: str, motif: str) -> list[int]:
    out, start = [], sequence.find(motif)
    while start != -1:
        out.append(start + 1)  # 1-based
        start = sequence.find(motif, start + 1)
    return out


def place(hit: RawHit, transcript_id: str, sequence: str,
          lengths: tuple[int, ...] = (19, 20, 21)) -> Placement:
    """Locate one hit on the transcript, or explain why it could not be placed."""
    sense = hit.sense
    if any(hit.extra.get(key) for key in (
            "source_guide_overhang", "source_passenger_overhang",
            "source_guide_modifications", "source_passenger_modifications")):
        return Placement(reason="unsupported_source_chemistry_requires_explicit_adapter")
    if not sense:
        return Placement(reason="empty_sequence")
    if len(sense) not in lengths:
        return Placement(reason=f"length_{len(sense)}_not_supported")

    positions = _positions(sequence, sense)
    def core_at(start: int) -> Core:
        return Core(transcript_id, start, len(sense), sense, **{
            key: to_dna(hit.extra.get(key, "")) for key in
            ("ordered_guide", "passenger", "guide_overhang", "passenger_overhang")})
    if hit.start is not None:
        expected = int(hit.start)
        if (1 <= expected <= len(sequence) - len(sense) + 1
                and sequence[expected - 1: expected - 1 + len(sense)] == sense):
            return Placement(core=core_at(expected))
        # The tool reported a coordinate that does not agree with its own
        # sequence. Trust the sequence, but only if it is unambiguous there.
        if len(positions) == 1:
            return Placement(core=core_at(positions[0]))
        if positions:
            return Placement(reason="ambiguous_repeat")
        return Placement(reason="not_found_in_transcript")

    if len(positions) == 1:
        return Placement(core=core_at(positions[0]))
    if not positions:
        return Placement(reason="not_found_in_transcript")
    return Placement(reason="ambiguous_repeat")


def normalize_placed(hits: list[RawHit], transcript_id: str, sequence: str,
                     lengths: tuple[int, ...] = (19, 20, 21)
                     ) -> tuple[list[PlacedHit], NormalizeStats]:
    """Place and de-duplicate hits, carrying each hit's score and evidence."""
    stats = NormalizeStats()
    seen: set[tuple[str, int, str]] = set()
    placed: list[PlacedHit] = []
    for hit in sorted(hits, key=lambda h: (h.tool, h.sense, h.start or 0,
                     -(h.score if h.score is not None else float('-inf')), repr(sorted(h.extra.items())))):
        stats.total += 1
        import math
        if hit.score is not None and not math.isfinite(hit.score):
            stats.reject("nonfinite_score")
            continue
        hit.normalise_sequence()
        placement = place(hit, transcript_id, sequence, lengths)
        if not placement.ok:
            stats.reject(placement.reason)
            continue
        core = placement.core
        assert core is not None
        dedup_key = (hit.tool, core.start, core.molecule_id)
        if dedup_key in seen:
            stats.duplicates += 1
            continue
        seen.add(dedup_key)
        placed.append(PlacedHit(tool=hit.tool, core=core, score=hit.score,
                                vote_eligible=hit.vote_eligible,
                                extra=dict(hit.extra)))
        stats.placed += 1
    return placed, stats


def normalize(hits: list[RawHit], transcript_id: str, sequence: str,
              lengths: tuple[int, ...] = (19, 20, 21)) -> tuple[list[Core], NormalizeStats]:
    """Compatibility wrapper: place and de-duplicate, returning only cores."""
    placed, stats = normalize_placed(hits, transcript_id, sequence, lengths)
    return [item.core for item in placed], stats
