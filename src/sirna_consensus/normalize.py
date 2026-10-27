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

from .types import Core, RawHit

OFFSET_TOLERANCE = 2


@dataclass
class Placement:
    core: Core | None = None
    reason: str = ""

    @property
    def ok(self) -> bool:
        return self.core is not None


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
    if not sense:
        return Placement(reason="empty_sequence")
    if len(sense) not in lengths:
        return Placement(reason=f"length_{len(sense)}_not_supported")

    positions = _positions(sequence, sense)
    if hit.start is not None:
        expected = int(hit.start)
        if 1 <= expected <= len(sequence) - len(sense) + 1:
            if sequence[expected - 1: expected - 1 + len(sense)] == sense:
                return Placement(core=Core(transcript_id, expected, len(sense), sense))
        # The tool reported a coordinate that does not agree with its own
        # sequence. Trust the sequence, but only if it is unambiguous there.
        if len(positions) == 1:
            return Placement(core=Core(transcript_id, positions[0], len(sense), sense))
        if positions:
            nearest = min(positions, key=lambda p: abs(p - expected))
            return Placement(core=Core(transcript_id, nearest, len(sense), sense))
        return Placement(reason="start_disagrees_with_sequence")

    if len(positions) == 1:
        return Placement(core=Core(transcript_id, positions[0], len(sense), sense))
    if not positions:
        return Placement(reason="not_found_in_transcript")
    return Placement(reason="ambiguous_repeat")


def normalize(hits: list[RawHit], transcript_id: str, sequence: str,
              lengths: tuple[int, ...] = (19, 20, 21)) -> tuple[list[Core], NormalizeStats]:
    """Place and de-duplicate a hit list. Keeps one core per (tool, position)."""
    stats = NormalizeStats()
    seen: set[tuple[str, int, int]] = set()
    cores: list[Core] = []
    for hit in hits:
        stats.total += 1
        hit.normalise_sequence()
        placement = place(hit, transcript_id, sequence, lengths)
        if not placement.ok:
            stats.reject(placement.reason)
            continue
        core = placement.core
        assert core is not None
        dedup_key = (hit.tool, core.start, core.length)
        if dedup_key in seen:
            stats.duplicates += 1
            continue
        seen.add(dedup_key)
        cores.append(core)
        stats.placed += 1
    return cores, stats
