"""Core data model shared by every stage.

Everything downstream of a tool adapter speaks in :class:`Core`: a stretch of the
target transcript on the *sense* strand, in DNA alphabet, with 1-based inclusive
coordinates.  A guide is always the reverse complement of its core, so no stage
ever invents a base.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

DNA_COMPLEMENT = str.maketrans("ACGT", "TGCA")
VALID_LENGTHS = (19, 20, 21)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(DNA_COMPLEMENT)[::-1]


def to_dna(sequence: str) -> str:
    return sequence.upper().replace("U", "T").replace("T", "T").strip()


@dataclass(frozen=True)
class Core:
    """A candidate guide/passenger pair anchored on the target transcript."""

    transcript_id: str
    start: int
    length: int
    sense: str

    def __post_init__(self) -> None:
        if self.length != len(self.sense):
            raise ValueError(
                f"length {self.length} disagrees with sequence {self.sense!r}")
        if self.start < 1:
            raise ValueError("start is 1-based and inclusive")

    @property
    def end(self) -> int:
        return self.start + self.length - 1

    @property
    def guide(self) -> str:
        return reverse_complement(self.sense)

    @property
    def key(self) -> str:
        return f"{self.transcript_id}:{self.start}:{self.length}"

    def interval(self) -> tuple[int, int]:
        return self.start, self.end


@dataclass
class RawHit:
    """One tool's claim about one candidate, before normalisation.

    Adapters decode their own tool's quirks (strand, overhang, alphabet) and hand
    back a sense-strand DNA core.  ``start`` may be ``None`` when the tool does
    not report coordinates; the normaliser then locates the sequence itself.
    """

    tool: str
    sense: str
    length: int
    start: int | None = None
    score: float | None = None
    score_kind: str = "tool_score"
    extra: dict[str, Any] = field(default_factory=dict)

    def normalise_sequence(self) -> "RawHit":
        self.sense = to_dna(self.sense)
        self.length = len(self.sense)
        return self


@dataclass
class ToolRun:
    """Bookkeeping for one adapter invocation."""

    tool: str
    status: str  # ok | unavailable | failed | empty
    hits: list[RawHit] = field(default_factory=list)
    detail: str = ""
    version: str = ""
    command: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.status == "ok"


@dataclass
class Candidate:
    """A consensus candidate: one representative core plus its evidence."""

    core: Core
    votes: int
    tools: list[str]
    tool_scores: dict[str, float]
    filter_flags: list[str] = field(default_factory=list)
    soft_flags: list[str] = field(default_factory=list)
    rank: int = 0

    @property
    def consensus_tier(self) -> str:
        return "unanimous" if self.votes == len(self.tools) else f"{self.votes}of_n"

    @property
    def passed(self) -> bool:
        return not self.filter_flags
