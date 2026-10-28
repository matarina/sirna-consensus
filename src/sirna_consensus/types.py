"""Core data model shared by every stage.

Everything downstream of a tool adapter speaks in :class:`Core`: a stretch of the
target transcript on the *sense* strand, in DNA alphabet, with 1-based inclusive
coordinates.  A guide retains its original sequence; its paired portion is the reverse
complement of the core. Documented overhangs and passenger identity are separate.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
import hashlib

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
    ordered_guide: str = ""
    passenger: str = ""
    guide_overhang: str = ""
    passenger_overhang: str = ""

    def __post_init__(self) -> None:
        if self.length <= 0:
            raise ValueError("core length must be positive")
        if self.guide_overhang and not self.ordered_guide:
            raise ValueError("guide overhang requires the original ordered guide")
        if self.passenger_overhang and not self.passenger:
            raise ValueError("passenger overhang requires the original passenger")
        if self.length != len(self.sense):
            raise ValueError(
                f"length {self.length} disagrees with sequence {self.sense!r}")
        if self.start < 1:
            raise ValueError("start is 1-based and inclusive")
        for sequence in (self.sense, self.ordered_guide, self.passenger,
                         self.guide_overhang, self.passenger_overhang):
            if set(sequence) - set("ACGT"):
                raise ValueError("molecular sequences must use uppercase DNA alphabet")
        if self.ordered_guide and self.ordered_guide != reverse_complement(self.sense) + self.guide_overhang:
            raise ValueError("guide does not match paired core plus documented 3-prime overhang")
        if self.passenger and self.passenger != self.sense + self.passenger_overhang:
            raise ValueError("passenger does not match paired core plus documented 3-prime overhang")

    @property
    def end(self) -> int:
        return self.start + self.length - 1

    @property
    def guide(self) -> str:
        return self.ordered_guide or reverse_complement(self.sense)

    @property
    def duplex_format(self) -> str:
        passenger = f"passenger{len(self.passenger)}-overhang{len(self.passenger_overhang)}" if self.passenger else "passenger_unknown"
        return f"binding{self.length}-guide{len(self.guide)}-overhang{len(self.guide_overhang)}-{passenger}"

    @property
    def molecule_id(self) -> str:
        text = "|".join((self.guide, self.passenger, self.guide_overhang,
                         self.passenger_overhang, self.duplex_format))
        return hashlib.sha256(text.encode()).hexdigest()

    @property
    def key(self) -> str:
        base = f"{self.transcript_id}:{self.start}:{self.length}"
        return base + (":" + self.molecule_id[:16] if self.guide_overhang or self.passenger else "")

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
    #: A generator hit is a design claim; an evaluator hit is evidence. A hit
    #: marked False must not contribute a design vote or be used as clean
    #: specificity evidence.
    vote_eligible: bool = True

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
    role: str = "generator"          # generator | evaluator
    acceptance: str = ""
    screen_state: str = ""           # transcriptome_screened | input_only | failed
    selection_limits: dict[str, Any] = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.status == "ok"


@dataclass
class Candidate:
    """A consensus candidate: one representative core plus its evidence."""

    core: Core
    votes: int                       # distinct generator tools
    tools: list[str]                 # generator tools
    tool_scores: dict[str, float]
    acting_tools: int = 0            # generator tools that actually ran
    evaluator_tools: list[str] = field(default_factory=list)
    specificity_status: str = "not_screened"
    off_target_min: int | None = None
    filter_flags: list[str] = field(default_factory=list)
    soft_flags: list[str] = field(default_factory=list)
    score: float = 0.0
    score_parts: dict[str, float] = field(default_factory=dict)
    rank: int = 0
    nearby_tools: list[str] = field(default_factory=list)

    @property
    def consensus_tier(self) -> str:
        if self.acting_tools <= 0:
            return "no-generator"
        if self.votes == self.acting_tools:
            return "unanimous"
        return f"{self.votes}of_{self.acting_tools}"

    @property
    def passed(self) -> bool:
        return not self.filter_flags
