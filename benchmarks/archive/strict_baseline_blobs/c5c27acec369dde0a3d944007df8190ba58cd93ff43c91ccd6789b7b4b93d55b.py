"""Sequence-level hard filters and soft risk flags.

These are the published design rules (GC window, homopolymer, internal repeats,
strand-end preferences, seed characteristics). They are applied *after* the
consensus vote, on purpose: the external tools must not be allowed to decide what
survives, otherwise "use your own thresholds" would be a fiction.

Hard flags reject. Soft flags only annotate.
"""
from __future__ import annotations

import re
from dataclasses import asdict, dataclass

from .types import Core, reverse_complement

GC_MIN, GC_MAX = 0.30, 0.55
GU_RICH_CUTOFF = 0.75


@dataclass(frozen=True)
class FilterConfig:
    """Configurable hard/soft filter thresholds.

    Defaults are the published rules; callers may override any threshold via the
    CLI or the ``Config`` API.
    """
    gc_min: float = GC_MIN
    gc_max: float = GC_MAX
    max_homopolymer: int = 4
    gu_rich_cutoff: float = GU_RICH_CUTOFF
    seed_gc_high: float = 0.75
    seed_gc_low: float = 0.25

    def __post_init__(self) -> None:
        if not 0.0 <= self.gc_min < self.gc_max <= 1.0:
            raise ValueError("gc_min must be < gc_max, both in [0, 1]")
        if self.max_homopolymer < 2:
            raise ValueError("max_homopolymer must be >= 2")
        if not 0.0 < self.gu_rich_cutoff <= 1.0:
            raise ValueError("gu_rich_cutoff must be in (0, 1]")
        if not 0.0 <= self.seed_gc_low < self.seed_gc_high <= 1.0:
            raise ValueError("seed_gc_low must be < seed_gc_high, both in [0, 1]")

    def as_dict(self) -> dict[str, float | int]:
        return asdict(self)


DEFAULT_FILTERS = FilterConfig()


def gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / max(1, len(sequence))


def _low_complexity(sequence: str) -> bool:
    return len({sequence[i:i + 2] for i in range(len(sequence) - 1)}) <= 3


def _internal_palindrome(guide: str) -> bool:
    return len(guide) >= 16 and guide[:8] in reverse_complement(guide[8:])


def hard_flags(core: Core, config: FilterConfig = DEFAULT_FILTERS) -> list[str]:
    """Reasons this core is rejected outright."""
    flags: list[str] = []
    gc = gc_fraction(core.sense)
    if not config.gc_min <= gc <= config.gc_max:
        flags.append("gc_outside_window")
    if re.search(rf"(.)\1{{{config.max_homopolymer - 1},}}", core.sense):
        flags.append("homopolymer")
    if _low_complexity(core.sense):
        flags.append("low_complexity")
    guide = core.guide
    if guide[0] not in "AT":
        flags.append("guide_5p_not_AU")
    if _internal_palindrome(guide):
        flags.append("internal_palindrome")
    return flags


def soft_flags(core: Core, config: FilterConfig = DEFAULT_FILTERS) -> list[str]:
    """Non-fatal annotations kept next to the candidate."""
    flags: list[str] = []
    guide = core.guide
    if (guide.count("G") + guide.count("T")) >= len(guide) * config.gu_rich_cutoff:
        flags.append("guide_gu_rich")
    seed = guide[1:8]
    gc_seed = (seed.count("G") + seed.count("C")) / len(seed)
    if gc_seed > config.seed_gc_high:
        flags.append("seed_gc_high")
    elif gc_seed < config.seed_gc_low:
        flags.append("seed_gc_low")
    if core.sense[-1] not in "AT":
        flags.append("sense_3p_not_AU")
    return flags


def seed_gc(core: Core) -> float:
    seed = core.guide[1:8]
    return (seed.count("G") + seed.count("C")) / len(seed)
