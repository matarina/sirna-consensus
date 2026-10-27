"""Sequence-level hard filters and soft risk flags.

These are the published design rules (GC window, homopolymer, internal repeats,
strand-end preferences, seed characteristics). They are applied *after* the
consensus vote, on purpose: the external tools must not be allowed to decide what
survives, otherwise "use your own thresholds" would be a fiction.

Hard flags reject. Soft flags only annotate.
"""
from __future__ import annotations

import re

from .types import Core, reverse_complement

GC_MIN, GC_MAX = 0.30, 0.55
HOMOPOLYMER = re.compile(r"(.)\1{3,}")
GU_RICH_CUTOFF = 0.75


def gc_fraction(sequence: str) -> float:
    return (sequence.count("G") + sequence.count("C")) / max(1, len(sequence))


def _low_complexity(sequence: str) -> bool:
    return len({sequence[i:i + 2] for i in range(len(sequence) - 1)}) <= 3


def _internal_palindrome(guide: str) -> bool:
    return len(guide) >= 16 and guide[:8] in reverse_complement(guide[8:])


def hard_flags(core: Core, gc_range: tuple[float, float] = (GC_MIN, GC_MAX),
               max_homopolymer: int = 4) -> list[str]:
    """Reasons this core is rejected outright."""
    flags: list[str] = []
    gc = gc_fraction(core.sense)
    if not gc_range[0] <= gc <= gc_range[1]:
        flags.append("gc_outside_window")
    if re.search(rf"(.)\1{{{max_homopolymer - 1},}}", core.sense):
        flags.append("homopolymer")
    if _low_complexity(core.sense):
        flags.append("low_complexity")
    guide = core.guide
    if guide[0] not in "AT":
        flags.append("guide_5p_not_AU")
    if _internal_palindrome(guide):
        flags.append("internal_palindrome")
    return flags


def soft_flags(core: Core) -> list[str]:
    """Non-fatal annotations kept next to the candidate."""
    flags: list[str] = []
    guide = core.guide
    if (guide.count("G") + guide.count("T")) >= len(guide) * GU_RICH_CUTOFF:
        flags.append("guide_gu_rich")
    seed = guide[1:8]
    gc_seed = (seed.count("G") + seed.count("C")) / len(seed)
    if gc_seed > 0.75:
        flags.append("seed_gc_high")
    elif gc_seed < 0.25:
        flags.append("seed_gc_low")
    if core.sense[-1] not in "AT":
        flags.append("sense_3p_not_AU")
    return flags


def seed_gc(core: Core) -> float:
    seed = core.guide[1:8]
    return (seed.count("G") + seed.count("C")) / len(seed)
