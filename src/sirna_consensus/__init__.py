"""sirna-consensus — consensus de novo siRNA design from several external tools.

Input: one human transcript in FASTA.
Output: ranked candidate guides with an explicit consensus tier.

The pipeline deliberately does not re-implement what the tools already do. It
adds the layer none of them provide: one coordinate axis, independent votes, and
filters that stay under the caller's control.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .consensus import Cluster, agree, cluster  # noqa: E402
from .filters import hard_flags, soft_flags  # noqa: E402
from .normalize import normalize, place  # noqa: E402
from .pipeline import Config, run  # noqa: E402
from .types import Candidate, Core, RawHit, ToolRun, reverse_complement  # noqa: E402

__all__ = [
    "__version__", "Config", "run", "Core", "RawHit", "Candidate", "ToolRun",
    "Cluster", "cluster", "agree", "normalize", "place", "hard_flags",
    "soft_flags", "reverse_complement",
]
