"""sirna-consensus — consensus de novo siRNA design from several external tools.

Input: one human transcript in FASTA.
Output: ranked candidate guides with an explicit consensus tier.

The pipeline deliberately does not re-implement what the tools already do. It
adds the layer none of them provide: one coordinate axis, independent votes, and
filters that stay under the caller's control.
"""
from __future__ import annotations

__version__ = "1.1.0"

from .consensus import Cluster, agree, bounded_cluster, cluster  # noqa: E402
from .filters import FilterConfig, hard_flags, soft_flags  # noqa: E402
from .normalize import PlacedHit, normalize, normalize_placed, place  # noqa: E402
from .pipeline import Config, run  # noqa: E402
from .types import Candidate, Core, RawHit, ToolRun, reverse_complement  # noqa: E402

__all__ = [
    "Candidate",
    "Cluster",
    "Config",
    "Core",
    "FilterConfig",
    "PlacedHit",
    "RawHit",
    "ToolRun",
    "__version__",
    "agree",
    "bounded_cluster",
    "cluster",
    "hard_flags",
    "normalize",
    "normalize_placed",
    "place",
    "reverse_complement",
    "run",
    "soft_flags",
]
