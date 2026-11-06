"""siRNA de novo duplex design with separate scoring and screening evidence."""
from __future__ import annotations

__version__ = "2.0.0"

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
