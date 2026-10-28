"""Site-level consensus.

An exact-sequence intersection between three tools is usually empty: they emit
different lengths (19/21), different offsets and different strands for the *same*
duplex. So agreement is defined on the **site**, not on the string:

    two cores agree when they lie on the same transcript, start within
    ``offset_tolerance`` of each other, and share at least ``length - 2`` bases.

Clusters are built with union-find. A cluster's vote count is the number of
*distinct tools* in it (a tool appearing twice still counts once), so a single
verbose tool cannot outvote the others.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .types import Core

OFFSET_TOLERANCE = 2
SHARED_SLACK = 2


def _overlap(a: Core, b: Core) -> int:
    return min(a.end, b.end) - max(a.start, b.start) + 1


def agree(a: Core, b: Core, offset_tolerance: int = OFFSET_TOLERANCE,
          shared_slack: int = SHARED_SLACK) -> bool:
    if a.transcript_id != b.transcript_id:
        return False
    if abs(a.start - b.start) > offset_tolerance:
        return False
    return _overlap(a, b) >= min(a.length, b.length) - shared_slack


@dataclass
class Cluster:
    members: list[Core] = field(default_factory=list)
    tools: dict[str, float | None] = field(default_factory=dict)

    @property
    def votes(self) -> int:
        return len(self.tools)

    def add(self, tool: str, core: Core, score: float | None) -> None:
        self.members.append(core)
        previous = self.tools.get(tool)
        self.tools[tool] = (max(previous, score) if previous is not None and score is not None
                            else score if score is not None else previous)

    def anchor(self) -> Core:
        """Leftmost geometry in the cluster, used as the bounded-site anchor."""
        return min(self.members, key=lambda c: (c.start, c.length))

    def representative(self) -> Core:
        """Most-shared geometry in the cluster; longest, then leftmost."""
        counts: dict[tuple[int, int], int] = {}
        for member in self.members:
            counts[(member.start, member.length)] = counts.get(
                (member.start, member.length), 0) + 1
        best = max(counts.items(),
                   key=lambda kv: (kv[1], kv[0][1], -kv[0][0]))[0]
        for member in self.members:
            if (member.start, member.length) == best:
                return member
        raise AssertionError("unreachable")


class _UnionFind:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))

    def find(self, node: int) -> int:
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left: int, right: int) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def cluster(items: list[tuple[str, Core, float | None]],
            offset_tolerance: int = OFFSET_TOLERANCE,
            shared_slack: int = SHARED_SLACK) -> list[Cluster]:
    """Fuzzy grouping: union-find over :func:`agree`.

    Position-indexed neighbour comparison: two cores can only agree when their
    starts differ by at most ``offset_tolerance``, so the sweep is O(n log n)
    instead of the old all-pairs O(n^2).

    Coarse by construction — a tool that proposes every window chains the whole
    transcript into one component. Use :func:`bounded_cluster` or :func:`bucket`
    unless you explicitly want that.
    """
    ordered = sorted(items, key=lambda x: (x[1].start, x[1].length))
    uf = _UnionFind(len(ordered))
    for i in range(len(ordered)):
        start_i = ordered[i][1].start
        for j in range(i + 1, len(ordered)):
            start_j = ordered[j][1].start
            if start_j - start_i > offset_tolerance:
                break
            if agree(ordered[i][1], ordered[j][1], offset_tolerance, shared_slack):
                uf.union(i, j)
    buckets: dict[int, Cluster] = {}
    for index, (tool, core, score) in enumerate(ordered):
        root = uf.find(index)
        buckets.setdefault(root, Cluster()).add(tool, core, score)
    return list(buckets.values())


def bounded_cluster(items: list[tuple[str, Core, float | None]],
                    offset_tolerance: int = OFFSET_TOLERANCE,
                    shared_slack: int = SHARED_SLACK) -> list[Cluster]:
    """Bounded site grouping without transitive chains.

    Each cluster is anchored at its leftmost member. A new hit may join only if
    its start is within ``offset_tolerance`` of that anchor and it agrees with
    the anchor. This keeps a chain A-B-C where A and C are farther apart than
    ``offset_tolerance`` from collapsing into one site. A hit that bridges two
    anchors is assigned to the nearest anchor, not to both.
    """
    ordered = sorted(items, key=lambda x: (x[1].start, x[1].length))
    clusters: list[Cluster] = []
    for tool, core, score in ordered:
        target: Cluster | None = None
        best_distance: int | None = None
        for candidate in clusters:
            anchor = candidate.anchor()
            distance = abs(core.start - anchor.start)
            if (distance <= offset_tolerance
                    and agree(anchor, core, offset_tolerance, shared_slack)
                    and (best_distance is None or distance < best_distance)):
                best_distance = distance
                target = candidate
        if target is None:
            target = Cluster()
            clusters.append(target)
        target.add(tool, core, score)
    return clusters


def bucket(items: list[tuple[str, Core, float | None]]) -> list[Cluster]:
    """Group by exact site identity ``(transcript, start, length)``.

    Why exact, when the tools use different conventions?  Because
    :mod:`sirna_consensus.normalize` has already resolved those conventions: a
    hit is converted to its sense core on the transcript, so two tools that
    really mean the same duplex land on the same key.  Single-linkage merging on
    top of that is too coarse — with ``offset_tolerance=2`` a tool that proposes
    every window dissolves the whole transcript into one cluster, because
    position *n* always agrees with position *n+1*.
    """
    buckets: dict[tuple[str, int, int], Cluster] = {}
    for tool, core, score in items:
        key = (core.transcript_id, core.start, core.molecule_id)
        buckets.setdefault(key, Cluster()).add(tool, core, score)
    return list(buckets.values())


def by_votes(clusters: list[Cluster]) -> dict[int, list[Cluster]]:
    tiers: dict[int, list[Cluster]] = {}
    for item in clusters:
        tiers.setdefault(item.votes, []).append(item)
    return dict(sorted(tiers.items(), reverse=True))
