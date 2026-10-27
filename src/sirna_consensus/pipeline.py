"""Orchestration.

    FASTA -> validate -> fan out to tools -> normalise -> consensus -> filter -> rank

The expensive per-candidate analysis (structure, specificity) is *not* re-run
here; it is left to whatever the tools already do. What this pipeline adds is the
part none of them do: putting independent opinions on one axis and refusing to
let a single one of them decide.

Design note: the tools are told to do their own window enumeration. There is no
duplicated enumeration stage — the only coordinate work is placing the candidates
that come back (see :mod:`sirna_consensus.normalize`).
"""
from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import __version__, report
from .consensus import Cluster, bucket
from .fasta import read_transcript, sha256
from .normalize import normalize
from .rank import build_candidates, diverse_top
from .tools import DEFAULT_TOOLS, ToolEnv, get_tool
from .types import ToolRun

SCHEMA_VERSION = 1


@dataclass
class Config:
    fasta: str
    out_dir: str
    tools: tuple[str, ...] = DEFAULT_TOOLS
    lengths: tuple[int, ...] = (19, 20, 21)
    min_votes: int = 2
    top_k: int = 10
    min_spacing: int = 30
    site_mode: str = "exact"
    include_singletons: bool = True
    offset_tolerance: int = 2
    jobs: int = 3
    test: bool = False
    transcriptome: str = ""
    extra: dict = field(default_factory=dict)

    def resolve(self) -> "Config":
        self.fasta = str(Path(self.fasta).expanduser().resolve())
        self.out_dir = str(Path(self.out_dir).expanduser().resolve())
        return self


def _select_tools(config: Config) -> tuple[list, list[ToolRun]]:
    env = ToolEnv()
    if config.transcriptome:
        env.transcriptome = str(Path(config.transcriptome).expanduser().resolve())
    usable, skipped = [], []
    for name in config.tools:
        tool = get_tool(name, env)
        ok, detail = tool.availability()
        if ok:
            usable.append(tool)
        else:
            skipped.append(tool.wrap("unavailable", detail))
    return usable, skipped


def run(config: Config) -> dict:
    config.resolve()
    out_dir = Path(config.out_dir)
    work_dir = out_dir / "work"
    work_dir.mkdir(parents=True, exist_ok=True)

    transcript_id, sequence = read_transcript(config.fasta)
    input_sha = sha256(sequence)

    tools, runs = _select_tools(config)
    if not tools:
        detail = "; ".join(f"{r.tool}: {r.detail}" for r in runs) or "no tools requested"
        raise RuntimeError("no external tool is usable: " + detail)

    def execute(tool) -> ToolRun:
        target = work_dir / tool.name
        target.mkdir(parents=True, exist_ok=True)
        try:
            return tool.run(Path(config.fasta), target, transcript_id)
        except Exception as exc:  # noqa: BLE001 - one tool must not kill the run
            return tool.wrap("failed", f"{type(exc).__name__}: {exc}")

    with ThreadPoolExecutor(max_workers=max(1, min(config.jobs, len(tools)))) as pool:
        runs.extend(pool.map(execute, tools))
    runs.sort(key=lambda r: (r.tool not in [t.name for t in tools], r.tool))

    placed: list[tuple[str, object, float | None]] = []
    normalization = {}
    for tool_run in runs:
        if not tool_run.ok or not tool_run.hits:
            continue
        scores = {(hit.start, len(hit.sense)): hit.score for hit in tool_run.hits}
        cores, stats = normalize(tool_run.hits, transcript_id, sequence, config.lengths)
        normalization[tool_run.tool] = stats
        for core in cores:
            placed.append((tool_run.tool, core,
                           scores.get((core.start, core.length))))

    if not placed:
        raise RuntimeError(
            "no candidate could be placed on the transcript; check the tool logs "
            f"under {work_dir}")

    # Tool scores live on wildly different scales (a rule score can run 0..10
    # while an off-target cleanliness is 0..1). Normalise per tool by min-max
    # within this run so the tie-break is not decided by unit choice.
    spread: dict[str, list[float]] = {}
    for tool, _core, score in placed:
        if score is not None:
            spread.setdefault(tool, []).append(score)

    def scaled(tool: str, score: float | None) -> float | None:
        if score is None:
            return None
        values = spread.get(tool, [])
        low, high = (min(values), max(values)) if values else (score, score)
        return 0.5 if high == low else (score - low) / (high - low)

    placed = [(tool, core, scaled(tool, score)) for tool, core, score in placed]

    if config.site_mode == "fuzzy":
        from .consensus import cluster as fuzzy_cluster
        clusters: list[Cluster] = fuzzy_cluster(placed, config.offset_tolerance)
    else:
        clusters = bucket(placed)
    candidates = build_candidates(clusters, min_votes=config.min_votes,
                                 include_singletons=config.include_singletons)
    passing = [c for c in candidates if c.passed]
    top = diverse_top(passing, config.top_k, config.min_spacing)

    acting = sorted({r.tool for r in runs if r.ok and r.hits})
    result = {
        "schema_version": SCHEMA_VERSION,
        "tool_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "research_only",
        "status_reason": ("Candidates are computational predictions; no knockdown is "
                          "established by this run."),
        "transcript_id": transcript_id,
        "transcript_length": len(sequence),
        "input_sha256": input_sha,
        "config": {
            "tools": list(config.tools),
            "acting_tools": acting,
            "lengths": list(config.lengths),
            "min_votes": config.min_votes,
            "top_k": config.top_k,
            "min_spacing": config.min_spacing,
            "offset_tolerance": config.offset_tolerance,
            "source": "external tools; no local enumeration",
        },
        "tools": [
            {"tool": r.tool, "status": r.status, "detail": r.detail,
             "version": r.version, "hits": len(r.hits)}
            for r in runs
        ],
        "normalization": {
            name: {"total": s.total, "placed": s.placed, "duplicates": s.duplicates,
                   "rejected": s.reasons}
            for name, s in normalization.items()
        },
        # `sites` is every agreeing site; `candidates` is what survived the vote
        # threshold. Keeping them apart stops "rejected by filters" from also
        # swallowing the sites that were never candidates.
        "counts": {
            "hits": sum(len(r.hits) for r in runs),
            "placed": len(placed),
            "sites": len(clusters),
            "below_min_votes": len(clusters) - len(candidates),
            "candidates": len(candidates),
            "rejected_by_filters": len(candidates) - len(passing),
            "consensus_unanimous": sum(1 for c in candidates
                                       if c.votes == len(acting)),
            "passing": len(passing),
            "top": len(top),
        },
        "candidates": candidates,
        "top": top,
        "limitations": [
            "Site-level consensus, not sequence-level intersection: tools emit "
            "different lengths and offsets for the same duplex.",
            "A tool that returns only its own Top-N cannot be distinguished from "
            "one that filtered a site out.",
            "No wet-lab confirmation; tool scores are not comparable across tools.",
            "Off-target scope is whatever each tool was given, unless a shared "
            "transcriptome was supplied.",
        ],
    }
    report.write_all(result, out_dir)
    return result


def load_result(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())
