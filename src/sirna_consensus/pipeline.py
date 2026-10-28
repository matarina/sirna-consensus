"""Orchestration.

    FASTA -> validate -> fan out to tools -> normalise -> consensus -> filter -> rank

The expensive per-candidate analysis (structure, specificity) is *not* re-run
here; it is left to whatever the tools already do. What this pipeline adds is the
part none of them do: putting independent opinions on one axis and refusing to
let a single one of them decide.

Design note: the tools are told to do their own window enumeration. There is no
duplicated enumeration stage — the only coordinate work is placing the candidates
that come back (see :mod:`sirna_consensus.normalize`).

Tool roles are explicit:

* **generators** contribute design votes (siRNAdesigner, siRNAforge);
* **evaluators** contribute specificity evidence, not design votes (SIREN).

A SIREN enumeration hit is therefore not a positive design vote. It can only
support a generated site when a transcriptome-backed screen completed and the
site's off-target burden is within the configured threshold.
"""
from __future__ import annotations

import json
import uuid
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from . import __version__, filters, report
from .consensus import Cluster, agree, bounded_cluster, bucket
from .consensus import cluster as fuzzy_cluster
from .fasta import read_transcript, sha256, sha256_file
from .normalize import PlacedHit, normalize_placed
from .rank import build_candidates, diverse_top
from .tools import DEFAULT_TOOLS, ToolEnv, get_tool
from .types import Core, ToolRun

SCHEMA_VERSION = 2


@dataclass
class Config:
    fasta: str
    out_dir: str
    tools: tuple[str, ...] = DEFAULT_TOOLS
    lengths: tuple[int, ...] = (19, 20, 21)
    min_votes: int = 2
    top_k: int = 10
    min_spacing: int = 30
    site_mode: str = "bounded"       # bounded | exact | fuzzy
    candidate_mode: str = "exact"    # exact molecules; legacy_site for baseline only
    include_singletons: bool = True
    offset_tolerance: int = 2
    jobs: int = 3
    test: bool = False
    transcriptome: str = ""
    strict_specificity: bool = False
    siren_max_offtargets: int = 0
    filter_config: filters.FilterConfig = field(
        default_factory=lambda: filters.DEFAULT_FILTERS)
    extra: dict = field(default_factory=dict)

    def resolve(self) -> "Config":
        self.fasta = str(Path(self.fasta).expanduser().resolve())
        self.out_dir = str(Path(self.out_dir).expanduser().resolve())
        if self.transcriptome:
            self.transcriptome = str(
                Path(self.transcriptome).expanduser().resolve())
        return self

    def validate(self) -> None:
        if not self.tools:
            raise ValueError("at least one tool is required")
        if self.min_votes < 1:
            raise ValueError("min_votes must be >= 1")
        if self.top_k < 1:
            raise ValueError("top_k must be >= 1")
        if self.min_spacing < 0:
            raise ValueError("min_spacing must be >= 0")
        if self.offset_tolerance < 0:
            raise ValueError("offset_tolerance must be >= 0")
        if self.jobs < 1:
            raise ValueError("jobs must be >= 1")
        if not self.lengths or any(length <= 0 for length in self.lengths):
            raise ValueError("lengths must be positive integers")
        if self.site_mode not in {"bounded", "exact", "fuzzy"}:
            raise ValueError("site_mode must be bounded, exact or fuzzy")
        if self.candidate_mode not in {"exact", "legacy_site"}:
            raise ValueError("candidate_mode must be exact or legacy_site")
        if self.siren_max_offtargets < 0:
            raise ValueError("siren_max_offtargets must be >= 0")
        # FilterConfig validates itself at construction.


def _select_tools(config: Config) -> tuple[list, list[ToolRun], ToolEnv]:
    env = ToolEnv()
    if config.transcriptome:
        env.transcriptome = config.transcriptome
    env.siren_max_offtargets = config.siren_max_offtargets
    usable, skipped = [], []
    for name in config.tools:
        tool = get_tool(name, env)
        ok, detail = tool.availability()
        if ok:
            usable.append(tool)
        else:
            skipped.append(tool.wrap("unavailable", detail))
    return usable, skipped, env


def _specificity_status(
        core: Core,
        evaluator_hits: list[PlacedHit],
        evaluator_runs: list[ToolRun],
        config: Config) -> tuple[str, list[str], int | None]:
    """Classify a generated site against available evaluator evidence."""
    matches = [hit for hit in evaluator_hits
               if hit.tool in {r.tool for r in evaluator_runs}
               and core.transcript_id == hit.core.transcript_id
               and core.guide == hit.core.guide]
    evaluator_tools = sorted({hit.tool for hit in matches})

    if not evaluator_runs:
        return "not_screened", evaluator_tools, None
    if not matches:
        return "not_screened", evaluator_tools, None

    verified = [hit for hit in matches if hit.extra.get("screen_state") == "transcriptome_screened"
                and hit.extra.get("transcriptome_backed") is True
                and isinstance(hit.extra.get("offtargets"), int)
                and hit.extra["offtargets"] >= 0]
    if verified:
        burden = max(hit.extra["offtargets"] for hit in verified)
        return ("clean" if burden <= config.siren_max_offtargets else "offtarget",
                evaluator_tools, burden)
    states = {hit.extra.get("screen_state", "") for hit in matches}
    if "transcriptome_screened" in states:
        return "screen_failed", evaluator_tools, None
    if "input_only" in states:
        return "input_only", evaluator_tools, None
    if "failed" in states:
        return "screen_failed", evaluator_tools, None
    return "not_screened", evaluator_tools, None


def run(config: Config) -> dict:
    config.resolve()
    config.validate()
    out_dir = Path(config.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_id = f"{stamp}-{uuid.uuid4().hex[:8]}"
    work_dir = out_dir / "work" / run_id
    work_dir.mkdir(parents=True, exist_ok=True)

    transcript_id, sequence = read_transcript(config.fasta)
    input_sha = sha256(sequence)

    transcriptome_sha = None
    if config.transcriptome:
        transcriptome_path = Path(config.transcriptome)
        if transcriptome_path.is_file():
            transcriptome_sha = sha256_file(transcriptome_path)

    tools, runs, tool_env = _select_tools(config)
    if not tools:
        detail = "; ".join(f"{r.tool}: {r.detail}" for r in runs) or "no tools requested"
        raise RuntimeError("no external tool is usable: " + detail)

    def execute(tool) -> ToolRun:
        # Each run gets a unique work directory. A reused output directory can no
        # longer feed stale *_csv / off_targets_summary.tsv files to a new run.
        target = work_dir / tool.name
        target.mkdir(parents=True, exist_ok=True)
        try:
            return tool.run(Path(config.fasta), target, transcript_id)
        except Exception as exc:  # noqa: BLE001 - one tool must not kill the run
            return tool.wrap("failed", f"{type(exc).__name__}: {exc}")

    with ThreadPoolExecutor(max_workers=max(1, min(config.jobs, len(tools)))) as pool:
        runs.extend(pool.map(execute, tools))
    runs.sort(key=lambda r: (r.tool not in [t.name for t in tools], r.tool))

    placed_generator: list[tuple[str, object, float | None]] = []
    evaluator_hits: list[PlacedHit] = []
    normalization = {}
    normalized_hits = []
    for tool_run in runs:
        if not tool_run.ok or not tool_run.hits:
            continue
        placed_hits, stats = normalize_placed(
            tool_run.hits, transcript_id, sequence, config.lengths)
        normalization[tool_run.tool] = stats
        normalized_hits.extend({**asdict(hit), "role": tool_run.role} for hit in placed_hits)
        if tool_run.role == "generator":
            for item in placed_hits:
                if item.vote_eligible:
                    placed_generator.append((tool_run.tool, item.core, item.score))
        else:
            evaluator_hits.extend(placed_hits)

    if not placed_generator and not evaluator_hits:
        raise RuntimeError(
            "no candidate could be placed on the transcript; check the tool logs "
            f"under {work_dir}")

    # Tool scores live on wildly different scales. Normalise per tool by min-max
    # within this run so the tie-break is not decided by unit choice. Scores are
    # carried through placement; they are never re-joined by raw coordinates.
    spread: dict[str, list[float]] = {}
    for tool, _core, score in placed_generator:
        if score is not None:
            spread.setdefault(tool, []).append(score)

    def scaled(tool: str, score: float | None) -> float | None:
        if score is None:
            return None
        values = spread.get(tool, [])
        low, high = (min(values), max(values)) if values else (score, score)
        return 0.5 if high == low else (score - low) / (high - low)

    placed_generator = [
        (tool, core, scaled(tool, score))
        for tool, core, score in placed_generator
    ]

    if config.candidate_mode == "exact":
        clusters = bucket(placed_generator)
    elif config.site_mode == "bounded":
        clusters: list[Cluster] = bounded_cluster(
            placed_generator, config.offset_tolerance)
    elif config.site_mode == "fuzzy":
        clusters = fuzzy_cluster(placed_generator, config.offset_tolerance)
    else:
        clusters = bucket(placed_generator)

    generator_runs = [
        r for r in runs if r.role == "generator" and r.ok and r.hits]
    evaluator_runs = [
        r for r in runs if r.role == "evaluator" and r.ok and r.hits]
    acting_generator = sorted({r.tool for r in generator_runs})

    candidates = build_candidates(
        clusters,
        min_votes=config.min_votes,
        include_singletons=config.include_singletons,
        acting_tools=len(acting_generator),
        filter_config=config.filter_config,
    )

    # Attach specificity evidence to generated sites. SIREN hits are evidence,
    # not design votes, and input-only/failed screens do not become
    # transcriptome-wide specificity claims.
    for candidate in candidates:
        candidate.nearby_tools = sorted({tool for tool, core, _ in placed_generator
                                        if agree(candidate.core, core, config.offset_tolerance)})
        status, evaluator_tools, off_min = _specificity_status(
            candidate.core, evaluator_hits, evaluator_runs, config)
        candidate.specificity_status = status
        candidate.evaluator_tools = evaluator_tools
        candidate.off_target_min = off_min
        if status == "offtarget":
            candidate.filter_flags.append("specificity_offtarget")
        if config.strict_specificity and status != "clean":
            candidate.filter_flags.append("specificity_unverified")

    passing = [c for c in candidates if c.passed]
    top = diverse_top(
        passing, config.top_k, config.min_spacing, min_votes=config.min_votes)

    if not placed_generator:
        status = "screening_only"
        status_reason = (
            "No candidate-generating tool produced placeable sites; evaluator "
            "evidence alone does not create design consensus.")
    elif not candidates:
        status = "no_consensus_candidates"
        status_reason = (
            "No site reached the vote threshold and singleton policy; relax "
            "min_votes or include singletons for inspection.")
    else:
        status = "research_only"
        status_reason = (
            "Candidates are computational predictions; no knockdown is "
            "established by this run.")

    consensus_candidates = sum(
        1 for c in candidates if c.votes >= config.min_votes)
    exploratory_singletons = sum(
        1 for c in candidates if c.votes == 1 and config.min_votes > 1)

    counts = {
        "hits": sum(len(r.hits) for r in runs),
        "generator_hits": sum(len(r.hits) for r in generator_runs),
        "evaluator_hits": sum(len(r.hits) for r in evaluator_runs),
        "placed": sum(stats.placed for stats in normalization.values()),
        "sites": len(clusters),
        "consensus_candidates": consensus_candidates,
        "exploratory_singletons": exploratory_singletons,
        "below_min_votes": len(clusters) - consensus_candidates,
        "candidates": len(candidates),
        "rejected_by_filters": len(candidates) - len(passing),
        "consensus_unanimous": sum(
            1 for c in candidates
            if c.votes >= config.min_votes
            and len(acting_generator) > 0
            and c.votes == len(acting_generator)),
        "passing": len(passing),
        "top": len(top),
        "evaluator_clean": sum(
            1 for c in candidates if c.specificity_status == "clean"),
    }

    result = {
        "schema_version": SCHEMA_VERSION,
        "implementation_hashes": {str(p.relative_to(Path(__file__).parent)): sha256_file(p)
                                  for p in sorted(Path(__file__).parent.rglob("*.py"))},
        "run_id": run_id,
        "tool_version": __version__,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": status,
        "status_reason": status_reason,
        "transcript_id": transcript_id,
        "transcript_length": len(sequence),
        "input_sha256": input_sha,
        "config": {
            "tools": list(config.tools),
            "generator_tools": acting_generator,
            "evaluator_tools": sorted({r.tool for r in evaluator_runs}),
            "acting_tools": acting_generator,
            "lengths": list(config.lengths),
            "min_votes": config.min_votes,
            "top_k": config.top_k,
            "min_spacing": config.min_spacing,
            "offset_tolerance": config.offset_tolerance,
            "site_mode": config.site_mode,
            "candidate_mode": config.candidate_mode,
            "include_singletons": config.include_singletons,
            "jobs": config.jobs,
            "strict_specificity": config.strict_specificity,
            "siren_max_offtargets": config.siren_max_offtargets,
            "filter_config": config.filter_config.as_dict(),
            "transcriptome": config.transcriptome,
            "transcriptome_sha256": transcriptome_sha,
            "siren_threads": tool_env.siren_threads,
            "tool_timeout_s": tool_env.timeout_s,
            "source": "external tools; no local enumeration",
        },
        "tools": [
            {
                "tool": r.tool,
                "role": r.role,
                "acceptance": r.acceptance,
                "status": r.status,
                "detail": r.detail,
                "version": r.version,
                "command": r.command,
                "hits": len(r.hits),
                "screen_state": r.screen_state,
                "selection_limits": r.selection_limits,
            }
            for r in runs
        ],
        "normalization": {
            name: {
                "total": s.total,
                "placed": s.placed,
                "duplicates": s.duplicates,
                "rejected": s.reasons,
            }
            for name, s in normalization.items()
        },
        "counts": counts,
        "normalized_hits": normalized_hits,
        "raw_output_hashes": {str(p.relative_to(out_dir)): sha256_file(p)
                              for p in sorted(work_dir.rglob("*")) if p.is_file()},
        "candidates": candidates,
        "top": top,
        "limitations": [
            "Consensus votes are generator-tool votes. SIREN is an evaluator; "
            "its enumeration hits and input-only fallback are not design votes.",
            "Specificity is separate evidence. A transcriptome-backed SIREN "
            "screen with off-targets is a hard eligibility failure; input-only "
            "or failed screens are not transcriptome-wide claims.",
            "Exact molecules are preserved by default; legacy site mode uses bounded clustering; "
            "exact identity and fuzzy single-linkage modes are explicit choices.",
            "Top-K enforces min_votes. Singleton sites remain in candidates.tsv "
            "and all_candidates.fa for inspection, not in the research-selected set.",
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
