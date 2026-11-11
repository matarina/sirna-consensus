"""Orchestration for candidates emitted by registered external tools."""
from __future__ import annotations

import fcntl
import hashlib
import json
import shutil
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field
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


def _native_hit_record(hit, run: ToolRun, raw_hashes: dict[str, str]) -> dict:
    """Stable, JSON-safe representation of one adapter-native hit."""
    extra = dict(hit.extra or {})
    guide = extra.get("ordered_guide", "")
    passenger = extra.get("passenger", "")
    guide_overhang = extra.get("guide_overhang", "")
    passenger_overhang = extra.get("passenger_overhang", "")
    duplex_format = (f"binding{hit.length}-guide{len(guide)}-overhang{len(guide_overhang)}-"
                     f"passenger{len(passenger)}-overhang{len(passenger_overhang)}")
    return {
        "tool": hit.tool,
        "role": run.role,
        "coordinate": hit.start,
        "length": hit.length,
        "sense_dna": hit.sense,
        "guide_dna": extra.get("ordered_guide", ""),
        "passenger_dna": extra.get("passenger", ""),
        "guide_overhang": guide_overhang,
        "passenger_overhang": passenger_overhang,
        "duplex_format": duplex_format,
        "chemistry": hit.chemistry or extra.get("chemistry") or "unknown",
        "native_score": hit.score,
        "native_score_kind": hit.score_kind,
        "vote_eligible": hit.vote_eligible,
        "extra": extra,
        "tool_version": run.version,
        "raw_output_hashes": raw_hashes,
    }


def _normalized_hit_record(item: PlacedHit, run: ToolRun,
                           raw_hashes: dict[str, str]) -> dict:
    core = item.core
    return {
        "tool": item.tool,
        "role": run.role,
        "transcript_id": core.transcript_id,
        "start": core.start,
        "end": core.end,
        "length": core.length,
        "sense_dna": core.sense,
        "guide_dna": core.guide,
        "passenger_dna": core.passenger,
        "guide_overhang": core.guide_overhang,
        "passenger_overhang": core.passenger_overhang,
        "chemistry": core.chemistry,
        "molecule_id": core.molecule_id,
        "native_score": item.score,
        "native_score_kind": item.extra.get("score_kind", "tool_score"),
        "vote_eligible": item.vote_eligible,
        "extra": dict(item.extra),
        "tool_version": run.version,
        "raw_output_hashes": raw_hashes,
    }


def _run_id(config, input_sha: str, transcript_id: str) -> str:
    """Content-derived id so run.json is byte-reproducible across regenerations.

    A timestamp/uuid id would make every regeneration change run.json (and every
    selection.json that records its hash), breaking freeze --verify for no semantic
    reason. Two identical runs therefore share a run id and work directory, so the
    digest must cover every field that can change the outputs (filters, transcriptome,
    specificity) -- anything affecting results that is omitted would collide.

    ``transcript_id`` (the FASTA header) is included because tool output filenames
    derive from the record id: same sequence, different header, same work dir would
    otherwise let a second run re-glob the first run's files.
    """
    fields = {k: v for k, v in asdict(config).items() if k not in ('out_dir', 'fasta', 'jobs')}
    fields['input'] = input_sha
    fields['transcript_id'] = transcript_id
    manifest = Path(__file__).resolve().parents[2] / "benchmarks" / "generator_manifest.json"
    if manifest.is_file():
        fields['generator_manifest_sha256'] = sha256_file(manifest)
    return hashlib.sha256(json.dumps(fields, sort_keys=True, default=str).encode()).hexdigest()[:16]


@dataclass
class Config:
    fasta: str
    out_dir: str
    tools: tuple[str, ...] = DEFAULT_TOOLS
    lengths: tuple[int, ...] = (19, 20, 21)
    min_votes: int | None = None
    top_k: int | None = None
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
    profile: str = 'tools'
    filter_policy: str = 'soft'
    scorer_manifest: str = ''
    screen_manifest: str = ''
    target_gene: str = ''
    bowtie: str = 'bowtie'
    strategy: str = 'auto'
    evidence_run: str = ''

    def resolve(self) -> "Config":
        if self.min_votes is None:
            self.min_votes = 2 if self.profile == 'legacy-v2' else 1
        if self.top_k is None:
            self.top_k = 10 if self.profile == 'legacy-v2' else 5
        self.fasta = str(Path(self.fasta).expanduser().resolve())
        self.out_dir = str(Path(self.out_dir).expanduser().resolve())
        for name in ('scorer_manifest', 'screen_manifest', 'evidence_run'):
            if getattr(self, name):
                setattr(self, name, str(Path(getattr(self, name)).expanduser().resolve()))
        if self.transcriptome:
            self.transcriptome = str(
                Path(self.transcriptome).expanduser().resolve())
        return self

    def validate(self) -> None:
        if self.profile not in ('tools', 'legacy-v2'):
            raise ValueError('unknown profile')
        if self.scorer_manifest or self.screen_manifest or self.evidence_run:
            raise ValueError('tool profiles do not accept v3 model/reference/replay options; use select')
        if self.filter_policy not in ('soft', 'hard', 'off'):
            raise ValueError('invalid filter policy')
        if self.strategy not in ('auto', 'quality', 'native', 'model'):
            raise ValueError('invalid strategy')
        if not self.tools:
            raise ValueError("at least one tool is required")
        if self.min_votes < 1:
            raise ValueError("min_votes is below the profile minimum")
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

    transcript_id, sequence = read_transcript(config.fasta)
    input_sha = sha256(sequence)
    run_id = _run_id(config, input_sha, transcript_id)
    # One lock per out_dir. Identical runs share a work dir (fixed run id); *different*
    # runs into one out_dir also write out_dir/run.json|report.html, so they must not
    # interleave either. try/finally so a raised run still releases the lock even if the
    # traceback (and thus the frame) is retained by a caller.
    lock_handle = (out_dir / ".run.lock").open("a")
    fcntl.flock(lock_handle, fcntl.LOCK_EX)
    try:
        return _run_locked(config, out_dir, run_id, transcript_id, sequence, input_sha)
    finally:
        fcntl.flock(lock_handle, fcntl.LOCK_UN)
        lock_handle.close()


def _run_locked(config, out_dir, run_id, transcript_id, sequence, input_sha) -> dict:
    work_dir = out_dir / "work" / run_id
    # Clear the whole run dir, not just configured tools: an unavailable tool would
    # otherwise leave a stale subdir that the raw-output glob still hashes, making
    # run.json history-dependent even at a fixed run id.
    shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    transcriptome_sha = None
    if config.transcriptome:
        transcriptome_path = Path(config.transcriptome)
        if transcriptome_path.is_file():
            transcriptome_sha = sha256_file(transcriptome_path)

    manifest_path = Path(__file__).resolve().parents[2] / "benchmarks" / "generator_manifest.json"
    generator_manifest_sha = sha256_file(manifest_path) if manifest_path.is_file() else None

    tools, runs, tool_env = _select_tools(config)
    if not tools:
        detail = "; ".join(f"{r.tool}: {r.detail}" for r in runs) or "no tools requested"
        raise RuntimeError("no external tool is usable: " + detail)

    def execute(tool) -> ToolRun:
        # The whole run dir was cleared before fan-out, so a tool can never re-glob a
        # previous run's *_csv / off_targets_summary.tsv files.
        target = work_dir / tool.name
        target.mkdir(parents=True, exist_ok=True)
        try:
            run = tool.run(Path(config.fasta), target, transcript_id)
            # Accept older third-party adapters that returned only a hit list;
            # they remain usable but receive an explicit unpinned version.
            if isinstance(run, list):
                detail = "legacy adapter returned hits without provenance"
                run = (tool.wrap("ok", detail, run) if hasattr(tool, "wrap") else
                       ToolRun(tool=tool.name, status="ok", detail=detail,
                               hits=run, version=tool.version(), role=tool.role,
                               acceptance=getattr(tool, "acceptance", "")))
            # Adapters may provide a narrower hash set (for example only the
            # parsed table); fill it from the complete work directory so every
            # native record has reproducibility evidence.
            if not run.raw_output_hashes:
                run.raw_output_hashes = (tool.hash_outputs(target)
                                         if hasattr(tool, "hash_outputs") else {})
            return run
        except Exception as exc:  # noqa: BLE001 - one tool must not kill the run
            detail = f"{type(exc).__name__}: {exc}"
            return (tool.wrap("failed", detail) if hasattr(tool, "wrap") else
                    ToolRun(tool=tool.name, status="failed", detail=detail,
                            version=tool.version(), role=getattr(tool, "role", "generator"),
                            acceptance=getattr(tool, "acceptance", "")))

    with ThreadPoolExecutor(max_workers=max(1, min(config.jobs, len(tools)))) as pool:
        runs.extend(pool.map(execute, tools))
    runs.sort(key=lambda r: (r.tool not in [t.name for t in tools], r.tool))

    placed_generator: list[tuple[str, object, float | None]] = []
    evaluator_hits: list[PlacedHit] = []
    normalization = {}
    normalized_hits = []
    native_candidates: dict[str, list[dict]] = {}
    evaluator_evidence: dict[str, list[dict]] = {}
    for tool_run in runs:
        raw_hashes = dict(tool_run.raw_output_hashes)
        if tool_run.hits:
            records = [
                _native_hit_record(hit, tool_run, raw_hashes)
                for hit in tool_run.hits
            ]
            (native_candidates if tool_run.role == "generator" else evaluator_evidence)[tool_run.tool] = records
        if not tool_run.ok or not tool_run.hits:
            continue
        placed_hits, stats = normalize_placed(
            tool_run.hits, transcript_id, sequence, config.lengths)
        normalization[tool_run.tool] = stats
        normalized_hits.extend(
            _normalized_hit_record(item, tool_run, raw_hashes)
            for item in placed_hits)
        if tool_run.role == "generator":
            for item in placed_hits:
                if item.vote_eligible:
                    placed_generator.append((tool_run.tool, item.core, item.score))
        else:
            evaluator_hits.extend(placed_hits)

    generator_formats = {
        tool: sorted({row["duplex_format"] for row in rows})
        for tool, rows in native_candidates.items()
        if any(r.get("role") == "generator" for r in rows)
    }
    format_sets = [set(values) for values in generator_formats.values() if values]
    shared_generator_formats = sorted(set.intersection(*format_sets)) if format_sets else []
    format_compatibility = {
        "generator_formats": generator_formats,
        "shared_formats": shared_generator_formats,
        "status": ("compatible" if len(format_sets) < 2 or shared_generator_formats
                    else "mismatch_requires_site_or_format_specific_comparison"),
    }

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

    # Always compute a separately labelled site-overlap diagnostic.  It never
    # changes exact-molecule candidates or invents a shared molecule, but it
    # makes format/offset agreement visible when tools emit different duplexes.
    site_diagnostic = bounded_cluster(placed_generator, config.offset_tolerance)
    site_agreement = {
        "sites": len(site_diagnostic),
        "multi_tool_sites": sum(1 for c in site_diagnostic if c.votes > 1),
        "tool_sets": [sorted(c.tools) for c in site_diagnostic if c.votes > 1],
        "note": "diagnostic site overlap; not exact-molecule consensus",
    }

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
        "exact_molecule_sites": sum(1 for c in bucket(placed_generator)),
        "site_diagnostic_multi_tool": site_agreement["multi_tool_sites"],
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

    # This is intentionally a separate view from ``native_candidates``: it is
    # the union after placement/clustering and may contain evidence from more
    # than one generator.  No row here is presented as a fallback emission by a
    # tool that did not produce it.
    candidate_union = [
        {
            "key": c.core.key,
            "molecule_id": c.core.molecule_id,
            "transcript_id": c.core.transcript_id,
            "start": c.core.start,
            "end": c.core.end,
            "length": c.core.length,
            "sense_dna": c.core.sense,
            "guide_dna": c.core.guide,
            "passenger_dna": c.core.passenger,
            "guide_overhang": c.core.guide_overhang,
            "passenger_overhang": c.core.passenger_overhang,
            "chemistry": c.core.chemistry,
            "source_tools": c.tools,
            "votes": c.votes,
            "score": c.score,
        }
        for c in candidates
    ]

    result = {
        "schema_version": SCHEMA_VERSION,
        "implementation_hashes": {str(p.relative_to(Path(__file__).parent)): sha256_file(p)
                                  for p in sorted(Path(__file__).parent.rglob("*.py"))},
        "run_id": run_id,
        "tool_version": __version__,
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
            "generator_manifest_sha256": generator_manifest_sha,
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
                "raw_output_hashes": r.raw_output_hashes,
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
        "native_candidates": native_candidates,
        "evaluator_evidence": evaluator_evidence,
        "format_compatibility": format_compatibility,
        "site_agreement": site_agreement,
        "candidate_union": candidate_union,
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
