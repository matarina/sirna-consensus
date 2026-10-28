"""Output writers: candidate table, research-selected FASTA, run record, HTML report.

Every number in the HTML comes from ``run.json``; nothing is typed by hand.
Coordinates are 1-based and inclusive on the input transcript, and the nucleic
acid alphabet is stated per column because ordering a duplex needs RNA while
aligning back needs DNA.

The exported ``top_guides.fa`` contains only the selected Top-K. Every candidate,
including filter-rejected and singleton exploratory sites, is written to
``all_candidates.fa`` with explicit status flags.
"""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from . import filters
from .types import Candidate

TSV_COLUMNS = [
    "rank", "consensus_tier", "votes", "acting_tools", "tools",
    "evaluator_tools", "transcript_id", "target_start", "target_end",
    "core_length", "target_sequence_dna", "guide_sequence_dna",
    "guide_sequence_rna", "passenger_sequence_rna", "score",
    "vote_component", "tool_score_component", "quality_component",
    "gc_fraction", "seed_gc", "specificity_status", "off_target_min",
    "hard_flags", "soft_flags", "passed", "tool_scores", "molecule_id", "duplex_format",
    "guide_overhang", "passenger_overhang", "nearby_tools",
]


def _dna_to_rna(sequence: str) -> str:
    return sequence.upper().replace("T", "U")


def candidate_row(candidate: Candidate) -> dict[str, str]:
    core = candidate.core
    return {
        "rank": str(candidate.rank),
        "consensus_tier": candidate.consensus_tier,
        "votes": str(candidate.votes),
        "acting_tools": str(candidate.acting_tools),
        "tools": "+".join(candidate.tools),
        "evaluator_tools": "+".join(candidate.evaluator_tools),
        "transcript_id": core.transcript_id,
        "target_start": str(core.start),
        "target_end": str(core.end),
        "core_length": str(core.length),
        "target_sequence_dna": core.sense,
        "guide_sequence_dna": core.guide,
        "guide_sequence_rna": _dna_to_rna(core.guide),
        "passenger_sequence_rna": _dna_to_rna(core.passenger),
        "molecule_id": core.molecule_id,
        "duplex_format": core.duplex_format,
        "guide_overhang": core.guide_overhang,
        "passenger_overhang": core.passenger_overhang,
        "nearby_tools": "+".join(candidate.nearby_tools),
        "score": f"{candidate.score:.4f}",
        "vote_component": f"{candidate.score_parts.get('vote_component', 0.0):.4f}",
        "tool_score_component": f"{candidate.score_parts.get('tool_score_component', 0.0):.4f}",
        "quality_component": f"{candidate.score_parts.get('quality_component', 0.0):.4f}",
        "gc_fraction": f"{filters.gc_fraction(core.sense):.4f}",
        "seed_gc": f"{filters.seed_gc(core):.4f}",
        "specificity_status": candidate.specificity_status,
        "off_target_min": "" if candidate.off_target_min is None else str(candidate.off_target_min),
        "hard_flags": ",".join(candidate.filter_flags),
        "soft_flags": ",".join(candidate.soft_flags),
        "passed": "true" if candidate.passed else "false",
        "tool_scores": ";".join(
            f"{tool}={score:.4f}"
            for tool, score in sorted(candidate.tool_scores.items())),
    }


def write_candidates(result: dict, out_dir: Path) -> Path:
    path = out_dir / "candidates.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=TSV_COLUMNS, delimiter="\t",
            lineterminator="\n")
        writer.writeheader()
        for candidate in result["candidates"]:
            writer.writerow(candidate_row(candidate))
    return path


def _fasta_header(candidate: Candidate, rank: int, include_status: bool = False) -> str:
    core = candidate.core
    parts = [
        f"rank={rank}",
        f"tier={candidate.consensus_tier}",
        f"site={core.transcript_id}:{core.start}-{core.end}",
        f"gc={filters.gc_fraction(core.sense):.3f}",
        f"generators={','.join(candidate.tools)}",
        f"specificity={candidate.specificity_status}",
    ]
    if include_status:
        parts.append(f"passed={'true' if candidate.passed else 'false'}")
        parts.append(f"hard_flags={','.join(candidate.filter_flags) or '-'}")
        parts.append(f"evaluators={','.join(candidate.evaluator_tools) or '-'}")
    return " ".join(parts)


def write_fasta(result: dict, out_dir: Path) -> Path:
    """Order-ready file: only the selected Top-K, never rejected candidates."""
    path = out_dir / "top_guides.fa"
    with path.open("w") as handle:
        for rank, candidate in enumerate(result["top"], start=1):
            handle.write(f">{_fasta_header(candidate, rank)}\n"
                         f"{_dna_to_rna(candidate.core.guide)}\n")
    return path


def write_all_candidates_fasta(result: dict, out_dir: Path) -> Path:
    """Inspection-only FASTA: every candidate with its filter status."""
    path = out_dir / "all_candidates.fa"
    with path.open("w") as handle:
        for candidate in result["candidates"]:
            handle.write(
                f">{_fasta_header(candidate, candidate.rank, include_status=True)}\n"
                f"{_dna_to_rna(candidate.core.guide)}\n")
    return path


def write_run_record(result: dict, out_dir: Path) -> Path:
    path = out_dir / "run.json"
    payload = dict(result)
    payload["candidates"] = [
        {
            "key": c.core.key,
            "molecule_id": c.core.molecule_id,
            "duplex_format": c.core.duplex_format,
            "ordered_guide": c.core.ordered_guide,
            "passenger": c.core.passenger,
            "guide_overhang": c.core.guide_overhang,
            "passenger_overhang": c.core.passenger_overhang,
            "nearby_tools": c.nearby_tools,
            "score": c.score,
            "rank": c.rank,
            "votes": c.votes,
            "acting_tools": c.acting_tools,
            "tools": c.tools,
            "evaluator_tools": c.evaluator_tools,
            "transcript_id": c.core.transcript_id,
            "start": c.core.start,
            "end": c.core.end,
            "length": c.core.length,
            "sense_dna": c.core.sense,
            "guide_dna": c.core.guide,
            "hard_flags": c.filter_flags,
            "soft_flags": c.soft_flags,
            "specificity_status": c.specificity_status,
            "off_target_min": c.off_target_min,
            "score_parts": c.score_parts,
            "tool_scores": c.tool_scores,
        }
        for c in result["candidates"]
    ]
    payload["top"] = [c.core.key for c in result["top"]]
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False))
    return path


CSS = """
:root{--ink:#111820;--muted:#5c6b74;--line:#d5dee2;--navy:#123b4f;--teal:#0f7b74}
*{box-sizing:border-box}
body{margin:0;background:#eef1f2;color:var(--ink);
font:15px/1.65 'Charter',Georgia,'Noto Serif SC',serif}
.page{max-width:1100px;margin:0 auto;background:#fff;padding:38px 44px 52px;
box-shadow:0 2px 20px rgba(20,35,43,.09)}
h1{font-size:23px;margin:0 0 6px}
h2{font-size:16px;margin:30px 0 10px;border-bottom:1px solid var(--line);padding-bottom:5px}
code{font:12px ui-monospace,Menlo,Consolas,monospace;background:#eff3f4;padding:1px 4px;border-radius:3px}
.sub{color:#3b4b55;font-size:13.5px;margin:0 0 16px}
.pills{display:flex;gap:8px;flex-wrap:wrap;margin:0 0 22px}
.pill{padding:3px 10px;border-radius:3px;font:700 11px ui-sans-serif,system-ui}
.pill.ok{background:#e2f1e7;color:#1d6b3a}.pill.bad{background:#fbe6e2;color:#a4301f}
.pill.info{background:#e4eef5;color:#123b4f}.pill.warn{background:#f8ecd6;color:#8a5c12}
table{width:100%;border-collapse:collapse;font:12.5px ui-sans-serif,system-ui;margin:10px 0}
th{text-align:left;font-size:10.5px;letter-spacing:.05em;text-transform:uppercase;color:var(--muted);
padding:6px;border-bottom:1.4px solid var(--ink)}
td{padding:6px;border-bottom:1px solid #e9eef0}
td.n{text-align:right;font-variant-numeric:tabular-nums}
.ok{color:#1d6b3a}.bad{color:#a4301f}.warn{color:#8a5c12}
.note{background:#f7faf9;border-left:3px solid var(--teal);padding:11px 14px;font-size:12.5px;margin:14px 0}
.note.warnbox{background:#fdf6e8;border-left:3px solid #8a5c12}
footer{margin-top:30px;padding-top:14px;border-top:1px solid var(--line);font-size:12px;color:var(--muted)}
@media print{body{background:#fff}.page{box-shadow:none;max-width:none;padding:0}}
"""


def write_html(result: dict, out_dir: Path) -> Path:
    acting = result["config"]["acting_tools"]
    counts = result["counts"]
    config = result["config"]

    tool_rows = "".join(
        f"<tr><td><code>{html.escape(t['tool'])}</code></td>"
        f"<td>{html.escape(t.get('role', 'generator'))}</td>"
        f"<td class='{'ok' if t['status'] == 'ok' else 'bad'}'>"
        f"{html.escape(t['status'])}</td>"
        f"<td class='n'>{t['hits']}</td>"
        f"<td>{html.escape(t.get('screen_state') or '-')}</td>"
        f"<td>{html.escape(t['detail'])}</td></tr>"
        for t in result["tools"])

    norm_rows = "".join(
        f"<tr><td><code>{html.escape(name)}</code></td><td class='n'>{s['total']}</td>"
        f"<td class='n'>{s['placed']}</td><td class='n'>{s['duplicates']}</td>"
        f"<td>{html.escape(json.dumps(s['rejected'], ensure_ascii=False))}</td></tr>"
        for name, s in result["normalization"].items())

    top_rows = "".join(
        f"<tr><td class='n'>{rank}</td>"
        f"<td><span class='pill {'ok' if c.votes == len(acting) else 'warn'}'>"
        f"{html.escape(c.consensus_tier)}</span></td>"
        f"<td>{html.escape(c.core.transcript_id)}:{c.core.start}-{c.core.end}</td>"
        f"<td><code>{html.escape(c.core.guide)}</code></td>"
        f"<td class='n'>{filters.gc_fraction(c.core.sense):.3f}</td>"
        f"<td>{html.escape('+'.join(c.tools))}</td>"
        f"<td>{html.escape(c.specificity_status)}</td>"
        f"<td>{html.escape(','.join(c.soft_flags) or '-')}</td></tr>"
        for rank, c in enumerate(result["top"], start=1))

    blockers = "".join(f"<li>{html.escape(x)}</li>" for x in result["limitations"])

    empty_notice = ""
    if not result["top"]:
        empty_notice = (
            '<div class="note warnbox"><b>No candidate passed into the Top-K.</b> '
            f'{counts["consensus_candidates"]} consensus candidate(s) and '
            f'{counts["exploratory_singletons"]} exploratory singleton(s) were '
            f'assembled from {counts["sites"]} site(s); '
            f'{counts["rejected_by_filters"]} failed sequence/specificity filters. '
            f'Top-K also enforces <code>min_votes={config["min_votes"]}</code>.</div>')

    config_rows = "".join(
        f"<tr><td>{html.escape(key)}</td>"
        f"<td>{html.escape('-' if value is None else str(value))}</td></tr>"
        for key, value in config.items()
        if key not in {"source"})

    body = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>sirna-consensus — {html.escape(result['transcript_id'])}</title>
<style>{CSS}</style></head><body><div class="page">
<h1>Consensus de novo siRNA candidates</h1>
<p class="sub">{html.escape(result['transcript_id'])} · {result['transcript_length']} nt ·
input SHA-256 <code>{result['input_sha256'][:16]}…</code> ·
generated {html.escape(result['created_utc'])}</p>
<div class="pills">
  <span class="pill info">generator tools: {len(acting)}</span>
  <span class="pill warn">status {html.escape(result['status'])}</span>
  <span class="pill ok">{counts['consensus_unanimous']} unanimous generator sites</span>
  <span class="pill info">{counts['consensus_candidates']} consensus candidates</span>
  <span class="pill warn">{counts['exploratory_singletons']} exploratory singletons</span>
  <span class="pill bad">{counts['rejected_by_filters']} rejected by filters/specificity</span>
  <span class="pill warn">{counts['below_min_votes']} below the vote threshold</span>
</div>

<h2>How to read this</h2>
<div class="note">Agreement is counted at the <b>site</b> level, not the string level:
the tools emit different lengths and offsets for the same duplex, so candidates
within <code>offset_tolerance</code> and sharing all but 2 bases are grouped.
The default candidate mode preserves <b>exact molecules</b>. Optional legacy site grouping is <b>bounded</b> (an anchored cluster, not a transitive
chain). Votes are <b>distinct generator tools</b>; SIREN is an <b>evaluator</b> and
does not add a design vote. Sequence/specificity filters are applied <b>after</b>
voting. The research-selected file is <code>top_guides.fa</code>; every candidate,
including rejected ones, is in <code>all_candidates.fa</code>.</div>

<h2>Effective configuration</h2>
<table><thead><tr><th>setting</th><th>value</th></tr></thead>
<tbody>{config_rows}</tbody></table>

<h2>External tools</h2>
<table><thead><tr><th>tool</th><th>role</th><th>status</th><th>hits</th>
<th>screen state</th><th>detail</th></tr></thead>
<tbody>{tool_rows}</tbody></table>

<h2>Normalisation onto the transcript</h2>
<table><thead><tr><th>tool</th><th>received</th><th>placed</th><th>duplicates</th><th>rejected</th></tr></thead>
<tbody>{norm_rows}</tbody></table>

{empty_notice}
<h2>Top {len(result['top'])} consensus candidates</h2>
<table><thead><tr><th>#</th><th>tier</th><th>site</th><th>guide 5′→3′ (RNA)</th>
<th>GC</th><th>generators</th><th>specificity</th><th>soft flags</th></tr></thead>
<tbody>{top_rows}</tbody></table>

<h2>Limitations</h2><ul>{blockers}</ul>
<footer>Generated by sirna-consensus {html.escape(result['tool_version'])} ·
full machine-readable record in <code>run.json</code> ·
candidates in <code>candidates.tsv</code> ·
research-selected guides in <code>top_guides.fa</code> ·
inspection-only candidates in <code>all_candidates.fa</code></footer>
</div></body></html>
"""
    path = out_dir / "report.html"
    path.write_text(body)
    return path


def write_all(result: dict, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    return [
        write_candidates(result, out_dir),
        write_fasta(result, out_dir),
        write_all_candidates_fasta(result, out_dir),
        write_run_record(result, out_dir),
        write_html(result, out_dir),
    ]
