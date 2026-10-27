"""Output writers: candidate table, order-ready FASTA, run record, HTML report.

Every number in the HTML comes from ``run.json``; nothing is typed by hand.
Coordinates are 1-based and inclusive on the input transcript, and the nucleic
acid alphabet is stated per column because ordering a duplex needs RNA while
aligning back needs DNA.
"""
from __future__ import annotations

import csv
import html
import json
from pathlib import Path

from . import filters
from .types import Candidate

TSV_COLUMNS = [
    "rank", "consensus_tier", "votes", "tools", "transcript_id",
    "target_start", "target_end", "core_length",
    "target_sequence_dna", "guide_sequence_dna",
    "guide_sequence_rna", "passenger_sequence_rna",
    "gc_fraction", "seed_gc", "hard_flags", "soft_flags", "passed", "tool_scores",
]


def _dna_to_rna(sequence: str) -> str:
    return sequence.upper().replace("T", "U")


def candidate_row(candidate: Candidate, acting: int) -> dict[str, str]:
    core = candidate.core
    tier = "unanimous" if candidate.votes == acting else f"{candidate.votes}-of-{acting}"
    return {
        "rank": str(candidate.rank),
        "consensus_tier": tier,
        "votes": str(candidate.votes),
        "tools": "+".join(candidate.tools),
        "transcript_id": core.transcript_id,
        "target_start": str(core.start),
        "target_end": str(core.end),
        "core_length": str(core.length),
        "target_sequence_dna": core.sense,
        "guide_sequence_dna": core.guide,
        "guide_sequence_rna": _dna_to_rna(core.guide),
        "passenger_sequence_rna": _dna_to_rna(core.sense),
        "gc_fraction": f"{filters.gc_fraction(core.sense):.4f}",
        "seed_gc": f"{filters.seed_gc(core):.4f}",
        "hard_flags": ",".join(candidate.filter_flags),
        "soft_flags": ",".join(candidate.soft_flags),
        "passed": "true" if candidate.passed else "false",
        "tool_scores": ";".join(f"{tool}={score:.4f}"
                                for tool, score in sorted(candidate.tool_scores.items())),
    }


def write_candidates(result: dict, out_dir: Path) -> Path:
    path = out_dir / "candidates.tsv"
    acting = len(result["config"]["acting_tools"])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TSV_COLUMNS, delimiter="\t")
        writer.writeheader()
        for candidate in result["candidates"]:
            writer.writerow(candidate_row(candidate, acting))
    return path


def write_fasta(result: dict, out_dir: Path) -> Path:
    path = out_dir / "top_guides.fa"
    tier_by_id = {c.core.key: c for c in result["candidates"]}
    with path.open("w") as handle:
        for rank, candidate in enumerate(result["top"], start=1):
            core = candidate.core
            header = (f"rank={rank} tier={candidate.votes}/{len(result['config']['acting_tools'])} "
                      f"site={core.transcript_id}:{core.start}-{core.end} "
                      f"gc={filters.gc_fraction(core.sense):.3f} "
                      f"tools={','.join(candidate.tools)}")
            handle.write(f">{header}\n{_dna_to_rna(core.guide)}\n")
        for candidate in tier_by_id.values():
            if candidate in result["top"]:
                continue
            core = candidate.core
            handle.write(f">site={core.transcript_id}:{core.start}-{core.end} "
                         f"tier={candidate.votes}/"
                         f"{len(result['config']['acting_tools'])}\n{_dna_to_rna(core.guide)}\n")
    return path


def write_run_record(result: dict, out_dir: Path) -> Path:
    path = out_dir / "run.json"
    payload = dict(result)
    payload["candidates"] = [
        {"rank": c.rank, "votes": c.votes, "tools": c.tools,
         "transcript_id": c.core.transcript_id, "start": c.core.start,
         "end": c.core.end, "length": c.core.length,
         "sense_dna": c.core.sense, "guide_dna": c.core.guide,
         "hard_flags": c.filter_flags, "soft_flags": c.soft_flags,
         "tool_scores": c.tool_scores}
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
.ok{color:#1d6b3a}.bad{color:#a4301f}
.note{background:#f7faf9;border-left:3px solid var(--teal);padding:11px 14px;font-size:12.5px;margin:14px 0}
footer{margin-top:30px;padding-top:14px;border-top:1px solid var(--line);font-size:12px;color:var(--muted)}
@media print{body{background:#fff}.page{box-shadow:none;max-width:none;padding:0}}
"""


def write_html(result: dict, out_dir: Path) -> Path:
    acting = result["config"]["acting_tools"]
    counts = result["counts"]

    tool_rows = "".join(
        f"<tr><td><code>{html.escape(t['tool'])}</code></td>"
        f"<td class='{'ok' if t['status'] == 'ok' else 'bad'}'>{html.escape(t['status'])}</td>"
        f"<td class='n'>{t['hits']}</td><td>{html.escape(t['detail'])}</td></tr>"
        for t in result["tools"])

    norm_rows = "".join(
        f"<tr><td><code>{html.escape(name)}</code></td><td class='n'>{s['total']}</td>"
        f"<td class='n'>{s['placed']}</td><td class='n'>{s['duplicates']}</td>"
        f"<td>{html.escape(json.dumps(s['rejected'], ensure_ascii=False))}</td></tr>"
        for name, s in result["normalization"].items())

    top_rows = "".join(
        f"<tr><td class='n'>{rank}</td>"
        f"<td><span class='pill {'ok' if c.votes == len(acting) else 'warn'}'>"
        f"{c.votes}/{len(acting)}</span></td>"
        f"<td>{html.escape(c.core.transcript_id)}:{c.core.start}-{c.core.end}</td>"
        f"<td><code>{html.escape(c.core.guide)}</code></td>"
        f"<td class='n'>{filters.gc_fraction(c.core.sense):.3f}</td>"
        f"<td>{html.escape('+'.join(c.tools))}</td>"
        f"<td>{html.escape(','.join(c.soft_flags) or '-')}</td></tr>"
        for rank, c in enumerate(result["top"], start=1))

    blockers = "".join(f"<li>{html.escape(x)}</li>" for x in result["limitations"])

    body = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>sirna-consensus — {html.escape(result['transcript_id'])}</title>
<style>{CSS}</style></head><body><div class="page">
<h1>Consensus de novo siRNA candidates</h1>
<p class="sub">{html.escape(result['transcript_id'])} · {result['transcript_length']} nt ·
input SHA-256 <code>{result['input_sha256'][:16]}…</code> ·
generated {html.escape(result['created_utc'])}</p>
<div class="pills">
  <span class="pill info">tools acting: {len(acting)}</span>
  <span class="pill warn">status {html.escape(result['status'])}</span>
  <span class="pill ok">consensus {counts['consensus_unanimous']} unanimous</span>
  <span class="pill bad">candidates rejected by filters: {counts['clusters'] - counts['passing']}</span>
</div>

<h2>How to read this</h2>
<div class="note">Agreement is counted at the <b>site</b> level, not the string level:
the tools emit different lengths and offsets for the same duplex, so candidates
within 2 nt and sharing all but 2 bases are treated as one site. A candidate's
vote count is the number of <i>distinct tools</i> that proposed it. Sequence
filters (GC window, homopolymer, complexity, guide 5′ end, internal palindrome)
are applied <b>after</b> voting, so no single tool decides what survives.</div>

<h2>External tools</h2>
<table><thead><tr><th>tool</th><th>status</th><th>hits</th><th>detail</th></tr></thead>
<tbody>{tool_rows}</tbody></table>

<h2>Normalisation onto the transcript</h2>
<table><thead><tr><th>tool</th><th>received</th><th>placed</th><th>duplicates</th><th>rejected</th></tr></thead>
<tbody>{norm_rows}</tbody></table>

<h2>Top {len(result['top'])} consensus candidates</h2>
<table><thead><tr><th>#</th><th>votes</th><th>site</th><th>guide 5′→3′ (RNA)</th>
<th>GC</th><th>tools</th><th>soft flags</th></tr></thead><tbody>{top_rows}</tbody></table>

<h2>Limitations</h2><ul>{blockers}</ul>
<footer>Generated by sirna-consensus {html.escape(result['tool_version'])} ·
full machine-readable record in <code>run.json</code> ·
candidates in <code>candidates.tsv</code> · order-ready guides in <code>top_guides.fa</code></footer>
</div></body></html>
"""
    path = out_dir / "report.html"
    path.write_text(body)
    return path


def write_all(result: dict, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    return [write_candidates(result, out_dir), write_fasta(result, out_dir),
            write_run_record(result, out_dir), write_html(result, out_dir)]
