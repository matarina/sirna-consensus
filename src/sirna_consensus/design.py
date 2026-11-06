"""Schema-v3 exhaustive duplex design; external opinions never define the universe."""

from __future__ import annotations

import csv
import html
import json
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from pathlib import Path

from . import __version__
from .fasta import sha256_file
from .normalize import normalize_placed
from .selection import attach_evidence, enumerate_duplexes, load_candidates, panels, rank_candidates

DEFAULT_POLICY_PATH = Path(__file__).resolve().parents[2] / 'models' / 'default_policy.json'


def run_design(config, out, run_id, tid, sequence, input_sha):
    from .pipeline import _select_tools

    work = out / "work" / run_id
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    tools, runs, _env = _select_tools(config)

    def execute(tool):
        target = work / tool.name
        target.mkdir()
        try:
            return tool.run(Path(config.fasta), target, tid)
        except Exception as exc:
            return tool.wrap("failed", f"{type(exc).__name__}: {exc}")

    if tools:
        with ThreadPoolExecutor(max_workers=min(config.jobs, len(tools))) as pool:
            runs.extend(pool.map(execute, tools))
    rows = enumerate_duplexes(tid, sequence)
    native, normalized, stats = [], [], {}
    for run in sorted(runs, key=lambda r: r.tool):
        if not run.ok:
            continue
        placed, stat = normalize_placed(run.hits, tid, sequence, config.lengths)
        stats[run.tool] = asdict(stat)
        for hit in placed:
            normalized.append({**asdict(hit), "role": run.role})
            if run.role == "generator" and hit.vote_eligible:
                native.append(
                    {
                        "core": hit.core,
                        "tools": [hit.tool],
                        "tool_scores": {hit.tool: hit.score} if hit.score is not None else {},
                    }
                )
    evidence_sha = None
    if config.evidence_run:
        replayed, source = load_candidates(Path(config.evidence_run), evidence_only=True)
        if source["transcript_id"] != tid or source["input_sha256"] != input_sha:
            raise ValueError("replayed evidence target mismatch")
        for hit in replayed:
            c = hit["core"]
            if sequence[c.start - 1 : c.end] != c.sense:
                raise ValueError("replayed evidence coordinate mismatch")
        native.extend(replayed)
        evidence_sha = sha256_file(Path(config.evidence_run))
    attach_evidence(rows, native)
    scores, scorer = {}, {"status": "not_requested"}
    manifest = config.scorer_manifest
    policy_path = DEFAULT_POLICY_PATH
    if not manifest and config.strategy == "auto" and policy_path.is_file():
        policy = json.loads(policy_path.read_text())
        if policy.get("promoted"):
            manifest = str(policy_path.parent / policy["scorer_manifest"])
            if sha256_file(Path(manifest)) != policy["scorer_manifest_sha256"]:
                raise ValueError("default model policy checksum mismatch")
    strategy = "quality" if config.strategy == "auto" else config.strategy
    if manifest:
        from .scoring import score_candidates

        scores, scorer = score_candidates(rows, sequence, Path(manifest).resolve(), work / "scorer")
        if len(scores) != len(rows):
            raise ValueError('the run scorer must support every RNA 19+2 candidate; use select for mixed-format diagnostics')
        if config.strategy == "auto":
            strategy = "model"
    if strategy == "model" and not manifest:
        raise ValueError("model strategy requires a scorer manifest")
    ranked = rank_candidates(rows, strategy, config.filter_policy, scores, config.filter_config)
    for row in ranked:
        row["screen"] = {"status": "unverified", "reason": "no reference screen requested"}
        if len(row["tools"]) < config.min_votes:
            row["eligible"] = False
            row["hard_flags"].append("below_exact_molecule_vote_threshold")
    screen_meta = {"status": "not_requested"}
    if config.screen_manifest:
        from .screening import screen

        try:
            screen_meta = screen(
                ranked,
                Path(config.screen_manifest),
                config.target_gene,
                config.bowtie,
                work / "screen",
            )
        except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as exc:
            for row in ranked:
                row["screen"] = {"status": "unverified", "reason": str(exc)}
            screen_meta = {"status": "failed", "reason": str(exc)}
    if config.strict_specificity:
        for row in ranked:
            if row["screen"]["status"] != "pass":
                row["eligible"] = False
                row["hard_flags"].append("specificity_unverified_or_rejected")
    selected = {
        "research_top": panels(ranked, config.top_k),
        "diverse_top": panels(ranked, config.top_k, config.min_spacing),
        "screened_top": panels(ranked, config.top_k, screened=True),
        "screened_diverse_top": panels(ranked, config.top_k, config.min_spacing, screened=True),
    }

    def count(name):
        return sum(map(len, selected[name].values()))

    record = {
        "schema_version": 3,
        "kind": "exhaustive_duplex_design",
        "run_id": run_id,
        "tool_version": __version__,
        "transcript_id": tid,
        "transcript_length": len(sequence),
        "input_sha256": input_sha,
        "status": "research_only",
        "config": {**asdict(config), "source": "exhaustive target-derived RNA 19+2 duplexes"},
        "provenance": {
            "scorer": scorer,
            "screen": screen_meta,
            "evidence_run_sha256": evidence_sha,
        },
        "tools": [{**asdict(r), "hits": len(r.hits)} for r in sorted(runs, key=lambda r: r.tool)],
        "normalization": stats,
        "normalized_hits": normalized,
        "candidates": ranked,
        **selected,
        "default_panel": "diverse_top",
        "top": [
            r["key"]
            for r in ranked
            if r["key"] in {k for v in selected["diverse_top"].values() for k in v}
        ],
        "counts": {
            "enumerated": len(rows),
            "candidates": len(rows),
            "hits": sum(len(r.hits) for r in runs),
            "eligible": sum(r["eligible"] for r in ranked),
            "with_site_support": sum(bool(r["site_evidence"]) for r in ranked),
            "with_exact_molecule_support": sum(bool(r["tools"]) for r in ranked),
            "top": count("diverse_top"),
            "research_top": count("research_top"),
            "screened_top": count("screened_top"),
            "screened_diverse_top": count("screened_diverse_top"),
        },
        "panel_fill_rates": {
            name: {group: len(keys) / config.top_k for group, keys in groups.items()}
            for name, groups in selected.items()
        },
        "independent_validation_passed": False,
        "limitations": [
            "Scores are relative rankings, not knockdown probabilities.",
            "Unmeasured molecules have unknown activity.",
            "A screen pass means only the recorded computational criterion.",
            "Site support and unknown chemistry are not exact molecule votes.",
        ],
        "implementation_hashes": {
            str(p.relative_to(Path(__file__).parent)): sha256_file(p)
            for p in sorted(Path(__file__).parent.rglob("*.py"))
        },
        "raw_output_hashes": {
            str(p.relative_to(out)): sha256_file(p) for p in sorted(work.rglob("*")) if p.is_file()
        },
    }
    (out / "run.json").write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    write_exports(record, out)
    return record


def write_exports(record, out):
    rows = record["candidates"]
    fields = [
        "key",
        "start",
        "end",
        "guide_rna",
        "passenger_rna",
        "chemistry",
        "score",
        "ranking",
        "eligible",
        "screen_status",
        "exact_tools",
        "site_tools",
        "flags",
    ]
    with (out / "candidates.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow(
                dict(
                    key=r["key"],
                    start=r["start"],
                    end=r["start"] + r["length"] - 1,
                    guide_rna=r["guide_dna"].replace("T", "U"),
                    passenger_rna=r["passenger"].replace("T", "U"),
                    chemistry=r["chemistry"],
                    score=r["score"],
                    ranking=r["ranking"],
                    eligible=r["eligible"],
                    screen_status=r["screen"]["status"],
                    exact_tools=",".join(r["tools"]),
                    site_tools=",".join(sorted({h["tool"] for h in r["site_evidence"]})),
                    flags=",".join(r["hard_flags"] + r["soft_flags"]),
                )
            )
    by_key = {r["key"]: r for r in rows}
    for name in ("research_top", "diverse_top", "screened_top", "screened_diverse_top"):
        keys = [key for values in record[name].values() for key in values]
        with (out / (name + ".fa")).open("w") as handle:
            for key in keys:
                r = by_key[key]
                handle.write(
                    f">{key} guide chemistry={r['chemistry']}\n{r['guide_dna'].replace('T', 'U')}\n"
                    f">{key} passenger chemistry={r['chemistry']}\n{r['passenger'].replace('T', 'U')}\n"
                )
    shutil.copyfile(out / "diverse_top.fa", out / "top_guides.fa")
    (out / "all_candidates.fa").write_text(
        "".join(f">{r['key']}\n{r['guide_dna'].replace('T', 'U')}\n" for r in rows)
    )
    top = [by_key[k] for k in record["top"]]
    table = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(v))}</td>"
            for v in (
                r["start"],
                r["guide_dna"].replace("T", "U"),
                r["passenger"].replace("T", "U"),
                round(r["score"], 5),
                r["screen"]["status"],
            )
        )
        + "</tr>"
        for r in top
    )
    (out / "report.html").write_text(
        '<!doctype html><html lang="zh"><meta charset="utf-8">'
        "<title>siRNA de novo design</title><style>body{font:16px sans-serif;max-width:1200px;margin:40px auto}"
        "td,th{padding:9px;border:1px solid #ccc}table{border-collapse:collapse}</style>"
        f"<h1>{html.escape(record['transcript_id'])} · RNA 19+2</h1>"
        "<p>完整枚举 → 工具位点证据 → 排序 → 参考筛查 → 间隔 Top5</p>"
        "<p>默认显示研究候选；筛查状态为 unverified 时不表示已通过。评分不是敲低概率。</p>"
        f"<p>候选 {len(rows)}；默认面板 {len(top)}；筛查且满足间隔 {record['counts']['screened_diverse_top']}。</p>"
        "<table><tr><th>位置</th><th>Guide 5′→3′</th><th>Passenger 5′→3′</th><th>评分</th><th>筛查</th></tr>"
        + table
        + '</table><p><a href="candidates.tsv">全部候选</a> · <a href="run.json">完整审计记录</a></p></html>'
    )
