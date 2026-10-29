#!/usr/bin/env python3
"""生成中文论文风格的真实运行报告 (report_zh.html)。

结构: 摘要 → 背景 → 方法 → 结果 → 图 → 讨论 → 结论 → 数据与复现。
所有数字与图都从本目录的 run.json (schema v2) 直接读出, 不做重新计算。
图以内联 SVG 绘制, 因此报告是单文件、离线可读。

用法:  python make_report_zh.py [run.json] [-o report_zh.html]
"""
from __future__ import annotations

import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent

TOOL_ROLE = {"generator": "生成器 (投设计票)", "evaluator": "评估器 (只提供证据)"}
STATUS_TEXT = {"ok": "成功", "failed": "失败", "unavailable": "不可用"}
LIMITS_ZH = (
    ("Consensus votes are generator-tool votes",
     "共识票是<b>生成器工具</b>的票。SIREN 是评估器, 其枚举命中与 input-only 回退都不算设计票。"),
    ("Specificity is separate evidence",
     "特异性是<b>独立证据</b>。转录组支持且脱靶数超阈才是硬性不合格; input-only 或失败的筛查"
     "都不能当作全转录组结论。"),
    ("Site-level consensus uses bounded offset clustering",
     "位点共识默认采用有界偏移聚类 (锚定簇, 非传递链); exact 与 fuzzy 是显式选项。"),
    ("Top-K enforces min_votes",
     "Top-K 强制满足投票阈值。单例位点仍留在 <code>candidates.tsv</code> 与 "
     "<code>all_candidates.fa</code> 中供查看, 但不进入订单集。"),
    ("A tool that returns only its own Top-N",
     "只返回自身 Top-N 的工具, 与「过滤掉了该位点」的工具无法区分。"),
    ("No wet-lab confirmation",
     "无湿实验确认; 不同工具的打分不可直接比较。"),
    ("Off-target scope is whatever each tool was given",
     "脱靶覆盖范围取决于各工具实际拿到的序列, 除非提供了共享转录组。"),
)
FLAG_TEXT = {
    "guide_5p_not_AU": "引导链 5′ 端非 A/U",
    "sense_3p_not_AU": "sense 链 3′ 端非 A/U",
    "homopolymer": "存在同聚物",
    "gc_outside_window": "GC 含量超窗",
    "seed_gc_low": "种子区 GC 偏低",
    "specificity_offtarget": "特异性脱靶",
    "specificity_unverified": "特异性未验证",
}

CSS = """
:root{--ink:#111820;--muted:#5c6b74;--line:#d5dee2;--navy:#123b4f;--teal:#0f7b74;
--red:#a4301f;--green:#1d6b3a;--amber:#8a5c12}
*{box-sizing:border-box}
body{margin:0;background:#eef1f2;color:var(--ink);
font:15px/1.8 'Noto Serif SC','Songti SC','Charter',Georgia,serif}
.page{max-width:960px;margin:0 auto;background:#fff;padding:46px 56px 60px;
box-shadow:0 2px 20px rgba(20,35,43,.09)}
h1{font-size:22px;line-height:1.5;margin:0 0 8px}
h2{font-size:16px;margin:34px 0 10px;color:var(--navy);
border-bottom:1px solid var(--line);padding-bottom:5px}
h3{font-size:14px;margin:22px 0 8px;color:var(--navy)}
h3 .num,h2 .num{color:var(--teal)}
p{margin:10px 0;text-align:justify}
code{font:12px ui-monospace,Menlo,Consolas,monospace;background:#eff3f4;padding:1px 4px;border-radius:3px}
.byline{color:#3b4b55;font-size:13px;margin:0 0 4px}
.meta{color:var(--muted);font-size:12.5px;margin:0 0 20px}
.abstract{background:#f7faf9;border:1px solid var(--line);border-left:3px solid var(--teal);
padding:14px 18px;font-size:13.5px;margin:18px 0}
.abstract h2{border:0;margin:0 0 6px;font-size:13px;letter-spacing:.1em;color:var(--teal)}
.kw{font-size:12.5px;color:#3b4b55;margin:10px 0 0}
table{width:100%;border-collapse:collapse;font:12.5px ui-sans-serif,system-ui;margin:12px 0}
caption{text-align:left;font-size:12px;color:var(--muted);padding:0 0 6px}
th{text-align:left;font-size:10.5px;letter-spacing:.05em;color:var(--muted);
padding:6px;border-top:1.4px solid var(--ink);border-bottom:1px solid var(--ink)}
td{padding:6px;border-bottom:1px solid #e9eef0;vertical-align:top}
tbody tr:last-child td{border-bottom:1.4px solid var(--ink)}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.seq{font:12px ui-monospace,Menlo,Consolas,monospace;white-space:nowrap}
.ok{color:var(--green)}.bad{color:var(--red)}.warn{color:var(--amber)}
.note{background:#f7faf9;border-left:3px solid var(--teal);padding:11px 14px;font-size:12.5px;margin:14px 0}
.note.warnbox{background:#fdf6e8;border-left:3px solid var(--amber)}
figure{margin:20px 0;padding:0}
figcaption{font-size:12px;color:#3b4b55;margin-top:8px;line-height:1.6}
figcaption b{color:var(--ink)}
svg{max-width:100%;height:auto;display:block}
ol.refs,ul.plain{font-size:12.5px;padding-left:22px}
ol.refs li,ul.plain li{margin:7px 0}
footer{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);font-size:12px;color:var(--muted)}
@media print{body{background:#fff}.page{box-shadow:none;max-width:none;padding:0}}
"""


def e(v) -> str:
    return html.escape(str(v), quote=True)


def key_of(c: dict) -> str:
    return f"{c['transcript_id']}:{c['start']}:{c['length']}"


def flag_text(codes) -> str:
    return "、".join(FLAG_TEXT.get(c, c) for c in codes) or "—"


# ---------------------------------------------------------------- figures
def fig_funnel(counts: dict) -> str:
    steps = [
        ("工具原始命中", counts["hits"], False),
        ("归一化并定位到转录本", counts["placed"], False),
        ("合并后的位点", counts["sites"], False),
        ("达成共识的候选", counts["consensus_candidates"], False),
        ("通过序列过滤", counts["passing"], False),
        ("最终 Top-K 输出", counts["top"], False),
        ("低于投票阈值 (丢弃)", counts["below_min_votes"], True),
        ("被过滤器淘汰", counts["rejected_by_filters"], True),
    ]
    w, row, pad, gap = 640, 30, 190, 9
    h = len(steps) * (row + gap) + 20
    top = max(v for _, v, _ in steps)
    bars = []
    for i, (label, value, drop) in enumerate(steps):
        y = 14 + i * (row + gap)
        bw = max(3, value / top * (w - pad - 90))
        color = "#a4301f" if drop else "#123b4f"
        opacity = ' opacity="0.65"' if drop else ""
        bars.append(
            f'<text x="0" y="{y+20}" font-size="12" fill="#3b4b55">{e(label)}</text>'
            f'<rect x="{pad}" y="{y+4}" width="{bw:.1f}" height="{row-6}" rx="2" fill="{color}"'
            f'{opacity}/>'
            f'<text x="{pad+bw+8:.1f}" y="{y+20}" font-size="12" font-weight="700"'
            f' fill="{color}">{value:,}</text>')
    return (f'<svg viewBox="0 0 {w} {h}" role="img" '
            f'aria-label="运行流程漏斗">{"".join(bars)}</svg>')


def fig_map(run: dict, top_keys: set[str]) -> str:
    length = run["transcript_length"]
    cands = run["candidates"]
    w, h, pad, track = 900, 168, 72, 44
    def x(pos):
        return pad + (pos - 1) / max(1, length - 1) * (w - 2 * pad)

    def track_row(y, items, color, radius, label):
        out = [f'<text x="0" y="{y+4}" font-size="11.5" fill="#3b4b55">{e(label)}</text>']
        for c in items:
            out.append(f'<circle cx="{x(c["start"]):.1f}" cy="{y}" r="{radius}" '
                       f'fill="{color}" opacity="0.85"/>')
        return "".join(out)

    top10 = [c for c in cands if key_of(c) in top_keys]
    passing = [c for c in cands if not c["hard_flags"]] if "hard_flags" in cands[0] else cands
    bottom = track + 92
    ticks = "".join(
        f'<line x1="{x(p):.1f}" y1="{track}" x2="{x(p):.1f}" y2="{bottom+6}" '
        f'stroke="#c9d4d8" stroke-width="1"/>'
        f'<text x="{x(p):.1f}" y="{bottom+22}" font-size="10.5" fill="#7a878e" '
        f'text-anchor="middle">{p}</text>'
        for p in range(200, length, 200))
    return f'''<svg viewBox="0 0 {w} {h}" role="img" aria-label="候选在转录本上的分布">
{ticks}
<line x1="{pad}" y1="{track}" x2="{w-pad}" y2="{track}" stroke="#123b4f" stroke-width="2.5"/>
{track_row(track, cands, "#9db4c0", 3.4, f"全部共识候选 (n={len(cands)})")}
{track_row(track+46, passing, "#0f7b74", 4.5, f"通过序列过滤 (n={len(passing)})")}
{track_row(bottom, top10, "#a4301f", 6, f"最终 Top-{len(top10)} (n={len(top10)})")}
<text x="{pad}" y="{track-20}" font-size="11" fill="#7a878e">5′</text>
<text x="{w-pad}" y="{track-20}" font-size="11" fill="#7a878e" text-anchor="end">3′ ({length} nt)</text>
</svg>'''


def fig_flags(hard: dict, soft: dict) -> str:
    items = [(FLAG_TEXT.get(k, k) + " (硬)", v, "#a4301f") for k, v in hard.items()]
    items += [(FLAG_TEXT.get(k, k) + " (软)", v, "#8a5c12") for k, v in soft.items()]
    items.sort(key=lambda t: t[1])
    w, row, pad = 640, 26, 210
    h = len(items) * row + 16
    top = max(v for _, v, _ in items)
    bars = []
    for i, (label, v, color) in enumerate(items):
        y = 10 + i * row
        bw = max(3, v / top * (w - pad - 70))
        bars.append(
            f'<text x="0" y="{y+15}" font-size="11.5" fill="#3b4b55">{e(label)}</text>'
            f'<rect x="{pad}" y="{y+4}" width="{bw:.1f}" height="{row-9}" rx="2" fill="{color}"/>'
            f'<text x="{pad+bw+8:.1f}" y="{y+15}" font-size="11.5" font-weight="700" '
            f'fill="{color}">{v}</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="过滤命中分布">{"".join(bars)}</svg>'


# ---------------------------------------------------------------- report
def build(run: dict) -> str:
    c = run["counts"]
    cfg = run["config"]
    tid = run["transcript_id"]

    top_keys = set(run["top"])
    top = sorted((x for x in run["candidates"] if key_of(x) in top_keys),
                 key=lambda x: x["rank"])

    hard: dict[str, int] = {}
    soft: dict[str, int] = {}
    for cand in run["candidates"]:
        for f in cand["hard_flags"]:
            hard[f] = hard.get(f, 0) + 1
        for f in cand["soft_flags"]:
            soft[f] = soft.get(f, 0) + 1

    acting = cfg["acting_tools"]
    n_gen = len(acting)
    rows_tools = "".join(
        f"<tr><td><b>{e(t['tool'])}</b></td><td>{e(TOOL_ROLE.get(t['role'], t['role']))}</td>"
        f"<td class='{'ok' if t['status'] == 'ok' else 'bad'}'>"
        f"{e(STATUS_TEXT.get(t['status'], t['status']))}</td>"
        f"<td class='n'>{t['hits']}</td><td class='seq'>{e(t.get('screen_state') or '—')}</td></tr>"
        for t in run["tools"])

    rows_top = "".join(
        f"<tr><td class='n'>{x['rank']}</td><td class='n'>{x['start']}–{x['end']}</td>"
        f"<td class='n'>{x['length']}</td>"
        f"<td>{e(' + '.join(x['tools']))} ({x['votes']}/{n_gen})</td>"
        f"<td class='n'>{x['score_parts']['tool_score_component']:.3f}</td>"
        f"<td class='seq'>{e(x['guide_dna'])}</td>"
        f"<td>{flag_text(x['hard_flags'] + x['soft_flags'])}</td></tr>"
        for x in top)

    rows_norm = "".join(
        f"<tr><td>{e(name)}</td><td class='n'>{s['total']}</td><td class='n'>{s['placed']}</td>"
        f"<td class='n'>{s['duplicates']}</td>"
        f"<td>{e(json.dumps(s['rejected'], ensure_ascii=False) if s['rejected'] else '—')}</td></tr>"
        for name, s in sorted(run["normalization"].items()))

    def flag_rows(counter):
        if not counter:
            return "<tr><td colspan='2'>—</td></tr>"
        return "".join(f"<tr><td>{e(FLAG_TEXT.get(k, k))}</td><td class='n'>{v}</td></tr>"
                       for k, v in sorted(counter.items(), key=lambda kv: -kv[1]))

    limits = "".join(
        f"<li>{zh}<br><span class='seq'>({e(code[:70])}…)</span></li>"
        for code in run["limitations"]
        for zh in [next((t for k, t in LIMITS_ZH if code.startswith(k)), e(code))])

    created = run["created_utc"][:19].replace("T", " ") + " UTC"

    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>基于多工具位点共识的 {e(tid)} 去新 siRNA 设计: 真实运行报告</title>
<style>{CSS}</style></head><body><div class="page">

<h1>基于多工具位点共识的 {e(tid)} (AREG) 去新 siRNA 候选设计:<br>一次真实运行的完整报告</h1>
<p class="byline">sirna-consensus v{e(run['tool_version'])} · 运行 ID <code>{e(run['run_id'])}</code></p>
<p class="meta">生成于 {e(created)} · schema v{run['schema_version']} ·
输入 {e(tid)} ({run['transcript_length']} nt) · 输入 SHA-256 <code>{e(run['input_sha256'][:16])}…</code></p>

<div class="abstract">
<h2>摘要</h2>
<p><b>背景.</b> 去新 siRNA 设计缺少可用的近缘同源参照, 单一工具的排序往往不可复核。
本文报告以外部独立设计工具的<b>位点级共识</b>替代单一打分器的做法。</p>
<p><b>方法.</b> 以单个转录本 FASTA 为唯一输入, 并行调用两个生成器工具
(siRNAdesigner、siRNAforge) 与一个评估器工具 (SIREN), 将其输出归一化并定位到转录本,
按起始位点与长度归并为位点, 要求 ≥{cfg['min_votes']} 个独立工具一致方可成为候选,
再施加序列过滤器并按间距选取 Top-K。评估器只提供筛查证据, 不投设计票。</p>
<p><b>结果.</b> 三个工具共产出 <b>{c['hits']:,}</b> 条原始命中, 全部 {c['placed']:,} 条成功定位
(零拒绝、零重复), 归并为 <b>{c['sites']}</b> 个位点。其中 <b>{c['consensus_unanimous']}</b> 个位点
由两个生成器<b>完全一致</b>支持, 构成共识候选集; <b>{c['below_min_votes']}</b> 个位点因未达投票阈值被丢弃。
经序列过滤淘汰 {c['rejected_by_filters']} 条后余 {c['passing']} 条, 最终输出 <b>{c['top']}</b> 条订单候选
(最短位点间距 {cfg['min_spacing']} nt)。</p>
<p><b>结论.</b> 在无转录组输入时, 全部 {c['consensus_candidates']} 条候选的特异性状态均为
<code>input_only</code>, 即<b>尚未获得全转录组特异性证据</b>; 本结果是计算预测, 不代表敲低效力已被验证。</p>
<p class="kw"><b>关键词:</b> siRNA 设计; 工具共识; 脱靶评估; 可复现性; 转录本定位</p>
</div>

<h2><span class="num">1.</span> 背景</h2>
<p>siRNA 的活性与脱靶行为由引导链的序列、位置与热力学共同决定, 现有工具对这些因素的加权各不相同。
当目标基因缺乏同源参照 (即 "去新" 设定) 时, 无法用同一套标注来调参, 因此单一工具的排序难以独立复核。</p>
<p>本工作采取的策略是: 不训练新模型, 而是<b>用一个转录本 FASTA 驱动既有的开源工具</b>,
把它们的<b>独立意见放在同一坐标轴上</b>, 并以共识而非单一分数决定取舍。
关键设计约束有三条: (i) 输入默认就是 FASTA, 而不是基因名; (ii) 本包不做窗口枚举, 枚举交给各工具;
(iii) 共识发生在<b>位点</b>层面而不是字符串层面 —— 不同工具对同一条双链会给出不同的长度与偏移。</p>

<h2><span class="num">2.</span> 方法</h2>
<h3>2.1 输入与工具</h3>
<p>输入为 {e(tid)} 的转录本序列 ({run['transcript_length']} nt),
运行前记录其 SHA-256 以确保可复现。三个外部工具各自运行在独立子进程中 (互不共享依赖):
生成器 {e('、'.join(acting))} 负责提出设计票; 评估器 SIREN 只产出证据。
工具全部命中最多保留, 由本包统一解析。</p>

<h3>2.2 归一化与位点共识</h3>
<p>工具输出被转换为 1-based 闭合区间, 并把 <b>结合核心</b>、原始引导链/过客链与
3′ 悬垂拆成独立字段。随后按 <code>offset_tolerance={cfg['offset_tolerance']}</code> 归并位点;
投票单位是<b>不同的生成器工具</b>, 同一工具多个命中不重复计票。
保留 <code>min_votes={cfg['min_votes']}</code> 票以上的位点为候选; 本次同时关闭了单例
(<code>include_singletons={cfg['include_singletons']}</code>)。</p>

<h3>2.3 过滤与排序</h3>
<p>过滤在投票<b>之后</b>执行, 因此不会把一个位点从其他工具的一致意见中移除。
硬过滤含 GC 窗口 ({cfg['filter_config']['gc_min']}–{cfg['filter_config']['gc_max']})、
最长同聚物 ({cfg['filter_config']['max_homopolymer']}) 等; 软标记 (如 5′ 端非 A/U) 仅告警。
排序按位点级分数, 并以 ≥{cfg['min_spacing']} nt 的最小间距选取 Top-{cfg['top_k']}。</p>

<h3>2.4 特异性证据</h3>
<p>SIREN 的命中被单独归类。只有转录组支持的筛查 (<code>transcriptome_screened</code>)
才构成特异性结论; 本次未提供转录组, 因此全部候选停留在 <code>input_only</code>。
在位点层面, 评估器结果只能<b>支持或淘汰</b>一条已生成的候选, 不能凭空产生设计票。</p>

<h3>2.5 复现</h3>
<p>运行写入了全部归一化命中、原始输出哈希与实现文件哈希, 命令为:<br>
<code>sirna-consensus run --fasta examples/AREG/{e(tid)}.fa --out examples/AREG
--tools {e(','.join(cfg['tools']))} --min-votes {cfg['min_votes']} --top-k {cfg['top_k']}
--exclude-singletons</code></p>

<h2><span class="num">3.</span> 结果</h2>
<h3>3.1 工具执行</h3>
<table><caption>表 1. 外部工具执行情况。角色与命中数均摘自运行记录。</caption>
<thead><tr><th>工具</th><th>角色</th><th>状态</th><th class="n">命中数</th><th>筛查状态</th></tr></thead>
<tbody>{rows_tools}</tbody></table>
<p>三个工具全部执行成功。SIREN 的 {c['evaluator_hits']} 条命中中绝大多数是<b>评估器枚举</b>,
其筛查状态为 <code>input_only</code>, 即只在输入序列内部比较, <b>不是</b>全转录组结论。</p>

<h3>3.2 归一化与定位</h3>
<table><caption>表 2. 归一化与转录本定位。零拒绝说明工具输出与本包的解析完全兼容。</caption>
<thead><tr><th>工具</th><th class="n">原始命中</th><th class="n">成功定位</th><th class="n">重复</th>
<th>拒绝原因</th></tr></thead><tbody>{rows_norm}</tbody></table>

<h3>3.3 共识与输出</h3>
<p>归并后共有 {c['sites']} 个位点, 其中 <b>{c['consensus_unanimous']}</b> 个位点被两个生成器
完全一致地支持 (即投票 {n_gen}/{n_gen}), {c['below_min_votes']} 个位点未达阈值。
过滤淘汰 {c['rejected_by_filters']} 条后, {c['passing']} 条进入排序, 输出 Top-{len(top)}:</p>
<table><caption>表 3. 最终 Top-{len(top)} 候选。"全局排名" 是全部 {c['passing']} 条通过过滤候选中的位次,
因最短间距约束而不连续; 括号内为达成一致的生成器工具数与总数。原始打分为工具内归一化值,
跨工具不可直接比较。</caption>
<thead><tr><th class="n">全局排名</th><th class="n">位置</th><th class="n">长度</th><th>投票</th>
<th class="n">原始打分</th><th>引导链 (5′→3′, DNA)</th><th>标记</th></tr></thead>
<tbody>{rows_top}</tbody></table>

<h3>3.4 过滤与特异性</h3>
<table><caption>表 4. 硬过滤命中分布 (全部 {c['consensus_candidates']} 条共识候选)。
单条候选可命中多个原因, 故合计 ({sum(hard.values())}) 不小于被淘汰数 ({c['rejected_by_filters']})。</caption>
<thead><tr><th>原因</th><th class="n">条数</th></tr></thead><tbody>{flag_rows(hard)}</tbody></table>
<table><caption>表 5. 软标记分布 (仅告警, 不淘汰)。</caption>
<thead><tr><th>原因</th><th class="n">条数</th></tr></thead><tbody>{flag_rows(soft)}</tbody></table>
<p>候选的特异性状态全部为 <code>input_only</code> ({c['candidates']} 条),
<code>clean</code> 计数为 {c['evaluator_clean']}。因此本报告<b>不主张</b>任何脱靶结论。</p>

<h2><span class="num">4.</span> 图</h2>
<figure>
{fig_funnel(c)}
<figcaption><b>图 1. 运行流程漏斗。</b>深色为进入下一阶段的量, 浅红色为丢弃量。
从 {c['hits']:,} 条原始命中收敛到 {c['top']} 条最终候选, 主要损失来自投票阈值
({c['below_min_votes']} 个位点) 与序列过滤。</figcaption>
</figure>
<figure>
{fig_map(run, top_keys)}
<figcaption><b>图 2. 候选在转录本上的位置分布。</b>横轴为 1-based 转录本坐标 (每 200 nt 一格)。
三个轨道自上而下依次为: 全部共识候选、通过序列过滤的候选、最终 Top-{len(top)}。
Top-K 因最短 {cfg['min_spacing']} nt 间距约束, 在候选密集区互相避让, 因而比全部候选稀疏得多;
两个生成器的一致性因此主要集中在 272–1066 nt 的编码区一带。</figcaption>
</figure>
<figure>
{fig_flags(hard, soft)}
<figcaption><b>图 3. 过滤命中的原因构成。</b>硬过滤 (红) 决定淘汰, 软标记 (琥珀) 仅记录。
最主要的硬过滤原因是引导链 5′ 端非 A/U。</figcaption>
</figure>

<h2><span class="num">5.</span> 讨论与局限</h2>
<ul class="plain">{limits}</ul>
<p>需要额外强调三点。其一, 共识票来自生成器工具, SIREN 的枚举命中与 input-only 回退
<b>不是</b>设计票, 本报告全程按此处理。其二, 特异性是<b>独立证据</b>:
只有转录组支持且脱靶数超阈时才是硬性不合格, 而本次既非 "合格" 也非 "不合格", 而是 "未验证"。
其三, Top-K 强制满足投票阈值, 单例位点仍保留在 <code>all_candidates.fa</code> 中供查看。</p>

<h2><span class="num">6.</span> 结论</h2>
<p>在一个 1234 nt 的靶转录本上, 两个独立生成器工具在 {c['consensus_unanimous']} 个位点上一致,
经序列过滤与间距约束后给出 {c['top']} 条订单候选, 全过程可复现且记录完整。
下一步最关键的补强是<b>提供转录组并重跑</b>以获得特异性结论, 以及用带有实测敲低数据的
数据集 (见 <code>benchmarks/</code>) 做效力对照。</p>

<h2><span class="num">7.</span> 数据与复现</h2>
<table><caption>表 6. 本报告关联的输出文件 (位于 <code>examples/AREG/</code>)。</caption>
<thead><tr><th>文件</th><th>内容</th></tr></thead><tbody>
<tr><td><code>run.json</code></td><td>完整运行记录: 全部归一化命中、配置、哈希、候选</td></tr>
<tr><td><code>candidates.tsv</code></td><td>全部 {c['consensus_candidates']} 条候选的表格</td></tr>
<tr><td><code>all_candidates.fa</code></td><td>全部候选序列, 带通过/淘汰标记</td></tr>
<tr><td><code>top_guides.fa</code></td><td>仅 {c['top']} 条最终候选</td></tr>
<tr><td><code>report.html</code></td><td>英文版自动报告</td></tr>
<tr><td><code>report_zh.html</code></td><td>本中文报告</td></tr>
</tbody></table>
<ol class="refs">
<li>siRNAdesigner —— 上游去新设计脚本 (本仓库 <code>vendor/</code> 内置回退副本)。</li>
<li>siRNAforge —— 容器化设计工具, 固定镜像 <code>ghcr.io/austin-s-h/sirnaforge:sha-4c05e4b</code>。</li>
<li>SIREN —— 枚举式脱靶评估器, 本流程中仅作证据来源。</li>
</ol>

<footer>本报告由 <code>examples/AREG/make_report_zh.py</code> 从 <code>run.json</code> 直接渲染,
所有数值与图均未重新计算。生成时间 {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')} UTC。
计算预测, 非湿实验验证。</footer>
</div></body></html>
"""


def main(argv: list[str]) -> int:
    src = Path(argv[1]) if len(argv) > 1 and not argv[1].startswith("-") else HERE / "run.json"
    out = HERE / "report_zh.html"
    if "-o" in argv:
        out = Path(argv[argv.index("-o") + 1])
    run = json.loads(src.read_text())
    if run.get("schema_version") not in (1, 2):
        raise SystemExit(f"未知的 schema 版本: {run.get('schema_version')!r}")
    out.write_text(build(run), encoding="utf-8")
    print(f"已写入 {out} ({out.stat().st_size} 字节)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
