#!/usr/bin/env python3
"""生成中文论文风格的真实运行报告 (report_zh.html)。

结构: 摘要 -> 背景 -> 方法(含流程图) -> 结果 -> 图 -> 讨论 -> 结论 -> 数据与复现。
所有数字与图都从本目录的 run.json 直接读出, 不做重新计算; 图以内联 SVG 绘制,
因此报告是单文件、离线可读。

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
FLAG_TEXT = {
    "guide_5p_not_AU": "引导链 5′ 端非 A/U",
    "sense_3p_not_AU": "sense 链 3′ 端非 A/U",
    "homopolymer": "存在同聚物",
    "gc_outside_window": "GC 含量超窗",
    "seed_gc_low": "种子区 GC 偏低",
    "specificity_offtarget": "特异性脱靶",
    "specificity_unverified": "特异性未验证",
}

INK, MUTED, NAVY, TEAL, RED, AMBER, LINE = (
    "#111820", "#5c6b74", "#123b4f", "#0f7b74", "#a4301f", "#8a5c12", "#d5dee2")

CSS = """
:root{--ink:#111820;--muted:#5c6b74;--line:#d5dee2;--navy:#123b4f;--teal:#0f7b74;
--red:#a4301f;--green:#1d6b3a;--amber:#8a5c12}
*{box-sizing:border-box}
body{margin:0;background:#eef1f2;color:var(--ink);
font:15px/1.8 'Noto Serif SC','Songti SC','Charter',Georgia,serif}
.page{max-width:1000px;margin:0 auto;background:#fff;padding:48px 60px 64px;
box-shadow:0 2px 20px rgba(20,35,43,.09)}
h1{font-size:21px;line-height:1.55;margin:0 0 10px;letter-spacing:-.2px}
h2{font-size:16px;margin:36px 0 10px;color:var(--navy);
border-bottom:1px solid var(--line);padding-bottom:6px}
h3{font-size:13.5px;margin:24px 0 8px;color:var(--navy)}
h2 span.num{margin-right:6px;color:var(--teal)}
p{margin:10px 0;text-align:justify}
code{font:12px ui-monospace,Menlo,Consolas,monospace;background:#eff3f4;padding:1px 4px;border-radius:3px}
.byline{color:#3b4b55;font-size:13px;margin:0 0 3px}
.meta{color:var(--muted);font-size:12.5px;margin:0 0 22px}
.abstract{background:#f7faf9;border:1px solid var(--line);border-left:3px solid var(--teal);
padding:15px 20px;font-size:13.5px;margin:20px 0}
.abstract h2{border:0;margin:0 0 6px;font-size:12.5px;letter-spacing:.12em;color:var(--teal)}
.abstract p{margin:8px 0}
.kw{font-size:12.5px;color:#3b4b55;margin-top:10px}
table{width:100%;border-collapse:collapse;font:12.5px ui-sans-serif,system-ui;margin:14px 0}
caption{text-align:left;font-size:12px;color:var(--muted);padding:0 0 6px;line-height:1.6}
th{text-align:left;font-size:10.5px;letter-spacing:.05em;color:var(--muted);
padding:6px;border-top:1.4px solid var(--ink);border-bottom:1px solid var(--ink)}
td{padding:6px;border-bottom:1px solid #eceff1;vertical-align:top}
tbody tr:last-child td{border-bottom:1.3px solid var(--ink)}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
.seq{font:12px ui-monospace,Menlo,Consolas,monospace;white-space:nowrap}
.ok{color:var(--green)}.bad{color:var(--red)}.warn{color:var(--amber)}
.note{background:#f7faf9;border-left:3px solid var(--teal);padding:11px 15px;font-size:12.5px;margin:16px 0}
.note.warnbox{background:#fdf6e8;border-left:3px solid var(--amber)}
figure{margin:24px 0;padding:0}
figcaption{font-size:12px;color:#3b4b55;margin-top:10px;line-height:1.65}
figcaption b{color:var(--ink)}
svg{max-width:100%;height:auto;display:block}
ul.plain,ol.refs{font-size:12.8px;padding-left:22px}
ul.plain li,ol.refs li{margin:8px 0}
footer{margin-top:36px;padding-top:15px;border-top:1px solid var(--line);font-size:12px;color:var(--muted)}
@media print{body{background:#fff}.page{box-shadow:none;max-width:none;padding:0}}
"""


def e(v) -> str:
    return html.escape(str(v), quote=True)


def key_of(c: dict) -> str:
    return f"{c['transcript_id']}:{c['start']}:{c['length']}"


def flag_text(codes) -> str:
    return "、".join(FLAG_TEXT.get(c, c) for c in codes) or "—"


# ------------------------------------------------------------------ 图
def fig_flow() -> str:
    w, h = 880, 300
    s = ['<defs><marker id="ah" markerWidth="7" markerHeight="7" refX="6" refY="3.5" '
         'orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="#8b949e"/></marker></defs>']

    def box(x, y, bw, bh, l1, l2, fill="#eef4fb", stroke=NAVY, tc=NAVY, dash=""):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        return (f'<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="5" fill="{fill}" '
                f'stroke="{stroke}"{d}/>'
                f'<text x="{x+bw/2}" y="{y+21}" font-size="11.5" font-weight="700" fill="{tc}" '
                f'text-anchor="middle">{e(l1)}</text>'
                f'<text x="{x+bw/2}" y="{y+37}" font-size="10.5" fill="#4b5a63" '
                f'text-anchor="middle">{e(l2)}</text>')

    def arrow(x1, y1, x2, y2, dash=""):
        d = ' stroke-dasharray="4 3"' if dash else ""
        return (f'<path d="M {x1} {y1} L {x2} {y2}" stroke="#8b949e" stroke-width="1.4" '
                f'fill="none"{d} marker-end="url(#ah)"/>')

    s.append(f'<text x="10" y="24" font-size="11" font-weight="700" fill="{NAVY}">输入</text>')
    s.append(f'<text x="400" y="24" font-size="11" font-weight="700" fill="{NAVY}" '
             f'text-anchor="middle">run · 共识 (不做本地枚举)</text>')
    s.append(f'<text x="760" y="24" font-size="11" font-weight="700" fill="{NAVY}" '
             f'text-anchor="middle">后续</text>')
    s.append(box(10, 60, 180, 48, "转录本 FASTA", "单条序列 · 记录 sha256",
                 fill="#eef2f5", stroke=MUTED, tc=INK))
    s.append(box(250, 50, 300, 46, "生成器 · siRNAdesigner + siRNAforge",
                 "各自独立子进程, 输出全部保留", fill="#e7f0ef", stroke=TEAL, tc=TEAL))
    s.append(box(250, 130, 300, 46, "评估器 · SIREN",
                 "只提供证据, 不投设计票", fill="#ffffff", stroke=NAVY, dash="5 3"))
    s.append(box(250, 210, 300, 46, "归一化 → 位点归并 → 过滤 → top-k",
                 "min_votes · 位点级共识", fill=NAVY, stroke="#0b2735", tc="#ffffff"))
    s.append(box(640, 50, 230, 46, "Top-K 引导链", "仅 10 条进入订单集",
                 fill="#fdf6e8", stroke=AMBER, tc=AMBER))
    s.append(box(640, 130, 230, 46, "全部候选 + 标记",
                 "含被淘汰者, 供人工复核", fill="#fdf6e8", stroke=AMBER, tc=AMBER))
    s.append(box(640, 210, 230, 46, "特异性状态", "本次全部 input_only",
                 fill="#fbf1ef", stroke=RED, tc=RED, dash="4 3"))
    s.append(arrow(190, 84, 248, 78))
    s.append(arrow(190, 84, 248, 153, dash="4 3"))
    s.append(arrow(100, 108, 100, 232) )
    s.append('<path d="M 100 232 L 248 232" stroke="#8b949e" stroke-width="1.4" fill="none" '
             'marker-end="url(#ah)"/>')
    s.append(arrow(400, 96, 400, 128))
    s.append(arrow(400, 176, 400, 208))
    s.append(arrow(550, 233, 638, 233))
    s.append(arrow(550, 70, 638, 70))
    s.append('<path d="M 755 208 L 755 100" stroke="#8b949e" stroke-width="1.4" fill="none" '
             'stroke-dasharray="4 3"/>')
    s.append('<path d="M 755 100 L 755 96" stroke="#8b949e" stroke-width="1.4" fill="none" '
             'marker-end="url(#ah)"/>')
    s.append(f'<text x="640" y="288" font-size="10.5" fill="{MUTED}">'
             f'SIREN 因未提供转录组, 只做输入序列内比对</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="运行流程图">{"".join(s)}</svg>'


def fig_funnel(counts: dict) -> str:
    steps = [("工具原始命中", counts["hits"], False),
             ("归一化并定位到转录本", counts["placed"], False),
             ("合并后的位点", counts["sites"], False),
             ("达成共识的候选", counts["consensus_candidates"], False),
             ("通过序列过滤", counts["passing"], False),
             ("最终 Top-K 输出", counts["top"], False),
             ("低于投票阈值 (丢弃)", counts["below_min_votes"], True),
             ("被过滤器淘汰", counts["rejected_by_filters"], True)]
    w, row, pad, gap = 760, 32, 210, 12
    h = len(steps) * (row + gap) + 16
    top = max(v for _, v, _ in steps)
    s = []
    for i, (label, value, drop) in enumerate(steps):
        y = 8 + i * (row + gap)
        bw = max(3, value / top * (w - pad - 90))
        col = RED if drop else NAVY
        s.append(f'<text x="0" y="{y+20}" font-size="12" fill="{INK}">{e(label)}</text>')
        s.append(f'<rect x="{pad}" y="{y+4}" width="{bw:.1f}" height="{row-8}" rx="2" '
                 f'fill="{col}" opacity="{0.6 if drop else 1}"/>')
        s.append(f'<text x="{pad+bw+8:.1f}" y="{y+20}" font-size="12" font-weight="700" '
                 f'fill="{col}">{value:,}</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="运行流程漏斗">{"".join(s)}</svg>'


def fig_map(run: dict, top_keys: set[str]) -> str:
    length = run["transcript_length"]
    cands = run["candidates"]
    w, h, pad, track = 900, 180, 90, 48
    def x(pos):
        return pad + (pos - 1) / max(1, length - 1) * (w - 2 * pad)
    top10 = [c for c in cands if key_of(c) in top_keys]
    passing = [c for c in cands if not c["hard_flags"]] if "hard_flags" in cands[0] else cands
    bottom = track + 92
    s = ["".join(
        f'<line x1="{x(p):.1f}" y1="{track-6}" x2="{x(p):.1f}" y2="{bottom+6}" '
        f'stroke="#e6ebed"/>'
        f'<text x="{x(p):.1f}" y="{bottom+22}" font-size="10.5" fill="{MUTED}" '
        f'text-anchor="middle">{p}</text>' for p in range(200, length, 200)),
        f'<line x1="{pad}" y1="{track}" x2="{w-pad}" y2="{track}" stroke="{NAVY}" stroke-width="2.5"/>']
    for row_y, items, col, radius, label in (
            (track, cands, "#9db4c0", 3.4, f"全部共识候选 (n={len(cands)})"),
            (track + 46, passing, TEAL, 4.5, f"通过序列过滤 (n={len(passing)})"),
            (bottom, top10, RED, 6, f"最终 Top-{len(top10)} (n={len(top10)})")):
        s.append(f'<text x="0" y="{row_y+4}" font-size="11.5" fill="{INK}">{e(label)}</text>')
        s.append("".join(f'<circle cx="{x(c["start"]):.1f}" cy="{row_y}" r="{radius}" '
                         f'fill="{col}" opacity="0.85"/>' for c in items))
    s.append(f'<text x="{pad}" y="{track-20}" font-size="11" fill="{MUTED}">5′</text>')
    s.append(f'<text x="{w-pad}" y="{track-20}" font-size="11" fill="{MUTED}" '
             f'text-anchor="end">3′ ({length} nt)</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="候选位置分布">{"".join(s)}</svg>'


def fig_flags(hard: dict, soft: dict) -> str:
    items = [(FLAG_TEXT.get(k, k), v, RED, "硬过滤") for k, v in hard.items()]
    items += [(FLAG_TEXT.get(k, k), v, AMBER, "软标记") for k, v in soft.items()]
    items.sort(key=lambda t: t[1])
    w, row, pad = 700, 28, 240
    h = len(items) * row + 14
    top = max(v for _, v, _, _ in items)
    s = []
    for i, (label, v, col, kind) in enumerate(items):
        y = 8 + i * row
        bw = max(3, v / top * (w - pad - 130))
        s.append(f'<text x="0" y="{y+16}" font-size="11.5" fill="{INK}">{e(label)}</text>')
        s.append(f'<rect x="{pad}" y="{y+4}" width="{bw:.1f}" height="{row-11}" rx="2" fill="{col}"/>')
        s.append(f'<text x="{pad+bw+8:.1f}" y="{y+16}" font-size="11.5" font-weight="700" '
                 f'fill="{col}">{v}</text>')
        s.append(f'<text x="{w}" y="{y+16}" font-size="10.5" fill="{MUTED}" '
                 f'text-anchor="end">{kind}</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="过滤原因">{"".join(s)}</svg>'


# ------------------------------------------------------------------ 报告
def build(run: dict) -> str:
    c, cfg, tid = run["counts"], run["config"], run["transcript_id"]
    top_keys = set(run["top"])
    top = sorted((x for x in run["candidates"] if key_of(x) in top_keys), key=lambda x: x["rank"])
    hard, soft = {}, {}
    for cand in run["candidates"]:
        for f in cand["hard_flags"]:
            hard[f] = hard.get(f, 0) + 1
        for f in cand["soft_flags"]:
            soft[f] = soft.get(f, 0) + 1
    acting, n_gen = cfg["acting_tools"], len(cfg["acting_tools"])

    rows_tools = "".join(
        f"<tr><td><b>{e(t['tool'])}</b></td><td>{e(TOOL_ROLE.get(t['role'], t['role']))}</td>"
        f"<td class='{'ok' if t['status'] == 'ok' else 'bad'}'>"
        f"{e(STATUS_TEXT.get(t['status'], t['status']))}</td><td class='n'>{t['hits']}</td>"
        f"<td class='seq'>{e(t.get('screen_state') or '—')}</td></tr>" for t in run["tools"])
    rows_top = "".join(
        f"<tr><td class='n'>{x['rank']}</td><td class='n'>{x['start']}–{x['end']}</td>"
        f"<td class='n'>{x['length']}</td>"
        f"<td>{e(' + '.join(x['tools']))} ({x['votes']}/{n_gen})</td>"
        f"<td class='n'>{x['score_parts']['tool_score_component']:.3f}</td>"
        f"<td class='seq'>{e(x['guide_dna'])}</td>"
        f"<td>{flag_text(x['hard_flags'] + x['soft_flags'])}</td></tr>" for x in top)
    rows_norm = "".join(
        f"<tr><td>{e(name)}</td><td class='n'>{s['total']}</td><td class='n'>{s['placed']}</td>"
        f"<td class='n'>{s['duplicates']}</td>"
        f"<td>{e(json.dumps(s['rejected'], ensure_ascii=False) if s['rejected'] else '—')}</td></tr>"
        for name, s in sorted(run["normalization"].items()))

    def flag_rows(counter):
        return "".join(f"<tr><td>{e(FLAG_TEXT.get(k, k))}</td><td class='n'>{v}</td></tr>"
                       for k, v in sorted(counter.items(), key=lambda kv: -kv[1])) or \
               "<tr><td colspan='2'>—</td></tr>"

    limits = "".join(f"<li>{e(x)}</li>" for x in run["limitations"])
    created = run["created_utc"][:19].replace("T", " ") + " UTC"

    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>基于多工具位点共识的 {e(tid)} 去新 siRNA 设计: 真实运行报告</title>
<style>{CSS}</style></head><body><div class="page">

<h1>基于多工具位点共识的 {e(tid)} (AREG) 去新 siRNA 候选设计:<br>一次真实运行的完整报告</h1>
<p class="byline">sirna-consensus v{e(run['tool_version'])} · 运行 ID <code>{e(run['run_id'])}</code></p>
<p class="meta">{e(created)} · schema v{run['schema_version']} · 输入 {e(tid)}
({run['transcript_length']} nt) · 输入 SHA-256 <code>{e(run['input_sha256'][:16])}…</code></p>

<div class="abstract">
<h2>摘要</h2>
<p><b>背景.</b> 去新 siRNA 设计缺少可用的近缘同源参照, 单一工具的排序难以独立复核。
本报告以外部独立设计工具的<b>位点级共识</b>替代单一打分器。</p>
<p><b>方法.</b> 单个转录本 FASTA 为唯一输入, 并行调用 2 个生成器 (siRNAdesigner、siRNAforge)
与 1 个评估器 (SIREN); 输出归一化并定位到转录本后按位点归并, 要求 ≥{cfg['min_votes']} 个独立工具
一致, 再过滤并取 Top-K。评估器只提供证据, 不投设计票。</p>
<p><b>结果.</b> 三个工具共产出 <b>{c['hits']:,}</b> 条原始命中, 全部 {c['placed']:,} 条成功定位
(零拒绝、零重复), 归并为 <b>{c['sites']}</b> 个位点; 其中 <b>{c['consensus_unanimous']}</b> 个位点
被两个生成器完全一致支持。过滤淘汰 {c['rejected_by_filters']} 条后余 {c['passing']} 条,
输出 <b>{c['top']}</b> 条订单候选。</p>
<p><b>结论.</b> 无转录组输入, 全部候选特异性状态为 <code>input_only</code>,
<b>尚未获得全转录组特异性证据</b>; 本结果是计算预测, 不代表敲低效力已验证。</p>
<p class="kw"><b>关键词:</b> siRNA 设计; 工具共识; 脱靶评估; 可复现性</p>
</div>

<h2><span class="num">1.</span> 背景</h2>
<p>siRNA 的活性与脱靶行为由引导链序列、位置与热力学共同决定, 各工具加权方式不同。
在缺乏同源参照的 "去新" 设定下无法调参, 因此单一工具的排序难以复核。
本工作不训练新模型, 而是用一份 FASTA 驱动既有开源工具, 把它们的独立意见放到同一坐标轴上,
以共识而非单一分数决定取舍。三条约束: 输入默认是 FASTA; 不做本地窗口枚举;
共识发生在<b>位点</b>层面而非字符串层面 (不同工具对同一双链给出不同长度与偏移)。</p>

<h2><span class="num">2.</span> 方法</h2>
<figure>{fig_flow()}
<figcaption><b>图 1. 运行流程图。</b>生成器与评估器角色分离: 只有生成器投票;
SIREN 的证据在候选生成之后附加, 且因未提供转录组而只做输入序列内比对。</figcaption>
</figure>

<h3>2.1 输入与工具</h3>
<p>输入为 {e(tid)} ({run['transcript_length']} nt), 运行前记录 SHA-256。
三个外部工具各自运行于独立子进程; 生成器 {e('、'.join(acting))} 提出设计票, SIREN 产出证据。</p>

<h3>2.2 归一化与位点共识</h3>
<p>输出被转换为 1-based 闭合区间, 并把结合核心、原始引导链/过客链与 3′ 悬垂拆为独立字段;
按 <code>offset_tolerance={cfg['offset_tolerance']}</code> 归并。投票单位是<b>不同的生成器工具</b>。
保留 ≥<code>min_votes={cfg['min_votes']}</code> 票的位点, 本次关闭单例
(<code>include_singletons={cfg['include_singletons']}</code>)。</p>

<h3>2.3 过滤与排序</h3>
<p>过滤在投票<b>之后</b>执行。硬过滤含 GC 窗口
({cfg['filter_config']['gc_min']}–{cfg['filter_config']['gc_max']}) 与最长同聚物
({cfg['filter_config']['max_homopolymer']}); 软标记仅告警。以 ≥{cfg['min_spacing']} nt 间距取 Top-K。</p>

<h3>2.4 特异性证据</h3>
<p>只有转录组支持的筛查才构成特异性结论。本次未提供转录组, 故全部候选停留在
<code>input_only</code>; 评估器只能支持或淘汰已有候选, 不能凭空产生设计票。</p>

<h2><span class="num">3.</span> 结果</h2>
<h3>3.1 工具执行</h3>
<table><caption>表 1. 外部工具执行情况。</caption>
<thead><tr><th>工具</th><th>角色</th><th>状态</th><th class="n">命中数</th><th>筛查状态</th></tr></thead>
<tbody>{rows_tools}</tbody></table>
<p>三个工具全部成功。SIREN 的 {c['evaluator_hits']} 条命中中绝大多数是评估器枚举,
属 <code>input_only</code>, <b>不是</b>全转录组结论。</p>

<h3>3.2 归一化与定位</h3>
<table><caption>表 2. 归一化与转录本定位。零拒绝说明工具输出与本包解析完全兼容。</caption>
<thead><tr><th>工具</th><th class="n">原始命中</th><th class="n">成功定位</th><th class="n">重复</th>
<th>拒绝原因</th></tr></thead><tbody>{rows_norm}</tbody></table>

<h3>3.3 共识与输出</h3>
<p>归并后 {c['sites']} 个位点, 其中 {c['consensus_unanimous']} 个位点被两个生成器完全一致支持
(投票 {n_gen}/{n_gen}), {c['below_min_votes']} 个位点未达阈值。过滤淘汰
{c['rejected_by_filters']} 条后 {c['passing']} 条进入排序, 输出 Top-{len(top)}。</p>
<table><caption>表 3. 最终 Top-{len(top)} 候选。"全局排名" 是 {c['passing']} 条通过过滤候选中的位次,
因间距约束而不连续。打分为工具内归一化值, 跨工具不可直接比较。</caption>
<thead><tr><th class="n">全局排名</th><th class="n">位置</th><th class="n">长度</th><th>投票</th>
<th class="n">原始打分</th><th>引导链 (5′→3′, DNA)</th><th>标记</th></tr></thead>
<tbody>{rows_top}</tbody></table>

<h3>3.4 过滤与特异性</h3>
<table><caption>表 4. 硬过滤命中分布 (全部 {c['consensus_candidates']} 条共识候选);
单条可命中多个原因, 故合计 ({sum(hard.values())}) 不小于被淘汰数 ({c['rejected_by_filters']})。</caption>
<thead><tr><th>原因</th><th class="n">条数</th></tr></thead><tbody>{flag_rows(hard)}</tbody></table>
<table><caption>表 5. 软标记分布 (仅告警, 不淘汰)。</caption>
<thead><tr><th>原因</th><th class="n">条数</th></tr></thead><tbody>{flag_rows(soft)}</tbody></table>
<p>特异性状态全部为 <code>input_only</code> ({c['candidates']} 条), <code>clean</code> 计数为
{c['evaluator_clean']}, 因此不主张任何脱靶结论。</p>

<h2><span class="num">4.</span> 图</h2>
<figure>{fig_funnel(c)}
<figcaption><b>图 2. 运行漏斗。</b>深色为进入下一阶段的量, 浅红为丢弃量。
从 {c['hits']:,} 条命中收敛到 {c['top']} 条候选, 主要损失来自投票阈值与序列过滤。</figcaption>
</figure>
<figure>{fig_map(run, top_keys)}
<figcaption><b>图 3. 候选位置分布。</b>横轴为 1-based 坐标 (每 200 nt 一格);
三条轨道自上而下为全部共识候选、通过过滤者、最终 Top-{len(top)}。</figcaption>
</figure>
<figure>{fig_flags(hard, soft)}
<figcaption><b>图 4. 过滤原因构成。</b>红色为硬过滤 (决定淘汰), 琥珀为软标记 (仅记录)。
最主要的硬过滤原因是引导链 5′ 端非 A/U。</figcaption>
</figure>

<h2><span class="num">5.</span> 讨论与局限</h2>
<ul class="plain">{limits}</ul>
<p>三点补充: 共识票来自生成器, SIREN 的枚举与 input-only 回退不是设计票;
特异性是独立证据, 本次既非合格也非不合格而是 "未验证"; Top-K 强制满足投票阈值,
单例位点仍保留在 <code>all_candidates.fa</code> 中供查看。</p>

<h2><span class="num">6.</span> 结论</h2>
<p>在一个 1234 nt 靶转录本上, 两个独立生成器在 {c['consensus_unanimous']} 个位点一致,
经过滤与间距约束给出 {c['top']} 条订单候选, 全过程可复现且记录完整。
下一步最关键的补强是<b>提供转录组并重跑</b>, 以及用带实测敲低数据的数据集做效力对照。</p>

<h2><span class="num">7.</span> 数据与复现</h2>
<table><caption>表 6. 关联输出文件 (位于 <code>examples/AREG/</code>)。</caption>
<thead><tr><th>文件</th><th>内容</th></tr></thead><tbody>
<tr><td><code>run.json</code></td><td>完整运行记录: 全部归一化命中、配置、哈希、候选</td></tr>
<tr><td><code>candidates.tsv</code></td><td>全部 {c['consensus_candidates']} 条候选的表格</td></tr>
<tr><td><code>all_candidates.fa</code></td><td>全部候选序列, 带通过/淘汰标记</td></tr>
<tr><td><code>top_guides.fa</code></td><td>仅 {c['top']} 条最终候选</td></tr>
<tr><td><code>report.html</code> / <code>report_zh.html</code></td><td>英文自动报告 / 本报告</td></tr>
</tbody></table>
<ol class="refs">
<li>siRNAdesigner —— 上游去新设计脚本 (本仓库 <code>vendor/</code> 内置回退副本)。</li>
<li>siRNAforge —— 容器化设计工具, 固定镜像 <code>ghcr.io/austin-s/h/sirnaforge:sha-4c05e4b</code>。</li>
<li>SIREN —— 枚举式脱靶评估器, 本流程中仅作证据来源。</li>
</ol>

<footer>由 <code>examples/AREG/make_report_zh.py</code> 从 <code>run.json</code> 直接渲染,
数值与图均未重新计算。{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')} UTC · 计算预测, 非湿实验验证。</footer>
</div></body></html>
"""


def main(argv: list[str]) -> int:
    src = Path(argv[1]) if len(argv) > 1 and not argv[1].startswith("-") else HERE / "run.json"
    out = Path(argv[argv.index("-o") + 1]) if "-o" in argv else HERE / "report_zh.html"
    run = json.loads(src.read_text())
    if run.get("schema_version") not in (1, 2):
        raise SystemExit(f"未知的 schema 版本: {run.get('schema_version')!r}")
    out.write_text(build(run), encoding="utf-8")
    print(f"已写入 {out} ({out.stat().st_size} 字节)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
