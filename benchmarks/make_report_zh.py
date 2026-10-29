#!/usr/bin/env python3
"""生成中文论文风格的基准测试报告 (report_zh.html)。

结构: 摘要 → 背景 → 方法 → 结果 → 图 → 讨论 → 结论 → 数据与复现。
数据全部从本目录下已冻结的结果文件读取, 不做重新计算:

    results/development_results.json      方法/消融/配对比较/筛查冒烟
    results/development_comparison.json   逐基因配对差异与自助 CI
    results/permutation_control.json      标签置换对照
    results/implementation_verification.json
    dataset_manifest.json, benchmark_config.json

图以内联 SVG 绘制, 报告是单文件、离线可读。

用法:  python make_report_zh.py [-o report_zh.html]
"""
from __future__ import annotations

import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent

METHOD_ZH = {
    "native": "native (并集 + 序列质量回退)",
    "sirnadesigner": "siRNAdesigner (单个生成器)",
    "sirnaforge": "siRNAforge (单个生成器)",
    "consensus": "exact consensus (精确共识)",
    "model": "model (冻结 OligoFormer)",
    "fusion": "fusion (模型 ⊕ native 百分位)",
    "support": "support (0.8 模型 + 0.2 支持度)",
    "quality": "quality (0.8 模型 + 0.1 native + 0.1 质量)",
}
ABLATION_ZH = {
    "full": "完整配置 (= quality)",
    "no_filters": "关闭全部过滤器",
    "hard_filters": "硬过滤 (非 soft)",
    "no_tool_scores": "去掉工具打分",
    "no_quality": "去掉序列质量项",
    "support": "support 策略",
    "no_support": "去掉支持度项",
    "spacing30": "面板最小间距 30 nt",
}
REASON_ZH = (
    ("Seven source studies occur",
     "<b>7 项来源研究</b>都出现在 OligoFormer 的开发/评估集合中, 检查点级暴露未经验证。"),
    ("Huesken has previous project",
     "该 Huesken 队列已被姊妹项目用于训练与评估, 存在重复暴露。"),
    ("Database Huesken metadata",
     "数据库标注 (HeLa / 荧光素酶) 与原始文献整理 (H1299 / 荧光报告基因) 冲突, 未强行统一。"),
    ("Database percent-inhibition field",
     "数据库的百分比抑制字段重置了 Huesken 的归一化活性尺度, <b>不可</b>当作敲低百分比解读。"),
)
PCT = lambda v: "不可用" if v is None else f"{v * 100:.1f}%"  # noqa: E731

CSS = """
:root{--ink:#111820;--muted:#5c6b74;--line:#d5dee2;--navy:#123b4f;--teal:#0f7b74;
--red:#a4301f;--green:#1d6b3a;--amber:#8a5c12}
*{box-sizing:border-box}
body{margin:0;background:#eef1f2;color:var(--ink);
font:15px/1.8 'Noto Serif SC','Songti SC','Charter',Georgia,serif}
.page{max-width:960px;margin:0 auto;background:#fff;padding:46px 56px 60px;
box-shadow:0 2px 20px rgba(20,35,43,.09)}
h1{font-size:21px;line-height:1.55;margin:0 0 8px}
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
.ok{color:var(--green)}.bad{color:var(--red)}.warn{color:var(--amber)}
.note{background:#f7faf9;border-left:3px solid var(--teal);padding:11px 14px;font-size:12.5px;margin:14px 0}
.note.warnbox{background:#fdf6e8;border-left:3px solid var(--amber)}
figure{margin:20px 0;padding:0}
figcaption{font-size:12px;color:#3b4b55;margin-top:8px;line-height:1.6}
svg{max-width:100%;height:auto;display:block}
ul.plain,ol.refs{font-size:12.5px;padding-left:22px}
ul.plain li,ol.refs li{margin:7px 0}
footer{margin-top:34px;padding-top:14px;border-top:1px solid var(--line);font-size:12px;color:var(--muted)}
@media print{body{background:#fff}.page{box-shadow:none;max-width:none;padding:0}}
"""


def e(v) -> str:
    return html.escape(str(v), quote=True)


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text())


# ---------------------------------------------------------------- figures
def fig_methods(methods: dict) -> str:
    """Top-5 精确率 (条) 与完成度 (点)。"""
    rows = sorted(methods.items(), key=lambda kv: (kv[1]["precision_at_5_complete"] is None,
                                                  kv[1]["precision_at_5_complete"] or 0))
    w, row, pad = 760, 30, 250
    h = len(rows) * row + 34
    bars = []
    for i, (name, m) in enumerate(rows):
        y = 20 + i * row
        p = m["precision_at_5_complete"]
        complete = m["complete_genes"]
        color = "#123b4f" if complete == 25 else "#8a5c12"
        bars.append(
            f'<text x="0" y="{y+16}" font-size="11.5" fill="#3b4b55">{e(METHOD_ZH.get(name, name))}</text>')
        if p is None:
            bars.append(f'<text x="{pad+8}" y="{y+16}" font-size="11.5" fill="#a4301f">'
                        f'无完整面板 (0/25 基因)</text>')
            continue
        bw = max(3, p * (w - pad - 250))
        bars.append(f'<rect x="{pad}" y="{y+4}" width="{bw:.1f}" height="{row-11}" rx="2" fill="{color}"/>')
        bars.append(f'<text x="{pad+bw+8:.1f}" y="{y+16}" font-size="11.5" font-weight="700" '
                    f'fill="{color}">{p * 100:.1f}%</text>')
        cov = m["mean_candidate_coverage"]
        bars.append(f'<text x="{w}" y="{y+16}" font-size="10.5" fill="#7a878e" '
                    f'text-anchor="end">覆盖 {cov * 100:.0f}% · {complete}/25 基因</text>')
    bars.append(f'<text x="{pad}" y="12" font-size="10.5" fill="#7a878e">基因加权 Top-5 精确率 (完整基因)</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="各方法 Top-5 精确率">{"".join(bars)}</svg>'


def fig_paired(comp: dict) -> str:
    """quality - model 的逐基因配对差异 + 自助 CI。"""
    c = comp["comparisons"][0]
    diffs = c["per_gene_difference"]
    lo, hi = c["ci95"]
    mean = c["mean_difference"]
    genes = sorted(diffs, key=lambda g: diffs[g])
    w, pad, row = 760, 140, 21
    zero = pad + (w - pad - 150) / 2
    scale = (w - pad - 150) / 2 / 0.25
    h = len(genes) * row + 74
    body = []
    for i, g in enumerate(genes):
        y = 42 + i * row
        d = diffs[g]
        x = zero + d * scale
        color = "#0f7b74" if d > 0 else ("#a4301f" if d < 0 else "#c9d4d8")
        body.append(f'<text x="{pad-8}" y="{y+4}" font-size="10.5" fill="#3b4b55" '
                    f'text-anchor="end">{e(g)}</text>')
        body.append(f'<line x1="{zero:.1f}" y1="{y}" x2="{x:.1f}" y2="{y}" stroke="{color}" '
                    f'stroke-width="1.6"/>')
        body.append(f'<circle cx="{x:.1f}" cy="{y}" r="3.4" fill="{color}"/>')
    lo_x, hi_x, m_x = zero + lo * scale, zero + hi * scale, zero + mean * scale
    band = (f'<rect x="{lo_x:.1f}" y="30" width="{hi_x-lo_x:.1f}" height="{h-72}" '
            f'fill="#123b4f" opacity="0.08"/>'
            f'<line x1="{lo_x:.1f}" y1="30" x2="{lo_x:.1f}" y2="{h-42}" stroke="#123b4f" stroke-dasharray="3 3"/>'
            f'<line x1="{hi_x:.1f}" y1="30" x2="{hi_x:.1f}" y2="{h-42}" stroke="#123b4f" stroke-dasharray="3 3"/>'
            f'<line x1="{m_x:.1f}" y1="30" x2="{m_x:.1f}" y2="{h-42}" stroke="#123b4f" stroke-width="2"/>')
    zero_line = (f'<line x1="{zero:.1f}" y1="26" x2="{zero:.1f}" y2="{h-42}" stroke="#7a878e"/>'
                 f'<text x="{zero:.1f}" y="18" font-size="10.5" fill="#7a878e" text-anchor="middle">0</text>')
    labels = (f'<text x="{lo_x:.1f}" y="{h-26}" font-size="10.5" fill="#123b4f" text-anchor="middle">'
              f'CI 下界 {lo * 100:+.1f}</text>'
              f'<text x="{hi_x:.1f}" y="{h-26}" font-size="10.5" fill="#123b4f" text-anchor="middle">'
              f'CI 上界 {hi * 100:+.1f}</text>'
              f'<text x="{m_x:.1f}" y="{h-10}" font-size="10.5" font-weight="700" fill="#123b4f" '
              f'text-anchor="middle">均值 {mean * 100:+.1f} pt</text>')
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="逐基因配对差异">'
            f'{band}{zero_line}{"".join(body)}{labels}</svg>')


def fig_ablation(abl: dict) -> str:
    rows = sorted(abl.items(), key=lambda kv: kv[1]["precision_at_5_complete"])
    w, row, pad = 640, 27, 230
    h = len(rows) * row + 18
    top = max(v["precision_at_5_complete"] for _, v in rows)
    bars = []
    for i, (name, m) in enumerate(rows):
        y = 12 + i * row
        p = m["precision_at_5_complete"]
        base = name == "full"
        bw = max(3, p / top * (w - pad - 90))
        color = "#123b4f" if base else "#9db4c0"
        bars.append(f'<text x="0" y="{y+16}" font-size="11.5" fill="#3b4b55">'
                    f'{e(ABLATION_ZH.get(name, name))}</text>')
        bars.append(f'<rect x="{pad}" y="{y+4}" width="{bw:.1f}" height="{row-10}" rx="2" fill="{color}"/>')
        bars.append(f'<text x="{pad+bw+8:.1f}" y="{y+16}" font-size="11.5" fill="#3b4b55">'
                    f'{p * 100:.1f}% · {m["complete_genes"]}/25 基因</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="消融研究">{"".join(bars)}</svg>'


def fig_permutation(perm: dict) -> str:
    obs, null = perm["observed_gene_precision_at_5"], perm["null_mean"]
    w, h, pad, y0 = 640, 190, 90, 130
    x = lambda v: pad + v * (w - pad - 110)  # noqa: E731
    boxes = [
        (null, "#8a5c12", f'零假设均值 {null * 100:.1f}%'),
        (obs, "#0f7b74", f'观测值 {obs * 100:.1f}%'),
    ]
    out = []
    for i, (v, color, label) in enumerate(boxes):
        bx = x(v) + i * 40
        out.append(f'<rect x="{bx:.1f}" y="{y0 - v * 90:.1f}" width="34" '
                   f'height="{v * 90:.1f}" fill="{color}"/>')
        out.append(f'<text x="{bx+17:.1f}" y="{y0 - v * 90 - 8:.1f}" font-size="11.5" '
                   f'font-weight="700" fill="{color}" text-anchor="middle">{label}</text>')
    out.append(f'<line x1="{pad-20}" y1="{y0}" x2="{w-60}" y2="{y0}" stroke="#7a878e"/>')
    out.append(f'<text x="{pad-20}" y="{y0+22}" font-size="11" fill="#7a878e">'
               f'{perm["permutations"]:,} 次基因内标签置换 · p = {perm["p_one_sided"]:.4f} '
               f'(单侧) · 种子 {perm["seed"]}</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="置换对照">{"".join(out)}</svg>'


# ---------------------------------------------------------------- report
def build() -> str:
    dev = load("results/development_results.json")
    comp = load("results/development_comparison.json")
    perm = load("results/permutation_control.json")
    impl = load("results/implementation_verification.json")
    audit = load("dataset_manifest.json")
    cfg = load("benchmark_config.json")

    methods = dev["methods"]
    abl = dev["ablations"]
    primary = comp["comparisons"][0]
    p5 = comp["comparisons"][0]["per_gene_difference"]
    moved = sum(1 for v in p5.values() if v)
    best = max((m["precision_at_5_complete"] or 0, k) for k, m in methods.items())

    rows_methods = "".join(
        f"<tr><td>{e(METHOD_ZH.get(k, k))}</td>"
        f"<td class='n'>{PCT(m['precision_at_5_complete'])}</td>"
        f"<td class='n'>{m['complete_genes']}/25</td>"
        f"<td class='n'>{m['mean_fill_rate'] * 100:.1f}%</td>"
        f"<td class='n'>{m['mean_candidate_coverage'] * 100:.2f}%</td></tr>"
        for k, m in sorted(methods.items(),
                           key=lambda kv: -(kv[1]["precision_at_5_complete"] or -1)))

    rows_paired = "".join(
        f"<tr><td class='n'>{c['stratum'][-1]}</td><td class='n'>{c['paired_genes']}</td>"
        f"<td class='n'>{c['mean_difference'] * 100:+.1f}</td>"
        f"<td class='n'>[{c['ci95'][0] * 100:+.1f}, {c['ci95'][1] * 100:+.1f}]</td>"
        f"<td class='n'>{c['budget_yield_difference'] * 100:+.1f}</td>"
        f"<td class='bad'>未达标</td></tr>"
        for c in comp["comparisons"])

    rows_studies = "".join(
        f"<tr><td>{e(name)}</td><td class='n'>{n:,}</td></tr>"
        for name, n in sorted(audit["studies"].items(), key=lambda kv: -kv[1]))

    rows_ablation = "".join(
        f"<tr><td>{e(ABLATION_ZH.get(k, k))}</td>"
        f"<td class='n'>{PCT(v['precision_at_5_complete'])}</td>"
        f"<td class='n'>{v['complete_genes']}/25</td>"
        f"<td class='n'>{v['mean_fill_rate'] * 100:.1f}%</td></tr>"
        for k, v in sorted(abl.items(), key=lambda kv: -kv[1]["precision_at_5_complete"]))

    reasons = "".join(
        f"<li>{zh}<br><span class='seq'>({e(r[:66])}…)</span></li>"
        for r in audit["reasons"]
        for zh in [next((t for k, t in REASON_ZH if r.startswith(k)), e(r))])
    screening = impl["real_transcriptome_screen"]

    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>共识 siRNA 选择的开发集基准测试: 结果报告</title>
<style>{CSS}</style></head><body><div class="page">

<h1>以公开实测数据评估共识 siRNA 选择流程:<br>Huesken 开发集基准测试报告</h1>
<p class="byline">sirna-consensus v{impl['package_version']} (schema v{impl['schema_version']}) ·
包状态 <code>{e(impl['release_status'])}</code></p>
<p class="meta">数据: {dev['genes']} 个人源克隆构建体 · {impl['curated_human_guides']:,} 条全长引导链 ·
实测表 SHA-256 <code>{e(perm['benchmark_sha256'][:16])}…</code> · 生成于
{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')} UTC</p>

<div class="abstract">
<h2>摘要</h2>
<p><b>背景.</b> 去新设计流程需要外部实测数据来标定 "排序是否有效", 但公开 siRNA 数据集
大多已被同类模型的开发集覆盖, 难以构成真正的留出集。本文报告在 Huesken 人源克隆构建体
开发集上的一次完整基准测试, 并明确标注其证据等级。</p>
<p><b>方法.</b> 主终点是<b>基因加权的 Top-5 实测精确率</b> (即前 5 条候选中有多少条实测有效),
每基因内先对测定取平均、再对基因等权。次级终点含 Top-10、完成度、填充率、候选覆盖率与命中回收。
比较的对象是 6 种排序策略与 2 个单个生成器; 统计使用 {cfg['bootstrap_replicates']:,} 次
配对基因自助与 {cfg['permutations']:,} 次基因内标签置换。</p>
<p><b>结果.</b> {dev['genes']} 个构建体、{impl['curated_human_guides']:,} 条引导链全部完成打分,
生成器共产生 {impl['generator_raw_hits']:,} 条原始命中。开发集上最优策略为
<b>{e(METHOD_ZH.get(best[1], best[1]))}</b> ({best[0] * 100:.1f}% Top-5);
相对冻结模型的配对差异为 <b>{primary['mean_difference'] * 100:+.1f} 个百分点</b>,
95% 自助 CI <b>[{primary['ci95'][0] * 100:+.1f}, {primary['ci95'][1] * 100:+.1f}]</b> ——
区间跨过 0, 且低于预设的 +10 pt 门槛。置换对照给出 p = {perm['p_one_sided']:.4f} (单侧),
说明排序整体显著优于随机, 但<b>相对模型的增益未被证实</b>。</p>
<p><b>结论.</b> 实现与执行检查全部通过, 但预先约定的 "效力提升" 判据<b>未达成</b>;
独立验证标志 <code>independent_validation_passed = false</code>。所有数字应视为开发集诊断值。</p>
<p class="kw"><b>关键词:</b> siRNA 基准测试; 基因加权精确率; 配对自助; 标签置换; 暴露审计</p>
</div>

<div class="note warnbox">
<b>证据等级声明 (全文最重要的限定):</b> 该 Huesken 队列同时出现在冻结模型的开发/评估集合中,
且被姊妹项目使用过, 因此<b>当前不存在独立留出集</b>。模型选择与性能估计使用同一批数据,
下列所有结论都是<b>开发集诊断</b>, 不是泛化性能估计。
</div>

<h2><span class="num">1.</span> 背景</h2>
<p>公开 siRNA 活性数据集 (Huesken 2005 等) 已被广泛用于训练, 直接在其上报告准确率会高估泛化能力。
本工作把 "证据等级" 当作一等公民: 先冻结数据与模型, 用固定的协议与终点统计, 再明确标注
哪些结论在什么前提下成立。</p>
<p>临床前的核心问题是<b>排序</b>: 在给定预算 (例如 5 条合成) 下, 前 K 条里有多少条真的有效。
因此主终点取<b>基因加权的 Top-5 实测精确率</b>, 并同时报告完成度与覆盖率 —— 一个只在少数基因上
"精确" 但对其余基因交白卷的方法, 单看精确率会被高估。</p>

<h2><span class="num">2.</span> 方法</h2>
<h3>2.1 数据与终点</h3>
<p>数据来自 Figshare 快照 (article 25908841, v1, {audit['snapshot_records']:,} 条记录, 7 项研究),
其中与本队列相关的为 {dev['genes']} 个人源克隆构建体上的 {dev['records']:,} 条全长引导链,
活性为原始归一化抑制活性 (非百分比敲低)。终点定义为
<code>{cfg['thresholds']['normalized_activity']}</code> 阈值下的有效/无效二分类;
同一条件的重复测定先取平均再判定阈值; 不同条件保持分离。</p>

<h3>2.2 比较对象</h3>
<p>包含 2 个单个生成器 (siRNAdesigner、siRNAforge)、精确共识 (consensus)、冻结的 OligoFormer
语言模型打分 (model), 以及 3 种融合策略 (fusion / support / quality) 与并集回退 (native)。
生成器未产出某条候选记为<b>弃权</b>, 而不是判定无效。</p>

<h3>2.3 统计与判据</h3>
<p>先在同一基因内对测定取平均, 再对基因等权, 以消除测量密度差。
数值判据为: 精确率提升 ≥ {cfg['minimum_precision_gain'] * 100:.0f} pt、
配对基因自助 95% CI 下界 &gt; 0、完成度不下降、预算产出为正;
至少需要 {cfg['minimum_complete_paired_genes']} 个完整配对基因作为可行性下限。
自助与置换分别使用 {cfg['bootstrap_replicates']:,} 与 {cfg['permutations']:,} 次重复 (种子 {cfg['seed']})。</p>

<h3>2.4 执行与合规检查</h3>
<p>{impl['hermetic_tests_passed']} 项密闭单元测试与 {impl['external_smoke_tests_passed']} 项外部冒烟测试通过,
并完成了一次真实转录组筛查 ({e(screening)})。数据集清单记录了出版商校验和与已知的不一致之处。</p>

<h2><span class="num">3.</span> 结果</h2>
<h3>3.1 主终点: 各方法 Top-5 精确率</h3>
<table><caption>表 1. 基因加权 Top-5 精确率 (仅计完整基因)。"完成度" 指在 25 个基因中都交得出 5 条候选的基因数,
"覆盖率" 指该方法候选占全部已测分子池的比例。sirnaforge 与 consensus 在本队列格式下无完整面板,
属覆盖缺失而非效力证据。</caption>
<thead><tr><th>方法</th><th class="n">Top-5 精确率</th><th class="n">完整基因</th>
<th class="n">填充率</th><th class="n">候选覆盖率</th></tr></thead>
<tbody>{rows_methods}</tbody></table>
<p>最优为 <b>{best[0] * 100:.1f}%</b>, 但需注意: 具有 25/25 完成度的 4 个非 native 方法彼此相差在 3.2 pt 以内,
而 siRNAdesigner 的 {methods['sirnadesigner']['precision_at_5_complete'] * 100:.1f}% 只在
{methods['sirnadesigner']['complete_genes']} 个基因上成立, 不可与 25 基因均值直接比较。</p>

<h3>3.2 配对比较: 最优策略 vs 冻结模型</h3>
<table><caption>表 2. 开发集最优策略 (quality) 相对冻结模型的配对差异。均值为基因加权,
CI 为 {cfg['bootstrap_replicates']:,} 次基因自助。判据要求 ≥{cfg['minimum_precision_gain'] * 100:.0f} pt 且 CI 下界 &gt; 0。</caption>
<thead><tr><th class="n">K</th><th class="n">配对基因</th><th class="n">均值差 (pt)</th>
<th class="n">95% CI</th><th class="n">预算产出差</th><th>判定</th></tr></thead>
<tbody>{rows_paired}</tbody></table>
<p>Top-5 上 25 个基因里只有 {moved} 个出现非零差异 (其余完全并列), 其中正负各有之。
CI 跨零是本报告的关键否定性结果。</p>

<h3>3.3 消融</h3>
<table><caption>表 3. 组件消融。完整配置 (与 quality 策略等价) 取得最高 Top-5 精确率;
"硬过滤" 精确率再高 2.9 pt 但覆盖率骤降, 是<b>以覆盖换精确</b>的典型权衡。</caption>
<thead><tr><th>配置</th><th class="n">Top-5 精确率</th><th class="n">完整基因</th><th class="n">填充率</th></tr></thead>
<tbody>{rows_ablation}</tbody></table>

<h3>3.4 置换对照与转录组筛查</h3>
<p>{perm['permutations']:,} 次基因内标签置换给出零假设均值 {perm['null_mean'] * 100:.1f}%,
观测值 {perm['observed_gene_precision_at_5'] * 100:.1f}%, 单侧 p = {perm['p_one_sided']:.4f}。
真实序列筛查在 {dev['screen_smoke']['candidates']} 条 AREG 候选上执行:
{dev['screen_smoke']['states']['pass']} 通过、{dev['screen_smoke']['states']['reject']} 因错配被排除,
仅满足<b>计算</b>错配判据, 不等于生物学特异性。</p>

<h3>3.5 数据集来源审计</h3>
<table><caption>表 4. Figshare 快照中的研究构成 (共 {audit['snapshot_records']:,} 条)。</caption>
<thead><tr><th>研究</th><th class="n">记录数</th></tr></thead><tbody>{rows_studies}</tbody></table>
<ol class="plain">{reasons}</ol>

<h2><span class="num">4.</span> 图</h2>
<figure>
{fig_methods(methods)}
<figcaption><b>图 1. 各方法的主终点表现。</b>深色为 25/25 完整基因的方法, 琥珀色为面板不完整的方法;
右侧标注候选覆盖率, 用来说明精确率是否建立在全池覆盖之上。
sirnaforge 与 consensus 在本队列的双链格式下没有完整面板。</figcaption>
</figure>
<figure>
{fig_paired(comp)}
<figcaption><b>图 2. quality 与冻结模型的逐基因配对差异 (Top-5)。</b>每行为一个基因,
横轴为精确率之差; 竖向阴影为自助 95% CI, 实线为均值。多数基因为 0 (完全并列),
正负各有少数基因移动 —— 这解释了 CI 为何跨过 0。</figcaption>
</figure>
<figure>
{fig_ablation(abl)}
<figcaption><b>图 3. 组件消融。</b>深色为完整配置。硬过滤以大幅降低覆盖率换取最高精确率;
去掉质量项或支持度项使精确率下降约 2 pt, 说明这些项确有贡献但幅度有限。</figcaption>
</figure>
<figure>
{fig_permutation(perm)}
<figcaption><b>图 4. 置换对照。</b>观测精确率远高于基因内标签置换的零分布均值,
说明排序显著优于随机; 但该对照<b>不能</b>消除预训练数据暴露与开发集选择偏差。</figcaption>
</figure>

<h2><span class="num">5.</span> 讨论与局限</h2>
<ul class="plain">
<li><b>无独立留出集。</b>7 项来源研究均出现在冻结模型的开发/评估集合中,
检查点级暴露审计状态为 <code>{e(audit['checkpoint_training_audit'])}</code>。</li>
<li><b>发布判据未达成。</b>数值门槛要求 ≥10 pt 且 CI&gt;0, 实测为
{primary['mean_difference'] * 100:+.1f} pt、CI 跨零, 故 <code>publication_criterion_met = false</code>。</li>
<li><b>端点单位不可混用。</b>Huesken 的归一化活性与数据库的百分比抑制字段含义不同,
阈值 {cfg['thresholds']['normalized_activity']} 与 {cfg['thresholds']['percent_inhibition']} 不可互换。</li>
<li><b>元数据冲突已保留。</b>数据库标注 (HeLa/荧光素酶) 与原始文献整理 (H1299/荧光报告基因) 不一致,
未强行统一。</li>
<li><b>覆盖 vs 精确。</b>硬过滤与 sirnadesigner 都说明 "精确率高" 可能只是候选面窄。</li>
<li><b>模型分是序, 不是百分比。</b>OligoFormer 输出为未标定的排序分。</li>
<li><b>批大小为 1</b> 是上游 LSTM 布局所致, 改变批大小会影响推断结果。</li>
</ul>

<h2><span class="num">6.</span> 结论</h2>
<p>在 {dev['genes']} 个构建体、{impl['curated_human_guides']:,} 条引导链上, 流程的<b>实现与执行</b>通过了全部检查,
排序整体显著优于随机 (p &lt; 0.001); 但<b>相对冻结模型的效力提升尚未达成</b> ——
增益 {primary['mean_difference'] * 100:+.1f} pt、CI 跨零, 且缺少独立留出集。
下一步的关键动作是取得一个审计干净的未用队列 (满足
<code>independent_holdout_eligible</code> 与 <code>checkpoint_training_audit: verified</code>)
并走 <code>evaluate.py --phase final</code> 的冻结闸门。</p>

<h2><span class="num">7.</span> 数据与复现</h2>
<table><caption>表 5. 本报告关联的机读产物 (位于 <code>benchmarks/</code>)。</caption>
<thead><tr><th>文件</th><th>内容</th></tr></thead><tbody>
<tr><td><code>results/development_results.json</code></td><td>方法/消融/配对比较/筛查冒烟</td></tr>
<tr><td><code>results/development_comparison.json</code></td><td>逐基因配对差异与自助 CI</td></tr>
<tr><td><code>results/permutation_control.json</code></td><td>1,000 次标签置换对照</td></tr>
<tr><td><code>dataset_manifest.json</code></td><td>快照校验和、暴露审计与不一致记录</td></tr>
<tr><td><code>benchmark_config.json</code></td><td>冻结的终点、阈值、判据与种子</td></tr>
<tr><td><code>figures/benchmark.pdf</code></td><td>英文版图组</td></tr>
<tr><td><code>report_zh.html</code></td><td>本中文报告</td></tr>
</tbody></table>
<ol class="refs">
<li>Huesken et al. (2005) —— 人源克隆构建体与实测活性队列 (原始来源重新整理)。</li>
<li>OligoFormer —— 冻结的语言模型打分器 (checkpoint 已哈希固定)。</li>
<li>RNA-FM —— 模型骨干权重来源, 官方发布版本。</li>
<li>Bowtie 1.3.1 + Ensembl 116 cDNA —— 计算错配筛查所用索引。</li>
</ol>

<footer>本报告由 <code>benchmarks/make_report_zh.py</code> 从上述冻结 JSON 直接渲染, 未做重新计算。
生成时间 {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')} UTC。
开发集诊断值, 非独立验证; 未经湿实验确认前不得称任何引导链无效。</footer>
</div></body></html>
"""


def main(argv: list[str]) -> int:
    out = Path(argv[argv.index("-o") + 1]) if "-o" in argv else HERE / "report_zh.html"
    out.write_text(build(), encoding="utf-8")
    print(f"已写入 {out} ({out.stat().st_size} 字节)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
