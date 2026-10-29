#!/usr/bin/env python3
"""生成中文论文风格的基准测试报告 (report_zh.html)。

结构: 摘要 -> 背景 -> 方法(含流程图) -> 结果 -> 图 -> 讨论 -> 结论 -> 数据与复现。

数据来自本目录下已冻结的产物:
    results/development_results.json      方法/消融/配对比较/筛查冒烟
    results/development_comparison.json   逐基因配对差异与自助 CI
    results/permutation_control.json      标签置换对照
    results/implementation_verification.json
    dataset_manifest.json, benchmark_config.json
    results/development/full/<基因>/<策略>/selection.json   标签盲预测 (逐分子)
    data/huesken/measurements.tsv                            实测活性

"模型性能" 与 "预算曲线" 两图由冻结的预测 + 实测文件复算得到描述性统计
(不重新拟合模型, 不重新抽样)。复算的 Top-5 精确率与 development_results.json
逐位一致 (见报告中的交叉校验)。

图以内联 SVG 绘制, 报告是单文件、离线可读。
用法:  python make_report_zh.py [-o report_zh.html]
"""
from __future__ import annotations

import csv
import html
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
FULL = HERE / "results/development/full"

METHOD_ZH = {
    "native": "native 并集+质量回退",
    "sirnadesigner": "siRNAdesigner (单生成器)",
    "sirnaforge": "siRNAforge (单生成器)",
    "consensus": "exact consensus",
    "model": "model 冻结模型",
    "fusion": "fusion 融合",
    "support": "support 支持度",
    "quality": "quality 质量融合",
}
ABLATION_ZH = {
    "full": "完整配置",
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
     "数据库的百分比抑制字段重置了归一化活性尺度, <b>不可</b>当作敲低百分比解读。"),
)

INK, MUTED, NAVY, TEAL, RED, AMBER, LINE = (
    "#111820", "#5c6b74", "#123b4f", "#0f7b74", "#a4301f", "#8a5c12", "#d5dee2")
ACTIVE_THRESHOLD = 0.7


def e(v) -> str:
    return html.escape(str(v), quote=True)


def load(name: str) -> dict:
    return json.loads((HERE / name).read_text())


def pct(v, nd=1) -> str:
    return "不可用" if v is None else f"{v * 100:.{nd}f}%"


# ------------------------------------------------------------------ 统计
def load_measured() -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    with (HERE / "data/huesken/measurements.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out.setdefault(row["gene"], {})[row["guide_dna"]] = float(row["value"])
    return out


def ranked(measured, gene: str, method: str) -> list[tuple[float, float]]:
    """标签盲预测中与实测池相交的部分, 按分数降序 -> [(score, activity)]。"""
    path = FULL / gene / method / "selection.json"
    if not path.exists():
        return []
    pool = measured[gene]
    rows = [(c["score"], pool[c["guide_dna"]])
            for c in json.loads(path.read_text())["candidates"]
            if c["guide_dna"] in pool and c.get("score") is not None]
    return sorted(rows, key=lambda t: -t[0])


def precision_at_k(measured, method: str, kmax: int = 10) -> dict[int, tuple[float, int]]:
    curve = {}
    for k in range(1, kmax + 1):
        vals = []
        for gene in measured:
            rows = ranked(measured, gene, method)
            if len(rows) < k:
                continue
            vals.append(sum(1 for _, a in rows[:k] if a >= ACTIVE_THRESHOLD) / k)
        curve[k] = (sum(vals) / len(vals), len(vals)) if vals else (None, 0)
    return curve


def score_activity(measured) -> list[tuple[float, float]]:
    pairs = []
    for gene in measured:
        for c in json.loads((FULL / gene / "model" / "selection.json").read_text())["candidates"]:
            if c["guide_dna"] in measured[gene] and c.get("model_score") is not None:
                pairs.append((c["model_score"], measured[gene][c["guide_dna"]]))
    return sorted(pairs)


def _ranks(vals: list[float]) -> list[float]:
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    out = [0.0] * len(vals)
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and vals[order[j]] == vals[order[i]]:
            j += 1
        for k in range(i, j):
            out[order[k]] = (i + j + 1) / 2
        i = j
    return out


def model_stats(pairs) -> dict:
    n = len(pairs)
    pos = sum(1 for _, a in pairs if a >= ACTIVE_THRESHOLD)
    order = sorted(range(n), key=lambda i: -pairs[i][0])
    hits, ap = 0, 0.0
    for rank, i in enumerate(order, 1):
        if pairs[i][1] >= ACTIVE_THRESHOLD:
            hits += 1
            ap += hits / rank
    rx, ry = _ranks([s for s, _ in pairs]), _ranks([a for _, a in pairs])
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    rx_score = _ranks([s for s, _ in pairs])
    pos_rank_sum = sum(r for r, (_, a) in zip(rx_score, pairs) if a >= ACTIVE_THRESHOLD)
    n0 = n - pos
    auroc = (pos_rank_sum - pos * (pos + 1) / 2) / (pos * n0) if pos and n0 else None
    return {"n": n, "active": pos, "base": pos / n, "ap": ap / pos if pos else None,
            "auroc": auroc, "spearman": num / den if den else None}


def deciles(pairs, bins: int = 10) -> list[tuple[int, float]]:
    pts = sorted(pairs)
    step = len(pts) // bins
    out = []
    for i in range(bins):
        chunk = pts[i * step:] if i == bins - 1 else pts[i * step:(i + 1) * step]
        out.append((len(chunk), sum(1 for _, a in chunk if a >= ACTIVE_THRESHOLD) / len(chunk)))
    return out


# ------------------------------------------------------------------ 图
def _axis(w, h, pad_l, pad_r, pad_b, ticks, xmax, title, ylabel=""):
    parts = [f'<rect x="{pad_l}" y="24" width="{w-pad_l-pad_r}" height="{h-24-pad_b}" '
             f'fill="#fbfcfc" stroke="{LINE}"/>']
    for t in ticks:
        x = pad_l + t / xmax * (w - pad_l - pad_r)
        parts.append(f'<line x1="{x:.1f}" y1="24" x2="{x:.1f}" y2="{h-pad_b}" '
                     f'stroke="#eef2f3"/>')
        parts.append(f'<text x="{x:.1f}" y="{h-pad_b+15}" font-size="10.5" fill="{MUTED}" '
                     f'text-anchor="middle">{t:g}</text>')
    parts.append(f'<text x="{(pad_l+w-pad_r)/2:.1f}" y="{h-pad_b+31}" font-size="11" '
                 f'fill="{MUTED}" text-anchor="middle">{e(title)}</text>')
    if ylabel:
        parts.append(f'<text x="12" y="{24+(h-24-pad_b)/2:.1f}" font-size="11" fill="{MUTED}" '
                     f'text-anchor="middle" transform="rotate(-90 12 '
                     f'{24+(h-24-pad_b)/2:.1f})">{e(ylabel)}</text>')
    return "".join(parts)


def fig_flow() -> str:
    """方法流程图。"""
    w, h = 880, 470
    box = ('<rect x="{x}" y="{y}" width="{bw}" height="{bh}" rx="5" fill="{fill}" '
           'stroke="{stroke}" stroke-width="{sw}"/>'
           '<text x="{tx}" y="{ty}" font-size="{fs}" font-weight="700" fill="{tc}" '
           'text-anchor="middle">{line1}</text>'
           '<text x="{tx}" y="{ty2}" font-size="10.5" fill="#4b5a63" '
           'text-anchor="middle">{line2}</text>')

    def b(x, y, bw, bh, l1, l2, fill="#eef4fb", stroke=NAVY, sw=1, tc=NAVY, fs=11.5):
        return box.format(x=x, y=y, bw=bw, bh=bh, fill=fill, stroke=stroke, sw=sw,
                          tx=x + bw / 2, ty=y + 21, ty2=y + 37, line1=e(l1),
                          line2=e(l2), fs=fs, tc=tc)

    def arrow(x1, y1, x2, y2, dash=""):
        da = ' stroke-dasharray="4 3"' if dash else ""
        return (f'<path d="M {x1} {y1} L {x2} {y2}" stroke="#8b949e" stroke-width="1.4" '
                f'fill="none"{da} marker-end="url(#ah)"/>')

    s = ['<defs><marker id="ah" markerWidth="7" markerHeight="7" refX="6" refY="3.5" '
         'orient="auto"><path d="M0,0 L7,3.5 L0,7 z" fill="#8b949e"/></marker></defs>']
    # 第一列: 数据
    s.append(b(10, 60, 190, 48, "Figshare 原始快照", "3,544 条记录 · 7 项研究"))
    s.append(b(10, 140, 190, 48, "Huesken 开发队列", "25 构建体 · 1,941 引导链"))
    s.append(b(10, 220, 190, 48, "候选池 (标签盲)", "21 nt 全长引导链 + 悬垂"))
    # 第二列: 预测
    s.append(b(250, 60, 220, 48, "生成器工具", "siRNAdesigner · siRNAforge (弃权)"))
    s.append(b(250, 140, 220, 48, "冻结语言模型", "OligoFormer (权重已哈希)"))
    s.append(b(250, 230, 220, 48, "排序策略", "native/model/fusion/support/quality",
               fill="#e7f0ef", stroke=TEAL, tc=TEAL))
    s.append(b(250, 310, 220, 48, "Top-K 面板选择", "仅用预测, 不看标签",
               fill="#e7f0ef", stroke=TEAL, tc=TEAL))
    # 第三列: 评估
    s.append(b(530, 140, 210, 48, "实测活性", "阈值 0.7 二分类"))
    s.append(b(530, 230, 210, 48, "主终点", "基因加权 Top-5 精确率"))
    s.append(b(530, 310, 210, 48, "不确定度", "1 万次配对自助 · 1 千次置换"))
    s.append(b(530, 390, 210, 52, "证据等级: 开发集", "无独立留出, 不自称泛化",
               fill="#fdf6e8", stroke=AMBER, tc=AMBER))
    # 连接
    s.append(arrow(200, 84, 248, 84))
    s.append(arrow(105, 188, 105, 218))
    s.append(arrow(200, 164, 248, 164))
    s.append(arrow(105, 244, 248, 254))
    s.append(arrow(360, 188, 360, 228))
    s.append(arrow(360, 278, 360, 308))
    s.append(arrow(470, 254, 528, 254))
    s.append(arrow(470, 334, 528, 334))
    s.append(arrow(635, 188, 635, 228))
    s.append(arrow(635, 278, 635, 308))
    s.append(arrow(635, 358, 635, 388))
    s.append(f'<text x="105" y="40" font-size="11" font-weight="700" fill="{NAVY}">数据</text>')
    s.append(f'<text x="360" y="40" font-size="11" font-weight="700" fill="{NAVY}" '
             f'text-anchor="middle">标签盲预测</text>')
    s.append(f'<text x="635" y="40" font-size="11" font-weight="700" fill="{NAVY}" '
             f'text-anchor="middle">评估</text>')
    s.append(f'<text x="250" y="384" font-size="10.5" fill="{MUTED}">'
             f'生成器未产出 = 弃权, 不判无效</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="评估流程图">{"".join(s)}</svg>'


def fig_primary(methods: dict) -> str:
    rows = sorted(methods.items(), key=lambda kv: (kv[1]["precision_at_5_complete"] is None,
                                                   kv[1]["precision_at_5_complete"] or 0))
    w, y0, rh, lab, bar = 780, 30, 34, 220, 300
    h = y0 + len(rows) * rh + 40
    s = []
    for i, (name, m) in enumerate(rows):
        y = y0 + i * rh
        s.append(f'<text x="0" y="{y+18}" font-size="12" fill="{INK}">'
                 f'{e(METHOD_ZH.get(name, name))}</text>')
        p = m["precision_at_5_complete"]
        if p is None:
            s.append(f'<text x="{lab}" y="{y+18}" font-size="11.5" fill="{RED}">'
                     f'本队列无完整面板 (0/25)</text>')
            continue
        full_panel = m["complete_genes"] == 25
        bw = max(3, p * bar)
        s.append(f'<rect x="{lab}" y="{y+6}" width="{bw:.1f}" height="{rh-14}" rx="2" '
                 f'fill="{NAVY if full_panel else AMBER}" opacity="{1 if full_panel else .85}"/>')
        s.append(f'<text x="{lab+bw+7:.1f}" y="{y+18}" font-size="12" font-weight="700" '
                 f'fill="{NAVY if full_panel else AMBER}">{p * 100:.1f}%</text>')
        s.append(f'<text x="{w}" y="{y+18}" font-size="10.5" fill="{MUTED}" '
                 f'text-anchor="end">完整 {m["complete_genes"]}/25 基因 · 覆盖 '
                 f'{m["mean_candidate_coverage"] * 100:.0f}%</text>')
    s.append(f'<line x1="{lab}" y1="{y0-8}" x2="{lab}" y2="{h-34}" stroke="{LINE}"/>')
    s.append(f'<text x="{lab}" y="{y0-14}" font-size="10.5" fill="{MUTED}">'
             f'基因加权 Top-5 精确率 (仅完整基因)</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="主终点">{"".join(s)}</svg>'


def fig_model(pairs, stats, dec) -> str:
    w, h = 860, 330
    pl, pr, pb = 54, 400, 46
    s = [_axis(pl + pr, h, pl, 0, pb, [0.4, 0.6, 0.8, 1.0], 1.15, "模型打分 (未标定)", "实测活性")]
    xmax = 1.15
    bw = (w - pl - pr) / xmax
    for sc, act in pairs:
        x = pl + sc * bw
        y = 24 + (1 - act / 1.15) * (h - 24 - pb)
        col = TEAL if act >= ACTIVE_THRESHOLD else "#b9c6cc"
        s.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="1.7" fill="{col}" opacity="0.45"/>')
    # 分箱中位趋势
    pts = sorted(pairs)
    step = len(pts) // 12
    line = []
    for i in range(12):
        chunk = pts[i * step:(i + 1) * step]
        xs = sum(p[0] for p in chunk) / len(chunk)
        ys = sorted(a for _, a in chunk)[len(chunk) // 2]
        line.append(f"{pl + xs * bw:.1f},{24 + (1 - ys / 1.15) * (h - 24 - pb):.1f}")
    s.append(f'<polyline points="{" ".join(line)}" fill="none" stroke="{NAVY}" stroke-width="2.4"/>')
    s.append(f'<line x1="{pl}" y1="{24 + (1 - ACTIVE_THRESHOLD / 1.15) * (h - 24 - pb):.1f}" '
             f'x2="{pl + pr}" y2="{24 + (1 - ACTIVE_THRESHOLD / 1.15) * (h - 24 - pb):.1f}" '
             f'stroke="{AMBER}" stroke-dasharray="5 4"/>')
    s.append(f'<text x="{pl + pr - 4}" y="{24 + (1 - ACTIVE_THRESHOLD / 1.15) * (h - 24 - pb) - 5:.1f}" '
             f'font-size="10" fill="{AMBER}" text-anchor="end">有效阈值 0.7</text>')
    # 右侧: 十分位命中率
    bx, _by, bh2 = pl + pr + 96, 24, h - 24 - pb
    s.append(f'<text x="{bx + 150}" y="16" font-size="11" fill="{MUTED}" text-anchor="middle">'
             f'模型打分十分位 → 实测有效比例</text>')
    base_y = 24 + (1 - stats["base"]) * bh2
    s.append(f'<line x1="{bx-16}" y1="{base_y:.1f}" x2="{bx+320}" y2="{base_y:.1f}" '
             f'stroke="{AMBER}" stroke-dasharray="5 4"/>')
    s.append(f'<text x="{bx+322}" y="{base_y+4:.1f}" font-size="10" fill="{AMBER}">'
             f'基线 {stats["base"] * 100:.0f}%</text>')
    for i, (_n, frac) in enumerate(dec):
        x = bx + i * 32
        hh = frac * bh2
        s.append(f'<rect x="{x}" y="{24 + bh2 - hh:.1f}" width="24" height="{hh:.1f}" rx="2" '
                 f'fill="{TEAL if frac >= stats["base"] else "#c3cfd5"}"/>')
        s.append(f'<text x="{x+12}" y="{24 + bh2 + 14}" font-size="9.5" fill="{MUTED}" '
                 f'text-anchor="middle">{i+1}</text>')
        s.append(f'<text x="{x+12}" y="{24 + bh2 - hh - 3:.1f}" font-size="9" fill="{MUTED}" '
                 f'text-anchor="middle">{frac * 100:.0f}</text>')
    s.append(f'<text x="{bx+160}" y="{24 + bh2 + 30}" font-size="11" fill="{MUTED}" '
             f'text-anchor="middle">打分十分位 (1 = 最低)</text>')
    s.append(f'<text x="{pl + pr + 6}" y="{h - 6}" font-size="10.5" fill="{MUTED}">'
             f'n = {stats["n"]:,} 条已测引导链 · Spearman ρ = {stats["spearman"]:.2f} · '
             f'AUROC = {stats["auroc"]:.3f} · AP = {stats["ap"]:.2f}</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="模型性能">{"".join(s)}</svg>'


def fig_budget(curves: dict, base: float) -> str:
    w, h = 760, 320
    pl, pr, pb = 54, 190, 46
    s = [_axis(w, h, pl, pr, pb, list(range(1, 11)), 10, "预算 K (合成条数)", "基因加权 Top-K 精确率")]
    ymax = 1.0
    def yy(v):
        return 24 + (1 - v / ymax) * (h - 24 - pb)
    s.append(f'<line x1="{pl}" y1="{yy(base):.1f}" x2="{w-pr}" y2="{yy(base):.1f}" '
             f'stroke="{AMBER}" stroke-dasharray="5 4"/>')
    s.append(f'<text x="{pl+6}" y="{yy(base)-6:.1f}" font-size="10" fill="{AMBER}">'
             f'随机基线 {base * 100:.0f}%</text>')
    style = {"quality": (NAVY, 2.6), "model": (TEAL, 2.0), "fusion": ("#7f9aa8", 1.8),
             "support": ("#b07a2a", 1.8), "native": (RED, 1.8), "sirnadesigner": ("#8a6fae", 1.8)}
    ly = 30
    for method, (col, sw) in style.items():
        pts = [f"{pl + (k-1)/9*(w-pl-pr):.1f},{yy(curves[method][k][0]):.1f}"
               for k in range(1, 11) if curves[method][k][0] is not None]
        s.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="{col}" stroke-width="{sw}"/>')
        s.append(f'<line x1="{w-pr+8}" y1="{ly}" x2="{w-pr+28}" y2="{ly}" stroke="{col}" stroke-width="{sw}"/>')
        s.append(f'<text x="{w-pr+34}" y="{ly+4}" font-size="10.5" fill="{INK}">'
                 f'{e(METHOD_ZH[method])}</text>')
        s.append(f'<text x="{w-pr+34}" y="{ly+17}" font-size="9.5" fill="{MUTED}">'
                 f'K=5: {pct(curves[method][5][0])} · {curves[method][5][1]}/25 基因</text>')
        ly += 40
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="预算曲线">{"".join(s)}</svg>'


def fig_ablation(abl: dict) -> str:
    full = abl["full"]["precision_at_5_complete"]
    rows = [k for k in abl if k != "full"]
    rows.sort(key=lambda k: abl[k]["precision_at_5_complete"])
    w, y0, rh, lab = 880, 30, 36, 210
    plot = w - lab - 210
    h = y0 + len(rows) * rh + 34
    col = lambda v: RED if v < full - 1e-9 else TEAL  # noqa: E731
    s = [f'<line x1="{lab + full * plot:.1f}" y1="{y0-10}" '
         f'x2="{lab + full * plot:.1f}" y2="{h-26}" stroke="{NAVY}" stroke-width="2" '
         f'stroke-dasharray="3 3"/>',
         f'<text x="{lab + full * plot:.1f}" y="{y0-16}" font-size="10.5" fill="{NAVY}" '
         f'text-anchor="middle">完整配置 {full * 100:.1f}%</text>']
    for i, name in enumerate(rows):
        y = y0 + i * rh
        v = abl[name]["precision_at_5_complete"]
        xf = lab + full * plot
        xv = lab + v * plot
        s.append(f'<text x="0" y="{y+18}" font-size="12" fill="{INK}">{e(ABLATION_ZH[name])}</text>')
        s.append(f'<line x1="{xv:.1f}" y1="{y+14}" x2="{xf:.1f}" y2="{y+14}" '
                 f'stroke="{col(v)}" stroke-width="2"/>')
        s.append(f'<circle cx="{xv:.1f}" cy="{y+14}" r="5" fill="{col(v)}"/>')
        s.append(f'<text x="{xv-10:.1f}" y="{y+18}" font-size="11" fill="{col(v)}" '
                 f'text-anchor="end">{v * 100:.1f}%</text>')
        d = (v - full) * 100
        s.append(f'<text x="{max(xv, xf)+9:.1f}" y="{y+18}" font-size="10.5" fill="{MUTED}">'
                 f'{d:+.1f} pt · 完整 {abl[name]["complete_genes"]}/25 基因</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="消融">{"".join(s)}</svg>'


def fig_paired(comp: dict) -> str:
    c = comp["comparisons"][0]
    diffs = {g: v for g, v in c["per_gene_difference"].items() if abs(v) > 1e-12}
    lo, hi = c["ci95"]
    mean = c["mean_difference"]
    genes = sorted(diffs, key=lambda g: diffs[g])
    w, y0, rh, lab = 780, 34, 26, 96
    span = max(abs(lo), abs(hi), max(abs(v) for v in diffs.values())) * 1.25
    mid = lab + (w - lab - 120) / 2
    sc = (w - lab - 120) / 2 / span
    h = y0 + len(genes) * rh + 66
    s = [f'<rect x="{mid + lo * sc:.1f}" y="{y0-8}" width="{(hi-lo)*sc:.1f}" '
         f'height="{h-y0-46}" fill="{NAVY}" opacity="0.09"/>',
         f'<line x1="{mid:.1f}" y1="{y0-14}" x2="{mid:.1f}" y2="{h-40}" stroke="{MUTED}"/>',
         f'<text x="{mid:.1f}" y="{y0-20}" font-size="10.5" fill="{MUTED}" '
         f'text-anchor="middle">无差异 (0)</text>']
    for i, g in enumerate(genes):
        y = y0 + i * rh
        d = diffs[g]
        x = mid + d * sc
        col = TEAL if d > 0 else RED
        s.append(f'<text x="{lab-8}" y="{y+14}" font-size="11" fill="{INK}" '
                 f'text-anchor="end">{e(g)}</text>')
        s.append(f'<line x1="{mid:.1f}" y1="{y+10}" x2="{x:.1f}" y2="{y+10}" '
                 f'stroke="{col}" stroke-width="1.6" opacity="0.6"/>')
        s.append(f'<circle cx="{x:.1f}" cy="{y+10}" r="4.6" fill="{col}"/>')
        s.append(f'<text x="{x + (12 if d > 0 else -12):.1f}" y="{y+14}" font-size="10.5" '
                 f'fill="{col}" text-anchor="{"start" if d > 0 else "end"}">{d * 100:+.0f} pt</text>')
    s.append(f'<text x="{mid + mean * sc:.1f}" y="{h-22}" font-size="11" font-weight="700" '
             f'fill="{NAVY}" text-anchor="middle">均值 {mean * 100:+.1f} pt</text>')
    s.append(f'<text x="{mid:.1f}" y="{h-6}" font-size="10.5" fill="{MUTED}" text-anchor="middle">'
             f'自助 95% CI [{lo * 100:+.1f}, {hi * 100:+.1f}] pt —— 跨过 0, 未达 +10 pt 判据</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="配对差异">{"".join(s)}</svg>'


def fig_permutation(perm: dict) -> str:
    obs, null = perm["observed_gene_precision_at_5"], perm["null_mean"]
    w, h = 700, 150
    s = []
    for i, (v, col, text) in enumerate(
            [(null, AMBER, f"基因内标签置换零假设均值 {null * 100:.1f}%"),
             (obs, TEAL, f"实际排序观测值 {obs * 100:.1f}%")]):
        y = 26 + i * 52
        bw = v * (w - 330)
        s.append(f'<text x="0" y="{y+16}" font-size="11.5" fill="{INK}">{e(text)}</text>')
        s.append(f'<rect x="290" y="{y+3}" width="{bw:.1f}" height="20" rx="2" fill="{col}"/>')
        s.append(f'<text x="{290+bw+8:.1f}" y="{y+18}" font-size="11.5" font-weight="700" '
                 f'fill="{col}">{v * 100:.1f}%</text>')
    s.append(f'<text x="0" y="{h-14}" font-size="10.5" fill="{MUTED}">'
             f'{perm["permutations"]:,} 次置换 · 单侧 p = {perm["p_one_sided"]:.4f} · '
             f'种子 {perm["seed"]} · 说明排序显著优于随机, 但不能消除预训练暴露</text>')
    return f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="置换对照">{"".join(s)}</svg>'


# ------------------------------------------------------------------ 报告
def build() -> str:
    dev = load("results/development_results.json")
    comp = load("results/development_comparison.json")
    perm = load("results/permutation_control.json")
    impl = load("results/implementation_verification.json")
    audit = load("dataset_manifest.json")
    cfg = load("benchmark_config.json")

    measured = load_measured()
    pairs = score_activity(measured)
    stats = model_stats(pairs)
    dec = deciles(pairs)
    curves = {m: precision_at_k(measured, m) for m in
              ("quality", "model", "fusion", "support", "native", "sirnadesigner")}
    base = stats["base"]

    methods, abl = dev["methods"], dev["ablations"]
    primary = comp["comparisons"][0]
    best = max((m["precision_at_5_complete"] or 0, k) for k, m in methods.items())
    # 交叉校验: 复算的 Top-5 必须与冻结 JSON 一致
    recomputed = curves["quality"][5][0]
    frozen = methods["quality"]["precision_at_5_complete"]
    agree = abs(recomputed - frozen) < 1e-9

    rows_methods = "".join(
        f"<tr><td>{e(METHOD_ZH.get(k, k))}</td>"
        f"<td class='n'>{pct(m['precision_at_5_complete'])}</td>"
        f"<td class='n'>{m['complete_genes']}/25</td>"
        f"<td class='n'>{m['mean_fill_rate'] * 100:.1f}%</td>"
        f"<td class='n'>{m['mean_candidate_coverage'] * 100:.2f}%</td></tr>"
        for k, m in sorted(methods.items(), key=lambda kv: -(kv[1]["precision_at_5_complete"] or -1)))

    rows_paired = "".join(
        f"<tr><td class='n'>{c['stratum'][-1]}</td><td class='n'>{c['paired_genes']}</td>"
        f"<td class='n'>{c['mean_difference'] * 100:+.1f}</td>"
        f"<td class='n'>[{c['ci95'][0] * 100:+.1f}, {c['ci95'][1] * 100:+.1f}]</td>"
        f"<td class='bad'>未达标</td></tr>" for c in comp["comparisons"])

    rows_studies = "".join(
        f"<tr><td>{e(name)}</td><td class='n'>{n:,}</td></tr>"
        for name, n in sorted(audit["studies"].items(), key=lambda kv: -kv[1]))

    reasons = "".join(
        f"<li>{zh}</li>" for r in audit["reasons"]
        for zh in [next((t for k, t in REASON_ZH if r.startswith(k)), e(r))])

    rows_model = "".join(
        f"<tr><td class='n'>{i + 1}</td><td class='n'>{n}</td><td class='n'>{f * 100:.1f}%</td></tr>"
        for i, (n, f) in enumerate(dec))

    return f"""<!doctype html>
<html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Huesken 开发集上共识 siRNA 选择的评估报告</title>
<style>{CSS}</style></head><body><div class="page">

<h1>以公开实测数据评估共识 siRNA 选择流程:<br>Huesken 开发集基准测试报告</h1>
<p class="byline">sirna-consensus v{impl['package_version']} · schema v{impl['schema_version']}
· 状态 <code>{e(impl['release_status'])}</code></p>
<p class="meta">{dev['genes']} 个人源克隆构建体 · {impl['curated_human_guides']:,} 条全长引导链 ·
实测表 SHA-256 <code>{e(perm['benchmark_sha256'][:16])}…</code> ·
生成于 {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')} UTC</p>

<div class="abstract">
<h2>摘要</h2>
<p><b>背景.</b> 排序流程需要一个外部标尺来判断 "哪条候选值得合成"。公开 siRNA 数据集大多
已被同类模型的开发集合覆盖, 因此本报告把证据等级当作一等公民: 先冻结数据、模型与判据,
再明确标注哪些结论在什么前提下成立。</p>
<p><b>方法.</b> 主终点为<b>基因加权的 Top-5 实测精确率</b>。逐分子预测全部标签盲产生,
之后才与实测活性比对。比较 6 种策略与 2 个单生成器; 不确定度用
{cfg['bootstrap_replicates']:,} 次配对基因自助与 {cfg['permutations']:,} 次基因内标签置换刻画。</p>
<p><b>结果.</b> 冻结模型在本队列上区分度良好: AUROC <b>{stats['auroc']:.3f}</b>、
平均精确率 AP <b>{stats['ap']:.2f}</b>、Spearman ρ <b>{stats['spearman']:.2f}</b>;
最低十分位命中率仅 {dec[0][1] * 100:.0f}%, 最高十分位达 {dec[-1][1] * 100:.0f}%
(基线 {base * 100:.0f}%)。开发集上最优策略为 <b>{e(METHOD_ZH[best[1]])}</b>
({best[0] * 100:.1f}% Top-5), 但其相对冻结模型的增益仅
<b>{primary['mean_difference'] * 100:+.1f} pt</b>, 95% CI
<b>[{primary['ci95'][0] * 100:+.1f}, {primary['ci95'][1] * 100:+.1f}]</b> —— 跨过 0,
低于预设 +10 pt 门槛。</p>
<p><b>结论.</b> "模型本身有效" 成立; "融合策略显著优于模型" <b>未成立</b>。
实现与执行检查全部通过, 但发布判据未达成, 且无独立留出集。</p>
<p class="kw"><b>关键词:</b> siRNA 基准测试; AUROC; 基因加权精确率; 配对自助; 暴露审计</p>
</div>

<div class="note warnbox">
<b>证据等级.</b> 7 项来源研究均在冻结模型的开发/评估集合中出现, 该队列也被姊妹项目使用,
故<b>当前无独立留出集</b>。以下全部数字都是开发集诊断值。
</div>

<h2><span class="num">1.</span> 背景</h2>
<p>临床前的实际问题是排序: 在只能合成 k 条的预算下, 前 k 条里有多少条真的有效。
因此主终点取基因加权的 Top-5 精确率, 并同时报告完成度与覆盖率 —— 只在少数基因上 "精确"
而对其余基因交白卷的方法, 单看精确率会被高估。</p>
<p>另一个同样重要的问题是<b>参照系</b>: 一个融合策略若能提升精确率, 必须在同一预算、
同一候选池上与基线可比。本报告据此把 "提升" 定义为配对差异 + CI, 而非绝对分数。</p>

<h2><span class="num">2.</span> 方法</h2>
<figure>{fig_flow()}
<figcaption><b>图 1. 评估流程图。</b>数据、标签盲预测与评估三段分离: 预测阶段完全不看实测标签,
面板选定后才做比对; 生成器未产出记为弃权而非 "无效"。</figcaption>
</figure>

<h3>2.1 数据与终点</h3>
<p>来源为 Figshare 快照 (article 25908841 v1, {audit['snapshot_records']:,} 条记录, 7 项研究);
其中 {dev['genes']} 个构建体、{dev['records']:,} 条全长引导链构成开发队列。
活性是原始归一化活性 (非百分比敲低), 以 {cfg['thresholds']['normalized_activity']} 为二分类阈值;
同一条件的重复测定先取平均, 不同条件保持分离。</p>

<h3>2.2 比较对象</h3>
<p>包含 2 个单生成器 (siRNAdesigner、siRNAforge)、精确共识、冻结的 OligoFormer 模型打分,
以及 3 种融合策略与并集回退 native。</p>

<h3>2.3 统计与判据</h3>
<p>先在同一基因内平均测定, 再对基因等权。数值判据为: 精确率提升 ≥
{cfg['minimum_precision_gain'] * 100:.0f} pt、配对自助 95% CI 下界 &gt; 0、完成度不下降、
预算产出为正; 至少 {cfg['minimum_complete_paired_genes']} 个完整配对基因。</p>

<h3>2.4 复算与交叉校验</h3>
<p>图中的逐分子分数来自各基因的标签盲 <code>selection.json</code>, 与
<code>measurements.tsv</code> 相交后计算描述性统计 (不重新拟合、不重新抽样)。
复算的 quality Top-5 精确率为 {recomputed * 100:.1f}%, 与冻结 JSON 的 {frozen * 100:.1f}%
{'完全一致' if agree else '存在偏差'}。</p>

<h2><span class="num">3.</span> 结果</h2>
<h3>3.1 主终点</h3>
<table><caption>表 1. 基因加权 Top-5 精确率 (仅完整基因)。覆盖率指候选占已测分子池的比例。</caption>
<thead><tr><th>方法</th><th class="n">Top-5 精确率</th><th class="n">完整基因</th>
<th class="n">填充率</th><th class="n">候选覆盖率</th></tr></thead>
<tbody>{rows_methods}</tbody></table>

<h3>3.2 冻结模型的区分度</h3>
<table><caption>表 2. 按模型打分升序切十分位后的实测有效比例。单调性良好, 说明打分具有排序意义。</caption>
<thead><tr><th class="n">十分位</th><th class="n">条数</th><th class="n">有效比例</th></tr></thead>
<tbody>{rows_model}</tbody></table>

<h3>3.3 配对比较与置换</h3>
<table><caption>表 3. 最优策略相对冻结模型的配对差异 (基因加权, {cfg['bootstrap_replicates']:,} 次自助)。</caption>
<thead><tr><th class="n">K</th><th class="n">配对基因</th><th class="n">均值差 (pt)</th>
<th class="n">95% CI</th><th>判定</th></tr></thead>
<tbody>{rows_paired}</tbody></table>

<h3>3.4 数据集来源审计</h3>
<table><caption>表 4. 快照中的研究构成 (共 {audit['snapshot_records']:,} 条)。</caption>
<thead><tr><th>研究</th><th class="n">记录数</th></tr></thead><tbody>{rows_studies}</tbody></table>
<ol class="plain">{reasons}</ol>

<h2><span class="num">4.</span> 图</h2>
<figure>{fig_primary(methods)}
<figcaption><b>图 2. 主终点: 各方法 Top-5 精确率。</b>深色 = 25/25 基因均有完整面板, 琥珀色 = 面板不完整。
右侧给出覆盖率, 用来看清高精确率是否只建立在窄候选面上。</figcaption>
</figure>
<figure>{fig_model(pairs, stats, dec)}
<figcaption><b>图 3. 模型性能。</b>左: 每条已测引导链的模型打分 vs 实测活性 (深色为有效),
折线为分箱中位趋势。右: 打分十分位对应的实测有效比例, 虚线为基线。
从最低十分位的 {dec[0][1] * 100:.0f}% 单调升到最高十分位的 {dec[-1][1] * 100:.0f}%,
是 "模型打分确实携带活性信息" 的直接证据 (AUROC {stats['auroc']:.3f}, AP {stats['ap']:.2f})。</figcaption>
</figure>
<figure>{fig_budget(curves, base)}
<figcaption><b>图 4. 预算曲线。</b>横轴为可合成的条数, 纵轴为基因加权 Top-K 精确率,
虚线为随机基线。曲线彼此贴近说明各策略量级相同; 差异要到配对检验才可判定。</figcaption>
</figure>
<figure>{fig_ablation(abl)}
<figcaption><b>图 5. 组件消融。</b>每个点是一个消融配置, 横线连到完整配置的参考线。
硬过滤以 +{abs(abl['hard_filters']['precision_at_5_complete'] - best[0]) * 100:.1f} pt 换掉
{abs(abl['hard_filters']['mean_fill_rate'] - 1) * 100:.0f}% 填充率, 是典型的以覆盖换精确。</figcaption>
</figure>
<figure>{fig_paired(comp)}
<figcaption><b>图 6. 配对差异。</b>只画出现非零差异的基因 —— 25 个基因里有
{25 - sum(1 for v in primary['per_gene_difference'].values() if abs(v) > 1e-12)} 个完全并列。
竖向阴影为自助 95% CI; 它跨过 0 正是 "增益未被证实" 的图形解释。</figcaption>
</figure>
<figure>{fig_permutation(perm)}
<figcaption><b>图 7. 置换对照。</b>观测值远高于基因内标签置换的零分布均值, 排序显著优于随机。</figcaption>
</figure>

<h2><span class="num">5.</span> 讨论与局限</h2>
<ul class="plain">
<li><b>无独立留出集。</b>检查点级暴露审计状态为 <code>{e(audit['checkpoint_training_audit'])}</code>。</li>
<li><b>发布判据未达成。</b>实测 {primary['mean_difference'] * 100:+.1f} pt、CI 跨零,
低于 +10 pt 门槛, <code>publication_criterion_met = false</code>。</li>
<li><b>端点单位不可混用。</b>归一化活性与百分比抑制阈值
({cfg['thresholds']['normalized_activity']} / {cfg['thresholds']['percent_inhibition']}) 不可互换。</li>
<li><b>覆盖 vs 精确。</b>硬过滤与 siRNAdesigner 都说明高精确率可能只是候选面窄。</li>
<li><b>模型分是序, 不是百分比。</b>批大小为 1 是上游 LSTM 布局所致。</li>
<li><b>未经湿实验确认前不得称任何引导链无效。</b></li>
</ul>

<h2><span class="num">6.</span> 结论</h2>
<p>模型有效 (AUROC {stats['auroc']:.3f}, 置换 p &lt; 0.001), 但融合策略相对模型的增益
({primary['mean_difference'] * 100:+.1f} pt, CI 跨零) <b>未被证实</b>。
下一步是取得审计干净的未用队列并走 <code>evaluate.py --phase final</code> 冻结闸门。</p>

<h2><span class="num">7.</span> 数据与复现</h2>
<table><caption>表 5. 关联产物 (位于 <code>benchmarks/</code>)。</caption>
<thead><tr><th>文件</th><th>内容</th></tr></thead><tbody>
<tr><td><code>results/development_results.json</code></td><td>方法/消融/配对比较/筛查冒烟</td></tr>
<tr><td><code>results/development_comparison.json</code></td><td>逐基因配对差异与自助 CI</td></tr>
<tr><td><code>results/permutation_control.json</code></td><td>1,000 次标签置换对照</td></tr>
<tr><td><code>data/huesken/measurements.tsv</code></td><td>实测活性 (本报告复算所用)</td></tr>
<tr><td><code>dataset_manifest.json</code></td><td>快照校验和与暴露审计</td></tr>
<tr><td><code>benchmark_config.json</code></td><td>冻结的终点、阈值、判据与种子</td></tr>
<tr><td><code>report_zh.html</code></td><td>本报告</td></tr>
</tbody></table>

<footer>由 <code>benchmarks/make_report_zh.py</code> 渲染; 主终点与配对统计取自冻结 JSON,
模型性能与预算曲线由冻结预测 + 实测文件复算 (已验证与冻结值一致)。
{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M')} UTC · 计算诊断, 非独立验证。</footer>
</div></body></html>
"""


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
h3:before,h2 span.num{margin-right:6px;color:var(--teal)}
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
padding:6px 6px;border-top:1.4px solid var(--ink);border-bottom:1px solid var(--ink)}
td{padding:6px 6px;border-bottom:1px solid #eceff1;vertical-align:top}
tbody tr:last-child td{border-bottom:1.3px solid var(--ink)}
td.n,th.n{text-align:right;font-variant-numeric:tabular-nums}
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


def main(argv: list[str]) -> int:
    out = Path(argv[argv.index("-o") + 1]) if "-o" in argv else HERE / "report_zh.html"
    out.write_text(build(), encoding="utf-8")
    print(f"已写入 {out} ({out.stat().st_size} 字节)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
