"""Current strict-round section; every performance number comes from frozen artifacts."""

from __future__ import annotations

import html
import json
import os
from pathlib import Path

from sirna_consensus.fasta import sha256_file


def section(package):
    package = Path(package)
    summary_path = package / "benchmarks/confirmation_20261106/summary.json"
    if not summary_path.exists():
        return ""

    def read(p):
        return json.loads(Path(p).read_text())

    summary = read(summary_path)
    verification = read(summary["verification"])
    evaluation = read(summary["evaluation"])
    decision = read(summary["decision"])
    model = Path(summary["model_dir"])
    development = read(model / "results.json")
    dev_verify = read(model / "independent_verification.json")
    if not (verification["all_passed"] and dev_verify["all_passed"]):
        raise ValueError("report requires passing arithmetic verification")
    if verification["evaluation_sha256"] != sha256_file(Path(summary["evaluation"])) or dev_verify[
        "results_sha256"
    ] != sha256_file(model / "results.json"):
        raise ValueError("stale report verification")
    if decision["promote"] or evaluation["candidate_metrics"] is not None:
        raise ValueError("zero-ID round cannot report confirmed performance")
    esc = html.escape

    def link(path, label):
        return (
            '<a href="'
            + esc(os.path.relpath(path, package / "benchmarks"), quote=True)
            + '">'
            + esc(label)
            + "</a>"
        )

    def table(headers, rows):
        return (
            '<div class="scroll"><table><tr>'
            + "".join("<th>" + esc(v) + "</th>" for v in headers)
            + "</tr>"
            + "".join(
                "<tr>" + "".join("<td>" + esc(str(v)) + "</td>" for v in row) + "</tr>"
                for row in rows
            )
            + "</table></div>"
        )

    out = [
        "<h1>siRNA de novo 平台：严格预注册与 Huesken-only 开发</h1>",
        "<p>2026-11-06 · 当前轮次 · 主指标 macro_gene_assay_NDCG@5</p>",
        '<div class="notice"><strong>当前结论：独立确认样本为 0，尚未证实性能提升；保留冻结的默认模型。</strong>'
        "<p>已暴露、训练重叠或来源不清楚的队列已移出本轮独立测试资格，原始数据与历史结果仍保留。算术核查通过只证明计算及阻断行为正确，不构成效力验证。</p></div>",
        "<h2>冻结、开发、确认与回退流程</h2>",
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 170" role="img" aria-label="基线与协议冻结，Huesken开发，候选冻结，独立资格检查，无合格数据，保留默认模型" style="width:100%;min-width:600px">'
        '<defs><marker id="arrowStrict" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0L10 5L0 10z" fill="#278073"/></marker></defs>',
    ]
    labels = [
        "冻结基线／协议",
        "Huesken 嵌套分组开发",
        "冻结候选与独立重算",
        "独立资格：0 条",
        "门槛未满足：保留默认",
    ]
    for i, label in enumerate(labels):
        x = 5 + 200 * i
        out.append(
            f'<rect x="{x}" y="35" width="185" height="70" rx="8" fill="#eaf3f3" stroke="#278073"/><text x="{x + 92}" y="75" text-anchor="middle" font-size="13">{label}</text>'
        )
        if i < 4:
            out.append(
                f'<path d="M{x + 185} 70h15" stroke="#278073" marker-end="url(#arrowStrict)"/>'
            )
    out.append("</svg><h2>当前独立确认结果</h2>")
    out.append(
        table(
            ["指标", "候选", "冻结基线"],
            [
                [m, "未确认（无合格独立样本）", "未确认"]
                for m in (
                    "macro_gene_assay_NDCG@5",
                    "precision@5",
                    "Top3 至少一个活性候选",
                    "Top5 活性比例",
                    "面板完成率",
                    "预测覆盖／弃权",
                    "真实生成 Top5 预算收益及其界限",
                )
            ],
        )
    )
    out.append(
        "<p>未将未测试序列标为无效；不按重叠区或相同 19mer 转移全双链活性。实际生成面板收益只接受同实验条件下完整双链、突出端与已知化学修饰的精确匹配；其余保留未知。两面板共享的同分子未知项可配对抵消，其余按最坏情况给出收益界限。</p>"
    )
    out.append("<h2>队列资格与排除原因</h2>")
    out.append(
        table(
            ["队列", "本地记录数", "本轮用途", "资格审计原因"],
            [
                [
                    c["name"],
                    c.get("records", 0),
                    c["role"],
                    "; ".join(c.get("exclusion_reasons", [])),
                ]
                for c in summary["cohorts"]
            ],
        )
    )
    out.append(
        "<p>各行可能有重叠，不能相加当作独立样本量。Harborth 的本地映射表为空；OligoFormer 文件是模型登记，不是实验测试队列。Huesken-only 训练无法执行真正的多研究 LOSO；外部逐研究评分不能冒称 LOSO。本轮没有下载新数据，也没有给已查看结果重新赋予“独立”资格。</p>"
    )
    out.append("<h2>Huesken 开发诊断：不能用于默认模型晋升</h2>")
    out.append(
        f"<p>{development['records']} 条记录，{development['genes']} 个基因，{development['components']} 个基因／近重复序列连通组；5 外折、3 内折、32 个预注册模型配置。特征预处理仅在训练折拟合，模型选择仅使用 Huesken。最终配置：<code>{esc(json.dumps(development['selected_params']))}</code>。</p>"
    )
    out.append(
        table(
            ["开发折外指标", "值"],
            [
                [k, "未定义" if v is None else f"{v:.6f}"]
                for k, v in development["metrics"].items()
                if isinstance(v, (int, float)) or v is None
            ],
        )
    )
    out.append(
        "<p>这些数值是开发诊断，不能证明新生成序列有效，也未与历史混合训练模型的不同样本划分作提升比较。</p>"
    )
    out.append(
        "<h2>晋升门槛与失败记录</h2><p>须同时满足：至少 2 个独立研究与 5 个独立连通组、NDCG@5 &gt; 0.80、相对冻结基线增益 &gt; 0、配对基因／序列组 bootstrap 95% CI 下界 &gt; 0、完成率与覆盖无损失、实际生成 Top5 保守预算收益增量 &gt; 0，并通过独立重算。本轮未满足，因此不替换默认。</p>"
    )
    out.append("<p>未通过门槛：<code>" + esc(", ".join(summary["failed_gates"])) + "</code>。</p>")
    out.append(
        "<p>基线提交：<code>"
        + summary["baseline_commit"]
        + "</code>；协议 SHA256：<code>"
        + summary["protocol_sha256"]
        + "</code>。预注册记录保存在本地 Git；未声称具备第三方时间戳，亦不能排除未记录的历史暴露。</p>"
    )
    for path, label in [
        (summary_path, "当前轮次摘要"),
        (package / "benchmarks/confirmation_20261106/protocol.json", "冻结协议"),
        (summary["evaluation"], "evaluate.py 确认结果"),
        (summary["verification"], "verify_platform_v3.py 独立核查"),
        (summary["decision"], "保留默认／回退记录"),
        (summary["development_verification"], "开发算术核查"),
        (model / "results.json", "开发完整指标"),
        (package / "docs/guarded_confirmation.md", "复现命令与限制"),
    ]:
        out.append("<p>" + link(path, label) + "</p>")
    out.append(
        "<hr><h2>以下为冻结基线的历史开发分析</h2><p><strong>下文“四研究”“新模型”和增益均指旧轮次的混合开发模型，不是当前 Huesken-only 候选的独立确认结论。</strong></p>"
    )
    return "".join(out)
