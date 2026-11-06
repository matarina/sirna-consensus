# 研发记录总报告

稳定入口：`../three_benchmark_report_zh.html`。离线浏览无需 CDN。原报告已经按内容 SHA256 归档，历史冻结文件和数值工件不改写。

## 内容

- Git / 实验时间线、Sciabola 独立纠错及未晋升记录。
- 原共识路线与 v3 完整双链候选生成的流程与模型谱系。
- 可筛选的 15 个队列/来源条目，区分训练重叠、历史暴露、无法评估和新外部研究。
- 同池指标比较、历史开发 OOF 逐基因热图、工程特异性筛查和生成覆盖缺口。
- Katoh ROC / PR（SVG / PNG）与 AUC / AP，基于冻结的原始分数，固定归一化活性阈值≥0.7。
- 失败原因、默认替换门槛、后续研发和复现链接。

ROC / PR 为报告整合时新增的事后分析，没有加入原预测前协议，没有用于拟合或默认替换；单靶标不给出伪造的行级泛化 CI。AUC 使用阈值分组 ROC 梯形积分，并以独立正负样本对比较复算；AUC 与 AP 再用 scorer 环境的 sklearn 独立核对。AP 不等同于插值 PR 梯形面积。

## 重建与校验

在仓库根目录运行，使用已有环境：

```bash
python benchmarks/scripts/plot_journal_curves.py
python benchmarks/scripts/build_development_journal_zh.py
python benchmarks/scripts/verify_development_journal.py
```

采集 Python 需要现有 matplotlib、numpy、beautifulsoup4；指标核对调用现有 `.scorer-venv/bin/python`。不重新训练或重跑基准。`sources.json` 记录输入数值工件、模板、生成器和输出哈希。图表数值保存在 `curve_metrics.json`，原始 Sciabola 纠错结果的哈希和摘要保存在 `sciabola_adaptation_snapshot.json`，避免历史报告构建必须依赖仓库外文件。

`historical_report_preservation.json` 显式指定四份旧冻结文件的历史字节映射。`verify_development_journal.py` 在映射下验证原有 1849 / 85 / 593 / 163 个条目，不能把更新后的展示 HTML 当成过去冻结的原报告。旧协议、模型和实验结果没有变化；当前报告与原报告的字节当然不同。

当前页面提供 SVG 指标图、原生 JS 搜索/筛选/热图、静态 PNG/SVG 曲线。所有核心数值同时保留表格，便于审查和打印。

## 在线展示

共享的用户 systemd 服务 `report-http-8788.service` 监听 `0.0.0.0:8788`，服务目录通过符号链接映射到报告源文件，因此报告重建后直接展示新内容。

- 主动学习： http://10.0.11.40:8788/active_learning_report.html
- siRNA 综合研发记录： http://10.0.11.40:8788/benchmarks/three_benchmark_report_zh.html

两个入口 HTTP 200 且内容与磁盘文件一致，49 个本地链接/资源均 HTTP 200；结果记录于 `http_verification.json`。此检查未重建主动学习报告或改变其模型。
