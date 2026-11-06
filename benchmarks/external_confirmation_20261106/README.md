# 新外部实验评估：Katoh 2007

此前 0.7868 是 Huesken 内部开发 CV。本轮另行采集 Katoh 原始 Table S1，先冻结模型、原始来源、全部 702 IDs 和评估代码，再首次执行两个已冻结模型的外部预测。没有拟合、选参或默认替换。

## 结果与适用范围

- 默认模型：实测配对核心 NDCG@5 0.8291167145；Huesken-only 候选 0.8696036315。
- 简单规则 0.6722038098；均匀随机解析期望 0.4441065160。
- 两模型 Precision@5=0.8，Top3 至少一条达归一化活性阈值=1。阈值为研究内归一化活性 0.7，不是绝对 70% knockdown。
- 实际完整双链生成：各 698 个候选，默认 diverse panel 各 5/5；严格完整分子实验标签覆盖 0/5，活性全部未知，不做核心重叠活性转移。
- 原表缺少完整双链身份，生产 full-molecule adapter 对已测核心的支持 0/702；实际核心输入的 native/quality fallback 对应 0.6722。0.8291/0.8696 只能作为配对核心模型组件结果。
- 一研究、一 reporter，不做行级 bootstrap 冒充独立基因 CI。预算产出未确认；保留原默认。

## 独立性披露

实际训练集不含 Katoh / EGFP。对两个模型实际训练记录、历史 benchmark 和 Hu/Mix 公开参考（10,928 条，含重复来源记录）逐条查 Hamming≤2 双向/RC 及位移 17-mer，共 0 匹配，702 条全保留。没有按标签或模型分数筛选。

Taka.csv 以前已下载作序列查重参考；EGFP 出现在既往 Martinelli 修饰分子测试。因此这不是“项目从未接触过的盲测靶标”。没有发现 Katoh 曾用于两个冻结模型的拟合或调参；本轮仅验证这两个已知训练谱系的冻结模型，不能声称对训练过 Taka 的上游 OligoFormer 同样独立。

今后若使用这些结果选择特征/超参数，该表成为开发暴露数据，必须另取未使用研究确认改进。Jagla 表未恢复、ACON 请求获取、Davis 修饰域未运行，均不虚构第二项结果。

## 冻结顺序

1. 采集前协议：commit `3dabca9`，`acquisition_protocol.json`。
2. 原始 XLS/DOC、ID、近重复审计、指标、评分/验证代码冻结：commit `1131d5a`。
3. Git 忽略的已测核心 TSV 补入版本控制：commit `0aa3855`；仍早于首次预测。
4. 本轮结果与报告为后续提交；不修改上述协议、生产模型或旧轮次报告。

`test_protocol.json` 含 54 个输入/代码/模型/历史参考 SHA256。原始补充文件 URL、内容哈希、字节数和失败文件隔离记录见 `sources.json`。Git 中保存原始输入；复现不需要再次下载。`source_verification.json` 是预测后增加的独立源数据 QA，并非预测前另一次性能选择。

## 复现

在仓库根目录执行；已有采集 Python 提供 xlrd，`.scorer-venv` 提供锁定的 sklearn/numpy/scipy。没有新安装依赖。

```bash
python benchmarks/scripts/prepare_external_20261106.py
PYTHONPATH=src .scorer-venv/bin/python benchmarks/scripts/evaluate_external_20261106.py
PYTHONPATH=src .scorer-venv/bin/python benchmarks/scripts/verify_external_20261106.py
python benchmarks/scripts/verify_external_source_20261106.py
python benchmarks/scripts/render_external_report_20261106.py
PYTHONPATH=src .scorer-venv/bin/python benchmarks/scripts/verify_freezes.py --historical-code benchmarks/confirmation_20261106/historical_snapshots.json
```

评估脚本在每次评分前检查所有 54 项冻结哈希；固定 ID、请求哈希、checkpoint 哈希与结果关联。验证脚本独立实现/读取既有独立 `verify_guarded.py` 的特征和指标算术，不调用生产特征函数或生产 NDCG 函数；源审计另用分块索引查近重复，区别于生产向量距离。查重只读序列和谱系。

旧 `verify_platform_v3.py --evaluation` / `evaluate.py --phase confirmation` 针对原先零-ID 确认协议，不能将本轮新增 ID 塞入旧冻结协议。本轮使用单独注册的 `evaluate_external_20261106.py` 与独立 `verify_external_20261106.py`，不声称旧零-ID 验证通过等于本轮验证通过。旧 Huesken 平台验证与旧冻结保持原样；未更改旧协议去追认新增数据。

## 公开原始来源

Katoh T, Suzuki T. *Specific residues at every third position of siRNA shape its efficient RNAi activity*. Nucleic Acids Research 35:e27 (2007), https://doi.org/10.1093/nar/gkl1120 。原始 S1 / 方法经 NCBI PMC OA 归档提供，来源链接保存在 sources.json。
