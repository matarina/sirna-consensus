"""Build the Chinese Wang/Filhol report from verified, frozen evaluation artifacts."""
from __future__ import annotations

import hashlib
import html
import json
import os
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[2]
BENCH = PACKAGE / 'benchmarks'
DATA = PACKAGE.parent / 'independent_validation/denovo_two_benchmarks_20261105'
DEST = BENCH / 'three_benchmark_report_zh.html'
F19 = 'binding19-guide21-overhang2-passenger21-overhang2'
esc = html.escape
percent = lambda v: f'{100*v:.1f}%'
points = lambda v: f'{100*v:+.1f} pp'
interval = lambda v: f'[{100*v[0]:+.1f}, {100*v[1]:+.1f}] pp'


def main():
    results_path = DATA / 'outputs/results.json'
    results = json.loads(results_path.read_text())
    studies = results['studies']
    sensitivity = json.loads((DATA / 'outputs/conservative_sensitivity.json').read_text())
    verification = json.loads((DATA / 'outputs/verification.json').read_text())
    digest = hashlib.sha256(results_path.read_bytes()).hexdigest()
    assert verification['all_pass'] and verification['results_sha256'] == digest
    order = ('Wang2009', 'Filhol2012')
    labels = {'Wang2009': 'Wang 2009', 'Filhol2012': 'Filhol 2012'}
    chunks = []

    def paragraph(text, kind=''):
        attr = f' class="{kind}"' if kind else ''
        chunks.append(f'<p{attr}>{text}</p>')

    def section(text, ident):
        chunks.append(f'<h2 id="{ident}">{esc(text)}</h2>')

    def table(caption, headers, rows, ident=''):
        attr = f' id="{ident}"' if ident else ''
        chunks.append(f'<div class="table-wrap"><table{attr}><caption>{esc(caption)}</caption><thead><tr>' +
                      ''.join(f'<th scope="col">{esc(h)}</th>' for h in headers) + '</tr></thead><tbody>' +
                      ''.join('<tr>' + ''.join(f'<td>{cell}</td>' for cell in row) + '</tr>' for row in rows) +
                      '</tbody></table></div>')

    def artifact(path, text):
        relative = os.path.relpath(DATA / path, DEST.parent)
        return f'<a href="{esc(relative, quote=True)}">{esc(text)}</a>'

    chunks.append('<header><p class="eyebrow">siRNA de novo design · 已冻结计算评估</p>'
                  '<h1>Wang 与 Filhol benchmark 的计算可辩护性报告</h1>'
                  '<p class="sub">Wang 2009 · Filhol 2012｜完整生成测试、已测候选排序与证据边界</p>'
                  '<p class="meta">更新：2026-11-05 · Wang / Filhol 版本 · 保留原文件名 three_benchmark_report_zh.html</p></header>')
    chunks.append('<nav aria-label="目录"><a href="#authority">实验依据</a><a href="#leakage">泄漏审计</a>'
                  '<a href="#deployment">实际生成结果</a><a href="#ranking">排序结果</a>'
                  '<a href="#sensitivity">敏感性分析</a><a href="#adjustments">根因与调整</a><a href="#artifacts">产物</a></nav>')
    chunks.append('<div class="abstract"><h2>当前可以支持的结论</h2>')
    w = studies['Wang2009']['ranking_diagnostic']['model']
    f = studies['Filhol2012']['ranking_diagnostic']['model']
    paragraph(f'<strong>推荐 Wang 作为主要回顾性排序 benchmark，Filhol 作为补充挑战集。</strong>'
              f'在已测候选池内，模型 Top5 命中率分别为 <strong>{percent(w["top5_active"])}</strong> 和 '
              f'<strong>{percent(f["top5_active"])}</strong>。完整生成流程已运行，但实际生成 Top5 的实测覆盖不足，'
              '目前无法证明新设计的效力或完整 pipeline 的优越性。')
    paragraph('两项研究均有同行评审的原始实验、多个内源性靶基因和可核验补充数据。'
              '<strong>这里的“权威”指实验来源可信、方法与序列可追溯；没有依据把它们称为官方标准或公认金标准。</strong>'
              '研究质量、与本 pipeline 的兼容性、是否发生数据泄漏，是需要分别核查的问题。')
    chunks.append('</div>')

    section('1. 为什么选择这两份原始实验数据', 'authority')
    table('来源的可信依据与用于本 pipeline 的理由', ['数据集', '原始实验与质量依据', '与本 pipeline 的关联', '必须保留的限制'], [
        ['<a href="https://doi.org/10.1093/nar/gkp864">Wang et al., 2009<br><em>Nucleic Acids Research</em> 37:e152</a>',
         '独立于该论文训练数据的 140 条验证记录，覆盖 14 个 HeLa 内源性基因；作者报告双盲、随机实验顺序。'
         '每基因包括新算法、旧算法及随机位置设计；同批合成，两个实验日、三次转染重复与重复 PCR，每条共 12 次 knockdown 测量。',
         'Supplementary Table S3 提供 guide、passenger 及 3 / 30 nM 剩余 mRNA；可核查双链身份，并在固定剂量下评估 TopK。',
         '作者也开发了被比较的算法，不能视为对本 pipeline 的第三方盲测。面板混合预选与随机设计；passenger 的 dTdT 化学身份需单独保留。'],
        ['<a href="https://doi.org/10.1371/journal.pone.0048057">Filhol et al., 2012<br><em>PLOS ONE</em> 7:e48057</a>',
         '10 个内源性癌症相关基因；HeLa 20 nM 条件下用 qRT-PCR 测定；两个内参、独立生物实验与重复测量，含阳性和阴性对照。',
         'S1 发布两条链，S3 发布逐条实验数据；21 nt、19 bp 配对区加 2 nt 悬垂，与当前 designer / scorer 的常规 siRNA 格式接近。',
         '候选按 DSIR 高分预选，不能代表转录本上所有位置。正文称 128 条，S1 / S3 实际均有 129 个唯一实验 ID；须按原表审计。']
    ], 'authority-table')
    paragraph('推荐依据是原始实验设计、可追溯序列与当前任务的适配程度。期刊名称只提供发表背景；'
              '两份论文均未验证本 pipeline，也未对其任何新设计提供效力背书。'
              '这两项实验支持常规细胞转染 siRNA 的回顾性评估，不能直接外推到体内递送或任意修饰化学体系。')
    paragraph('原论文：<a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC2794195/">Wang 全文与实验方法</a>；'
              '<a href="https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0048057">Filhol 全文与补充材料</a>。'
              'Wang 使用 S3 验证面板，不使用该论文开发模型时分析的 Huesken 数据。')

    section('2. 样本整理、端点与泄漏审计', 'leakage')
    table('原始行、排除原因与最终候选池', ['研究', '原始记录', '既有 Sciabola 近重复排除', '无法映射', '重复分子合并', '最终唯一 duplex / 基因'], [
        [labels[name], str(studies[name]['cohort']['raw_records']), str(studies[name]['cohort']['exposure_excluded']),
         str(studies[name]['cohort']['mapping'].get('unmapped', 0)),
         str(studies[name]['cohort']['eligible_records'] - studies[name]['cohort']['eligible_unique_molecules']),
         f'{studies[name]["cohort"]["eligible_unique_molecules"]} / {studies[name]["cohort"]["genes"]}'] for name in order
    ], 'cohort-table')
    paragraph('<strong>Wang 主端点：</strong>3 nM，knockdown = 100 − percent mRNA remaining；30 nM 仅作预设敏感性分析，'
              '不跨剂量平均。相同 duplex 的重复观察在各剂量内分别取算术均值。'
              '<strong>Filhol 主端点：</strong>S3 最后一列 “%silencing (from statistical model)”，直接按 knockdown 读取；'
              '不把 DSIR 预测分数当作实测效力，也不按性能选择列。另一 dilution-series 列的文字说明存在歧义，未用于主结论。')
    paragraph('Filhol 使用与原始表中长度一致、全部 129 个配对核心均可唯一映射的历史 RefSeq；'
              '当前 RefSeq 的变更曾造成三个 CSNK2B 位点缺失。Wang 按预先固定的 RefSeq 映射，保留六条未映射记录及原因。'
              '全部原始记录、排除项、转录本和文件 SHA256 可从产物核查。')
    paragraph('对原始 269 条记录的 guide 配对 19-core，在正向与反向互补两种方向上，'
              '与公开 Hu2361 / Mix472 / Taka702 的 Hamming ≤ 2 匹配数均为 <strong>0</strong>。'
              '与此前已使用 Sciabola 核心匹配的五条记录已排除；两份新面板之间无 Hamming ≤ 2 匹配。'
              '更保守的连续 17-mer 位移重叠排除及整基因 HIF1A 排除另行报告。')
    paragraph('<strong>这项审计限定于已固定的参考序列。</strong>完整上游 checkpoint 训练与模型选择谱系尚未独立证明，'
              '因此不能写“保证完全无泄漏”或“从未见过这些基因”。本次未拟合或调参，使用原始冻结 checkpoint，'
              '未使用 Sciabola 适配模型或外部结果选出的 selector policy。已查看这些结果后，再据此调整的流程需要新的封存测试或前瞻实验。', 'notice')

    section('3. 本次实际执行的两类评估', 'protocol')
    table('两类评估回答不同问题', ['评估', '进入排序的候选', 'Top5 如何产生', '可解释范围'], [
        ['完整生成测试', '完整转录本 → siRNAdesigner / siRNAforge 原始生成候选',
         '冻结 native / model selector，soft 标注；每种 duplex 格式各取 Top5，之后再对接实测数据。',
         '生成召回、实际 Top5 实测覆盖及已知结果；未测分子的效力未知。'],
        ['已测候选排序诊断', '仅输入符合排除规则的已发表 measured duplexes',
         '在同一已测池内分别按原始冻结模型与 native quality fallback 取 Top5。',
         '已测候选间的排序能力；不能替代真实 de novo 选择验证。']
    ])
    paragraph('有效定义：实测 knockdown ≥ 70%。每个基因的 precision@5 = 有效数 / 5；各基因等权平均。'
              '随机期望直接取该基因同一候选池的有效比例，不通过 Monte Carlo 近似。'
              '95% CI 对基因作配对 bootstrap，10,000 次，seed 20261105；Spearman 使用真实秩相关。'
              '所有主分析基因均至少有五个候选。24 个 study–gene 任务对应 23 个不同基因，HIF1A 在两项研究均出现。')

    section('4. 完整生成结果：效力目前不可识别', 'deployment')
    rows = []
    for name in order:
        s = studies[name]; p = s['generated_pool']; n = s['cohort']['eligible_unique_molecules']
        rows.append([labels[name], f'{p["generated"]:,}', f'{p["exact_pool_matches"]}/{n}',
                     f'{p["guide_core_pool_matches"]}/{n}', percent(p['frozen_model_supported_generated']/p['generated']),
                     str(p['multi_tool_exact_molecules']), str(p['default_min_votes2_top_n'])])
    table('完整生成候选池；严格匹配要求两条链碱基和悬垂均一致',
          ['研究', '生成分子', '已测 duplex 召回', '19-core 召回（诊断）', '模型可打分比例', '双工具精确分子', '默认 ≥2 票输出数'], rows, 'generation-table')
    rows = []
    for name in order:
        for strategy in ('native', 'model'):
            for group, v in studies[name]['deployment'][strategy].items():
                fmt = '19+2' if group.startswith(F19) else '21 blunt'
                rows.append([labels[name], strategy + ('（native fallback）' if group.endswith('/native_fallback') else ''),
                             fmt, f'{v["exact_measured"]}/{v["selected"]}（{percent(v["coverage"])}）',
                             str(v['guide_core_matches']), str(v['genes_with_complete_measured_top5']), '未知：覆盖不足'])
    table('实际保存的 Top5 panel；选择之后才查询实验记录',
          ['研究', 'selector', '分子格式', '完整 duplex 实测覆盖', '19-core 匹配（诊断）', '完整测得 Top5 的基因数', '新设计效力'], rows, 'coverage-table')
    paragraph('19+2 与 21 blunt 各自是一个五条预算的 panel。19-core 相同但 passenger、悬垂或配对长度变化的分子，'
              '不能继承原实验标签。21 blunt 不受冻结 19-core 模型支持，model arm 对该格式走 native fallback。'
              'Wang 的小写 tt 表示 dTdT；当前软件只保存 DNA 字母碱基身份，不编码核糖 / 脱氧核糖差异，'
              '因此即便完整碱基匹配，也仍需单独核对化学身份。')
    bounds = studies['Filhol2012']['deployment']['model'][F19+'/model']['active_fraction_identification_bounds']
    paragraph(f'Filhol 的 model 19+2 panel 中，已测 5/50，有效 3 条；其全部 50 条的有效比例只能界定在 '
              f'{percent(bounds[0])}–{percent(bounds[1])}。Wang 对应 panel 为 0/70，效力界限为 0–100%。'
              '<strong>低覆盖表示缺少实验信息，不能写成有效率 0%，也不能直接用 scorer-only 命中率替代。</strong>')
    table('原始 research run Top10（min_votes=1，30 nt spacing）',
          ['研究', '保存的设计数', '完整 duplex 实测匹配', '19+2 / 21 blunt 数量'], [
              [labels[name], str(studies[name]['original_pipeline_top10']['original_pipeline_min_votes1_top_n']),
               str(studies[name]['original_pipeline_top10']['original_pipeline_top10_exact_measured']),
               f'{studies[name]["original_pipeline_top10"]["format_counts"].get(F19,0)} / '
               f'{studies[name]["original_pipeline_top10"]["format_counts"].get("binding21-guide21-overhang0-passenger21-overhang0",0)}'] for name in order
          ])

    section('5. 已测候选排序结果', 'ranking')
    rows = []
    for name in order:
        for strategy in ('native', 'model'):
            v = studies[name]['ranking_diagnostic'][strategy]
            rows.append([labels[name], strategy, percent(v['top5_active']), percent(v['random_expected']),
                         points(v['gain']), interval(v['gain_ci95']), f'{v["enrichment_ratio"]:.3f}×', f'{v["mean_gene_spearman"]:.3f}'])
    table('各基因等权；native 在该已测池中是无生成器证据的规则 quality fallback',
          ['研究', '排序', 'Top5 命中率', '随机期望', '富集差值', '差值 95% CI', '富集比', '平均基因 Spearman'], rows, 'ranking-table')
    table('模型相对规则比较器的改进', ['研究', 'model − native', '95% CI', '判读'], [
        [labels[name], points(studies[name]['model_minus_native']['gain']),
         interval(studies[name]['model_minus_native']['ci95']), '区间跨 0，尚未证明优于规则比较器'] for name in order
    ])
    paragraph('Wang 显示较强的已测候选排序能力；Filhol 的模型随机富集为正，但其平均基因 Spearman 接近零，'
              '且删除 HIF1A 后置信区间跨零，说明结论对候选与基因范围敏感。'
              '两者对随机的结果均不能推出模型优于规则比较器，更不能推出完整 de novo pipeline 的效力优势。'
              '这些 CI 只反映固定实验面板下的基因间抽样，不包含实验测量误差、预选偏倚或多个比较的校正。')

    section('6. 保守排除与剂量敏感性', 'sensitivity')
    rows = []
    modes = {'exclude_prior_shared_gene_HIF1A': '整基因移除 HIF1A',
             'exclude_shifted_17mer_matches': '移除连续 17-mer 重叠', 'exclude_both': '同时执行两项排除'}
    for name in order:
        for mode, label in modes.items():
            v = sensitivity['studies'][name]['evaluations'][mode]['summary']['model']
            rows.append([labels[name], label, f'{v["molecules"]} / {v["genes"]}', percent(v['top5']),
                         percent(v['random']), points(v['gain']), interval(v['gain_ci95'])])
    table('附加排除规则在计算研究级性能前记录；复用冻结分数，不拟合',
          ['研究', '附加排除', '分子 / 基因', 'Model Top5', '随机期望', '富集差值', '95% CI'], rows, 'sensitivity-table')
    paragraph('连续 17-mer 审计覆盖同寄存 Hamming 距离可能漏掉的位移重叠；在主队列外另外识别 Filhol 两条、Wang 一条 HIF1A 记录。'
              '整基因 HIF1A 排除同时避免此前 Sciabola 使用过该基因及两项研究共享该基因的影响。'
              'Wang 的富集在这些保守分析中保留；Filhol 整基因排除后的区间跨 0。不得依据哪个分析更好而切换主队列。')
    wd = studies['Wang2009']['ranking_diagnostic']
    table('Wang 30 nM 次端点：保持原来选中的设计不变', ['native Top5', 'model Top5', '同池随机期望'], [
        [percent(wd['native']['sensitivity_top5']), percent(wd['model']['sensitivity_top5']), percent(wd['model']['sensitivity_random'])]
    ])

    section('7. 每个基因的结果', 'genes')
    rows = []
    for name in order:
        for r in studies[name]['per_gene']:
            rows.append([labels[name], r['gene'], str(r['n']), percent(r['random_top5']),
                         percent(r['measured_panel']['native']['top5_active']), percent(r['measured_panel']['model']['top5_active']),
                         str(r['exact_pool_matches']), str(r['guide_core_pool_matches']), str(r['hard_filter_rejected_active'])])
    table('已测候选排序与完整生成召回并列展示；不混用两类分母',
          ['研究', '基因', '已测候选', '随机有效率', 'Native Top5', 'Model Top5', '完整生成池 duplex 匹配', 'Core 匹配', '被硬规则排除的有效 duplex'], rows)

    section('8. 根因与优先调整方法', 'adjustments')
    table('当前证据指向的阻塞点', ['问题', '本次证据', '建议调整'], [
        ['分子格式与投票不兼容', 'designer 输出 19+2，forge 输出 blunt21；两项研究均无 ≥2 工具支持的精确分子，默认两票输出为空。',
         '先定义共同的分子格式、悬垂来源与化学身份；或把 site agreement 单独作为证据特征，保留原分子身份。不得把异构分子算作同一 duplex 票。'],
        ['生成 / 模型支持覆盖不足', '冻结模型仅能打分约 11–12% 的生成分子；Wang 已测完整 duplex 召回仅 2/127。',
         '增加明确的 19+2 候选生成模式，考虑覆盖所有合规 target windows，再把生成器评分作为特征。改变 blunt21 的配对区或悬垂会产生新分子。'],
        ['硬过滤损失有效分子',
         '在已测池内，Filhol 与 Wang 分别有 '+str(studies['Filhol2012']['generated_pool']['hard_filter_rejected_active'])+' 与 '+
         str(studies['Wang2009']['generated_pool']['hard_filter_rejected_active'])+' 条有效 duplex 不满足当前硬规则；其中部分原本也未被生成。',
         '在独立开发集比较 hard 与 soft，并约束 TopK 完成率；统一原始 run 的硬过滤与 select 的 soft 标注口径。'],
        ['排序增益和跨实验适用性仍有限', 'Wang 已测排序有用，Filhol 较脆弱；两项 model − native 的 CI 都跨 0。',
         '保留冻结模型与简单规则基线；如需训练，用独立开发数据并按基因 / study / 序列近邻分组嵌套验证。若修正上游 LSTM 布局，须重新训练验证。'],
        ['缺少新设计的实验标签', '没有任何基因的实际生成 Top5 被完整测量。',
         '冻结调整后流程，在多个新基因上合成实际 Top5；固定剂量、化学与转染条件，对比同预算基线和随机设计。历史未测候选不能补成阴性。']
    ])
    paragraph('本次没有完成可验证的转录组脱靶筛查，因此效力结果不提供脱靶特异性结论。'
              '上述调整是后续开发方向；当前报告没有在这两个面板上调参，也不把看过标签后的流程修改称为新的独立测试。')

    section('9. 可复核产物与历史版本', 'artifacts')
    table('本页数值由固定 JSON 生成；未重新拟合模型', ['产物', '用途'], [
        [artifact('outputs/results.json', 'results.json'), '完整生成、覆盖、Top5、随机期望、CI 与逐基因结果'],
        [artifact('outputs/exposure_mapping_audit.json', 'exposure_mapping_audit.json'), '参考序列 SHA256、排除、映射与转录本选择'],
        [artifact('outputs/conservative_sensitivity.json', 'conservative_sensitivity.json'), '整基因与位移近重复排除结果'],
        [artifact('outputs/verification.json', 'verification.json'), '24 个任务、144 个五条 panel、原始端点与六条 fresh checkpoint 重现均通过'],
        [artifact('protocol.md', 'protocol.md')+' / '+artifact('preprocessing_notes.md', 'preprocessing_notes.md'), '预先记录的规则与原表对账说明'],
        [artifact('outputs/source_hashes.json', 'source_hashes.json'), '原始补充文件及转录本的校验值'],
        [artifact('report.html', '完整英文实验报告'), '实验范围与方法详情'],
        ['<a href="archive/three_benchmark_report_zh.before_wang_filhol_20261105.html">更新前中文页面（仅历史记录）</a>',
         '保留此前三个 benchmark 版本；历史指标与判定仅供追溯，不作为本次独立验证证据']
    ])
    paragraph('本页保留旧文件名，但当前仅评估 Wang / Filhol 两个面板。此前 siRNAEfficacyDB 子集有训练序列重叠，'
              '5-study 存在上游评估暴露，Sciabola 已用于适配，均不作为本页新的封存外部测试。')
    chunks.append('<pre>bash independent_validation/denovo_two_benchmarks_20261105/reproduce.sh\n'
                  'python sirna-consensus/benchmarks/scripts/build_wang_filhol_report_zh.py</pre>')
    paragraph('复现命令从 <code>/data1/maxw/proj/cargo_denovo</code> 执行。'
              '本页构建脚本只读取已验证结果；不会运行训练或重新生成 benchmark 分数。', 'note')
    chunks.append('<footer>结果 SHA256：<code>'+digest+'</code><br>'
                  '结论范围：Wang / Filhol 回顾性已测候选排序；完整 de novo 生成的实际选择仍缺少充分实测覆盖。'
                  '同行评审来源、序列去重与模型无泄漏是三项独立证据要求。</footer>')

    css = '''
    :root{--ink:#14242d;--muted:#52636e;--navy:#153f53;--teal:#08756b;--line:#d8e1e4;--bg:#edf2f4}
    *{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.85 "Noto Serif CJK SC","Source Han Serif SC","Songti SC",serif}
    main{max-width:1200px;margin:0 auto;background:#fff;padding:40px 48px 56px;box-shadow:0 2px 20px #10203010}
    h1{font-size:28px;line-height:1.4;margin:0 0 10px}h2{font-size:19px;color:var(--navy);margin:34px 0 12px;border-bottom:1px solid var(--line);padding-bottom:7px}
    p{margin:10px 0}.eyebrow{color:var(--teal);font:12px/1.5 system-ui,sans-serif;letter-spacing:.08em}.sub{color:var(--muted);font-size:15px}.meta,.note,footer{color:var(--muted);font-size:12px}
    nav{display:flex;flex-wrap:wrap;gap:8px 18px;padding:14px 0;border-block:1px solid var(--line);font:13px/1.6 system-ui,sans-serif;margin:22px 0}
    a{color:#146281;text-underline-offset:3px;overflow-wrap:anywhere}.abstract{background:#f2f8f6;border:1px solid #cbded8;border-left:4px solid var(--teal);padding:16px 22px;margin:24px 0}.abstract h2{border:0;margin:0 0 8px;font-size:16px}
    .notice{background:#fff7e8;border-left:3px solid #bc862a;padding:12px 16px}.table-wrap{overflow-x:auto;margin:14px 0 20px}table{width:100%;border-collapse:collapse;font:13px/1.65 system-ui,"Noto Sans CJK SC",sans-serif}caption{text-align:left;color:var(--muted);font-size:12px;padding-bottom:8px}
    th{text-align:left;background:#ecf2f5;border-top:1.5px solid var(--navy);border-bottom:1px solid var(--line);font-size:12px}th,td{padding:9px 10px;vertical-align:top}td{border-bottom:1px solid var(--line)}tbody tr:last-child td{border-bottom:1px solid var(--navy)}tbody tr:nth-child(even){background:#fafcfd}
    code,pre{font:12px/1.6 ui-monospace,Menlo,Consolas,monospace;overflow-wrap:anywhere}code{background:#f0f4f6;padding:2px 4px}pre{background:#f3f6f7;padding:15px;overflow:auto}footer{border-top:1px solid var(--line);padding-top:16px;margin-top:28px}
    @media(max-width:720px){main{padding:24px 16px}h1{font-size:23px}h2{font-size:17px}table{min-width:640px}.abstract{padding:12px 14px}}
    @media print{body{background:white}main{max-width:none;padding:0;box-shadow:none}nav{display:none}h2{break-after:avoid}.table-wrap{overflow:visible}tr{break-inside:avoid}a{color:inherit}}
    '''
    document = '<!doctype html>\n<html lang="zh-CN"><head><meta charset="utf-8">\n' \
               '<meta name="viewport" content="width=device-width,initial-scale=1">\n' \
               '<title>Wang 与 Filhol benchmark 的计算可辩护性报告</title>\n<style>'+css+'</style></head>\n<body><main>\n' + \
               '\n'.join(chunks) + '\n</main></body></html>\n'
    DEST.write_text(document)
    print(json.dumps({'output':str(DEST),'results_sha256':digest,'studies':list(order)},ensure_ascii=False))


if __name__ == '__main__':
    main()
