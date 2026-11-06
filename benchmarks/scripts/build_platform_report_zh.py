"""Generate current report solely from versioned development/evaluation artifacts."""
from __future__ import annotations
import hashlib
import html
import json
import math
from collections import defaultdict
import os
import shutil
from pathlib import Path

PKG=Path(__file__).resolve().parents[2]
BENCH=PKG/'benchmarks';MODEL=PKG/'models/development_20261106';OUT=BENCH/'results/platform_v3'
DEST=BENCH/'three_benchmark_report_zh.html'
HIST=BENCH/'archive/three_benchmark_report_zh.before_platform_v3_20261106.html'
esc=html.escape


def read(p):return json.loads(p.read_text())
def link(p,label):return f'<a href="{esc(os.path.relpath(p,DEST.parent),quote=True)}">{esc(label)}</a>'
def pct(x):return f'{100*x:.1f}%'
def num(x):return '未定义' if x is None else f'{x:.3f}'


def main():
    results=read(MODEL/'results.json');registry=read(MODEL/'data_registry.json');evaluation=read(OUT/'evaluation.json')
    verify=read(MODEL/'verification.json');independent=read(MODEL/'independent_verification.json')
    assert verify['all_passed'] and independent['all_passed']
    assert hashlib.sha256((MODEL/'results.json').read_bytes()).hexdigest()==independent['results_sha256']==evaluation['model_results_sha256']
    if not HIST.exists():
        HIST.parent.mkdir(exist_ok=True)
        shutil.copyfile(DEST,HIST)
    data=read(MODEL/'development_records.json')
    analytic_random={}
    for study in sorted({r['study'] for r in data}):
        units=defaultdict(list)
        for r in data:
            if r['study']==study:units[r['unit']].append(r)
        p5=[];ndcg=[]
        for rows in units.values():
            p5.append(sum(r['value']>=r['active_threshold'] for r in rows)/len(rows))
            gains=[2**r['target_rank']-1 for r in rows]
            discounts=[1/math.log2(i+2) for i in range(min(5,len(rows)))]
            ideal=sum(v*d for v,d in zip(sorted(gains,reverse=True),discounts))
            if ideal:ndcg.append(sum(gains)/len(gains)*sum(discounts)/ideal)
        analytic_random[study]={'precision5':sum(p5)/len(p5),'ndcg5':sum(ndcg)/len(ndcg)}
    (OUT/'random_expectation.json').write_text(json.dumps(analytic_random,indent=2)+'\n')
    from build_confirmation_report_zh import section
    chunks=[section(PKG)]
    def paragraph(text):chunks.append('<p>'+text+'</p>')
    def heading(text,ident):chunks.append(f'<h2 id="{ident}">{esc(text)}</h2>')
    def table(headers,rows):chunks.append('<div class="scroll"><table><thead><tr>'+''.join('<th>'+esc(h)+'</th>' for h in headers)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(cell)+'</td>' for cell in row)+'</tr>' for row in rows)+'</tbody></table></div>')
    chunks.append('<h2>历史开发：Wang／Filhol 回归与四研究分析</h2><p>2026-11-06 · schema v3 · 保留历史报告文件名</p>')
    gain=results['selected_vs_quality'];m=results['metrics'];split=read(MODEL/'split_manifest.json')
    chunks.append('<div class="notice"><strong>结论边界：已改善开发期排序与候选覆盖；尚无新设计效力的独立实验证明。</strong>')
    paragraph(f'整合 {registry["records"]} 条唯一记录、{registry["genes"]} 个基因，按基因和近重复序列形成 {split["groups"]} 个验证组。'
        f'嵌套分组 OOF 的 NDCG@5：新模型 {num(m["selected"]["ndcg5"])}，质量规则 {num(m["quality"]["ndcg5"])}；'
        f'配对增益 {gain["delta"]:+.3f}，95% 区间 [{gain["ci95"][0]:+.3f}, {gain["ci95"][1]:+.3f}]。'
        f'预先指定的默认模型门槛{"通过" if results["promoted"] else "未通过"}。')
    lower=[s for s,v in evaluation['study_results'].items() if v['metrics']['selected']['precision5']<v['metrics']['quality']['precision5']]
    uncertain=[s for s,v in evaluation['study_results'].items() if v['selected_vs_quality']['ci95'][0]<=0<=v['selected_vs_quality']['ci95'][1]]
    paragraph('Top5活性比例低于质量规则的研究：'+esc('、'.join(lower) or '无')+'；单研究NDCG增益区间跨零：'+esc('、'.join(uncertain) or '无')+'。'
        '不能把综合NDCG提升解释成已证明每个靶点更有效。四份数据均已接触并用于开发，今后不再称独立最终测试集。')
    chunks.append('</div>')
    heading('1. 当前默认流程与分子身份','workflow')
    chunks.append('''<svg viewBox="0 0 1120 390" role="img" aria-label="候选枚举与工具证据独立进入排序，再分别输出研究和筛查面板" style="width:100%">
    <defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="9" refY="3" orient="auto"><path d="M0,0 L0,6 L9,3 z" fill="#476575"/></marker></defs>
    <g fill="#edf5f3" stroke="#476575" stroke-width="2"><rect x="15" y="145" width="150" height="65" rx="12"/><rect x="205" y="45" width="240" height="70" rx="12"/><rect x="205" y="240" width="240" height="70" rx="12"/><rect x="495" y="145" width="200" height="65" rx="12"/><rect x="755" y="30" width="340" height="70" rx="12"/><rect x="755" y="165" width="340" height="70" rx="12"/><rect x="755" y="300" width="340" height="70" rx="12"/></g>
    <g fill="none" stroke="#476575" stroke-width="2" marker-end="url(#arrow)"><path d="M165 165 L185 165 L185 80 L205 80"/><path d="M165 190 L185 190 L185 275 L205 275"/><path d="M445 80 L470 80 L470 165 L495 165"/><path d="M445 275 L470 275 L470 190 L495 190"/><path d="M695 165 L725 165 L725 65 L755 65"/><path d="M695 190 L725 190 L725 200 L755 200"/><path d="M925 235 L925 300"/></g>
    <g font-size="17" fill="#183a48" text-anchor="middle"><text x="90" y="182">转录本 FASTA</text><text x="325" y="75">完整枚举 RNA 19+2</text><text x="325" y="99">21nt双链 · 明确overhang</text><text x="325" y="270">外部工具 · 独立证据</text><text x="325" y="294">坐标 / 链方向 / 分子身份</text><text x="595" y="173">质量规则 / 分组新模型</text><text x="595" y="198">经验过滤默认作提示</text><text x="925" y="60">研究 Top5 + ≥30nt 间隔 Top5</text><text x="925" y="85">默认导出完整 guide / passenger</text><text x="925" y="195">固定参考 · 双链全配对区 + 真正3′UTR</text><text x="925" y="221">失败或不完整 = unverified</text><text x="925" y="329">筛查通过 Top5 + 间隔 Top5</text><text x="925" y="354">不足五条如实报告，不补未验证序列</text></g></svg>''')
    paragraph('枚举不算工具投票；工具21bp配对区的guide前19nt对应sense后19nt，按坐标与方向附加位点证据。'
        '只有完整两条链和已知化学身份一致，才算同一分子的投票。未知passenger或修饰不推定为全RNA。'
        '边界不足或含N的窗口不生成；研究面板和筛查面板分开，支持legacy-v2回放。')
    heading('2. 为什么 Wang 与 Filhol 仍适合作为回归基准','authority')
    table(['研究','可追溯性与适配理由','使用限制'],[
        ['<a href="https://doi.org/10.1093/nar/gkp864">Wang 2009 · NAR</a>','S3提供双链及3/30nM实验；14个HeLa内源性基因，随机与算法设计并存，实验重复可核查。','作者开发算法时使用Huesken；当前验证表有重复分子、6条不可唯一映射，passenger dTdT与全RNA不同。'],
        ['<a href="https://doi.org/10.1371/journal.pone.0048057">Filhol 2012 · PLOS ONE</a>','S1/S3提供序列与逐条qRT-PCR；10个HeLa内源性基因，19+2格式接近当前设计。','按DSIR高分预选；原表129条，与正文128不一致；主端点固定为统计模型silencing列。']])
    paragraph('“权威”指原始实验、同行评审及可核查补充材料，不代表官方金标准，也不等于对本平台的第三方盲测。'
        '目前两者用于开发回归和挑战分析；新的独立泛化证据须来自未接触数据或前瞻实验。'
        +link(HIST,'升级前 Wang／Filhol 详细报告（含原实验与旧结果）')+'。')
    heading('3. 数据登记、训练与防泄漏','data')
    table(['开发数据','唯一记录','当前状态'],[[esc(s),n,'已接触；非独立最终测试'] for s,n in registry['studies'].items()])
    paragraph(f'原始总记录 {registry["raw_records"]}；隔离 {len(registry["quarantine"])} 条；同分子同条件技术重复均值合并。'
        '不同研究、剂量与实验条件不混合数值尺度，目标为条件内活性中秩。HGNC固定文件解析旧符号，保留原名与解析路径。')
    paragraph('同HGNC基因、19core正反向Hamming≤2、共享连续17-mer构成连通验证组。外层五折、内层三折，'
        '折内拟合加权标准化、权重和超参数；训练依次平衡研究、基因、条件、记录。主指标对基因—实验条件单位宏平均，'
        '不按siRNA行数加权；跨研究同基因仍在同一验证组。区间是固定OOF预测的连通组bootstrap，不包含实验误差、重新分折或训练不确定性，也不能纠正开发者反复观察数据的选择偏差。最终全量拟合只用于部署，不作为回顾性性能证据。')
    paragraph('从零训练Ridge和浅树；特征仅来自配对19nt序列、组成及末端不对称，不输入旧模型输出。'
        'OligoFormer冻结预测只作描述性对照：它见过Huesken且上游谱系不完整。Huesken原上游工作区缺失，当前整理表已固定哈希，不能声称原始来源已完整恢复。')
    heading('4. 各研究 OOF 表现：先看差异，再看总体','ranking')
    table(['数据','n / 基因','新模型 NDCG@5','质量规则','增益 [95%区间]','新模型 P@5','质量规则 P@5','真实 Spearman'],[
        [esc(s),f'{v["records"]} / {v["genes"]}',num(v['metrics']['selected']['ndcg5']),num(v['metrics']['quality']['ndcg5']),
         f'{v["selected_vs_quality"]["delta"]:+.3f} [{v["selected_vs_quality"]["ci95"][0]:+.3f}, {v["selected_vs_quality"]["ci95"][1]:+.3f}]',
         pct(v['metrics']['selected']['precision5']),pct(v['metrics']['quality']['precision5']),num(v['metrics']['selected']['spearman'])]
        for s,v in evaluation['study_results'].items()])
    paragraph('P@5活性阈值：Huesken normalized_activity≥0.7，其余主端点敲低≥70%。NDCG不是成功率。'
        '不同实验尺度的相关性和活性比例分别计算再宏平均，不池化原始实验值。单次随机排列只是可重放对照；不能把一次随机Top5的波动当成随机期望。')
    table(['同池随机期望（精确计算）','NDCG@5','P@5'],[[esc(study),num(v['ndcg5']),pct(v['precision5'])] for study,v in analytic_random.items()])
    table(['总体方法','NDCG@5','P@5','Spearman'],[[esc(name),num(v['ndcg5']),pct(v['precision5']),num(v['spearman'])] for name,v in m.items()])
    table([f'共同候选池对照（n={evaluation["common_comparator_records"]}）','NDCG@5','P@5'],[
        [esc(name),num(v['ndcg5']),pct(v['precision5'])] for name,v in evaluation['common_comparator'].items()])
    paragraph('共同候选池表仅使用所有方法均有预测的同一批记录，重新计算该池的真实相关性与NDCG增益。'
        '外部工具原生分数来自固定旧运行，按修正后的位点规则重新附加；这些旧工具分数也可能存在上游数据接触。')
    heading('5. 留一研究验证：同时移除相关基因／序列组','loso')
    def avg(rows,key):
        vals=[r[key] for r in rows if r[key] is not None];return sum(vals)/len(vals) if vals else None
    table(['留出研究','训练 / 测试','额外清除','NDCG@5','质量规则','P@5'],[
        [esc(r['study']),f'{r["train_records"]} / {r["test_records"]}',r['purged_records'],
         num(avg(r.get('metrics',[]),'ndcg5')),num(avg(r.get('quality_metrics',[]),'ndcg5')),pct(avg(r.get('metrics',[]),'precision5')) if r.get('metrics') else '不可计算'] for r in results['loso']])
    heading('6. de novo 覆盖与工程回归','deployment')
    rows=[]
    for study in ('Wang2009','Filhol2012'):
        ds=[r for r in evaluation['deployment'] if r['study']==study]
        for arm in ('quality_spaced','native_spaced','model_spaced','model_unconstrained','model_hard_spaced'):
            vs=[r['arms'][arm] for r in ds]
            rows.append([study,arm,sum(v['selected'] for v in vs),sum(v['paired_core_matches'] for v in vs),sum(v['full_strand_base_matches'] for v in vs),sum(v['known_chemistry_matches'] for v in vs),'未知'])
    table(['研究','同一枚举池 / 面板','选出数','19core实测匹配','完整双链碱基匹配','已知化学匹配','新设计活性'],rows)
    paragraph('以上de novo排序使用对应基因的外层留出checkpoint；先选序列，再查实测覆盖。仅19core相同不能移植实验活性。'
        'passenger缺失、overhang或dTdT不同都不能称化学身份相同。新生成且未实验测量的序列活性保持未知。')
    pipeline=read(OUT/'pipeline_runs.json')
    paragraph(f'另完成 Wang／Filhol 共 {len(pipeline["runs"])} 个 schema-v3 完整输出回归：复用固定的外部工具原始结果，'
        '重新枚举、匹配证据、评分并写出双链与JSON。该工程回归使用全量开发模型，不是独立效力测试；'
        '下表的排序诊断未运行全人类参考筛查，不能将其面板称为筛查通过。')
    heading('7. 根因、已修复问题与下一轮调整','causes')
    rows=[]
    for study in ('Wang2009','Filhol2012'):
        aa=[r for r in evaluation['ablations'] if r['study']==study]
        rows.append([study,sum(r['measured'] for r in aa),sum(r['native_projected_core_recall'] for r in aa),sum(r['enumeration_19core_recall'] for r in aa),sum(r['hard_filter_rejected'] for r in aa),sum(r['hard_filter_rejected_active'] for r in aa),sum(r['exact_two_tool_enumerated'] for r in aa)])
    table(['研究','实测分子','工具位点召回（含长guide投影）','完整枚举召回','硬过滤删除','其中活性分子','确证双工具分子票'],rows)
    paragraph('<strong>已证实的工程问题：</strong>工具预筛限制候选覆盖；不同分子格式无法靠投票合并；旧prefix连接可错配位点；'
        '硬规则会删除实验活性分子。已实施完整枚举、坐标/几何证据、化学身份和默认软提示，分别显示过滤与间隔损失。')
    paragraph('<strong>尚未证明的生物学解释：</strong>Filhol的DSIR预选、跨研究剂量/检测体系差异、化学格式未知及序列特征不足，'
        '都可能限制泛化；本次结果不能区分其因果贡献。排序指标与70%阈值命中不是同一目标，需要分别评估。')
    paragraph('<strong>下一轮候选方法：</strong>若产品目标改为Top5活性比例，先预注册该指标和模型升级门槛，再在嵌套训练内比较分类/成对排序；'
        '实验条件与化学体系足够时加入条件分层或多任务头；另做折内比较的靶RNA可及性与热力学特征。'
        '不得根据当前Filhol或Sciabola结果反复挑模型后再称其独立测试。优先收集新基因的前瞻Top5实验。')
    heading('8. 参考筛查、验证与可复现产物','artifacts')
    screening_path=OUT/'screening_integration.json'
    if screening_path.exists():
        screening=read(screening_path)
        paragraph('人类参考集成状态：'+esc(screening['status'])+'。'+link(screening_path,'参考与筛查实测记录')+'。')
        if screening.get('reference'):
            rr=screening['reference'];paragraph(f'{esc(rr["release"])}：参考转录本 {rr["transcripts"]}、基因身份 {rr["genes"]}；3′UTR无法投影 {rr["utr3_unprojectable"]}。'
                '两条链均筛查；exact/1 mismatch非目标命中拒绝，2 mismatch和seed负担为提示，缺失参考或映射不完整不算通过。')
        table(['完整参考集成样例','枚举候选','通过','拒绝','未验证','筛查且满足间隔的Top5'],[
            [esc(r['gene']),r['counts']['enumerated'],r['candidate_screen_status'].get('pass',0),
             r['candidate_screen_status'].get('reject',0),r['candidate_screen_status'].get('unverified',0),
             r['counts']['screened_diverse_top']] for r in screening.get('runs',[])])
        for r in screening.get('runs',[]):
            full_run=read(Path(r['run']))
            absent=sum(c['screen']['status']=='unverified' and not c['screen'].get('intended_gene_exact_hit',False)
                       for c in full_run['candidates'])
            if absent:
                paragraph(f'{esc(r["gene"])} 有 {absent} 个候选在固定参考中未找到目标基因完全匹配，保持unverified；'
                    '不能据此推定实验无效，也不进入筛查通过面板。')
    else:
        paragraph('完整GENCODE参考正在构建；尚未产生全参考筛查结果，当前面板不应称筛查通过。')
    paragraph(f'独立算术核查通过：折外ID覆盖、基因/17mer/正反向Hamming≤2隔离、等权重重算、Ridge独立重拟合、真实Spearman、NDCG和P@5。'
        f'Ridge最大误差 {independent["independent_ridge_max_error"]:.2e}；checkpoint重放最大误差 {verify["checkpoint_prediction_max_error"]:.1e}。')
    for p,label in [(MODEL/'results.json','完整训练结果'),(MODEL/'data_registry.json','数据登记与来源哈希'),(MODEL/'split_manifest.json','分组、权重与折内模型'),
                    (MODEL/'independent_verification.json','独立算术核查'),(OUT/'evaluation.json','逐基因评估与面板消融'),(OUT/'pipeline_runs.json','24个完整输出回归'),
                    (PKG/'docs/prospective_validation.md','前瞻实验预注册模板'),(HIST,'旧版Wang／Filhol结果')]:paragraph(link(p,label))
    paragraph('结果SHA256：<code>'+independent['results_sha256']+'</code>。旧freeze未改写；旧代码以原哈希字节归档，新实现单独冻结。')
    # Per-gene evidence is generated, not an editorial explanation detached from numbers.
    heading('9. 逐基因／实验条件结果','pergene')
    for study,v in evaluation['study_results'].items():
        chunks.append('<h3>'+esc(study)+'</h3>')
        q={r['unit']:r for r in v['metrics']['quality']['per_gene_assay']}
        table(['基因 / 条件','n','新模型NDCG','规则NDCG','新模型P@5','规则P@5','Spearman'],[
            [esc(r['gene'])+'<br><small>'+esc(r['unit'].split('|')[-1])+'</small>',r['n'],num(r['ndcg5']),num(q[r['unit']]['ndcg5']),pct(r['precision5']),pct(q[r['unit']]['precision5']),num(r['spearman'])]
            for r in v['metrics']['selected']['per_gene_assay']])
    DEST.write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">'
        '<title>siRNA平台升级与开发验证</title><style>body{font:16px/1.7 system-ui,sans-serif;max-width:1200px;margin:40px auto;padding:0 24px;color:#173541;background:#fcfdfd}'
        'h1{font-size:30px}h2{margin-top:42px;border-top:1px solid #cad9de;padding-top:20px}table{border-collapse:collapse;width:100%;font-size:14px}td,th{border:1px solid #d1dfe2;padding:9px;text-align:left}'
        'th{background:#e7f0ef}.scroll{overflow:auto}.notice{padding:20px;background:#eaf3f3;border-left:5px solid #278073}a{color:#086d80}small{color:#536b76}code{overflow-wrap:anywhere}</style>'
        +''.join(chunks)+'</html>')
    print(DEST)


if __name__=='__main__':main()
