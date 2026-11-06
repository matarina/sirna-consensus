"""Build the living R&D report from frozen outcomes; no fitting or evaluation."""
import hashlib
import html
import json
import os
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
B=ROOT/'benchmarks'
J=B/'development_journal'
SOURCES={}


def read(path):
    p=ROOT/path
    SOURCES[str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text())


def table(head, rows):
    return '<div class="scroll"><table><thead><tr>'+''.join('<th>'+html.escape(x)+'</th>' for x in head)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+str(x)+'</td>' for x in row)+'</tr>' for row in rows)+'</tbody></table></div>'


def bars(rows, percent=False):
    text=''
    for i,(name,value,color) in enumerate(rows):
        val=f'{value:.1%}' if percent else f'{value:.4f}'
        y=25+i*43
        text+=f'<text x="0" y="{y}" font-size="13">{html.escape(name)}</text><rect x="205" y="{y-17}" width="{value*430}" height="24" rx="4" fill="{color}"><title>{html.escape(name)}: {val}</title></rect><text x="{212+value*430}" y="{y}" font-size="13">{val}</text>'
    return f'<svg viewBox="0 0 720 {len(rows)*43+12}" role="img" aria-label="同一数据集内的方法表现对比">{text}</svg>'


def main():
    platform=read('benchmarks/platform_v3_summary.json')
    hu=read('models/huesken_only_20261106/results.json')
    mixed=read('models/development_20261106/results.json')
    cohort=read('benchmarks/confirmation_20261106/cohorts.json')
    strict=read('benchmarks/confirmation_20261106/summary.json')
    external=read('benchmarks/external_confirmation_20261106/outputs/results.json')
    verify=read('benchmarks/external_confirmation_20261106/outputs/independent_verification.json')
    sv=read('benchmarks/external_confirmation_20261106/outputs/source_verification.json')
    sci_path=ROOT.parent/'independent_validation/sciabola/adaptation_20261105/outputs/results.json'
    sci=read('benchmarks/development_journal/sciabola_adaptation_snapshot.json')
    if sci_path.exists() and hashlib.sha256(sci_path.read_bytes()).hexdigest()!=sci['original_sha256']:
        raise ValueError('Historical Sciabola snapshot differs from original')
    preserve=read('benchmarks/development_journal/historical_report_preservation.json')
    curve=read('benchmarks/development_journal/curve_metrics.json')
    names={'current_default':'当前默认模型','huesken_only_challenger':'Huesken-only 候选','quality_rules':'既有简单规则','uniform_random_expectation':'随机选择（期望）'}
    ext_rows=[]
    for k,label in names.items():
        m=external['metrics'][k]
        ext_rows.append([label,f'{m["ndcg5"]:.4f}',f'{m["precision5"]:.1%}',f'{m["top3_at_least_one_active"]:.1%}',f'{m["spearman"]:.4f}',f'{m["fill_rate"]:.0%}',f'{m["prediction_coverage"]:.0%}'])
    roles={
     'huesken2005':('开发','训练与嵌套分组开发；不是独立测试'),
     'sciabola2013':('已暴露','先做适配，随后纳入混合训练；不能再次作为新模型的独立确认'),
     'legacy_holdout':('已暴露','五研究历史 holdout；已反复评估，且有上游模型选择及 N→A 参考历史'),
     'Harborth':('不可评估','单独映射文件 0 条，不等于历史 Harborth 没有实验数据'),
     'Martinelli':('已暴露','历史评估；修饰双链不属于当前 all-RNA 主评估域'),
     'Shabalina':('已暴露','历史评估，且缺少完整双链身份'),
     'siRNAEfficacyDB':('训练重叠','此前子集 470/470 核心出现在公开 Hu 训练参考，不能声称独立'),
     'UBE2H_snapshot':('已暴露','历史单靶标展示，含训练基因 UBE2H'),
     'OligoFormer':('不是队列','模型清单不是实验测试集；上游训练谱系需要单独审计')}
    cohorts=[]
    for r in cohort['cohorts']:
        role,reason=roles[r['name']]
        cohorts.append(dict(name=r['name'],n=r.get('records',0),role=role,reason=reason))
    cohorts += [dict(name='Wang2009',n=129,role='开发',reason='历史外部比较后纳入混合模型训练；现在归入开发数据'),dict(name='Filhol2012',n=129,role='开发',reason='同上；化学身份不完整，需要保持实验域说明'),dict(name='Katoh2007 / EGFP',n=702,role='新外部评估',reason='训练集外的新研究；曾作为序列查重参考，非从未接触的盲测；仅 1 个 reporter'),dict(name='Jagla2005',n=None,role='未获取',reason='未恢复可认证的 601 行序列/活性表，不虚构成绩'),dict(name='siRNADiscovery / ACON',n=None,role='未获取',reason='102 条新实验需向作者请求；公开旧合集不能替代新实验数据'),dict(name='Davis2025',n=None,role='未运行',reason='化学修饰域，本轮没有执行 stress test')]
    cohort_table=table(['队列 / 研究','本地评估记录','当前证据身份','原因'],[[html.escape(r['name']),r['n'] if r['n'] is not None else '—','<span class="tag">'+r['role']+'</span>',r['reason']] for r in cohorts])
    timeline=[
      ('2026-10-28','d5065d0','v1 输出与共识语义修复','完善输出正确性、来源追踪和 benchmark 框架。','工程'),
      ('2026-10-28','31d518e','完整分子身份与实验评估','保留 guide / passenger / overhang，加入可复现的实验标签评估；序列相似不等于同一个实验分子。','工程'),
      ('2026-11-02','751ad88','历史 holdout 与冻结','引入权重独立 holdout、复现冻结与 pipeline 加固。后来暴露审计进一步收紧独立资格，不能沿用早期“独立”标签。','历史'),
      ('2026-11-04','95d42b9','选择策略与报告迭代','增加 benchmark selector policy 和 guide21 评分路径；对应早期三 benchmark 报告。','历史'),
      ('2026-11-05','d0de4bb','Sciabola 适配与独立纠错','实验记录日期为 11-05，工件于 11-06 基线一并纳入。修复 Ridge 基因权重和 Pearson 误标 Spearman；重跑后所有适配方案均未达晋升门槛。','纠错'),
      ('2026-11-05 → 11-06','d0de4bb','Wang / Filhol 与 v3 架构整合','从工具交集改为目标派生完整双链穷举；GENCODE v50 筛查，四研究混合开发模型。既有数据全部按已暴露开发证据解释。','开发'),
      ('2026-11-06 14:00','d0de4bb','冻结可运行基线','固定实现、模型、数据与旧轮次结果。此后所有替换需通过独立预注册门槛。','冻结'),
      ('2026-11-06 14:02–14:31','f8adbf9','Huesken-only 严格研发轮次','先注册，再仅用 Huesken 训练。HGNC 大小写解析失败被保留并修复；候选 CV 为开发诊断。旧队列独立资格为 0，保留默认。','开发'),
      ('2026-11-06 15:04–15:27','1131d5a','新外部数据采集与测试冻结','采集协议 3dabca9 → 原始 Katoh S1 → 固定全部 702 IDs / 代码 / 哈希 → 0aa3855 补入被 Git 忽略的核心表，均在首次评分前。','冻结'),
      ('2026-11-06 15:37','dc12833','完成 Katoh 外部测试','两个模型冻结评分、实际生成、23 项独立数值检查与 5 项源数据核查；单基因证据不足，默认不替换。','外部')]
    timeline_html=''.join(f'<article class="event"><div class="event-date">{date}<code>{commit}</code></div><div><span class="tag">{kind}</span><h3>{title}</h3><p>{body}</p></div></article>' for date,commit,title,body,kind in timeline)
    arms={'frozen_baseline':'冻结 OligoFormer 基线','ridge':'Ridge head','ranking':'Pairwise head','head_finetune':'受限 head 微调','nested_selected':'嵌套选择'}
    sci_rows=[(label,sci['summary_top5'][k]['mean_precision'],'#8395a7' if k=='frozen_baseline' else '#9678c5') for k,label in arms.items()]
    study_rows=[]
    for name,r in platform['study_results'].items():
        a,b=r['metrics']['selected'],r['metrics']['quality']
        study_rows.append([name,r['records'],r['genes'],f'{a["ndcg5"]:.4f}',f'{b["ndcg5"]:.4f}',f'{a["precision5"]:.1%}','开发 OOF'])
    genes=[dict(gene=r.get('gene',r['gene_id']),study=r['unit'].split('|')[0],ndcg=r['ndcg5'],p5=r['precision5'],rho=r['spearman'],n=r['n']) for r in mixed['metrics']['selected']['per_gene_assay']]
    screen_rows=''
    for r in platform['human_screen_runs']:
        s=r['candidate_screen_status'];n=r['counts']['enumerated']
        segs=''.join(f'<span style="width:{100*s.get(k,0)/n}%" class="{k}" title="{k}: {s.get(k,0)}"></span>' for k in ['pass','reject','unverified'])
        screen_rows+=f'<div class="screenrow"><strong>{r["gene"]}</strong><div class="stack">{segs}</div><small>总候选 {n} · pass {s.get("pass",0)} / reject {s.get("reject",0)} / unverified {s.get("unverified",0)} · screened Top5 5/5</small></div>'
    commits=subprocess.check_output(['git','log','--date=iso-strict','--pretty=format:%h\t%ad\t%s'],cwd=ROOT,text=True).splitlines()
    data=dict(external={k:external['metrics'][k] for k in names},names=names,cohorts=cohorts,genes=genes)
    links=[('external_confirmation_20261106/report_zh.html','Katoh 新外部评估（详细）'),(preserve['archived_report'],'改版前严格轮次 + v3 完整报告'),('archive/three_benchmark_report_zh.before_platform_v3_20261106.html','Wang / Filhol 历史版本'),('archive/three_benchmark_report_zh.before_wang_filhol_20261105.html','早期三 benchmark 版本'),('computational_defensibility_report_zh.html','早期计算可辩护性报告（历史口径）'),('confirmation_20261106/protocol.json','Huesken-only 预注册协议'),('external_confirmation_20261106/test_protocol.json','Katoh 预测前冻结协议'),('../docs/prospective_validation.md','前瞻实验验证模板'),('development_journal/sources.json','本页数据来源及哈希'),('development_journal/historical_report_preservation.json','历史报告冻结映射')]
    evidence_links=''.join(f'<a class="evidence-link" href="{html.escape(url,quote=True)}"><span>↗</span>{label}</a>' for url,label in links)
    replacements={'DATA':json.dumps(data,ensure_ascii=False).replace('</','<\\/'),'TIMELINE':timeline_html,'COHORT_TABLE':cohort_table,
      'EXTERNAL_TABLE':table(['方法','NDCG@5','P@5 / Top5 活性比例','Top3≥1','Spearman','填满率','核心评分覆盖'],ext_rows),
      'EXTERNAL_CHART':bars([(names[k],external['metrics'][k]['ndcg5'],c) for k,c in zip(names,['#197d83','#9bd6bb','#92a0b4','#c8d1dc'])]),
      'SCI_CHART':bars(sci_rows,True),
      'STUDY_TABLE':table(['研究','记录','基因','混合模型 NDCG','规则 NDCG','混合模型 P@5','身份'],study_rows),
      'SCREEN':screen_rows,'LINKS':evidence_links,
      'GENE_TABLE':table(['基因','研究','n','NDCG@5','P@5','Spearman'],[[r['gene'],r['study'],r['n'],f'{r["ndcg"]:.3f}',f'{r["p5"]:.0%}',f'{r["rho"]:.3f}' if r['rho'] is not None else '—'] for r in genes]),
      'COMMITS':table(['Commit','提交时间（+08:00）','提交说明'],[[html.escape(x) for x in line.split('\t',2)] for line in commits]),
      'CURVE_TABLE':table(['方法','ROC AUC','Average precision (AP)','原始阳性 / 阴性'],[[names[k],f'{v["roc_auc"]:.4f}',f'{v["average_precision"]:.4f}',f'{curve["positives"]} / {curve["negatives"]}'] for k,v in curve['metrics'].items()]),
      'VERIFIED':str(sum(c['pass_'] for c in verify['checks'])+sum(sv['checks'].values()))}
    template=B/'templates/development_journal_zh.html'
    text=template.read_text()
    for key,value in replacements.items(): text=text.replace('@@'+key+'@@',value)
    if '@@' in text: raise ValueError('unresolved template placeholder')
    (B/'three_benchmark_report_zh.html').write_text(text)
    J.mkdir(exist_ok=True)
    manifest=dict(description='Living report provenance; rebuilding does not rerun training or benchmarks.',metrics_sources=SOURCES,
      template_sha256=hashlib.sha256(template.read_bytes()).hexdigest(),generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      report_sha256=hashlib.sha256(text.encode()).hexdigest(),commit_history=commits)
    (J/'sources.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')


if __name__=='__main__': main()
