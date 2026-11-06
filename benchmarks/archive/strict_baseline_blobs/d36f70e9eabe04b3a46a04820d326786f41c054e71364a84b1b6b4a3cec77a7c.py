"""Independent arithmetic checks for grouped ranker outputs and weighted ridge fits."""
import argparse
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.fasta import sha256_file
from sirna_consensus.learning import features
from sirna_consensus.types import reverse_complement


def main():
    p=argparse.ArgumentParser();p.add_argument('--model-dir',required=True);a=p.parse_args()
    out=Path(a.model_dir)
    data=json.loads((out/'development_records.json').read_text())
    split=json.loads((out/'split_manifest.json').read_text())
    pred=json.loads((out/'oof_predictions.json').read_text())
    results=json.loads((out/'results.json').read_text())
    ids={r['record_id']:i for i,r in enumerate(data)}
    oof={r['record_id']:r for r in pred}
    x=features([r['guide'] for r in data]);y=np.array([r['target_rank'] for r in data])
    errors=[];weight_errors=[];covered=[];overlap=[]
    for fold in split['folds']:
        tr=np.array([ids[i] for i in fold['train_ids']]);te=np.array([ids[i] for i in fold['test_ids']]);covered+=fold['test_ids']
        train=[data[i] for i in tr]
        studies={r['study'] for r in train}
        genes={s:{r['gene_id'] for r in train if r['study']==s} for s in studies}
        strata=defaultdict(list)
        for j,r in enumerate(train):strata[r['unit']].append(j)
        conds=defaultdict(set)
        for r in train:conds[(r['study'],r['gene_id'])].add(r['unit'])
        w=np.array([1/(len(studies)*len(genes[r['study']])*len(conds[(r['study'],r['gene_id'])])*len(strata[r['unit']])) for r in train]);w/=w.mean()
        weight_errors.append(float(np.max(np.abs(w-np.array(fold['train_weights'])))))
        # Independent mean-one weighted StandardScaler and ridge normal equations.
        mean=np.average(x[tr],axis=0,weights=w)
        scale=np.sqrt(np.average((x[tr]-mean)**2,axis=0,weights=w));scale[scale<1e-14]=1
        z=(x[tr]-mean)/scale;zt=(x[te]-mean)/scale
        ym=np.average(y[tr],weights=w)
        alpha=fold['family_params']['ridge']['alpha']
        coef=np.linalg.solve(z.T@(w[:,None]*z)+alpha*np.eye(z.shape[1]),z.T@(w*(y[tr]-ym)))
        errors.append(float(np.max(np.abs(zt@coef+ym-np.array([oof[data[i]['record_id']]['ridge'] for i in te])))))
        tg={data[i]['gene_id'] for i in tr};vg={data[i]['gene_id'] for i in te}
        overlap.append(len(tg&vg))
    groups=split['records'];bad=[];seen={}
    for r in data:
        for s in (r['core'],reverse_complement(r['core'])):
            for i in range(3):
                mer=s[i:i+17]
                if mer in seen and groups[seen[mer]]!=groups[r['record_id']]:bad.append((seen[mer],r['record_id']))
                seen[mer]=r['record_id']
    # Independent Hamming<=2 audit using two-position deletion signatures,
    # rather than the training code's vectorised Hamming distance implementation.
    signatures={};near_bad=[]
    for r in data:
        for s in (r['core'],reverse_complement(r['core'])):
            for left,right in itertools.combinations(range(19),2):
                key=(left,right,s[:left]+s[left+1:right]+s[right+1:])
                group=groups[r['record_id']]
                if key in signatures and signatures[key]!=group:
                    near_bad.append(r['record_id'])
                signatures[key]=group
    # Independently recompute metrics, including deterministic score ties.
    units=defaultdict(list)
    for r in data:units[r['unit']].append(r)
    ndcgs=[];rhos=[];p5=[]
    for unit,rows in units.items():
        ordered=sorted(rows,key=lambda r:(-oof[r['record_id']]['selected'],r['guide'],r['record_id']))
        k=min(5,len(rows));discount=1/np.log2(np.arange(k)+2)
        observed=np.array([r['target_rank'] for r in ordered[:k]])
        ideal=np.array(sorted([r['target_rank'] for r in rows],reverse=True)[:k])
        denom=np.sum(np.expm1(np.log(2)*ideal)*discount)
        if denom:ndcgs.append(float(np.sum(np.expm1(np.log(2)*observed)*discount)/denom))
        v=[r['value'] for r in rows];s=[oof[r['record_id']]['selected'] for r in rows]
        if np.ptp(v)>0 and np.ptp(s)>0:rhos.append(float(spearmanr(v,s).statistic))
        p5.append(np.mean([r['value']>=r['active_threshold'] for r in ordered[:k]]))
    expected=results['metrics']['selected']
    checks={'exact_oof_ids':len(covered)==len(set(covered))==len(data) and set(covered)==set(ids)==set(oof),
        'gene_disjoint':max(overlap)==0,'shared17mer_grouped':not bad,'Hamming2_grouped':not near_bad,
        'independent_weight_max_error':max(weight_errors),'independent_ridge_max_error':max(errors),
        'weighted_ridge_refits_pass':max(errors)<1e-10 and max(weight_errors)<1e-12,
        'ndcg5_error':abs(np.mean(ndcgs)-expected['ndcg5']),
        'true_spearman_error':abs(np.mean(rhos)-expected['spearman']),
        'precision5_error':abs(np.mean(p5)-expected['precision5']),
        'results_sha256':sha256_file(out/'results.json')}
    checks['all_passed']=all(checks[k] for k in ('exact_oof_ids','gene_disjoint','shared17mer_grouped','Hamming2_grouped','weighted_ridge_refits_pass')) and max(checks[k] for k in ('ndcg5_error','true_spearman_error','precision5_error'))<1e-12
    checks['all_passed']=bool(checks['all_passed'])
    (out/'independent_verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    print(json.dumps(checks,indent=2));assert checks['all_passed']


if __name__=='__main__':main()
