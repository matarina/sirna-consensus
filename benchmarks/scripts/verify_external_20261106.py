"""Separate checkpoint/metric arithmetic for the frozen external round.
Does not import either production feature encoder or production metric function.
"""
import hashlib
import importlib.util
import json
import math
import pickle
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'benchmarks/external_confirmation_20261106'
spec=importlib.util.spec_from_file_location('independent_math',ROOT/'benchmarks/scripts/verify_guarded.py')
independent=importlib.util.module_from_spec(spec)
spec.loader.exec_module(independent)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=json.loads((OUT/'test_protocol.json').read_text())
    records=json.loads((OUT/'records.json').read_text())
    predictions=json.loads((OUT/'outputs/predictions.json').read_text())
    results=json.loads((OUT/'outputs/results.json').read_text())
    checks=[]
    def check(name,ok,detail=''):
        checks.append(dict(name=name,pass_=bool(ok),detail=detail))
    check('all preregistered artifact hashes',all(sha(ROOT/f)==h for f,h in p['artifacts'].items()))
    check('exact frozen record IDs',p['test_ids']==predictions['record_ids']==[r['record_id'] for r in records])
    check('protocol and records linked',results['protocol_sha256']==sha(OUT/'test_protocol.json') and results['records_sha256']==sha(OUT/'records.json'))
    for name,path in p['models'].items():
        manifest_path=ROOT/path
        m=json.loads(manifest_path.read_text())
        bundle=pickle.loads((manifest_path.parent/m['checkpoint']).read_bytes())
        version='sequence_v1' if name=='current_default' else bundle['params']['feature']
        x=independent.encode([r['guide'] for r in records],version)
        if name=='current_default':
            model=bundle
        else:
            if bundle['scaler'] is not None:
                x=bundle['scaler'].transform(x)
            model=bundle['model']
        with threadpool_limits(limits=4):
            v=model.predict(x)
        check(name+' independent checkpoint features/predictions',np.allclose(v,predictions['scores'][name],rtol=0,atol=1e-12))
    for name,scores in predictions['scores'].items():
        actual=independent.metrics(records,scores)
        expected=results['metrics'][name]
        check(name+' independent metrics',all((v is None and expected[k] is None) or
              (v is not None and expected[k] is not None and abs(v-expected[k])<1e-12) for k,v in actual.items()))
        selected=sorted(zip(records,scores),key=lambda t:(-t[1],t[0]['guide'],t[0].get('molecule_id',''),t[0]['record_id']))[:5]
        check(name+' exact Top5 IDs',[r['record_id'] for r,s in selected]==expected['per_gene_assay'][0]['top_ids'])
    # Analytic random arithmetic independently from midrank expectation.
    vals=[r['value'] for r in records]; n=len(vals)
    rel=(rankdata(vals)-1)/(n-1)
    gains=np.expm1(np.log(2)*rel)
    ideal=np.sort(gains)[::-1][:5]@np.reciprocal(np.log2(np.arange(2,7)))
    random_ndcg=float(np.mean(gains)*np.reciprocal(np.log2(np.arange(2,7))).sum()/ideal)
    check('independent uniform random expectation',abs(random_ndcg-results['metrics']['uniform_random_expectation']['ndcg5'])<1e-12)
    for name,g in results['generation'].items():
        run=json.loads((ROOT/g['run_path']).read_text())
        selected=[r for r in run['candidates'] if r['key'] in run['top']]
        target=''.join((OUT/'egfp.fasta').read_text().splitlines()[1:])
        check(name+' complete generation universe',len(run['candidates'])==len(target)-22==g['candidates'])
        check(name+' target complementation and RNA19p2',all(target[r['start']-1:r['start']+18]==r['sense'] and
              r['ordered_guide'][:19].translate(str.maketrans('ACGT','TGCA'))[::-1]==r['sense'] and
              len(r['ordered_guide'])==len(r['passenger'])==21 for r in run['candidates']))
        check(name+' default diverse panel',len(selected)==5 and all(abs(a['start']-b['start'])>=30 for i,a in enumerate(selected) for b in selected[i+1:]))
        check(name+' missing activity retained unknown',g['known_active']==g['known_inactive']==0 and g['unknown']==5 and g['exact_full_molecule_measurement_coverage']==0)
        check(name+' true full-molecule adapter abstains',results['metrics'][name]['full_molecule_adapter']['provenance']['supported_count']==0)
    check('no invalid row bootstrap/no promotion',results['paired_gene_bootstrap_95CI'] is None and results['promotion'] is False)
    result=dict(passed=all(c['pass_'] for c in checks),checks=checks,
                results_sha256=sha(OUT/'outputs/results.json'),verifier_sha256=sha(Path(__file__)))
    (OUT/'outputs/independent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    if not result['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
