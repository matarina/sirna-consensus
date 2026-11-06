"""Same-universe measured ranking, de novo OOF coverage and engineering regressions."""
from __future__ import annotations
import argparse
import csv
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits
from scipy.stats import rankdata

PKG=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PKG/'src'))
from sirna_consensus.fasta import read_transcript, sha256_file
from sirna_consensus.learning import features, unit_metrics, macro, confidence
from sirna_consensus.selection import from_record, enumerate_duplexes, load_candidates, attach_evidence, rank_candidates, panels
from sirna_consensus.filters import hard_flags
from sirna_consensus.pipeline import Config, run
from sirna_consensus.types import Core


def read(path):return json.loads(Path(path).read_text())
def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--model-dir',required=True);ap.add_argument('--out',required=True);a=ap.parse_args()
    model_dir,out=Path(a.model_dir).resolve(),Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    records=read(model_dir/'development_records.json');oof=read(model_dir/'oof_predictions.json');splits=read(model_dir/'split_manifest.json')
    by_id={r['record_id']:r for r in oof};root=PKG.parent
    data=root/'independent_validation/denovo_two_benchmarks_20261105'
    source_pins={}
    # Reuse frozen single-candidate predictions, pinning the selection artifacts.
    paths=list((PKG/'benchmarks/results/development/full').glob('*/model/selection.json'))
    paths+=list((root/'independent_validation/sciabola/canonical/selections').glob('*/model/selection.json'))
    for study in ('Wang2009','Filhol2012'):
        paths+=list((data/study/'selections/measured').glob('*/model/selection.json'))
    frozen={}
    for path in sorted(paths):
        d=read(path);source_pins[str(path)]=sha256_file(path)
        for r in d['candidates']:
            if r.get('model_score') is not None:
                frozen[(r['transcript_id'],r['guide_dna'])]=r['model_score']
    candidate_sources={
        'huesken2005':PKG/'benchmarks/data/huesken/candidates.tsv',
        'sciabola2013':root/'independent_validation/sciabola/canonical/candidates.tsv',
        'Wang2009':data/'Wang2009/candidates.tsv','Filhol2012':data/'Filhol2012/candidates.tsv'}
    positions={}
    for study,path in candidate_sources.items():
        source_pins[str(path)]=sha256_file(path)
        with path.open() as h:
            for r in csv.DictReader(h,delimiter='\t'):
                c=from_record(r);positions[(study,c.transcript_id,c.guide)]=c.start
    targets=defaultdict(list)
    for i,r in enumerate(records):targets[(r['study'],r['gene_original'],r['transcript_id'])].append(i)
    native_predictions=np.full(len(records),np.nan)
    deployment=[];ablations=[]
    # Each gene uses the selected model from its held-out outer fold, never the pooled fit.
    fold_cache={}
    with threadpool_limits(limits=4):
        for (study,gene,tid),indices in sorted(targets.items()):
            if study=='huesken2005':
                fasta=PKG/f'benchmarks/data/huesken/transcripts/{gene}.fa'
                run_path=PKG/f'benchmarks/results/runs/huesken/{gene}/run.json'
            elif study=='sciabola2013':
                fasta=root/f'independent_validation/sciabola/canonical/transcripts/{tid}.fa'
                run_path=root/f'independent_validation/sciabola/runs/{tid}/run.json'
            else:
                fasta=data/study/'transcripts'/f'{gene}.fa'
                run_path=data/study/'runs_final'/gene/'run.json'
            _,sequence=read_transcript(fasta)
            source_pins[str(fasta)]=sha256_file(fasta)
            native=[]
            if run_path.exists():
                native,oldrun=load_candidates(run_path);source_pins[str(run_path)]=sha256_file(run_path)
            measured=[]
            for i in indices:
                r=records[i];start=r.get('start') or positions.get((study,tid,r['guide']))
                if not start:
                    # Records excluded from the historical evaluation can still map exactly.
                    found=[p+1 for p in range(len(sequence)-18) if sequence[p:p+19]==r['core']]
                    if len(found)!=1:raise ValueError('ambiguous development target position: '+r['record_id'])
                    start=found[0]
                c=Core(tid,start,19,r['core'],r['guide'],r['passenger'],r['guide'][19:],r['passenger'][19:] if r['passenger'] else '',r['chemistry'])
                assert sequence[c.start-1:c.end]==c.sense
                measured.append({'core':c,'tools':[],'tool_scores':{},'record_index':i})
            attach_evidence(measured,native)
            ranked=rank_candidates(measured,'native','soft')
            lookup={r['key']:r['score'] for r in ranked}
            for m in measured:native_predictions[m['record_index']]=lookup[m['core'].key]
            # Engineering ablation: explicit hard filtering and tool-universe losses.
            all_enum=enumerate_duplexes(tid,sequence);attach_evidence(all_enum,native)
            generated_cores={r['core'].sense for r in all_enum}
            native_cores={r['core'].sense for r in native if r['core'].length==19}
            native_projected={r['core'].sense[-19:] for r in native if r['core'].length>=19}
            rejected=[m for m in measured if hard_flags(m['core'])]
            fold=by_id[records[indices[0]]['record_id']]['fold']
            assert all(by_id[records[i]['record_id']]['fold']==fold for i in indices)
            ck=splits['folds'][fold]['checkpoints'][splits['folds'][fold]['selected_params']['family']]
            if fold not in fold_cache:fold_cache[fold]=pickle.loads((model_dir/ck['path']).read_bytes())
            if all_enum:
                values=fold_cache[fold].predict(features([r['core'].guide for r in all_enum]))
                model_scores={r['core'].key:float(v) for r,v in zip(all_enum,values)}
            else:model_scores={}
            def coverage(selected):
                pool=[bykey[k] for group in selected.values() for k in group]
                base_matches=[];core_matches=[];chem_matches=[]
                for r in pool:
                    cm=[records[i] for i in indices if records[i]['core']==r['sense']]
                    bm=[m for m in cm if m['guide']==r['guide_dna'] and m['passenger'] and m['passenger']==r['passenger']]
                    ch=[m for m in bm if m['chemistry']=='all_RNA']
                    core_matches+=bool(cm),;base_matches+=bool(bm),;chem_matches+=bool(ch),
                return {'selected':len(pool),'paired_core_matches':sum(core_matches),'full_strand_base_matches':sum(base_matches),
                    'known_chemistry_matches':sum(chem_matches),'full_base_coverage':sum(base_matches)/len(pool) if pool else 0,
                    'active_fraction':None,'unmeasured_activity':'unknown','selected_keys':[r['key'] for r in pool]}
            arms={}
            for strategy in ('quality','native','model'):
                rr=rank_candidates(all_enum,strategy,'soft',model_scores)
                bykey={r['key']:r for r in rr}
                for spaced in (False,True):
                    arms[strategy+('_spaced' if spaced else '_unconstrained')]=coverage(panels(rr,5,30 if spaced else 0))
                if strategy=='model':
                    hard=rank_candidates(all_enum,strategy,'hard',model_scores);bykey={r['key']:r for r in hard}
                    arms['model_hard_spaced']=coverage(panels(hard,5,30))
            deployment.append({'study':study,'gene':records[indices[0]]['gene'],'transcript_id':tid,'fold':fold,
                'candidates':len(all_enum),'site_supported':sum(bool(r['site_evidence']) for r in all_enum),
                'arms':arms,'scorer':'outer-fold held-out checkpoint','reference_screen':'not_run_in_this_ranking_diagnostic'})
            ablations.append({'study':study,'gene':records[indices[0]]['gene'],'measured':len(indices),
                'native_19core_recall':sum(records[i]['core'] in native_cores for i in indices),
                'native_projected_core_recall':sum(records[i]['core'] in native_projected for i in indices),
                'enumeration_19core_recall':sum(records[i]['core'] in generated_cores for i in indices),
                'hard_filter_rejected':len(rejected),'hard_filter_rejected_active':sum(records[m['record_index']]['value']>=records[m['record_index']]['active_threshold'] for m in rejected),
                'exact_two_tool_enumerated':sum(len(r['tools'])>=2 for r in all_enum)})
            print(json.dumps({'study':study,'gene':gene,'enumerated':len(all_enum)}),flush=True)
    predictions={name:np.asarray([r[name] for r in oof]) for name in ('selected','ridge','tree','quality','random')}
    predictions['native']=native_predictions
    predictions['frozen_oligoformer']=np.asarray([frozen.get((r['transcript_id'],r['guide']),np.nan) for r in records])
    common=np.isfinite(predictions['frozen_oligoformer'])
    def summarize(data,pred):
        m=unit_metrics(data,pred)
        return {'ndcg5':macro(m),'spearman':macro(m,'spearman'),'precision5':macro(m,'precision5'),'per_gene_assay':m}
    common_records=[dict(r) for r,c in zip(records,common) if c]
    common_units=defaultdict(list)
    for i,r in enumerate(common_records):common_units[r['unit']].append(i)
    for ids in common_units.values():
        ranks=rankdata([common_records[i]['value'] for i in ids])
        for i,rank in zip(ids,ranks):common_records[i]['target_rank']=float((rank-1)/(len(ids)-1)) if len(ids)>1 else .5
    study_results={}
    for study in sorted({r['study'] for r in records}):
        ix=np.asarray([i for i,r in enumerate(records) if r['study']==study])
        study_records=[records[int(i)] for i in ix]
        groups=np.asarray([splits['records'][r['record_id']] for r in study_records])
        study_results[study]={'records':len(ix),'genes':len({r['gene_id'] for r in study_records}),
            'metrics':{name:summarize(study_records,p[ix]) for name,p in predictions.items() if name!='frozen_oligoformer'},
            'selected_vs_quality':confidence(study_records,groups,predictions['selected'][ix],predictions['quality'][ix]),
            'frozen_available':int(common[ix].sum())}
    record={'schema_version':1,'status':'development_only','independent_validation_passed':False,
        'model_results_sha256':sha256_file(model_dir/'results.json'),'sources':source_pins,
        'study_results':study_results,'common_comparator_records':len(common_records),
        'common_comparator':{name:summarize(common_records,p[common]) for name,p in predictions.items()},
        'frozen_comparator_limitation':'Frozen OligoFormer saw Huesken; its upstream training lineage is incomplete. Cached raw predictions are descriptive comparators only.',
        'deployment':deployment,'ablations':ablations}
    write(out/'evaluation.json',record)
    # Actual schema-v3 pipeline replay for both requested benchmarks, including exports.
    runs=[]
    for study in ('Filhol2012','Wang2009'):
        for fasta in sorted((data/study/'transcripts').glob('*.fa')):
            directory=out/'runs'/study/fasta.stem
            result=run(Config(str(fasta),str(directory),tools=(),
                evidence_run=str(data/study/'runs_final'/fasta.stem/'run.json'),
                scorer_manifest=str(model_dir/'scorer_manifest.json')))
            runs.append({'study':study,'gene':fasta.stem,'counts':result['counts'],
                         'run_sha256':sha256_file(directory/'run.json'),'run':str(directory/'run.json')})
    write(out/'pipeline_runs.json',{'status':'engineering_regression; pooled development model, not independent efficacy test','runs':runs})
    print(json.dumps({'output':str(out),'runs':len(runs),'common_comparator_records':len(common_records)}))


if __name__=='__main__':main()
