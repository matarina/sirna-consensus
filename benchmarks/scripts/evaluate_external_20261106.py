"""One-shot frozen external evaluation. Never fits or selects models.
Measured-core ranking and full-molecule generation are separate endpoints.
"""
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from sirna_consensus.confirmation import assay_metrics, summarize
from sirna_consensus.rank import _quality_term
from sirna_consensus.scoring import score_candidates
from sirna_consensus.types import Core

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'benchmarks/external_confirmation_20261106'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def main():
    protocol = json.loads((OUT/'test_protocol.json').read_text())
    for p, h in protocol['artifacts'].items():
        assert sha(ROOT/p) == h, p
    records = json.loads((OUT/'records.json').read_text())
    assert [r['record_id'] for r in records] == protocol['test_ids']
    target = ''.join((OUT/'egfp.fasta').read_text().splitlines()[1:])
    cores = [dict(core=Core(r['transcript_id'],r['start'],19,r['core'],chemistry='all_RNA'),
                  tools=[], tool_scores={}) for r in records]
    pred, result, generation = {}, {}, {}
    work = OUT/'outputs'
    for name, relative in protocol['models'].items():
        manifest_path = ROOT/relative
        manifest = json.loads(manifest_path.read_text())
        # Explicit component-level request: only paired cores are observed.
        # This is NOT the full-molecule scorer adapter and has no invented overhang.
        request = dict(schema_version=1, scope='measured_paired_core_component_only',
                       scorer_manifest_sha256=sha(manifest_path),
                       candidates=[dict(key=r['record_id'],paired_guide=r['guide']) for r in records])
        source, dest = work/name/'component_input.json', work/name/'component_output.json'
        save(source,request)
        command = [s.replace('{input}',str(source)).replace('{output}',str(dest)) for s in manifest['command']]
        subprocess.run(command,check=True, cwd=ROOT)
        response = json.loads(dest.read_text())
        assert response['input_sha256'] == sha(source)
        assert set(response['scores']) == set(protocol['test_ids'])
        scores = [response['scores'][r['record_id']] for r in records]
        pred[name] = scores
        result[name] = summarize(assay_metrics(records,scores))
        # Run actual production adapter on observed cores: preserve unsupported status.
        available, provenance = score_candidates(cores,target,manifest_path,work/name/'observed_adapter')
        result[name]['full_molecule_adapter'] = dict(provenance=provenance,
                macro=summarize(assay_metrics(records,[available.get(c['core'].key) for c in cores])),
                note='Missing documented overhangs: this abstention is distinct from paired-core component ranking.')
        cmd = [sys.executable,'-m','sirna_consensus.cli','run','--fasta',str(OUT/'egfp.fasta'),
               '--out',str(work/name/'generation'),'--tools','','--scorer-manifest',str(manifest_path),
               '--strategy','model','--profile','rna19p2','--top-k','5','--min-spacing','30']
        subprocess.run(cmd,check=True,cwd=ROOT)
        run_path = work/name/'generation/run.json'
        run = json.loads(run_path.read_text())
        selected = [r for r in run['candidates'] if r['key'] in run['top']]
        # No complete experimental duplex sequence is available in S1. Therefore
        # none of the generated molecules has a matchable experimental identity.
        assert not any(r['full_molecule_known'] for r in records)
        generation[name] = dict(run_path=str(run_path.relative_to(ROOT)), run_sha256=sha(run_path),
            candidates=run['counts']['enumerated'], selected=len(selected), panel_fill_rate=len(selected)/5,
            candidates_with_experimental_19core=sum(any(r['core']==g['sense'] for r in records) for g in run['candidates']),
            exact_full_molecule_measurement_coverage=0.0,
            known_active=0, known_inactive=0, unknown=len(selected),
            active_budget_fraction_bounds=[0.0,len(selected)/5],
            specificity_verified_panel_fill_rate=0.0,
            specificity='not assessed; reporter is not a human reference target',
            panel=[{k:r[k] for k in ['key','start','sense','guide_dna','ordered_guide','passenger','molecule_id','score']} for r in selected])
    quality = [_quality_term(r['core']) for r in cores]
    pred['quality_rules'] = quality
    result['quality_rules'] = summarize(assay_metrics(records,quality))
    values = np.array([r['value'] for r in records])
    relevance = (rankdata(values)-1)/(len(values)-1)
    ideal = sum((2**v-1)/math.log2(i+2) for i,v in enumerate(sorted(relevance,reverse=True)[:5]))
    a = int(sum(r['value']>=r['active_threshold'] for r in records)); n=len(records)
    result['uniform_random_expectation'] = dict(ndcg5=float(np.mean(2**relevance-1)*sum(1/math.log2(i+2) for i in range(5))/ideal),
          precision5=a/n, top5_active_fraction=a/n, top3_at_least_one_active=1-math.comb(n-a,3)/math.comb(n,3),
          fill_rate=1.0,prediction_coverage=1.0, measurement_coverage=1.0, abstentions=0, budget_yield=a/n,
          spearman=0.0, active_records=a, records=n)
    # Actual full-molecule pipeline on experimental cores uses the documented
    # native/quality fallback rather than assigning a model score to unsupported input.
    result['observed_core_pipeline_fallback'] = result['quality_rules']
    ga,gb=generation['huesken_only_challenger'],generation['current_default']
    ka={r['molecule_id'] for r in ga['panel']}; kb={r['molecule_id'] for r in gb['panel']}
    budget_bound=[-len(kb-ka)/5,len(ka-kb)/5]
    save(work/'predictions.json',dict(record_ids=protocol['test_ids'],scores=pred))
    save(work/'results.json',dict(protocol_sha256=sha(OUT/'test_protocol.json'),
         records_sha256=sha(OUT/'records.json'), metrics=result, generation=generation,
         paired_ndcg_delta=result['huesken_only_challenger']['ndcg5']-result['current_default']['ndcg5'],
         paired_gene_bootstrap_95CI=None, CI_reason='Only one independent gene/study; row bootstrap prohibited.',
         actual_generated_paired_budget_yield_bounds=budget_bound,
         promotion=False, decision='Keep frozen default: one study/gene, no valid paired-gene CI or positive measured generation budget-yield evidence.',
         reused_for_training=False, independent_confirmation_scope='new-study measured-core retrospective ranking; not untouched blind, multi-gene or prospective validation'))
    print(json.dumps({k:{x:v.get(x) for x in ['ndcg5','precision5','spearman']} for k,v in result.items()},indent=2))


if __name__ == '__main__':
    main()
