"""Freeze v3 separately from historical experiments; never rewrite their locks."""
import json
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

PKG=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PKG/'src'))
from sirna_consensus.fasta import sha256_file


def read(p):return json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')


def main():
    bench=PKG/'benchmarks';out=bench/'results/platform_v3';model=PKG/'models/development_20261106'
    results=read(model/'results.json');evaluation=read(out/'evaluation.json');registry=read(model/'data_registry.json')
    verification=read(model/'verification.json');independent=read(model/'independent_verification.json')
    screen=read(out/'screening_integration.json')
    assert verification['all_passed'] and independent['all_passed'] and screen['status']=='passed'
    assert independent['results_sha256']==sha256_file(model/'results.json')==evaluation['model_results_sha256']
    tests={}
    for name in ('tests_unit_ml.xml','tests_external.xml'):
        tree=ET.parse(out/name)
        suites=tree.findall('.//testsuite')
        tests[name]={k:sum(int(s.attrib.get(k,0)) for s in suites) for k in ('tests','failures','errors','skipped')}
        assert tests[name]['failures']==tests[name]['errors']==0
    summary={'schema_version':1,'status':'development_only','independent_efficacy_validation_passed':False,
        'records':registry['records'],'genes':registry['genes'],'studies':registry['studies'],
        'model_results_sha256':sha256_file(model/'results.json'),'evaluation_sha256':sha256_file(out/'evaluation.json'),
        'metrics':{k:{n:v for n,v in m.items() if n!='per_gene_assay'} for k,m in results['metrics'].items()},
        'model_gate':results['selected_vs_quality'],'promoted':results['promoted'],
        'study_results':{s:{'records':v['records'],'genes':v['genes'],'selected_vs_quality':v['selected_vs_quality'],
            'metrics':{k:{n:val for n,val in m.items() if n!='per_gene_assay'} for k,m in v['metrics'].items()}}
            for s,v in evaluation['study_results'].items()},
        'independent_arithmetic_checks':independent,'checkpoint_oof_checks':verification,
        'reference':screen['reference'],'human_screen_runs':[{k:v for k,v in r.items() if k!='screen_metadata'} for r in screen['runs']],
        'tests':tests,'limitations':registry['provenance_limitations'],
        'primary_metric':'macro gene/assay NDCG@5; CI conditional on fixed OOF predictions',
        'prospective_protocol':'../docs/prospective_validation.md'}
    write(bench/'platform_v3_summary.json',summary)
    files=set(PKG.glob('src/**/*.py'))|set(PKG.glob('tests/**/*.py'))
    files.update(PKG.glob('docs/*.md'))
    files.update(PKG.glob('.github/workflows/*.yml'))
    files.update([PKG/'README.md',PKG/'pyproject.toml',PKG/'.gitignore',bench/'platform_v3_summary.json',
                  bench/'platform_v3_environment.json',bench/'historical_code_snapshots.json',bench/'three_benchmark_report_zh.html'])
    files.update(bench.glob('archive/legacy_source_blobs/*'))
    files.add(bench/'archive/three_benchmark_report_zh.before_platform_v3_20261106.html')
    for name in ('verify_freezes.py','run_generators.py','freeze_platform_v3.py','replay_legacy_source.py',
                 'verify_platform_v3.py','evaluate_platform_v3.py','build_platform_report_zh.py',
                 'validate_human_screen.py','reproduce_platform_v3.sh'):
        files.add(bench/'scripts'/name)
    files.update(p for p in model.iterdir() if p.is_file())
    files.add(PKG/'models/default_policy.json')
    files.update(p for p in out.rglob('*') if p.is_file() and 'work' not in p.parts)
    files.update(Path(p) for p in registry['sources'])
    files.update(Path(p) for p in evaluation['sources'])
    reference_path=PKG/'references/gencode_v50/built/screen_manifest.json'
    reference=read(reference_path);files.add(reference_path)
    for key in ('reference_fasta','utr3_fasta','seed_index'):
        files.add(reference_path.parent/reference[key])
    files.update(reference_path.parent/name for name in reference['bowtie_index_artifacts'])
    files.update(Path(p) for p in reference['source_artifacts'])
    tools=PKG/'tools/bowtie-1.3.1'
    files.add(tools/'installation.json')
    for name in ('bowtie','bowtie-build','bowtie-align-s','bowtie-align-l','bowtie-build-s','bowtie-build-l'):
        files.update(tools.rglob(name))
    def key(p):
        import os
        return os.path.relpath(p,bench)
    artifacts={key(p):sha256_file(p) for p in sorted(files)}
    lock={'schema_version':1,'status':'development_only','version':'2.0.0',
        'independent_validation_passed':False,'historical_locks_unchanged':True,
        'artifacts':artifacts}
    write(bench/'platform_v3_freeze.json',lock)
    print(json.dumps({'artifacts':len(artifacts),'lock':str(bench/'platform_v3_freeze.json')}))


if __name__=='__main__':main()
