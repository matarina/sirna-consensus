"""Two complete designs screened against the pinned ALL human reference."""
import argparse
import json
import sys
from collections import Counter
from pathlib import Path

PKG=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(PKG/'src'))
from sirna_consensus.pipeline import Config,run
from sirna_consensus.fasta import sha256_file


def main():
    p=argparse.ArgumentParser();p.add_argument('--manifest',required=True);p.add_argument('--bowtie',required=True)
    p.add_argument('--out',required=True);p.add_argument('--model-dir',required=True);a=p.parse_args()
    out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    manifest_path=Path(a.manifest).resolve();manifest=json.loads(manifest_path.read_text())
    source=PKG.parent/'independent_validation/denovo_two_benchmarks_20261105'
    record={'status':'running','manifest_sha256':sha256_file(manifest_path),'reference':{
        'release':manifest['release'],'transcripts':len(manifest['transcript_genes']),
        'genes':len(set(manifest['transcript_genes'].values())),
        'utr3_unprojectable':len(manifest.get('utr3_unprojectable',[])),
        'utr3_complete':manifest['utr3_complete'],'gene_identity':manifest.get('gene_identity','GENCODE_ID'),
        'reference_sha256':manifest['reference_sha256'],'utr3_sha256':manifest['utr3_sha256'],
        'source_artifacts':manifest['source_artifacts'],'build_command':manifest['build_command'],
        'build_version':manifest['build_version']},'runs':[]}
    for study,gene in (('Filhol2012','BCL2L1'),('Wang2009','HIF1A')):
        directory=out/'human_screen'/study/gene
        result=run(Config(str(source/study/'transcripts'/f'{gene}.fa'),str(directory),tools=(),
            scorer_manifest=str(Path(a.model_dir).resolve()/'scorer_manifest.json'),
            screen_manifest=str(manifest_path),target_gene=gene,bowtie=str(Path(a.bowtie).resolve())))
        item={'study':study,'gene':gene,'counts':result['counts'],
            'candidate_screen_status':dict(Counter(r['screen']['status'] for r in result['candidates'])),
            'screen_metadata':result['provenance']['screen'],'run':str(directory/'run.json'),
            'run_sha256':sha256_file(directory/'run.json')}
        record['runs'].append(item)
        (out/'screening_integration.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps({'gene':gene,'status':item['candidate_screen_status']}),flush=True)
    passed=all(r['screen_metadata'].get('status')=='complete' and r['screen_metadata'].get('seed_scope')=='3UTR'
               and r['counts']['screened_diverse_top']>0 for r in record['runs']) and manifest['utr3_complete']
    record['status']='passed' if passed else 'completed_with_unverified_or_empty_panels'
    record['independent_efficacy_validation']=False
    (out/'screening_integration.json').write_text(json.dumps(record,indent=2)+'\n')
    if not passed:raise RuntimeError('reference integration needs review; see screening_integration.json')


if __name__=='__main__':main()
