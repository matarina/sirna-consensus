"""Reuse source-curated Huesken constructs; preserve full guides and raw activity."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.types import Core, to_dna
from sirna_consensus.fasta import sha256_file, sha256


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,required=True,help='sibling sirnaforge-platform/data')
    p.add_argument('--out',type=Path,required=True)
    args=p.parse_args(); out=args.out;out.mkdir(parents=True,exist_ok=True)
    inserts=json.loads((args.source/'huesken_inserts.json').read_text())
    inserts={r['gene']:r for r in inserts if r['species']=='Homo sapiens'}
    with (args.source/'curated_measurements.tsv').open() as f:raw=list(csv.DictReader(f,delimiter='\t'))
    records,candidates,excluded=[],[],[]
    transcripts=out/'transcripts';transcripts.mkdir(exist_ok=True)
    for r in raw:
        insert=inserts[r['gene']];seq=to_dna(insert['sequence'])
        tid=r['gene']+'|'+r['transcript_id']+'|cloned_insert'
        guide=to_dna(r['guide_original']);sense=to_dna(r['target_sequence']);start=int(r['insert_start'])+1
        if set(seq)-set('ACGT'):
            excluded.append({'record_id':r['record_id'],'reason':'ambiguous_construct_sequence'});continue
        if seq[start-1:start-1+19]!=sense or len(guide)!=21:
            raise ValueError('source molecular mapping mismatch: '+r['record_id'])
        core=Core(tid,start,19,sense,ordered_guide=guide,guide_overhang=to_dna(r['guide_overhang']))
        target_sha=sha256(seq)
        (transcripts/(r['gene']+'.fa')).write_text('>'+tid+'\n'+seq+'\n')
        records.append({'record_id':r['record_id'],'study':'huesken2005','gene':r['gene'],
                        'assay':'H1299_reporter_unspecified_dose_time','transcript_id':tid,
                        'target_sha256':target_sha,'guide_dna':guide,'duplex_format':core.duplex_format,
                        'endpoint':'normalized_activity','value':r['normalized_inhibitory_activity'],
                        'replicate_group':'','source_url':r['source_url'],'split':'development',
                        'passenger':'','guide_overhang':core.guide_overhang,
                        'mapping_status':'original_construct_exact','source_record':r['record_id']})
        candidates.append({'transcript_id':tid,'start':start,'length':19,'sense':sense,
                           'ordered_guide':guide,'guide_overhang':core.guide_overhang,
                           'passenger':'','passenger_overhang':''})
    # Duplicate molecules are not independent observations; require source review.
    seen=set();unique_records=[];unique_candidates=[]
    for r,c in zip(records,candidates):
        key=(r['transcript_id'],r['guide_dna'])
        if key in seen:
            excluded.append({'record_id':r['record_id'],'reason':'duplicate_molecule_requires_replicate_review'});continue
        seen.add(key);unique_records.append(r);unique_candidates.append(c)
    for name,rows in [('measurements.tsv',unique_records),('candidates.tsv',unique_candidates)]:
        with (out/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    (out/'exclusions.json').write_text(json.dumps(excluded,indent=2)+'\n')
    manifest={'dataset':'Huesken original human cloned constructs','split':'development',
              'previous_project_exposure':True,'independent_holdout':False,
              'endpoint':'raw control-normalized inhibitory activity, not percent knockdown',
              'source_hashes':{name:sha256_file(args.source/name) for name in
                               ('curated_measurements.tsv','huesken_inserts.json','source_manifest.json')},
              'records':len(unique_records),'genes':len({r['gene'] for r in unique_records}),
              'exclusions':len(excluded),'artifacts':{str(p.relative_to(out)):sha256_file(p)
              for p in sorted(out.rglob('*')) if p.is_file() and p.name!='dataset_manifest.json'}}
    (out/'dataset_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ('records','genes','exclusions')}))

if __name__=='__main__':main()
