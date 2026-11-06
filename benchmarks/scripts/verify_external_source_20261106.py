"""Independent original-XLS and near-sequence audit; no model predictions.
Use acquisition Python with xlrd, independent of scorer environment.
"""
import csv
import hashlib
import itertools
import json
from collections import defaultdict
from pathlib import Path

import xlrd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'benchmarks/external_confirmation_20261106'

def rc(s):
    return ''.join({'A':'T','C':'G','G':'C','T':'A'}[c] for c in reversed(s))


def main():
    frozen=json.loads((OUT/'test_protocol.json').read_text())
    data=json.loads((OUT/'records.json').read_text())
    audit=json.loads((OUT/'sequence_audit.json').read_text())
    s=xlrd.open_workbook(OUT/'raw/katoh_original_table_s1.xls').sheet_by_index(0)
    by_pos={int(r['start']):r for r in data}
    source_ok=True
    for row in range(2,s.nrows):
        pos=int(s.cell_value(row,0)); sense=s.cell_value(row,1).replace('U','T')
        r=by_pos[pos]
        source_ok &= (r['source_row']==row+1 and r['core']==sense and r['guide']==rc(sense) and
                      abs(r['value']-s.cell_value(row,5)/100)<1e-15)
    reference=[]
    training_gene_set=set()
    for path in audit['reference_hashes']:
        p=Path(path)
        if p.name=='cohorts.json':
            continue
        if p.suffix=='.json':
            rows=json.loads(p.read_text())
            training_gene_set.update(r['gene_id'] for r in rows)
            reference.extend(r['core'] for r in rows)
        elif p.suffix=='.tsv':
            reference.extend(rc(r['guide_dna'].replace('U','T')[:19]) for r in csv.DictReader(p.open(),delimiter='\t'))
        elif p.suffix=='.csv':
            reference.extend(rc(r['siRNA'].replace('U','T')[:19]) for r in csv.DictReader(p.open()))
    # Pigeonhole indexing: <=2 substitutions must preserve one of 3 disjoint
    # blocks. Different algorithm from the curation's full vector distance scan.
    oriented={s for core in reference for s in [core,rc(core)]}
    blocks=defaultdict(set);mers=set()
    for core in oriented:
        for left,right in [(0,6),(6,12),(12,19)]:
            blocks[(left,right,core[left:right])].add(core)
        mers.update(core[j:j+17] for j in range(3))
    bad=[]
    for r in data:
        core=r['core'];candidates=set()
        for left,right in [(0,6),(6,12),(12,19)]:
            candidates.update(blocks[(left,right,core[left:right])])
        near=any(sum(a!=b for a,b in zip(core,t))<=2 for t in candidates)
        shared=any(core[j:j+17] in mers for j in range(3))
        if near or shared: bad.append(r['record_id'])
    checks=dict(original_source_rows_and_activity=source_ok and len(data)==702,
                no_training_gene_overlap=not ({r['gene_id'] for r in data}&training_gene_set),
                independent_no_forward_RC_hamming2_or_shift17=not bad,
                all_original_rows_preserved=len(json.loads((OUT/'all_records.json').read_text()))==s.nrows-2==len(data),
                original_source_hash=hashlib.sha256((OUT/'raw/katoh_original_table_s1.xls').read_bytes()).hexdigest()==frozen['artifacts']['benchmarks/external_confirmation_20261106/raw/katoh_original_table_s1.xls'])
    result=dict(passed=all(checks.values()),checks=checks, reference_records=len(reference), matched_ids=bad,
                note='Original normalized activity /100, not mirror rescaling; no full-duplex identity inferred.')
    (OUT/'outputs/source_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    assert result['passed']


if __name__=='__main__':
    main()
