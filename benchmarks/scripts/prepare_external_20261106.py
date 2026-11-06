"""Curate ORIGINAL Katoh S1 and freeze eligibility before any model predictions.
Run with the acquisition Python (xlrd/numpy); never trains or scores a model.
"""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import xlrd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'benchmarks/external_confirmation_20261106'
REF = ROOT.parent / 'independent_validation/sciabola/adaptation_20261105/external/refs'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')


def rc(s):
    return s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def main():
    protocol = json.loads((OUT / 'acquisition_protocol.json').read_text())
    for p, sha in protocol['artifacts'].items():
        assert digest(p) == sha, p
    sheet = xlrd.open_workbook(OUT / 'raw/katoh_original_table_s1.xls').sheet_by_index(0)
    assert sheet.nrows == 704 and sheet.ncols == 7
    rows = [sheet.row_values(i) for i in range(2, sheet.nrows)]
    target = ''.join(r[1][0] for r in rows) + rows[-1][1][1:]
    target = target.replace('U', 'T')
    assert len(target) == 720
    historical = []
    source_paths = set()
    for p in ['models/development_20261106/development_records.json',
              'models/huesken_only_20261106/development_records.json']:
        path = ROOT / p
        source_paths.add(path)
        for r in json.loads(path.read_text()):
            historical.append((r['record_id'], r['core'], p))
    registry_path = ROOT / 'benchmarks/confirmation_20261106/cohorts.json'
    source_paths.add(registry_path)
    for cohort in json.loads(registry_path.read_text())['cohorts']:
        if not cohort.get('records'):
            continue
        path = Path(cohort['path'])
        source_paths.add(path)
        assert digest(path) == cohort['sha256']
        for r in csv.DictReader(path.open(), delimiter='\t'):
            historical.append((r['record_id'], rc(r['guide_dna'].replace('U', 'T')[:19]), cohort['name']))
    # Extra conservative public upstream references, sequence-only, no labels.
    for name in ['Hu', 'Mix']:
        path = REF / (name + '.csv')
        source_paths.add(path)
        for i, r in enumerate(csv.DictReader(path.open())):
            historical.append((f'{name}:{i}', rc(r['siRNA'].replace('U', 'T')[:19]), name))
    oriented = defaultdict(list)
    for ident, core, origin in historical:
        assert len(core) == 19 and not set(core) - set('ACGT'), (ident, core)
        for s in (core, rc(core)):
            oriented[s].append([origin, ident])
    sequences = sorted(oriented)
    matrix = np.asarray([list(s.encode()) for s in sequences], dtype=np.uint8)
    kmers = defaultdict(set)
    for s in sequences:
        for i in range(3):
            kmers[s[i:i+17]].add(s)
    eligible, excluded, audit, all_records = [], [], [], []
    for i, row in enumerate(rows):
        pos, sense, activity = int(row[0]), row[1].replace('U', 'T'), float(row[5])
        assert pos == i+1 and target[i:i+19] == sense
        guide = rc(sense)
        record = dict(record_id=f'Katoh2007:EGFP:{pos:03d}', study='Katoh2007',
                      gene='EGFP', gene_id='REPORTER:EGFP',
                      unit='Katoh2007|REPORTER:EGFP|HeLa_EGFP_10nM_72h_flow',
                      condition='HeLa_EGFP_10nM_72h_flow', endpoint='normalized_activity',
                      value=activity/100, raw_activity=activity, active_threshold=0.7,
                      active_definition='normalized assay activity >=0.7; NOT absolute 70% knockdown',
                      start=pos, guide=guide, core=sense, passenger_core=sense,
                      passenger='', guide_overhang='', passenger_overhang='', chemistry='all_RNA',
                      full_molecule_known=False, transcript_id='Katoh2007_EGFP_ORF',
                      target_sha256=hashlib.sha256(target.encode()).hexdigest(),
                      source='10.1093/nar/gkl1120 Table S1', source_row=i+3,
                      mapping_status='original_19core_exact',
                      scope='measured paired-core ranker component; not full-molecule efficacy')
        all_records.append(record)
        distance = (matrix != np.frombuffer(sense.encode(), dtype=np.uint8)).sum(axis=1)
        near = [sequences[j] for j in np.flatnonzero(distance <= 2)]
        shared = sorted(set().union(*(kmers.get(sense[j:j+17], set()) for j in range(3))))
        reasons = ([] if not near else ['forward_or_RC_Hamming_le2']) + ([] if not shared else ['shared_shifted_17mer'])
        entry = dict(record_id=record['record_id'], min_hamming=int(distance.min()),
                     near_matches=[dict(sequence=s, provenance=oriented[s]) for s in near],
                     shared17_matches=[dict(sequence=s, provenance=oriented[s]) for s in shared],
                     eligible=not reasons, reasons=reasons)
        audit.append(entry)
        if reasons:
            excluded.append(entry)
        else:
            eligible.append(record)
    assert len(eligible)+len(excluded) == 702
    write('all_records.json', all_records)
    write('records.json', eligible)
    write('exclusions.json', excluded)
    write('sequence_audit.json', dict(original_records=702, eligible=len(eligible), excluded=len(excluded),
          exclusion_counts=dict(Counter(reason for r in excluded for reason in r['reasons'])),
          reference_rows=len(historical), reference_unique_oriented_cores=len(sequences),
          reference_hashes={str(p): digest(p) for p in sorted(source_paths)}, audit=audit,
          caveats=['Taka reference previously accessed only for sequence audit; not a never-accessed blind set.',
                   'EGFP occurred in prior modified Martinelli evaluation; this is new-study, not never-seen-target evidence.',
                   'Adjacent tiling windows are dependent: one gene/one study, no row bootstrap CI.',
                   'No upstream neural checkpoint is used; cannot claim Katoh independence for OligoFormer.']))
    (OUT / 'egfp.fasta').write_text('>Katoh2007_EGFP_ORF\n'+target+'\n')
    # These are observed paired cores, not fabricated complete 21mer duplexes.
    with (OUT / 'observed_cores.tsv').open('w') as f:
        w = csv.DictWriter(f, delimiter='\t', fieldnames=['transcript_id','start','length','sense','chemistry'])
        w.writeheader()
        for r in eligible:
            w.writerow(dict(transcript_id=r['transcript_id'],start=r['start'],length=19,sense=r['core'],chemistry='all_RNA'))
    print(json.dumps(dict(raw=702, eligible=len(eligible), excluded=len(excluded))))


if __name__ == '__main__':
    main()
