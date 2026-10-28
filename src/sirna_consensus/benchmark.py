"""Experiment-aware exact-guide benchmarking; no fuzzy efficacy transfer."""
from __future__ import annotations

import csv
import json
import math
import random
from collections import defaultdict
from pathlib import Path
from statistics import mean

from .types import to_dna

REQUIRED = ('record_id', 'study', 'gene', 'assay', 'transcript_id', 'target_sha256',
            'guide_dna', 'duplex_format', 'endpoint', 'value')
THRESHOLDS = {'percent_inhibition': 70., 'normalized_activity': .7}


def load_measurements(path: Path, sensitivity: float = .7) -> list[dict]:
    if sensitivity not in (.6, .7, .8):
        raise ValueError('only prespecified 0.6/0.7/0.8 thresholds supported')
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not set(REQUIRED) <= set(reader.fieldnames or []):
            raise ValueError('measurement table requires: ' + ', '.join(REQUIRED))
        rows = list(reader)
    seen, groups = set(), defaultdict(list)
    for row in rows:
        if any(not row.get(k, '').strip() for k in REQUIRED):
            raise ValueError('missing measurement field: ' + str(row.get('record_id')))
        if row['record_id'] in seen:
            raise ValueError('duplicate measurement record_id')
        seen.add(row['record_id'])
        row['guide_dna'] = to_dna(row['guide_dna'])
        if len(row['guide_dna']) not in (19,20,21) or set(row['guide_dna']) - set('ACGT'):
            raise ValueError('invalid measured guide')
        if row['endpoint'] not in THRESHOLDS:
            raise ValueError('endpoint must have an explicit validated transform')
        value = float(row['value'])
        if not math.isfinite(value):
            raise ValueError('nonfinite measurement')
        # Normalized reporter values may exceed one; never clip them.
        row['value'] = value
        key = tuple(row[k] for k in ('study','gene','assay','transcript_id','target_sha256',
                                      'guide_dna','duplex_format','endpoint'))
        groups[key].append(row)
    out = []
    for records in groups.values():
        if len(records) > 1:
            replicate = {r.get('replicate_group', '') for r in records}
            if len(replicate) != 1 or not next(iter(replicate)):
                raise ValueError('repeated guide/condition needs explicit technical replicate_group')
        row = dict(records[0]); row['value'] = mean(r['value'] for r in records)
        row['record_ids'] = [r['record_id'] for r in records]
        row['effective'] = row['value'] >= sensitivity * (100 if row['endpoint']=='percent_inhibition' else 1)
        out.append(row)
    return out


def metric(top: list[dict], measured: list[dict], k: int) -> dict:
    if k < 1:
        raise ValueError('k must be positive')
    identity = lambda r: (r['guide_dna'], r['duplex_format'])
    table = {identity(r): r for r in measured}
    unique, seen = [], set()
    for row in top:
        key = identity(row)
        if key not in seen:
            unique.append(row); seen.add(key)
        if len(unique) == k:
            break
    matched = [table[identity(r)] for r in unique if identity(r) in table]
    effective = sum(r['effective'] for r in matched)
    positives = sum(r['effective'] for r in measured)
    return {'requested_k': k, 'returned': len(unique), 'measured': len(matched),
            'effective': effective, 'available_measured': len(measured), 'available_effective': positives,
            'complete_panel': len(unique)==k,
            'complete_measured_panel': len(matched)==k,
            'measurement_coverage': len(matched)/len(unique) if unique else 0.,
            'fill_rate': len(unique)/k,
            'precision_at_k': effective/len(matched) if matched else None,
            'recovery_at_k': effective/positives if positives else None,
            # A slot utility, NOT a label assigned to an untested prediction.
            'observed_hits_per_budget': effective/k}


def discrimination(rows: list[dict], scores: dict) -> dict:
    pairs = [(scores[(r['guide_dna'],r['duplex_format'])], int(r['effective'])) for r in rows
             if (r['guide_dna'],r['duplex_format']) in scores]
    positives = sum(y for _,y in pairs); negatives = len(pairs)-positives
    continuous = [(scores[(r['guide_dna'],r['duplex_format'])], r.get('value')) for r in rows
                  if (r['guide_dna'],r['duplex_format']) in scores and r.get('value') is not None]
    def correlation(xs, ys):
        if len(xs)<2:return None
        mx,my=mean(xs),mean(ys)
        denominator=(sum((x-mx)**2 for x in xs)*sum((y-my)**2 for y in ys))**.5
        return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/denominator if denominator else None
    from .selection import percentiles
    rx=percentiles({str(i):x for i,(x,y) in enumerate(continuous)})
    ry=percentiles({str(i):y for i,(x,y) in enumerate(continuous)})
    result = {'scored': len(pairs), 'auroc': None, 'average_precision': None,
              'spearman_raw_activity':correlation(list(rx.values()),[ry[k] for k in rx])}
    if not positives or not negatives:
        return result
    pos = [s for s,y in pairs if y]; neg = [s for s,y in pairs if not y]
    result['auroc'] = sum((a>b)+.5*(a==b) for a in pos for b in neg)/(positives*negatives)
    buckets = defaultdict(list)
    for s,y in pairs:
        buckets[s].append(y)
    tp = total = 0; ap = 0.
    for score in sorted(buckets, reverse=True):
        labels = buckets[score]; added = sum(labels)
        tp += added; total += len(labels)
        ap += (added/positives)*(tp/total)
    result['average_precision'] = ap
    return result


def evaluate_record(record: dict, measurements: list[dict], k=5, track='measured',
                    panel='research_top', random_repeats=1000, seed=2026) -> list[dict]:
    if panel not in ('research_top','screened_top','diverse_top'):
        raise ValueError('unknown panel')
    if record.get('schema_version') != 2:
        raise ValueError('legacy runs must be rerun before efficacy evaluation')
    if track not in ('measured', 'deployment'):
        raise ValueError('unknown evaluation track')
    if track == 'deployment' and k > record.get('config', {}).get('top_k', 0):
        raise ValueError('requested deployment K exceeds the saved selection budget; regenerate panels')
    candidates = record['candidates']
    if record.get('kind') != 'exact_guide_selection':
        raise ValueError('use select/replay to create selection.json before evaluation')
    target_sha = record['provenance']['target_sha256']
    tid = record['provenance']['transcript_id']
    experiments = defaultdict(list)
    for row in measurements:
        if row['target_sha256'] == target_sha and row['transcript_id'] == tid:
            experiments[(row['study'],row['gene'],row['assay'],row['endpoint'],row['duplex_format'])].append(row)
    output = []
    by_key = {row['key']:row for row in candidates}
    for unit, measured in sorted(experiments.items()):
        table = {(r['guide_dna'],r['duplex_format']) for r in measured}
        groups = sorted({row['group'] for row in candidates if row['duplex_format']==unit[-1]})
        if not groups:
            groups = [unit[-1]+'/unavailable']
        for group in groups:
            ranked = [r for r in candidates if r['group']==group and r['eligible']]
            if panel == 'screened_top':
                ranked = [r for r in ranked if r.get('screen',{}).get('status')=='pass']
            scored = {}
            for r in ranked:
                scored.setdefault((r['guide_dna'], r['duplex_format']), r['score'])
            if track == 'measured':
                pool = [r for r in ranked if (r['guide_dna'],r['duplex_format']) in table]
                # Match BEFORE top-k for measured-candidate ranking; never for deployment.
                selected = pool
                if panel == 'diverse_top':
                    selected = []
                    for r in pool:
                        if all(abs(r['start']-other['start'])>=30 for other in selected): selected.append(r)
                        if len(selected)==k: break
            else:
                selected = [by_key[key] for key in record.get(panel,{}).get(group,[])][:k]
                pool = ranked
            m = metric(selected, measured, k)
            rng = random.Random(seed)
            random_values = []
            pool_by_guide = {}
            for r in pool:
                pool_by_guide.setdefault((r['guide_dna'],r['duplex_format']), r)
            unique_pool = list(pool_by_guide.values())
            for _ in range(random_repeats):
                if panel == 'diverse_top':
                    sample=[]
                    for r in rng.sample(unique_pool,len(unique_pool)):
                        if all(abs(r['start']-other['start'])>=30 for other in sample):sample.append(r)
                        if len(sample)==k:break
                else:
                    sample = rng.sample(unique_pool, min(k,len(unique_pool)))
                value = metric(sample, measured, k)['precision_at_k']
                if value is not None:
                    random_values.append(value)
            base = mean(random_values) if random_values else None
            output.append({'study':unit[0], 'gene':unit[1], 'assay':unit[2], 'endpoint':unit[3],
                           'duplex_format':unit[4], 'transcript_id':tid, 'target_sha256':target_sha, 'group':group, 'track':track, **m,
                           'candidate_coverage': len(set(scored)&table)/len(table) if table else 0.,
                           'random_precision':base,
                           'enrichment': m['precision_at_k']/base if base and m['precision_at_k'] is not None else None,
                           **discrimination(measured,scored)})
    return output


def paired_summary(left: list[dict], right: list[dict], replicates=10000, seed=2026) -> dict:
    """Pair identical assays, then average within genes; never pool guide rows."""
    identity = lambda r: tuple(r[k] for k in ('study','gene','assay','endpoint','duplex_format','track','requested_k','transcript_id','target_sha256'))
    a,b = {identity(r):r for r in left},{identity(r):r for r in right}
    if len(a)!=len(left) or len(b)!=len(right):
        raise ValueError('comparison must contain one method/group per assay unit')
    if set(a)!=set(b):
        raise ValueError('paired comparison requires the same prespecified assay units; include abstentions')
    differences, yields = defaultdict(list), defaultdict(list)
    for unit in sorted(a):
        x,y = a[unit], b[unit]
        gene = unit[1]
        yields[gene].append(x['observed_hits_per_budget']-y['observed_hits_per_budget'])
        if x['complete_measured_panel'] and y['complete_measured_panel']:
            differences[gene].append(x['precision_at_k']-y['precision_at_k'])
    values = [mean(v) for _,v in sorted(differences.items())]
    rng = random.Random(seed)
    boot = sorted(mean(rng.choices(values,k=len(values))) for _ in range(replicates)) if values else []
    lo = boot[int(.025*(len(boot)-1))] if boot else None
    hi = boot[int(.975*(len(boot)-1))] if boot else None
    completion = lambda d: mean(mean(r['complete_panel'] for r in d.values() if r['gene']==g)
                              for g in sorted({r['gene'] for r in d.values()})) if d else 0.
    delta = mean(values) if values else None
    gain = mean(mean(v) for v in yields.values()) if yields else None
    completion_delta = completion(a)-completion(b)
    return {'minimum_genes':5, 'evidence_status':'evaluable' if len(values)>=5 else 'insufficient_genes', 'paired_genes':len(values), 'eligible_genes':len(yields), 'mean_difference':delta,
            'ci95':[lo,hi], 'completion_difference':completion_delta,
            'budget_yield_difference':gain, 'per_gene_difference':{g:mean(v) for g,v in differences.items()},
            'numeric_criterion_met': bool(len(values) >= 5 and delta is not None and delta>=.1 and lo>0 and
                                          completion_delta>=0 and gain is not None and gain>0),
            'independent_validation_passed':False,
            'note':'numeric criteria alone do not establish holdout independence'}


def permutation_control(record, measurements, k=5, permutations=1000, seed=2026):
    """Shuffle aggregated molecular measurements within study/gene/assay/format."""
    strata = {(r['study'],r['endpoint'],r['duplex_format']) for r in measurements}
    if len(strata) > 1:
        raise ValueError('permutation control requires one study/endpoint/format stratum')
    rng = random.Random(seed)
    units = defaultdict(list)
    for i,r in enumerate(measurements):
        units[tuple(r[x] for x in ('study','gene','assay','endpoint','duplex_format'))].append(i)
    def statistic(rows):
        genes = defaultdict(list)
        for r in rows:
            if r['complete_measured_panel']:
                genes[r['gene']].append(r['precision_at_k'])
        return mean(mean(v) for v in genes.values()) if genes else None
    real = statistic(evaluate_record(record,measurements,k,random_repeats=0))
    null=[]
    for _ in range(permutations):
        shuffled=[dict(r) for r in measurements]
        for indexes in units.values():
            labels=[measurements[i]['effective'] for i in indexes];rng.shuffle(labels)
            for i,label in zip(indexes,labels):shuffled[i]['effective']=label
        value=statistic(evaluate_record(record,shuffled,k,random_repeats=0))
        if value is not None:null.append(value)
    return {'observed':real,'permutations':len(null),
            'p_one_sided':(1+sum(v>=real for v in null))/(1+len(null)) if real is not None and null else None,
            'null_mean':mean(null) if null else None}
