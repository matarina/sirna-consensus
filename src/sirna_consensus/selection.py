"""Exact-molecule selection and replay, independent of generator execution."""
from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict
from pathlib import Path

from . import filters
from .fasta import sha256_file
from .rank import _quality_term
from .types import Core, to_dna

STRATEGIES = ('native', 'model', 'fusion', 'support', 'quality', 'legacy')


def from_record(row: dict) -> Core:
    row = dict(row)
    row.setdefault("ordered_guide", row.get("guide_sequence_dna", row.get("guide_sequence_rna", "")))
    row.setdefault("passenger", row.get("passenger_sequence_rna", ""))
    return Core(row['transcript_id'], int(row.get('start', row.get('target_start', 0))),
                int(row.get('length', row.get('core_length', 0))),
                to_dna(row.get('sense', row.get('sense_dna', row.get('target_sequence_dna', '')))),
                **{k: to_dna(row.get(k, '')) for k in
                   ('ordered_guide', 'passenger', 'guide_overhang', 'passenger_overhang')})


def load_candidates(run_path: Path) -> tuple[list[dict], dict]:
    run = json.loads(run_path.read_text())
    if run.get('schema_version') != 2:
        raise ValueError('schema-v1 runs are inspection-only; rerun with molecular identity fixes')
    by_key = {}
    for hit in run['normalized_hits']:
        if hit['role'] != 'generator' or not hit['vote_eligible']:
            continue
        core = Core(**hit['core'])
        row = by_key.setdefault(core.key, {'core': core, 'tool_scores': {}, 'tools': []})
        tool = hit['tool']
        if tool not in row['tools']:
            row['tools'].append(tool)
        score = hit['score']
        if score is not None:
            if not math.isfinite(score):
                raise ValueError('nonfinite native score')
            row['tool_scores'][tool] = max(score, row['tool_scores'].get(tool, score))
    return list(by_key.values()), run


def supplied_candidates(path: Path, transcript_id: str, sequence: str) -> list[dict]:
    """Explicit molecules only; overhangs must be stated, never inferred."""
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    out, seen = [], set()
    for row in rows:
        core = from_record(row)
        if core.transcript_id != transcript_id:
            continue
        if sequence[core.start-1:core.end] != core.sense:
            raise ValueError(f'candidate does not match target: {core.key}')
        if core.key in seen:
            raise ValueError(f'duplicate supplied candidate: {core.key}')
        seen.add(core.key)
        out.append({'core': core, 'tools': [], 'tool_scores': {}})
    return out


def enumerate_guides(transcript_id: str, sequence: str, lengths=(19, 20, 21)) -> list[dict]:
    return [{'core': Core(transcript_id, i+1, n, sequence[i:i+n]),
             'tools': [], 'tool_scores': {}} for n in lengths
            for i in range(len(sequence)-n+1)]


def percentiles(values: dict[str, float]) -> dict[str, float]:
    """Midrank percentiles; singleton/constant scores have percentile 0.5."""
    ordered = sorted(values, key=lambda key: (values[key], key))
    result, i = {}, 0
    while i < len(ordered):
        j = i+1
        while j < len(ordered) and values[ordered[j]] == values[ordered[i]]:
            j += 1
        for key in ordered[i:j]:
            result[key] = ((i+j-1)/2 / (len(ordered)-1)) if len(ordered) > 1 else .5
        i = j
    return result


def rank_candidates(rows: list[dict], strategy='native', filter_policy='soft',
                    model_scores: dict | None = None, config=filters.DEFAULT_FILTERS,
                    components=('tools', 'support', 'quality')) -> list[dict]:
    if strategy not in STRATEGIES or filter_policy not in ('hard', 'soft', 'off'):
        raise ValueError('invalid ranking/filter policy')
    model_scores = model_scores or {}
    grouped = {}
    for row in rows:
        grouped.setdefault(row['core'].duplex_format, []).append(row)
    result = []
    for fmt, group in sorted(grouped.items()):
        tools = sorted({t for row in group for t in row['tools']})
        tool_p = {t: percentiles({r['core'].key: r['tool_scores'][t] for r in group
                                  if t in r['tool_scores']}) for t in tools}
        scores = {r['core'].key: model_scores[r['core'].key] for r in group
                  if r['core'].key in model_scores}
        if any(not math.isfinite(v) for v in scores.values()):
            raise ValueError('nonfinite model score')
        model_p = percentiles(scores)
        for row in group:
            core = row['core']; key = core.key
            flags = filters.hard_flags(core, config) if filter_policy != 'off' else []
            native = [tool_p[t][key] for t in tools if key in tool_p[t]]
            tool_score = sum(native)/len(native) if native and 'tools' in components else 0.
            support = len(row['tools']) / max(1, len(tools)) if 'support' in components else 0.
            quality = _quality_term(core) if 'quality' in components else 0.
            m = model_p.get(key)
            # Unsupported candidates are never mixed with model-scored candidates.
            method = strategy if m is not None or strategy in ('native', 'legacy') else 'native_fallback'
            if method in ('native', 'native_fallback'):
                score = tool_score if native else quality
            elif method == 'legacy':
                score = (2*len(row['tools']) if 'support' in components else 0) + .5*tool_score + .5*quality
            elif method == 'model':
                score = m
            elif method == 'fusion':
                terms = [m] + (native if 'tools' in components else [])
                score = sum(terms)/len(terms)
            elif method == 'support':
                score = .8*m + .2*support
            else:
                score = .8*m + .1*tool_score + .1*quality
            result.append({**asdict(core), 'key': key, 'guide_dna': core.guide,
                           'molecule_id': core.molecule_id, 'duplex_format': fmt,
                           'tools': sorted(row['tools']), 'tool_scores': row['tool_scores'],
                           'score': score, 'model_score': scores.get(key),
                           'model_status': 'scored' if m is not None else 'unsupported_or_unavailable',
                           'ranking': method, 'group': fmt + '/' + method,
                           'eligible': not flags or filter_policy != 'hard',
                           'hard_flags': flags if filter_policy == 'hard' else [],
                           'soft_flags': (flags if filter_policy == 'soft' else []) +
                               (filters.soft_flags(core, config) if filter_policy != 'off' else []),
                           'components': {'model': m, 'native': tool_score,
                                          'support': support, 'quality': quality}})
    return sorted(result, key=lambda r: (r['group'], -r['score'], r['molecule_id'], r['start']))


def panels(rows: list[dict], k=5, spacing=0, screened=False) -> dict[str, list[str]]:
    if k < 1 or spacing < 0:
        raise ValueError('k must be positive and spacing nonnegative')
    groups = {row['group']: [] for row in rows}
    seen = {group: set() for group in groups}
    for row in rows:
        kept = groups[row['group']]
        if not row['eligible'] or len(kept) >= k or row['guide_dna'] in seen[row['group']]:
            continue
        if screened and row.get('screen', {}).get('status') != 'pass':
            continue
        if all(abs(row['start']-other['start']) >= spacing for other in kept):
            kept.append(row); seen[row['group']].add(row['guide_dna'])
    return {group: [row['key'] for row in kept] for group, kept in groups.items()}


def write_selection(rows, out: Path, provenance: dict, k=5):
    out.mkdir(parents=True, exist_ok=True)
    record = {'schema_version': 2, 'kind': 'exact_guide_selection', 'config': {'top_k': k},
              'provenance': provenance, 'candidates': rows,
              'implementation_hashes': {str(p.relative_to(Path(__file__).parent)):sha256_file(p)
                                        for p in sorted(Path(__file__).parent.rglob('*.py'))},
              'research_top': panels(rows, k), 'screened_top': panels(rows, k, screened=True),
              'diverse_top': panels(rows, k, spacing=30),
              'independent_validation_passed': False}
    (out/'selection.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    by_key = {row['key']: row for row in rows}
    for name in ('research_top', 'screened_top', 'diverse_top'):
        with (out/(name+'.fa')).open('w') as handle:
            for group, keys in record[name].items():
                for rank, key in enumerate(keys, 1):
                    row = by_key[key]
                    handle.write(f'>{key} rank={rank} group={group}\n{row["guide_dna"].replace("T", "U")}\n')
    return record
