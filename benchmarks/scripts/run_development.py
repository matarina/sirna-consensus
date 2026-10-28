"""Reproducible development replay; never makes an independent efficacy claim."""
import argparse,json,sys
from pathlib import Path
from statistics import mean
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.selection import supplied_candidates,rank_candidates,write_selection,load_candidates
from sirna_consensus.fasta import read_transcript,sha256,sha256_file
from sirna_consensus.benchmark import load_measurements,evaluate_record,paired_summary
from sirna_consensus.scoring import score_candidates

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--data',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
p.add_argument('--tool-runs',type=Path,help='per-gene schema-v2 runs containing native generator evidence')
p.add_argument('--scorer-manifest',type=Path);p.add_argument('--genes',help='optional comma-separated development smoke subset')
p.add_argument('--random-repeats',type=int,default=1000)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
measured=load_measurements(a.data/'measurements.tsv');all_results={}
strategies=('native', 'sirnadesigner', 'sirnaforge', 'consensus') if a.tool_runs else ('native',)
if a.scorer_manifest: strategies += ('model','fusion','support','quality') if a.tool_runs else ('model','quality')
chosen=set(a.genes.split(',')) if a.genes else None
for fasta in sorted((a.data/'transcripts').glob('*.fa')):
 if chosen and fasta.stem not in chosen:continue
 tid,sequence=read_transcript(fasta);rows=supplied_candidates(a.data/'candidates.tsv',tid,sequence)
 if a.tool_runs:
  native,native_run=load_candidates(a.tool_runs/fasta.stem/'run.json')
  if native_run['input_sha256']!=sha256(sequence):raise ValueError('native evidence target mismatch')
  for row in rows:
   core=row['core']
   for hit in native:
    h=hit['core']
    if h.guide==core.guide and h.sense==core.sense and h.guide_overhang==core.guide_overhang:
     row['tools']=sorted(set(row['tools'])|set(hit['tools']))
     row['tool_scores'].update(hit['tool_scores'])
 scores={};scorer={}
 if a.scorer_manifest:
  scores,scorer=score_candidates(rows,sequence,a.scorer_manifest.resolve(),(a.out/fasta.stem/'scorer').resolve())
 for strategy in strategies:
  method_rows=rows
  ranking=strategy
  if strategy in ('sirnadesigner','sirnaforge'):
   method_rows=[{**r,'tools':[strategy],'tool_scores':{strategy:r['tool_scores'][strategy]} if strategy in r['tool_scores'] else {}} for r in rows if strategy in r['tools']]
   ranking='native'
  elif strategy=='consensus':
   method_rows=[r for r in rows if len(r['tools'])>=2];ranking='legacy'
  ranked=rank_candidates(method_rows,ranking,filter_policy='soft',model_scores=scores if strategy not in ('native','sirnadesigner','sirnaforge','consensus') else {})
  record=write_selection(ranked,a.out/fasta.stem/strategy,
                         {'transcript_id':tid,'target_sha256':sha256(sequence),'strategy':strategy,
                          'scorer':scorer,'data_split':'development','source_sha256':sha256_file(a.data/'candidates.tsv')},10)
  all_results.setdefault(strategy,[]).extend(evaluate_record(record,measured,k,'measured',random_repeats=a.random_repeats)
                                           for k in (5,10))
 print(fasta.stem, len(rows),'candidates',flush=True)
for strategy,batches in all_results.items():
 rows=[r for batch in batches for r in batch]
 (a.out/(strategy+'.json')).write_text(json.dumps({'results':rows,'split':'development','benchmark_sha256':sha256_file(a.data/'measurements.tsv'),
   'config':{'threshold':.7,'panel':'research_top','track':'measured'},'independent_validation_passed':False},indent=2)+'\n')
# Choose only on this explicitly exposed development cohort. No inference about final improvement.
summary={}
for strategy in all_results:
 rows=json.loads((a.out/(strategy+'.json')).read_text())['results']
 per_gene={}
 for r in rows:
  if r['requested_k']==5 and r['complete_measured_panel']:
   per_gene.setdefault(r['gene'],[]).append(r['precision_at_k'])
 summary[strategy]={'mean_gene_precision_at_5':mean(mean(v) for v in per_gene.values()) if per_gene else None,
                    'complete_genes':len(per_gene)}
valid=[s for s in strategies if summary[s]['mean_gene_precision_at_5'] is not None]
winner=max(valid,key=lambda s:(summary[s]['complete_genes'],summary[s]['mean_gene_precision_at_5'],-strategies.index(s))) if valid else None
(a.out/'development_summary.json').write_text(json.dumps({'strategies':summary,'development_winner':winner,
 'candidate_universe':'supplied assayed molecules; no generator emission is abstention evidence',
 'evaluation_scope':'native evidence replay plus frozen scorer' if a.tool_runs else 'sequence-quality and scorer smoke; not a generator comparison',
 'native_definition':'union with sequence-quality fallback for non-emitted guides' if a.tool_runs else 'sequence-quality-only baseline',
 'independent_validation_passed':False},indent=2)+'\n')

if winner and a.scorer_manifest:
 formats=sorted({r['duplex_format'] for r in measured})
 (a.out/'selection_policy.json').write_text(json.dumps({'schema_version':1,'status':'development_only',
  'strategies':{fmt:winner for fmt in formats},
  'scorer_manifest_sha256':sha256_file(a.scorer_manifest),
  'measurement_sha256':sha256_file(a.data/'measurements.tsv'),
  'independent_validation_passed':False},indent=2)+'\n')
