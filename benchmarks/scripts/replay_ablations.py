"""Replay fixed model/native evidence for all planned component removals."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.selection import from_record,rank_candidates,write_selection
from sirna_consensus.benchmark import load_measurements,evaluate_record
from sirna_consensus.fasta import sha256_file
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--development',required=True,type=Path);p.add_argument('--benchmark',required=True,type=Path)
p.add_argument('--out',required=True,type=Path);p.add_argument('--random-repeats',type=int,default=1000)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
measurements=load_measurements(a.benchmark)
conditions={'full':('quality','soft',('tools','support','quality'),'research_top'),
 'no_filters':('quality','off',('tools','support','quality'),'research_top'),
 'hard_filters':('quality','hard',('tools','support','quality'),'research_top'),
 'no_tool_scores':('quality','soft',('support','quality'),'research_top'),
 'no_quality':('quality','soft',('tools','support'),'research_top'),
 'support':('support','soft',('tools','support','quality'),'research_top'),
 'no_support':('support','soft',('tools','quality'),'research_top'),
 'spacing30':('quality','soft',('tools','support','quality'),'diverse_top')}
results={name:[] for name in conditions}
for path in sorted(a.development.glob('*/quality/selection.json')):
 frozen=json.loads(path.read_text());rows=[{'core':from_record(r),'tools':r['tools'],'tool_scores':r['tool_scores']} for r in frozen['candidates']]
 scores={r['key']:r['model_score'] for r in frozen['candidates'] if r['model_score'] is not None}
 for name,(strategy,policy,components,panel) in conditions.items():
  ranked=rank_candidates(rows,strategy,policy,scores,components=components)
  selected=write_selection(ranked,a.out/'runs'/path.parents[1].name/name,
                           {**frozen['provenance'],'ablation':name,'source_selection_sha256':sha256_file(path)},10)
  for k in (5,10):results[name]+=evaluate_record(selected,measurements,k,'measured',panel,a.random_repeats)
for name,rows in results.items():
 (a.out/(name+'.json')).write_text(json.dumps({'results':rows,'benchmark_sha256':sha256_file(a.benchmark),
  'config':{'threshold':.7,'track':'measured','panel':conditions[name][3]},'independent_validation_passed':False},indent=2)+'\n')
print(json.dumps({'ablations':len(results),'targets':len(results['full'])//2}))
