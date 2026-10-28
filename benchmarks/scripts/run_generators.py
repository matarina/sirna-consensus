"""Freeze both native generators on original constructs; no efficacy labels read."""
import argparse,json,sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.pipeline import Config,run
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--transcripts',required=True,type=Path);p.add_argument('--out',required=True,type=Path)
p.add_argument('--workers',type=int,default=2)
a=p.parse_args()
def execute(fasta):
 result=run(Config(str(fasta),str(a.out/fasta.stem),tools=('sirnadesigner','sirnaforge'),min_votes=1,top_k=10,site_mode='exact'))
 return {'gene':fasta.stem,'tools':result['tools'],'counts':result['counts']}
with ThreadPoolExecutor(max_workers=a.workers) as pool:
 summary=list(pool.map(execute,sorted(a.transcripts.glob('*.fa'))))
a.out.mkdir(parents=True,exist_ok=True)
(a.out/'generator_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'targets':len(summary),'status':'label-blind generator outputs frozen'}))
