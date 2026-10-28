"""Within-gene/assay label permutations for one prespecified benchmark stratum."""
import argparse,hashlib,json,random,sys
from collections import defaultdict
from pathlib import Path
from statistics import mean
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.benchmark import load_measurements,evaluate_record
from sirna_consensus.fasta import sha256_file
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--benchmark',required=True,type=Path);p.add_argument('--development',required=True,type=Path)
p.add_argument('--strategy',default='quality');p.add_argument('--permutations',type=int,default=1000)
p.add_argument('--out',required=True,type=Path);p.add_argument('--seed',type=int,default=2026)
a=p.parse_args();all_measured=load_measurements(a.benchmark)
strata={(r['study'],r['endpoint'],r['duplex_format']) for r in all_measured}
if len(strata)!=1:raise ValueError('run permutation controls separately by study/endpoint/format')
observed=[];null_sums=[0.]*a.permutations;genes=0
for path in sorted(a.development.glob('*/'+a.strategy+'/selection.json')):
 record=json.loads(path.read_text());tid=record['provenance']['transcript_id']
 measured=[r for r in all_measured if r['transcript_id']==tid and r['target_sha256']==record['provenance']['target_sha256']]
 real=evaluate_record(record,measured,5,random_repeats=0)
 values=[r['precision_at_k'] for r in real if r['complete_measured_panel']]
 if not values:continue
 observed.append(mean(values));genes+=1
 units=defaultdict(list)
 for i,r in enumerate(measured):units[(r['gene'],r['assay'])].append(i)
 rng=random.Random(a.seed+int(hashlib.sha256(tid.encode()).hexdigest()[:12],16))
 for j in range(a.permutations):
  shuffled=[dict(r) for r in measured]
  for indexes in units.values():
   labels=[measured[i]['effective'] for i in indexes];rng.shuffle(labels)
   for i,label in zip(indexes,labels):shuffled[i]['effective']=label
  metric=evaluate_record(record,shuffled,5,random_repeats=0)
  null_sums[j]+=mean(r['precision_at_k'] for r in metric if r['complete_measured_panel'])
null=[v/genes for v in null_sums] if genes else []
real=mean(observed) if observed else None
result={'stratum':next(iter(strata)),'genes':genes,'permutations':a.permutations,'seed':a.seed,
        'observed_gene_precision_at_5':real,'null_mean':mean(null) if null else None,
        'p_one_sided':(1+sum(v>=real for v in null))/(1+len(null)) if null else None,
        'benchmark_sha256':sha256_file(a.benchmark),'independent_validation_passed':False}
a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
