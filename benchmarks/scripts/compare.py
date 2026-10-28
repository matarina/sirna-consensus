"""Paired study/endpoint/format comparisons; never pool incompatible panels."""
import argparse,json,sys
from pathlib import Path
from collections import defaultdict
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.benchmark import paired_summary
from sirna_consensus.fasta import sha256_file
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--candidate',required=True);p.add_argument('--baseline',required=True)
p.add_argument('--out',required=True);p.add_argument('--seed',type=int,default=2026)
p.add_argument('--replicates',type=int,default=10000)
a=p.parse_args()
candidate=json.loads(Path(a.candidate).read_text());baseline=json.loads(Path(a.baseline).read_text())
if not candidate.get('benchmark_sha256') or candidate.get('benchmark_sha256') != baseline.get('benchmark_sha256'):
 raise ValueError('comparison requires the same pinned measurement snapshot')
for key in ('threshold','panel','track'):
 if key not in candidate.get('config',{}) or candidate['config'][key] != baseline.get('config',{}).get(key):
  raise ValueError('incompatible evaluation protocol: '+key)
left=candidate['results'];right=baseline['results']
def group(rows):
 d=defaultdict(list)
 for r in rows:d[tuple(r[k] for k in ('study','endpoint','duplex_format','track','requested_k'))].append(r)
 return d
l,r=group(left),group(right)
if set(l)!=set(r):raise ValueError('candidate and baseline have different evaluation strata')
output={'candidate_sha256':sha256_file(Path(a.candidate)),'baseline_sha256':sha256_file(Path(a.baseline)),
        'comparisons':[{'stratum':key,**paired_summary(l[key],r[key],a.replicates,a.seed)} for key in sorted(l)],
        'independent_validation_passed':False}
Path(a.out).write_text(json.dumps(output,indent=2)+'\n')
