"""Standalone scientific plots for benchmark tables; development status explicit."""
import argparse,json
from pathlib import Path
from collections import defaultdict
from statistics import mean
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--result',type=Path,action='append',required=True)
p.add_argument('--comparison',type=Path);p.add_argument('--selection',type=Path)
p.add_argument('--measurements',type=Path)
p.add_argument('--study');p.add_argument('--endpoint');p.add_argument('--format')
p.add_argument('--out',type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
fig,axes=plt.subplots(1,3,figsize=(15,4),constrained_layout=True)
if a.selection:
 d=json.loads(a.selection.read_text());by={r['key']:r for r in d['candidates']}
 activity={}
 if a.measurements:
  import csv
  with a.measurements.open() as handle:
   for row in csv.DictReader(handle,delimiter='\t'):
    if row['transcript_id']==d['provenance']['transcript_id'] and row['target_sha256']==d['provenance']['target_sha256']:
     activity[row['guide_dna'].replace('U','T')]=float(row['value'])
 measured=[r for r in d['candidates'] if r['guide_dna'] in activity]
 if measured:
  scatter=axes[0].scatter([r['start'] for r in measured],[0]*len(measured),c=[activity[r['guide_dna']] for r in measured],cmap='viridis',s=18)
  fig.colorbar(scatter,ax=axes[0],label='Raw measured activity')
 for i,(g,keys) in enumerate(d['research_top'].items(),1):
  keys=keys[:5]
  axes[0].scatter([by[k]['start'] for k in keys],[i]*len(keys),marker='|',s=80,label='Selected Top-5')
 axes[0].set(xlabel='Target position (1-based)',ylabel='Measured / selected',title='A  '+d['provenance']['transcript_id'].split('|')[0]+' measured guides')
 axes[0].legend(fontsize=5)
else:
 axes[0].text(.5,.5,'Supply --selection for a candidate map',ha='center',transform=axes[0].transAxes)
 axes[0].set_title('A  Candidate map')
for path in a.result:
 d=json.loads(path.read_text());values=defaultdict(lambda:defaultdict(list))
 for r in d['results']:
  if a.study and r['study']!=a.study:continue
  if a.endpoint and r['endpoint']!=a.endpoint:continue
  if a.format and r['duplex_format']!=a.format:continue
  if r['complete_measured_panel']:
   stratum=(r['study'],r['endpoint'],r['duplex_format'],r['group'])
   values[(stratum,r['requested_k'])][r['gene']].append(r['precision_at_k'])
 for stratum in sorted({key[0] for key in values}):
  xs=sorted(k for s,k in values if s==stratum)
  ys=[mean(mean(assays) for assays in values[(stratum,k)].values()) for k in xs]
  axes[1].plot(xs,ys,marker='o',label={'model':'Frozen model','quality':'Model + native/quality','native':'Native union + fallback'}.get(path.stem,path.stem))
axes[1].set(xlabel='K',ylabel='Measured precision',ylim=(0,1),title='B  Complete measured panels');axes[1].legend(fontsize=7)
if a.comparison:
 d=json.loads(a.comparison.read_text())
 for i,c in enumerate(d['comparisons']):
  values=c['per_gene_difference'];ys=[i+(j-(len(values)-1)/2)*.008 for j in range(len(values))];axes[2].scatter(list(values.values()),ys,alpha=.6,s=22)
  if c['mean_difference'] is not None:
   lo,hi=c['ci95'];m=c['mean_difference'];axes[2].errorbar(m,i,xerr=[[m-lo],[hi-m]],fmt='ks')
 axes[2].axvline(0,color='grey',lw=1)
 axes[2].set_yticks(range(len(d['comparisons'])),['Top-'+str(c['stratum'][-1]) for c in d['comparisons']])
else:axes[2].text(.5,.5,'Paired comparison unavailable',ha='center',transform=axes[2].transAxes)
axes[2].set(xlabel='Paired per-gene precision difference',title='C  Paired gene differences')
fig.suptitle('Retrospective development benchmark — independent validation not established')
for ext in ('png','pdf','svg'):fig.savefig(a.out/('benchmark.'+ext),dpi=180)
