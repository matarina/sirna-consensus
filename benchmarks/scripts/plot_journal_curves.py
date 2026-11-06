"""Post-hoc ROC/PR from frozen Katoh predictions. No model fitting."""
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
B=ROOT/'benchmarks'
J=B/'development_journal'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    source=B/'external_confirmation_20261106'
    records=json.loads((source/'records.json').read_text())
    saved=json.loads((source/'outputs/predictions.json').read_text())
    assert saved['record_ids']==[r['record_id'] for r in records]
    y=np.array([r['value']>=r['active_threshold'] for r in records],dtype=int)
    assert len(set(y))==2
    positive=int(y.sum());negative=len(y)-positive
    names={'current_default':'Current default','huesken_only_challenger':'Huesken-only candidate','quality_rules':'Existing quality rules'}
    fig,axes=plt.subplots(1,2,figsize=(12,5.2),layout='constrained')
    plt.rcParams.update({'svg.fonttype':'none','font.size':10})
    result={}
    for (key,name),color in zip(names.items(),['#197d83','#b68535','#8395a7']):
        scores=np.asarray(saved['scores'][key]);order=np.argsort(-scores,kind='stable');s=scores[order];labels=y[order]
        ends=np.r_[np.flatnonzero(np.diff(s)!=0),len(s)-1]
        tp=np.r_[0,np.cumsum(labels)[ends]]; fp=np.r_[0,np.cumsum(1-labels)[ends]]
        recall=tp/positive;fpr=fp/negative
        precision=np.r_[1,tp[1:]/(tp[1:]+fp[1:])]
        auc=float(np.trapezoid(recall,fpr));ap=float(np.sum(np.diff(recall)*precision[1:]))
        # Independent AUC as positive-negative pairwise concordance, including ties.
        a=scores[y==1];b=scores[y==0]
        independent_auc=float(np.mean((a[:,None]>b[None,:])+.5*(a[:,None]==b[None,:])))
        assert abs(auc-independent_auc)<1e-12
        result[key]=dict(roc_auc=auc,average_precision=ap,independent_pairwise_auc=independent_auc,
                        roc=dict(fpr=fpr.tolist(),tpr=recall.tolist()),
                        pr=dict(recall=recall.tolist(),precision=precision.tolist()))
        axes[0].plot(fpr,recall,color=color,lw=2,label=f'{name} (AUC {auc:.3f})')
        axes[1].step(recall,precision,where='pre',color=color,lw=2,label=f'{name} (AP {ap:.3f})')
    axes[0].plot([0,1],[0,1],ls='--',color='#aab5bb',label='Random (AUC 0.500)')
    axes[1].axhline(positive/len(y),ls='--',color='#aab5bb',label=f'Prevalence ({positive/len(y):.3f})')
    for ax,title,xlabel,ylabel in zip(axes,['ROC curve','Precision–Recall curve'],['False positive rate','Recall'],['True positive rate','Precision']):
        ax.set(xlim=(0,1),ylim=(0,1.02),title=title,xlabel=xlabel,ylabel=ylabel)
        ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.16);ax.legend(loc='lower right' if ax==axes[0] else 'lower left',fontsize=8)
    fig.suptitle(f'Katoh 2007 | measured paired-core component | n={len(y)} ({positive} active / {negative} inactive)\nActive: normalized assay activity >=0.7; post-hoc diagnostic, no gene-level CI',fontsize=11)
    figures=J/'figures';figures.mkdir(parents=True,exist_ok=True)
    fig.savefig(figures/'katoh_roc_pr.svg',metadata={'Date':None});fig.savefig(figures/'katoh_roc_pr.png',dpi=200)
    plt.close(fig)
    report=dict(scope='post-hoc diagnostic on frozen Katoh measured paired-core scores; not preregistered primary or promotion criterion',
                active_definition='original normalized activity >=0.7; not absolute 70% knockdown',
                records=len(y),positives=positive,negatives=negative,prevalence=positive/len(y),
                ci=None,ci_reason='One gene/study; tiling rows dependent; no row bootstrap CI.',
                metrics=result,sources={str(p.relative_to(ROOT)):sha(p) for p in [source/'records.json',source/'outputs/predictions.json',source/'test_protocol.json']})
    (J/'curve_metrics.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:{m:v[m] for m in ['roc_auc','average_precision']} for k,v in result.items()},indent=2))


if __name__=='__main__':main()
