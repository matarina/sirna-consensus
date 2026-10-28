"""Create compact, source-linked development results without a validation claim."""
import argparse,json,sys
from pathlib import Path
from collections import defaultdict,Counter
from statistics import mean
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.fasta import sha256_file
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--development',type=Path,required=True);p.add_argument('--comparison',type=Path,required=True)
p.add_argument('--ablations',type=Path,required=True);p.add_argument('--screen',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
summary=json.loads((a.development/'development_summary.json').read_text())
comparison=json.loads(a.comparison.read_text());screen=json.loads(a.screen.read_text())
def summarize(rows):
 rows=[r for r in rows if r['requested_k']==5]
 genes=defaultdict(list)
 for r in rows:genes[r['gene']].append(r)
 complete={g:[r['precision_at_k'] for r in rs if r['complete_measured_panel']] for g,rs in genes.items()}
 values=[mean(v) for v in complete.values() if v]
 return {'complete_genes':sum(bool(v) for v in complete.values()),'eligible_genes':len(genes),
         'precision_at_5_complete':mean(values) if values else None,
         'mean_fill_rate':mean(mean(r['fill_rate'] for r in rs) for rs in genes.values()) if genes else 0,
         'mean_candidate_coverage':mean(mean(r['candidate_coverage'] for r in rs) for rs in genes.values()) if genes else 0}
methods={name:summarize(json.loads((a.development/(name+'.json')).read_text())['results']) for name in summary['strategies']}
ablations={p.stem:summarize(json.loads(p.read_text())['results']) for p in a.ablations.glob('*.json')}
result={'status':'development_only_not_independently_validated','records':1941,'genes':25,
        'methods':methods,'ablations':ablations,'paired_comparison':comparison,
        'screen_smoke':{'candidates':len(screen['candidates']),
                        'states':dict(Counter(r.get('screen',{}).get('status','unverified') for r in screen['candidates'])),
                        'screened_panels':{g:len(keys) for g,keys in screen['screened_top'].items()},
                        'reference_sha256':screen['provenance']['screen'].get('reference_sha256')},
        'source_hashes':{str(p):sha256_file(p) for p in [a.comparison,a.screen,a.development/'development_summary.json']},
        'independent_validation_passed':False}
(a.out/'development_results.json').write_text(json.dumps(result,indent=2)+'\n')
lines=['# Development benchmark results','',
 'Implementation and execution checks passed; the agreed efficacy-improvement criterion did not.',
 '', 'This is exposed Huesken development data: 1,941 full guides across 25 human cloned constructs. '
 'It is not independent validation or endogenous knockdown evidence. Model selection and these '
 'performance estimates use the same development cohort.', '',
 '| Method | Mean Top-5 precision, complete genes | Complete genes / 25 | Mean fill rate |',
 '|---|---:|---:|---:|']
for name,r in methods.items():
 precision=f"{r['precision_at_5_complete']:.1%}" if r['precision_at_5_complete'] is not None else 'Unavailable'
 lines.append(f"| {name} | {precision} | {r['complete_genes']} | {r['mean_fill_rate']:.1%} |")
primary=next(c for c in comparison['comparisons'] if c['stratum'][-1]==5)
lo,hi=primary['ci95']
lines += ['',f"Development-selected quality fusion versus frozen model: **{primary['mean_difference']*100:+.1f} percentage points**, "
 f"paired 95% gene-bootstrap CI **[{lo*100:+.1f}, {hi*100:+.1f}] points**. "
 'The interval includes zero and the point estimate is below the +10-point target.', '',
 '`native` is the union with sequence-quality fallback; it is not an individual tool. '
 'siRNAdesigner precision uses only its 15 complete genes and should not be compared directly '
 'with the 25-gene means. siRNAforge and exact consensus have no complete measured panels '
 'in this duplex-format cohort; this is missing coverage, not evidence of inefficacy.', '',
 'Eight component ablations and 1,000 within-gene/assay label permutations were run. '
 'They are development diagnostics; pretrained exposure and development selection prevent '
 'interpreting their statistics as independent generalization evidence.', '',
 'The real AREG sequence screen completed against the pinned Ensembl 116 cDNA reference. '
 'Its separately screened panels satisfy the computational mismatch criterion only.', '',
 'See [figures](../figures/benchmark.pdf), [machine-readable results](development_results.json), '
 '[source audit](../dataset_manifest.json), and [reproduction instructions](../README.md).']
(a.out/'REVIEW.md').write_text('\n'.join(lines)+'\n')
print(a.out/'REVIEW.md')
