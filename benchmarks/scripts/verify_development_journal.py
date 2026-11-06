"""Verify journal provenance, historical archived bytes and frozen AUC/AP."""
import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

from bs4 import BeautifulSoup

ROOT=Path(__file__).resolve().parents[2]
B=ROOT/'benchmarks';J=B/'development_journal'
sys.path.insert(0,str(ROOT/'src'))


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    source=json.loads((J/'sources.json').read_text())
    for p,digest in source['metrics_sources'].items():assert sha(ROOT/p)==digest,p
    assert sha(B/'three_benchmark_report_zh.html')==source['report_sha256']
    assert sha(B/'templates/development_journal_zh.html')==source['template_sha256']
    assert sha(B/'scripts/build_development_journal_zh.py')==source['generator_sha256']
    preserve=json.loads((J/'historical_report_preservation.json').read_text())
    assert sha(B/preserve['archived_report'])==preserve['archived_sha256']
    spec=importlib.util.spec_from_file_location('freeze',B/'scripts/verify_freezes.py');freeze=importlib.util.module_from_spec(spec);spec.loader.exec_module(freeze)
    locks={}
    for name,mapping in preserve['locks'].items():
        present,bad,missing=freeze.check(B/name,mapping)
        assert not bad,(name,bad)
        assert missing==0,(name,missing)
        locks[name]=present
    soup=BeautifulSoup((B/'three_benchmark_report_zh.html').read_text(),'html.parser')
    ids={t['id'] for t in soup.find_all(id=True)};broken=[]
    for tag,attr in [('a','href'),('img','src')]:
        for t in soup.find_all(tag):
            link=t.get(attr,'')
            if not link or link.startswith(('https:','http:')):continue
            if link.startswith('#'):
                if link[1:] not in ids:broken.append(link)
            elif not (B/link.split('#')[0]).exists():broken.append(link)
    assert not broken,broken
    assert '@@' not in soup.get_text()
    embedded=json.loads(soup.find('script',id='journal-data').string)
    assert len(embedded['cohorts'])==15
    assert len(embedded['genes'])==59
    # Scikit-learn arithmetic independent of the custom plotting implementation.
    code='''import json,sys
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
from pathlib import Path
b=Path(sys.argv[1]);r=json.loads((b/'external_confirmation_20261106/records.json').read_text());p=json.loads((b/'external_confirmation_20261106/outputs/predictions.json').read_text());c=json.loads((b/'development_journal/curve_metrics.json').read_text());y=[x['value']>=x['active_threshold'] for x in r]
for name,m in c['metrics'].items():
 s=p['scores'][name]
 assert abs(roc_auc_score(y,s)-m['roc_auc'])<1e-12
 assert abs(average_precision_score(y,s)-m['average_precision'])<1e-12
print('independent sklearn AUC/AP: all 3 methods pass')
'''
    checked=subprocess.run([str(ROOT/'.scorer-venv/bin/python'),'-c',code,str(B)],check=True,capture_output=True,text=True)
    output=dict(passed=True,source_files=len(source['metrics_sources']),historical_locks=locks,
                broken_report_links=broken,cohorts=len(embedded['cohorts']),gene_assays=len(embedded['genes']),
                independent_AUC_AP=checked.stdout.strip(),report_sha256=source['report_sha256'])
    (J/'verification.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
