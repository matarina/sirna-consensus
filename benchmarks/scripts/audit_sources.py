"""Audit the publisher-linked snapshot without using it as an independent test."""
import collections,hashlib,json,xml.etree.ElementTree as ET,zipfile
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
raw=ROOT/'data/raw';ns={'x':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
meta=json.loads((raw/'figshare_25908841_v1.json').read_text())
files={'46578241':'target_gene_sequences.fa','46578250':'sirna_data.xlsx'}
artifacts={}
for f in meta['files']:
 p=raw/files[str(f['id'])];data=p.read_bytes()
 if hashlib.md5(data).hexdigest()!=f['computed_md5']:raise ValueError('publisher MD5 mismatch')
 artifacts[p.name]={'sha256':hashlib.sha256(data).hexdigest(),'publisher_md5':f['computed_md5'],
                    'source_url':f['download_url'],'bytes':len(data)}
with zipfile.ZipFile(raw/'sirna_data.xlsx') as z:
 strings=[''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('x:si',ns)]
 rows=[]
 for r in ET.fromstring(z.read('xl/worksheets/sheet1.xml')).findall('x:sheetData/x:row',ns):
  values={}
  for c in r:
   v=c.findtext('x:v','',ns)
   values[''.join(x for x in c.attrib['r'] if x.isalpha())]=strings[int(v)] if c.attrib.get('t')=='s' else v
  rows.append(values)
 headers=rows.pop(0)
 records=[{headers[k]:v for k,v in r.items()} for r in rows]
studies=collections.Counter(r['Authors'] for r in records)
# All seven source studies occur in OligoFormer's published development/evaluation datasets.
report={'created_utc':datetime.now(timezone.utc).isoformat(),'figshare_article':25908841,'version':meta['version'],
        'snapshot_records':len(records),'studies':dict(studies),'artifacts':artifacts,
        'figshare_license':meta['license'],'paper_license':'CC BY-NC 2.0 UK; provenance discrepancy retained',
        'redistribution_policy':'raw snapshot local only pending license reconciliation',
        'independent_holdout_eligible':False,
        'reasons':['Seven source studies occur in OligoFormer development/evaluation collections; checkpoint-specific exposure unverified',
                   'Huesken has previous project training and evaluation exposure',
                   'Database Huesken metadata (HeLa/luciferase) conflicts with original-source H1299/fluorescent-reporter curation',
                   'Database percent-inhibition field rescales Huesken normalized activity; do not interpret it as percent knockdown'],
        'checkpoint_training_audit':'unverified; no independent claim permitted',
        'external_holdout_records':0,'independent_validation_passed':False}
(ROOT/'dataset_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'snapshot_records':len(records),'studies':dict(studies),'independent_holdout_eligible':False}))
