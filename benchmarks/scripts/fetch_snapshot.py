"""Download the publisher-linked Figshare v1 snapshot and verify publisher MD5s."""
import hashlib,json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parents[1]/'data/raw';root.mkdir(parents=True,exist_ok=True)
meta=root/'figshare_25908841_v1.json'
def fetch(url,path):
 subprocess.run(['curl','--fail','--location','--retry','2','--max-time','180','--output',str(path),url],check=True)
fetch('https://api.figshare.com/v2/articles/25908841/versions/1',meta)
data=json.loads(meta.read_text());names={46578241:'target_gene_sequences.fa',46578250:'sirna_data.xlsx'}
for f in data['files']:
 if f['id'] not in names:continue
 path=root/names[f['id']]
 if not path.is_file() or hashlib.md5(path.read_bytes()).hexdigest()!=f['computed_md5']:
  fetch(f['download_url'],path)
 if hashlib.md5(path.read_bytes()).hexdigest()!=f['computed_md5']:raise ValueError('snapshot hash mismatch')
print(root)
