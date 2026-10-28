"""Create a checksum lock before final labels are opened, or verify an existing lock."""
import argparse,json,sys
from datetime import datetime,timezone
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.fasta import sha256_file
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--lock',type=Path,required=True);p.add_argument('--verify',action='store_true')
p.add_argument('files',nargs='*',type=Path)
a=p.parse_args()
if a.verify:
 d=json.loads(a.lock.read_text())
 for name,digest in d['artifacts'].items():
  if sha256_file(Path(name))!=digest:raise ValueError('freeze mismatch: '+name)
 print('freeze verified')
else:
 if a.lock.exists():raise ValueError('refusing to overwrite an existing freeze')
 if not a.files:raise ValueError('provide code, config, splits, model manifests and label-blind predictions')
 a.lock.write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
                             'artifacts':{str(f.resolve()):sha256_file(f) for f in a.files}},indent=2)+'\n')
