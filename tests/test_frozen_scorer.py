import hashlib,json,sys
from pathlib import Path
import pytest
from sirna_consensus.scoring import score_candidates
from sirna_consensus.types import Core


def test_frozen_scorer_checks_identity_formats_artifacts_and_replays(tmp_path):
    worker=tmp_path/'worker.py'
    worker.write_text('import sys,json,hashlib\nfrom pathlib import Path\np=Path(sys.argv[1]);d=json.loads(p.read_text())\nPath(sys.argv[2]).write_text(json.dumps({"input_sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"scores":{r["key"]:.8 for r in d["candidates"]}}))\n')
    checkpoint=tmp_path/'weights';checkpoint.write_text('fixture')
    digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    a=Core('TX',1,19,'AGCTTACGATCGTACCTGA');b=Core('TX',2,21,'GCTTACGATCGTACCTGATAC')
    manifest=tmp_path/'manifest.json'
    manifest.write_text(json.dumps({'name':'OligoFormer','input_convention':'paired_antisense_19',
     'supported_formats':[a.duplex_format],'artifacts':{str(worker):digest(worker),str(checkpoint):digest(checkpoint)},
     'checkpoint_sha256':digest(checkpoint),'command':[sys.executable,str(worker),'{input}','{output}']}))
    rows=[{'core':c} for c in (a,b)]
    scores,meta=score_candidates(rows,'AGCTTACGATCGTACCTGATACG',manifest,tmp_path/'work')
    assert scores=={a.key:.8} and meta['unsupported_count']==1
    again,meta=score_candidates(rows,'AGCTTACGATCGTACCTGATACG',manifest,tmp_path/'work')
    assert scores==again and meta['replayed']
    checkpoint.write_text('changed')
    with pytest.raises(ValueError,match='checksum mismatch'):
        score_candidates(rows,'AGCTTACGATCGTACCTGATACG',manifest,tmp_path/'work')
