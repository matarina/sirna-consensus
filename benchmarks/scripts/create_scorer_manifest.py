"""Pin the local frozen scorer and dependency source files; never assert training independence."""
import argparse,hashlib,json,importlib.metadata,sys
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--upstream',type=Path,required=True);p.add_argument('--rna-fm-checkpoint',type=Path,required=True)
p.add_argument('--out',type=Path,required=True);p.add_argument('--device',default='cpu')
a=p.parse_args();root=Path(__file__).resolve().parents[2];worker=root/'benchmarks/scripts/oligoformer_worker.py'
model=a.upstream/'scripts/model.py';features=a.upstream/'scripts/infer.py';checkpoint=a.upstream/'model/best_model.pth'
import fm
artifacts=[worker,model,features,checkpoint,a.rna_fm_checkpoint]+list(Path(fm.__file__).parent.rglob('*.py'))
digest=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
d={'name':'OligoFormer','input_convention':'paired_antisense_19',
   'supported_formats':['binding19-guide19-overhang0-passenger_unknown',
                        'binding19-guide21-overhang2-passenger_unknown',
                        'binding19-guide21-overhang2-passenger21-overhang2'],
   'checkpoint_sha256':digest(checkpoint),'artifacts':{str(p.resolve()):digest(p) for p in artifacts},
   'command':[sys.executable,str(worker),'--input','{input}','--output','{output}',
              '--model-source',str(model.resolve()),'--features-source',str(features.resolve()),
              '--checkpoint',str(checkpoint.resolve()),'--rna-fm-checkpoint',str(a.rna_fm_checkpoint.resolve()),
              '--device',a.device,'--threads','4'],
   'packages':{name:importlib.metadata.version(name) for name in ['torch','numpy','pandas','rna-fm']},
   'inference_policy':'eval mode, fixed seed 2026, batch size one, no fitting, local checkpoint only',
   'training_audit':'unverified checkpoint split; all siRNAEfficacyDB studies exposed in published model development/evaluation',
   'timeout_s':3600,'independent_validation_passed':False}
a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(d,indent=2)+'\n')
print(a.out)
