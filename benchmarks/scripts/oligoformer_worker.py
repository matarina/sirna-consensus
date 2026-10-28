"""Frozen OligoFormer inference with explicit context, eval mode, and local RNA-FM.

Run in the scorer environment. Uses audited upstream model and thermodynamic
functions; no fitting, downloads, arbitrary window selection or guide truncation.
"""
import argparse
import ast
import hashlib
import importlib.util
import itertools
import json
import os
import re
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True,type=Path);p.add_argument('--output',required=True,type=Path)
    p.add_argument('--model-source',required=True,type=Path);p.add_argument('--features-source',required=True,type=Path)
    p.add_argument('--checkpoint',required=True,type=Path);p.add_argument('--rna-fm-checkpoint',required=True,type=Path)
    p.add_argument('--device',default='cpu');p.add_argument('--threads',type=int,default=4)
    a=p.parse_args()
    # Upstream's positional embeddings bind to its module-global device.
    if a.device=='cpu':os.environ['CUDA_VISIBLE_DEVICES']=''
    import numpy as np
    import pandas as pd
    import torch
    import fm
    torch.set_num_threads(a.threads);torch.manual_seed(2026);np.random.seed(2026)
    if not a.rna_fm_checkpoint.is_file():raise ValueError('local RNA-FM checkpoint required; downloads disabled')
    spec=importlib.util.spec_from_file_location('frozen_oligo_model',a.model_source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.device=torch.device(a.device)
    model=module.Oligo().to(a.device)
    model.load_state_dict(torch.load(a.checkpoint,map_location=a.device,weights_only=True));model.eval()
    # RNA-FM's legacy checkpoint stores argparse.Namespace; read only pinned local artifacts.
    with torch.serialization.safe_globals([argparse.Namespace]):
        language,alphabet=fm.pretrained.rna_fm_t12(str(a.rna_fm_checkpoint))
    language=language.to(a.device).eval();convert=alphabet.get_batch_converter()
    allowed={'antiRNA','Calculate_DGH','Calculate_end_diff','calculate_td'}
    tree=ast.parse(a.features_source.read_text())
    nodes=[n for n in tree.body if (isinstance(n,ast.FunctionDef) and n.name in allowed) or
           (isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id in {'DeltaG','DeltaH'} for t in n.targets))]
    namespace={'np':np,'pd':pd,'itertools':itertools,'re':re}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(a.features_source),'exec'),namespace)
    request=json.loads(a.input.read_text());transcript=request['transcript'].replace('T','U')
    def embedding(sequence):
        _,_,tokens=convert([('rna',sequence)])
        with torch.inference_mode():
            return language(tokens.to(a.device),repr_layers=[12])['representations'][12][:,1:len(sequence)+1]
    def onehot(sequence):
        # sklearn OneHotEncoder used upstream sorts categories lexicographically.
        categories='ACGUX';x=torch.zeros((1,1,len(sequence),5),device=a.device)
        for i,c in enumerate(sequence):x[0,0,i,categories.index(c)]=1
        return x
    scores={}
    # Batch size one is deliberate: upstream LSTM uses sequence-first layout.
    # Changing batch size changes predictions for this published architecture.
    with torch.inference_mode():
        for row in request['candidates']:
            core=row['core'];guide=row['paired_guide'].replace('T','U');start=core['start']-1
            if len(guide)!=19 or core['length']!=19:raise ValueError('unsupported molecule')
            target=core['sense'].replace('T','U')
            if transcript[start:start+19]!=target:raise ValueError('target context mismatch')
            context='X'*max(0,19-start)+transcript[max(0,start-19):start+38]+'X'*max(0,start+38-len(transcript))
            if len(context)!=57:raise ValueError('invalid target context length')
            frame=namespace['calculate_td'](pd.DataFrame({'siRNA':[guide],'mRNA':[context]}))
            td=torch.tensor([frame.iloc[0]['td']],dtype=torch.float32,device=a.device)
            prediction,_,_=model(onehot(guide),onehot(context),embedding(guide),embedding(context),td)
            scores[row['key']]=float(prediction[0,1])*1.341
    a.output.write_text(json.dumps({'input_sha256':hashlib.sha256(a.input.read_bytes()).hexdigest(),
                                   'scores':scores,'eval_mode':True,'batch_size':1},indent=2)+'\n')

if __name__=='__main__':main()
