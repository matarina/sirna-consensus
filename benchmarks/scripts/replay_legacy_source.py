"""Materialise the exact hash-checked historical package before a legacy CLI run."""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))
from sirna_consensus.fasta import sha256_file


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--lock',choices=['development_freeze.json','final_freeze.json'],default='development_freeze.json')
    p.add_argument('arguments',nargs=argparse.REMAINDER)
    args=p.parse_args();bench=Path(__file__).resolve().parents[1]
    lock=json.loads((bench/args.lock).read_text())['artifacts']
    sources=json.loads((bench/'historical_code_snapshots.json').read_text())['locks'][args.lock]
    with tempfile.TemporaryDirectory(prefix='sirna-legacy-') as tmp:
        for key,path in sources.items():
            if not key.startswith('../src/'):
                continue
            artifact=bench/path
            if sha256_file(artifact)!=lock[key]:
                raise ValueError('historical source hash mismatch: '+key)
            target=Path(tmp)/key[len('../src/'):]
            target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(artifact.read_bytes())
        env=dict(os.environ);env['PYTHONPATH']=tmp
        command=args.arguments[1:] if args.arguments[:1]==['--'] else args.arguments
        if not command:p.error('supply legacy CLI arguments after --')
        return subprocess.run([sys.executable,'-m','sirna_consensus.cli',*command],env=env).returncode


if __name__=='__main__':raise SystemExit(main())
