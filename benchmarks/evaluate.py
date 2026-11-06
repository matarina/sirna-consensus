"""Evaluate schema-v2 selections using exact guides and explicit assay endpoints."""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from sirna_consensus.benchmark import evaluate_record, load_measurements, permutation_control
from sirna_consensus.fasta import sha256_file


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if '--phase' in arguments and arguments[arguments.index('--phase') + 1:][:1] == ['confirmation']:
        from sirna_consensus.strict_evaluation import main as confirmation_main
        return confirmation_main(arguments)
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--benchmark',required=True)
    p.add_argument('--run',required=True,action='append',help='selection.json; repeat for multiple targets/methods')
    p.add_argument('--track',choices=['measured','deployment'],default='measured')
    p.add_argument('--panel',choices=['research_top','screened_top','diverse_top'],default='research_top')
    p.add_argument('--top-k',default='5,10')
    p.add_argument('--threshold',type=float,choices=[.6,.7,.8],default=.7)
    p.add_argument('--random-repeats',type=int,default=1000)
    p.add_argument('--permutations',type=int,default=0)
    p.add_argument('--seed',type=int,default=2026)
    p.add_argument('--phase',choices=['development','final'],default='development')
    p.add_argument('--freeze-lock')
    p.add_argument('--cohort-manifest')
    p.add_argument('--protocol')
    p.add_argument('--output',required=True)
    args=p.parse_args(argv)
    if args.random_repeats<0 or args.permutations<0:
        p.error('repetitions cannot be negative')
    if args.phase == 'final':
        if not all((args.freeze_lock,args.cohort_manifest,args.protocol)):
            p.error('final evaluation requires --freeze-lock, --cohort-manifest and --protocol')
        import os
        lock_path=Path(args.freeze_lock)
        lock=json.loads(lock_path.read_text())['artifacts']
        def key(f): return os.path.relpath(Path(f).resolve(), lock_path.parent.resolve())
        def resolve(k): return Path(k) if Path(k).is_absolute() else lock_path.parent/k
        required=[args.benchmark, args.cohort_manifest, args.protocol, *args.run]
        if any(key(f) not in lock for f in required):
            p.error('final inputs are missing from the protocol freeze')
        for f in required:
            if not Path(f).is_file():
                p.error('final input not present: '+str(f))
        for k,digest in lock.items():
            target=resolve(k)
            # Absent pinned entries (e.g. gitignored run.json on a bare clone) are tolerated,
            # like verify_freezes; a present file must match. Dirs hash as absent, not crash.
            if target.is_file() and sha256_file(target)!=digest:
                p.error('freeze mismatch: '+k)
        cohort=json.loads(Path(args.cohort_manifest).read_text())
        if cohort.get('independent_holdout_eligible') is not True or cohort.get('checkpoint_training_audit')!='verified':
            p.error('independent source/checkpoint exposure audit has not passed')
        holdout=set(cohort.get('holdout_studies') or [])
        if not holdout:
            p.error('cohort manifest lists no holdout studies')
        protocol=json.loads(Path(args.protocol).read_text())
        if args.threshold != .7 or args.track != 'measured' or args.panel != 'research_top':
            p.error('final primary analysis requires the frozen 0.7 measured research protocol')
        if args.seed != protocol['seed']:
            p.error('seed differs from frozen protocol')
    measured=load_measurements(Path(args.benchmark),args.threshold)
    if args.phase == 'final':
        # The measured benchmark must be the audited holdout, not some other table.
        outside=sorted({r['study'] for r in measured}-holdout)
        if outside:
            p.error('benchmark contains studies outside the cohort holdout: '+', '.join(outside))
    output={'schema_version':2,'benchmark_sha256':sha256_file(Path(args.benchmark)),
            'config':vars(args),'results':[],'controls':[],'independent_validation_passed':False}
    for path in args.run:
        record=json.loads(Path(path).read_text())
        for k in map(int,args.top_k.split(',')):
            rows=evaluate_record(record,measured,k,args.track,args.panel,args.random_repeats,args.seed)
            output['results'].extend({**r,'run':path,'run_sha256':sha256_file(Path(path))} for r in rows)
            if args.permutations:
                if args.track!='measured' or args.panel!='research_top':
                    p.error('permutation control currently applies to measured research ranking')
                output['controls'].append({'run':path,'k':k,**permutation_control(record,measured,k,args.permutations,args.seed)})
    if args.phase=='final' and not output['results']:
        p.error('final evaluation produced no results (empty or unmatchable benchmark/run)')
    out=Path(args.output)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'assay_results':len(output['results']),'output':str(out)}))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
