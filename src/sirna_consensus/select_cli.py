"""CLI for label-blind exact-candidate selection and frozen-output replay."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from .fasta import read_transcript, sha256, sha256_file
from .selection import (
    STRATEGIES,
    enumerate_guides,
    load_candidates,
    rank_candidates,
    supplied_candidates,
    write_selection,
)


def add_parser(sub):
    p=sub.add_parser('select',help='rank exact guides from a v2 run or explicit transcript candidates')
    p.add_argument('--fasta',required=True)
    sources=p.add_mutually_exclusive_group(required=True)
    sources.add_argument('--run',help='schema-v2 run.json; replays all exact normalized hits')
    sources.add_argument('--candidates',help='label-free TSV of explicit molecules')
    sources.add_argument('--enumerate',action='store_true',help='enumerate guide-only 19/20/21-nt windows')
    # Native evidence may be attached onto a supplied candidate universe (the dev path):
    # candidates stay the assayed molecules; the run only donates exact-molecule votes.
    p.add_argument('--run-evidence',help='schema-v2 run.json whose exact hits annotate --candidates')
    p.add_argument('--out',required=True)
    p.add_argument('--strategy',choices=STRATEGIES,default='native')
    p.add_argument('--filter-policy',choices=['hard','soft','off'],default='soft')
    p.add_argument('--disable-components',default='',help='comma-separated tools,support,quality ablations')
    p.add_argument('--top-k',type=int,default=5)
    p.add_argument('--scorer-manifest')
    p.add_argument('--policy',help='frozen development-selected format/strategy policy')
    p.add_argument('--screen-manifest')
    p.add_argument('--target-gene')
    p.add_argument('--bowtie',default='bowtie')
    return p


def command(args):
    tid,sequence=read_transcript(args.fasta)
    out=Path(args.out).resolve()
    # Deterministic work dir: the recorded scorer command must not embed a random uuid,
    # or every rerun would change selection.json (and the frozen lock) for no reason.
    # Also scope it by the inputs, so distinct selections sharing one --out cannot reuse
    # each other's work/scorer or work/screen files; identical inputs still share it.
    key=sha256('|'.join([tid,sha256(sequence),args.candidates or args.run or 'enumeration',
                         args.run_evidence or '',args.scorer_manifest or '',args.screen_manifest or '',
                         args.target_gene or '',args.bowtie or '']))
    work=out/'work'/key[:16]
    provenance={'transcript_id':tid,'target_sha256':sha256(sequence),'strategy':args.strategy,
                'filter_policy':args.filter_policy,'source':args.run or args.candidates or 'enumeration'}
    if args.run:
        rows,run=load_candidates(Path(args.run))
        if run['input_sha256']!=sha256(sequence) or run['transcript_id']!=tid:
            raise ValueError('replay input target mismatch')
        for row in rows:
            core = row['core']
            if core.transcript_id != tid or sequence[core.start-1:core.end] != core.sense:
                raise ValueError('replayed molecule does not match the input target')
        provenance['source_sha256']=sha256_file(Path(args.run))
    elif args.candidates:
        rows=supplied_candidates(Path(args.candidates),tid,sequence)
        provenance['source_sha256']=sha256_file(Path(args.candidates))
    else:
        rows=enumerate_guides(tid,sequence)
    if args.run_evidence:
        if not args.candidates:
            raise ValueError('--run-evidence annotates --candidates')
        if not Path(args.run_evidence).is_file():
            raise ValueError('--run-evidence not found: '+args.run_evidence)
        native,evidence=load_candidates(Path(args.run_evidence))
        if evidence['input_sha256']!=sha256(sequence) or evidence['transcript_id']!=tid:
            raise ValueError('native evidence target mismatch')
        # Generators emit 19/20/21-nt molecules; the assayed candidates are 19-nt cores.
        # Join on the shared 19-nt binding site, aggregating every tool that emitted it.
        # Convention: the site is sense[:19] (5' prefix), matching siRNAdesigner's own
        # passenger[:19] core rule. A 3' suffix would be an equally defensible choice and
        # changes which hits attach; the position (start) is deliberately not part of the
        # key, so a repeated 19-mer legitimately attaches to every occurrence.
        by_site={}
        skipped=0
        for hit in native:
            site=hit['core'].sense[:19]
            if len(site)!=19:
                skipped+=1
                continue
            slot=by_site.setdefault(site,{'tools':[],'tool_scores':{}})
            slot['tools']=sorted(set(slot['tools'])|set(hit['tools']))
            for tool,value in hit['tool_scores'].items():
                slot['tool_scores'][tool]=max(value,slot['tool_scores'].get(tool,value))
        attached=0
        for row in rows:
            site=row['core'].sense[:19]
            hit=by_site.get(site) if len(site)==19 else None
            if hit:
                row['tools']=sorted(set(row['tools'])|set(hit['tools']))
                row['tool_scores'].update(hit['tool_scores'])
                attached+=1
        provenance['evidence_sha256']=sha256_file(Path(args.run_evidence))
        provenance['evidence_attached']=attached
        provenance['evidence_skipped_short']=skipped
    scores={}
    if args.scorer_manifest:
        from .scoring import score_candidates
        scores,provenance['scorer']=score_candidates(rows,sequence,Path(args.scorer_manifest).resolve(),work/'scorer')
    elif args.strategy not in ('native','legacy'):
        raise ValueError('model strategies require --scorer-manifest')
    disabled=set(filter(None,args.disable_components.split(',')))
    if disabled-{'tools','support','quality'}:
        raise ValueError('unknown ablation component')
    components=tuple(c for c in ('tools','support','quality') if c not in disabled)
    provenance['components']=components
    if args.policy:
        policy=json.loads(Path(args.policy).read_text())
        if not args.scorer_manifest or policy.get('scorer_manifest_sha256')!=sha256_file(Path(args.scorer_manifest)):
            raise ValueError('policy requires its pinned scorer manifest')
        ranked=[]
        for fmt in sorted({r['core'].duplex_format for r in rows}):
            method=policy['strategies'].get(fmt,'model')
            ranked.extend(rank_candidates([r for r in rows if r['core'].duplex_format==fmt],method,args.filter_policy,scores,components=components))
        provenance['policy_sha256']=sha256_file(Path(args.policy))
        provenance['policy_status']=policy.get('status','development_only')
    else:
        ranked=rank_candidates(rows,args.strategy,args.filter_policy,scores,components=components)
    if args.screen_manifest:
        from .screening import screen
        try:
            provenance['screen']=screen(ranked,Path(args.screen_manifest),args.target_gene,args.bowtie,work/'screen')
        except (ValueError,OSError,RuntimeError,subprocess.SubprocessError) as exc:
            for row in ranked:
                row['screen']={'status':'unverified'}
            provenance['screen']={'status':'failed','reason':str(exc)}
    record=write_selection(ranked,out,provenance,args.top_k)
    print(json.dumps({'candidates':len(ranked),'groups':len(record['research_top']),
                      'research_selected':sum(map(len,record['research_top'].values())),
                      'screened_selected':sum(map(len,record['screened_top'].values())),'output':str(out)}))
    return 0
