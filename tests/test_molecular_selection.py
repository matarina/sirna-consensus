import json
from dataclasses import asdict

import pytest

from sirna_consensus.types import Core, RawHit
from sirna_consensus.normalize import place, PlacedHit
from sirna_consensus.pipeline import Config, _specificity_status
from sirna_consensus.tools import get_tool
from sirna_consensus.types import ToolRun
from sirna_consensus.selection import rank_candidates, panels, load_candidates, write_selection


def test_staggered_duplex_survives_real_parser_normalizer_and_export(tmp_path):
    sequence='AGCTTACGATCGTACCTGATACGA'
    path=tmp_path/'_tx.csv'
    path.write_text('1,23,'+sequence+',7,TATCAGGTACGATCGTAAGCT,CTTACGATCGTACCTGATACG\n')
    hit=get_tool('sirnadesigner')._parse(path)[0]
    core=place(hit,'TX',sequence).core
    assert core.guide=='TATCAGGTACGATCGTAAGCT'
    assert core.passenger=='CTTACGATCGTACCTGATACG'
    assert core.guide_overhang=='CT' and core.passenger_overhang=='CG'
    ranked=rank_candidates([{'core':core,'tools':['designer'],'tool_scores':{'designer':7}}])
    result=write_selection(ranked,tmp_path/'out',{'transcript_id':'TX','target_sha256':'hash'})
    assert result['candidates'][0]['guide_dna']==core.guide
    assert core.guide.replace('T','U') in (tmp_path/'out/research_top.fa').read_text()


def test_wrong_coordinate_cannot_choose_between_repeats():
    motif='AGCTTACGATCGTACCTGATA'
    assert place(RawHit('a',motif,20,start=999),'TX',motif+'TTT'+motif).reason=='ambiguous_repeat'


def test_exact_offtarget_wins_over_clean_neighbor():
    seq='AGCTTACGATCGTACCTGATACGA'
    core=Core('TX',1,21,seq[:21]);neighbor=Core('TX',2,21,seq[1:22])
    def hit(c,n):return PlacedHit('siren',c,None,n==0,{'offtargets':n,'screen_state':'transcriptome_screened','transcriptome_backed':True})
    result=_specificity_status(core,[hit(core,10),hit(neighbor,0)],[ToolRun('siren','ok')],Config('x','y'))
    assert result[0]=='offtarget' and result[2]==10
    assert _specificity_status(core,[hit(neighbor,0)],[ToolRun('siren','ok')],Config('x','y'))[0]=='not_screened'


def test_replay_preserves_neighbor_scores_and_is_order_independent(tmp_path):
    seq='AGCTTACGATCGTACCTGATACGA';hits=[]
    for start,score in [(1,.9),(2,.1)]:
        c=Core('TX',start,21,seq[start-1:start+20])
        hits.append({'core':asdict(c),'tool':'a','role':'generator','score':score,'vote_eligible':True})
    p=tmp_path/'run.json';p.write_text(json.dumps({'schema_version':2,'normalized_hits':hits}))
    rows,_=load_candidates(p)
    one=rank_candidates(rows);two=rank_candidates(list(reversed(rows)))
    assert one==two
    assert one[0]['start']==1 and one[0]['tool_scores']['a']==.9
    assert len(one)==2


def test_unsupported_formats_never_share_model_panel():
    a=Core('TX',1,19,'AGCTTACGATCGTACCTGA');b=Core('TX',2,21,'GCTTACGATCGTACCTGATAC')
    rows=[{'core':c,'tools':['a'],'tool_scores':{'a':1}} for c in (a,b)]
    ranked=rank_candidates(rows,'model',model_scores={a.key:.9})
    assert len(panels(ranked))==2
    assert {r['ranking'] for r in ranked}=={'model','native_fallback'}
    assert ranked[1]['guide_dna']==b.guide


def test_screened_panel_backfills_without_relaxing_screen():
    seq='AGCTTACGATCGTACCTGATACGA'
    rows=rank_candidates([{'core':Core('TX',i+1,19,seq[i:i+19]),'tools':['a'],'tool_scores':{'a':10-i}} for i in range(5)])
    for i,r in enumerate(rows):r['screen']={'status':'reject' if i==0 else 'pass'}
    p=panels(rows,2,screened=True)
    assert list(p.values())[0]==[rows[1]['key'],rows[2]['key']]


def test_v1_is_rejected_for_scientific_replay(tmp_path):
    p=tmp_path/'old.json';p.write_text('{"schema_version":1}')
    with pytest.raises(ValueError,match='inspection-only'):load_candidates(p)


def test_legacy_support_ablation_changes_score():
    core=Core('TX',1,19,'AGCTTACGATCGTACCTGA')
    rows=[{'core':core,'tools':['a','b'],'tool_scores':{'a':1.,'b':1.}}]
    full=rank_candidates(rows,'legacy')[0]
    ablated=rank_candidates(rows,'legacy',components=('tools','quality'))[0]
    assert full['score']-ablated['score']==pytest.approx(4.)
