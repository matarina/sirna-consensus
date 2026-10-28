import csv
import json

import pytest

from sirna_consensus.benchmark import (load_measurements,metric,evaluate_record,
                                       paired_summary,discrimination)

A='AGCTTACGATCGTACCTGATA';B='GCTTACGATCGTACCTGATAC';FMT='guide21-passenger_unknown'

def measured(guide=A,effective=True,gene='G'):
    return {'guide_dna':guide,'duplex_format':FMT,'effective':effective,'gene':gene,
            'study':'S','assay':'a','endpoint':'percent_inhibition','transcript_id':'TX','target_sha256':'H'}


def record():
    candidates=[{'key':str(i),'guide_dna':guide,'duplex_format':FMT,'eligible':True,
                 'group':FMT+'/native','score':2-i} for i,guide in enumerate([B,A])]
    return {'schema_version':2,'kind':'exact_guide_selection','config':{'top_k':2},'candidates':candidates,
            'provenance':{'transcript_id':'TX','target_sha256':'H'},
            'research_top':{FMT+'/native':['0','1']}}


def test_shifted_guide_cannot_receive_experimental_credit():
    m=metric([measured(B)],[measured(A)],5)
    assert m['measured']==0 and m['precision_at_k'] is None
    assert m['recovery_at_k']==0


def test_unmeasured_candidates_change_coverage_not_precision():
    m=metric([measured(B),measured(A)],[measured(A)],5)
    assert m['precision_at_k']==1 and m['measurement_coverage']==.5
    assert m['fill_rate']==.4 and m['observed_hits_per_budget']==.2


def test_measured_projection_and_deployment_have_different_topk():
    r=record();m=[measured()]
    a=evaluate_record(r,m,k=1,track='measured',random_repeats=0)[0]
    b=evaluate_record(r,m,k=1,track='deployment',random_repeats=0)[0]
    assert a['precision_at_k']==1 and b['precision_at_k'] is None


def test_other_transcript_does_not_change_recovery():
    a=measured();b={**measured(B),'transcript_id':'OTHER','gene':'H'}
    one=evaluate_record(record(),[a],random_repeats=0)
    two=evaluate_record(record(),[a,b],random_repeats=0)
    assert one==two


def test_label_validation_and_replicates(tmp_path):
    row={**measured(),'record_id':'r1','value':'NA'};row.pop('effective')
    p=tmp_path/'measurements.tsv'
    def write(rows):
        with p.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    write([row])
    with pytest.raises(ValueError):load_measurements(p)
    row['value']='75';write([row]);assert load_measurements(p)[0]['effective']
    write([row,{**row,'record_id':'r2'}])
    with pytest.raises(ValueError,match='replicate_group'):load_measurements(p)
    row['replicate_group']='technical1';write([row,{**row,'record_id':'r2','value':'65'}])
    assert load_measurements(p)[0]['value']==70


def test_normalized_activity_not_clipped(tmp_path):
    row={**measured(),'record_id':'r1','value':'1.252','endpoint':'normalized_activity'};row.pop('effective')
    p=tmp_path/'m.tsv'
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(row),delimiter='\t');w.writeheader();w.writerow(row)
    assert load_measurements(p)[0]['value']==1.252


def test_ties_and_single_class_are_explicit():
    assert discrimination([measured()],{(A,FMT):1})['auroc'] is None
    scores={(A,FMT):.5,(B,FMT):.5}
    result=discrimination([measured(),measured(B,False)],scores)
    assert result['auroc']==.5 and result['average_precision']==.5


def test_paired_comparison_refuses_missing_assays():
    rows=evaluate_record(record(),[measured()],random_repeats=0)
    with pytest.raises(ValueError,match='prespecified'):paired_summary(rows,[])
    out=paired_summary(rows,rows,replicates=100)
    assert not out['numeric_criterion_met']
    assert out['mean_difference'] is None


def test_deployment_budget_must_have_been_generated():
    with pytest.raises(ValueError,match='saved selection budget'):
        evaluate_record(record(),[measured()],k=5,track='deployment')


def test_paired_comparison_checks_construct_identity():
    rows=evaluate_record(record(),[measured()],k=1,random_repeats=0)
    changed=[{**r,'target_sha256':'different_construct'} for r in rows]
    with pytest.raises(ValueError,match='prespecified'):paired_summary(rows,changed)


def test_one_gene_cannot_pass_improvement_gate():
    rows=evaluate_record(record(),[measured()],k=1,random_repeats=0)
    baseline=[{**r,'precision_at_k':0.,'observed_hits_per_budget':0.} for r in rows]
    result=paired_summary(rows,baseline,replicates=100)
    assert result['mean_difference']==1.
    assert result['evidence_status']=='insufficient_genes'
    assert not result['numeric_criterion_met']
