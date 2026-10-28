import json
import hashlib
import sys
from pathlib import Path
import pytest
from sirna_consensus.screening import screen
from sirna_consensus.types import Core
from sirna_consensus.selection import from_record


def test_overhang_without_its_strand_is_invalid():
    with pytest.raises(ValueError,match='ordered guide'):
        Core('TX',1,19,'A'*19,guide_overhang='TT')


def test_explicit_and_implicit_same_guide_have_same_key():
    core=Core('TX',1,19,'A'*19)
    explicit=Core('TX',1,19,'A'*19,ordered_guide='T'*19)
    assert core.key==explicit.key and core.molecule_id==explicit.molecule_id


def test_report_tsv_aliases_preserve_ordered_strands():
    c=from_record({'transcript_id':'TX','target_start':'1','core_length':'19',
                   'target_sequence_dna':'A'*19,'guide_sequence_rna':'U'*19+'AC',
                   'passenger_sequence_rna':'A'*19+'UU','guide_overhang':'AC','passenger_overhang':'UU'})
    assert c.guide=='T'*19+'AC' and c.passenger=='A'*19+'TT'


def test_passenger_geometries_are_distinct():
    a=Core('TX',1,19,'A'*19)
    b=Core('TX',1,19,'A'*19,passenger='A'*19)
    c=Core('TX',1,19,'A'*19,passenger='A'*19+'TT',passenger_overhang='TT')
    assert len({r.duplex_format for r in (a,b,c)})==3


def test_hashes_for_a_different_index_are_rejected(tmp_path):
    ref=tmp_path/'ref.fa';ref.write_text('>TX\n'+'A'*30+'\n')
    manifest={'reference_fasta':str(ref),'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),
              'bowtie_index':str(tmp_path/'actual'),'transcript_genes':{'TX':'G'},
              'bowtie_index_artifacts':{f'wrong.{i}.ebwt':'fake' for i in range(6)}}
    p=tmp_path/'manifest.json';p.write_text(json.dumps(manifest))
    with pytest.raises(ValueError,match='six pinned'):
        screen([],p,'G','not-invoked',tmp_path/'work')
