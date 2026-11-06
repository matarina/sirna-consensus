"""Real Bowtie integration on synthetic transcripts with known strand-specific hits."""

import json
from pathlib import Path

import pytest

from sirna_consensus.reference import build_reference
from sirna_consensus.screening import screen
from sirna_consensus.selection import enumerate_duplexes, rank_candidates
from sirna_consensus.types import reverse_complement


@pytest.mark.external
def test_both_strands_true_utr_seed_cache_and_corrupt_index(tmp_path, random_dna):
    tool = Path(__file__).resolve().parents[1] / "tools/bowtie-1.3.1/bowtie-1.3.1-linux-x86_64"
    if not (tool / "bowtie").is_file():
        pytest.skip("requires pinned local Bowtie 1.3.1")
    a = random_dna(150, seed=889)
    # Passenger has a full unintended target in B. Guide does not.
    core = a[30:49]
    b = random_dna(50, seed=890) + reverse_complement(core) + random_dna(81, seed=891)
    fasta = tmp_path / "tx.fa"
    fasta.write_text(f">TA\n{a}\n>TB\n{b}\n")
    gtf = tmp_path / "tx.gtf"
    gtf.write_text(
        "".join(
            f'chr1\ttest\texon\t1\t150\t.\t+\t.\tgene_id "G{g}"; transcript_id "T{g}"; gene_name "Gene{g}";\n'
            f'chr1\ttest\tCDS\t1\t90\t.\t+\t0\tgene_id "G{g}"; transcript_id "T{g}"; gene_name "Gene{g}";\n'
            for g in ("A", "B")
        )
    )
    out = tmp_path / "reference"
    manifest = build_reference(fasta, gtf, out, str(tool / "bowtie-build"))
    assert manifest["utr3_complete"]
    assert (out / "utr3.fa").read_text() == f">TA\n{a[90:]}\n>TB\n{b[90:]}\n"
    rows = enumerate_duplexes("TA", a)
    chosen = [r for r in rows if r["core"].start in (31, 95)]
    ranked = rank_candidates(chosen, "quality")
    metadata = screen(
        ranked, out / "screen_manifest.json", "GeneA", str(tool / "bowtie"), tmp_path / "screen"
    )
    assert metadata["seed_scope"] == "3UTR" and metadata["unique_queries"] == 4
    rejected = next(r for r in ranked if r["start"] == 31)
    passing = next(r for r in ranked if r["start"] == 95)
    assert rejected["screen"]["status"] == "reject"
    assert not rejected["screen"]["unintended_genes_by_mismatches"]["0"]
    assert rejected["screen"]["passenger_unintended_genes_by_mismatches"]["0"] == ["GB"]
    assert passing["screen"]["status"] == "pass"
    assert passing["screen"]["intended_isoform_coverage"] == 1
    # The cached counts must equal a direct search of other genes' true UTRs.
    for row in ranked:
        seed = reverse_complement(row["guide_dna"][1:8])
        assert row["screen"]["seed_unintended_gene_count"] == int(seed in b[90:])
    manifest["utr3_complete"] = False
    (out / "screen_manifest.json").write_text(json.dumps(manifest))
    screen(
        ranked, out / "screen_manifest.json", "GA", str(tool / "bowtie"), tmp_path / "incomplete"
    )
    assert all(r["screen"]["status"] == "unverified" for r in ranked)
    (out / "transcripts.1.ebwt").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="checksum"):
        screen(ranked, out / "screen_manifest.json", "GA", str(tool / "bowtie"), tmp_path / "bad")
