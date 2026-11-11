"""Gene-aware guide/passenger screening. Pass means the pinned criterion only."""

from __future__ import annotations

import gzip
import json
import subprocess
from pathlib import Path

from .fasta import sha256_file
from .types import reverse_complement


def fasta_records(path):
    name, parts = None, []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts).upper().replace("U", "T")
                name, parts = line[1:].strip().split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, "".join(parts).upper().replace("U", "T")


def _resolve_gene(gene, mapping, names, aliases=None):
    genes = set(mapping.values())
    if gene in genes:
        return gene
    if aliases and gene in aliases:
        candidates = set(aliases[gene]) & genes
        if len(candidates) == 1:
            return candidates.pop()
    matches = {g for g in genes if g.split(".")[0] == gene or names.get(g) == gene}
    # GENCODE may have a PAR_Y copy of a gene: do not silently choose one.
    if len(matches) != 1:
        raise ValueError("target gene is absent or ambiguous in the reference: " + str(gene))
    return matches.pop()


def screen(
    rows: list[dict],
    manifest_path: Path,
    gene: str,
    executable: str,
    work: Path,
    timeout=1800,
    threads=4,
) -> dict:
    for row in rows:
        row["screen"] = {"status": "unverified"}
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    mapping = manifest.get("transcript_genes", {})

    def path(value):
        p = Path(value)
        return p if p.is_absolute() else manifest_path.parent / p

    reference = path(manifest["reference_fasta"])
    index = path(manifest.get("bowtie_index") or "")
    artifacts = manifest.get("bowtie_index_artifacts", {})
    small = {
        index.name + s
        for s in (".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")
    }
    large = {name + "l" for name in small}
    if not gene or not mapping or set(artifacts) not in (small, large):
        raise ValueError(
            "screen requires target gene, transcript-gene mapping and all six pinned index files"
        )
    gene = _resolve_gene(gene, mapping, manifest.get("gene_names", {}), manifest.get('gene_aliases'))
    if sha256_file(reference) != manifest["reference_sha256"]:
        raise ValueError("reference checksum mismatch")
    for name, digest in artifacts.items():
        if Path(name).name != name or sha256_file(index.parent / name) != digest:
            raise ValueError("index checksum mismatch: " + name)
    utr = path(manifest["utr3_fasta"]) if manifest.get("utr3_fasta") else reference
    if manifest.get("utr3_fasta") and sha256_file(utr) != manifest["utr3_sha256"]:
        raise ValueError("3UTR reference checksum mismatch")
    work.mkdir(parents=True, exist_ok=True)
    # Only documented paired regions enter full-length mismatch screening.
    sequences = sorted(
        {s for row in rows for s in (row["sense"], reverse_complement(row["sense"]))}
    )
    query_ids = {s: i for i, s in enumerate(sequences)}
    queries = work / "screen_queries.fa"
    queries.write_text("".join(f">q{i}\n{s}\n" for i, s in enumerate(sequences)))
    command = [
        executable,
        "-f",
        "-v",
        "2",
        "-a",
        "--norc",
        "--quiet",
        "-p",
        str(threads),
        str(index),
        str(queries),
    ]
    alignments = work / "alignments.tsv"
    with alignments.open("w") as handle:
        subprocess.run(
            command, stdout=handle, stderr=subprocess.PIPE, text=True, check=True, timeout=timeout
        )
    burdens = [{0: set(), 1: set(), 2: set()} for _ in sequences]
    self_hit, unknown, intended = set(), set(), {}
    with alignments.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7 or not fields[0].startswith("q"):
                raise ValueError("malformed Bowtie alignment")
            i = int(fields[0][1:])
            if not 0 <= i < len(sequences) or fields[1] != "+":
                raise ValueError("unexpected Bowtie query/strand")
            other = mapping.get(fields[2])
            if not other:
                unknown.add(i)
                continue
            n = len([x for x in fields[7].split(",") if x]) if len(fields) > 7 else 0
            if n > 2:
                raise ValueError("unexpected mismatch count")
            if other == gene:
                if n == 0:
                    self_hit.add(i)
                    intended.setdefault(i, set()).add(fields[2])
            else:
                burdens[i][n].add(other)
    seeds = {
        reverse_complement(s[1:8])
        for row in rows
        for s in (row["guide_dna"], row.get("passenger") or row["sense"])
    }
    seed_genes = {seed: set() for seed in seeds}
    seed_counts = {}
    unmapped_seed = set()
    if manifest.get("seed_index"):
        cache = path(manifest["seed_index"])
        if sha256_file(cache) != manifest["seed_index_sha256"]:
            raise ValueError("seed index checksum mismatch")
        with gzip.open(cache, "rt") as handle:
            seed_data = json.load(handle)
        if seed_data["utr3_sha256"] != sha256_file(utr):
            raise ValueError("seed index does not match UTR reference")
        own = bytes.fromhex(seed_data["gene_seed_masks"].get(gene, "00" * 2048))
        for seed in seeds:
            code = 0
            for base in seed:
                code = code * 4 + "ACGT".index(base)
            seed_counts[seed] = seed_data["counts"][seed] - bool(own[code // 8] & (1 << (code % 8)))
    for transcript, sequence in () if seed_counts else fasta_records(utr):
        other = mapping.get(transcript)
        if other == gene:
            continue
        found = {sequence[i : i + 7] for i in range(len(sequence) - 6)} & seeds
        for seed in found:
            if other:
                seed_genes[seed].add(other)
            else:
                unmapped_seed.add(seed)
    utr_complete = manifest.get("utr3_complete", True)
    intended_total = sum(g == gene for g in mapping.values())
    for row in rows:
        gi, pi = query_ids[row["sense"]], query_ids[reverse_complement(row["sense"])]
        gs = reverse_complement(row["guide_dna"][1:8])
        ps = reverse_complement((row.get("passenger") or row["sense"])[1:8])
        near = burdens[gi][0] | burdens[gi][1] | burdens[pi][0] | burdens[pi][1]
        complete = (
            gi in self_hit
            and gi not in unknown
            and pi not in unknown
            and utr_complete
            and gs not in unmapped_seed
            and ps not in unmapped_seed
        )
        # Historical unknown passenger is not a verified two-strand duplex screen.
        if manifest.get("schema_version", 1) >= 2 and not row.get("passenger"):
            complete = False
        row["screen"] = {
            "status": ("reject" if near else "pass") if complete else "unverified",
            "query": row["sense"],
            "query_scope": "both documented paired strands",
            "unintended_genes_by_mismatches": {str(n): sorted(burdens[gi][n]) for n in range(3)},
            "passenger_unintended_genes_by_mismatches": {
                str(n): sorted(burdens[pi][n]) for n in range(3)
            },
            "seed_unintended_gene_count": seed_counts.get(gs, len(seed_genes[gs])),
            "passenger_seed_unintended_gene_count": seed_counts.get(ps, len(seed_genes[ps])),
            "seed_scope": "3UTR" if manifest.get("utr3_fasta") else "whole_cDNA",
            "seed_mapping_complete": gs not in unmapped_seed
            and ps not in unmapped_seed
            and utr_complete,
            "intended_gene_exact_hit": gi in self_hit,
            "intended_transcripts_exact": sorted(intended.get(gi, set())),
            "intended_isoform_coverage": len(intended.get(gi, set())) / intended_total,
        }
    version = subprocess.run(
        [executable, "--version"], capture_output=True, text=True, check=True, timeout=15
    )
    return {
        "status": "complete",
        "command": command,
        "version": version.stdout.splitlines()[0],
        "target_gene": gene,
        "manifest_sha256": sha256_file(manifest_path),
        "reference_sha256": manifest["reference_sha256"],
        "unique_queries": len(sequences),
        "alignment_sha256": sha256_file(alignments),
        "criterion": "no exact or one-mismatch site in unintended genes on either paired strand; 2 mismatches and seed burden advisory",
        "seed_scope": "3UTR" if manifest.get("utr3_fasta") else "whole_cDNA",
    }
