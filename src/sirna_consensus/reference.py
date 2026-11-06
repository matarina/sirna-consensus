"""Pinned reference construction, transcript/gene crosswalk and spliced 3'UTRs."""

from __future__ import annotations

import csv
import gzip
import io
import itertools
import json
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from .fasta import sha256_file
from .screening import fasta_records


def canonicalize_genes(mapping, names, hgnc):
    """Collapse annotated alternate-locus copies only through pinned HGNC identity."""
    approved, aliases = {}, defaultdict(set)
    with Path(hgnc).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            approved[row["symbol"]] = row["hgnc_id"]
            for field in ("prev_symbol", "alias_symbol"):
                for name in row.get(field, "").split("|"):
                    if name:
                        aliases[name].add(row["hgnc_id"])
    symbols = {ident: symbol for symbol, ident in approved.items()}
    canonical, new_names, crosswalk = {}, {}, defaultdict(set)
    for original, name in names.items():
        matches = {approved[name]} if name in approved else aliases.get(name, set())
        resolved = next(iter(matches)) if len(matches) == 1 else original
        canonical[original] = resolved
        new_names[resolved] = symbols.get(resolved, name)
        for alias in (original, original.split(".")[0], name, new_names[resolved]):
            if alias:
                crosswalk[alias].add(resolved)
    return (
        {tid: canonical[g] for tid, g in mapping.items()},
        new_names,
        {alias: sorted(ids) for alias, ids in crosswalk.items()},
    )


def build_seed_index(utr_path, mapping, destination):
    seeds = ["".join(s) for s in itertools.product("ACGT", repeat=7)]
    codes = {s: i for i, s in enumerate(seeds)}
    by_gene = defaultdict(set)
    for tid, sequence in fasta_records(utr_path):
        gene = mapping[tid]
        by_gene[gene].update(
            codes[s] for s in {sequence[i : i + 7] for i in range(len(sequence) - 6)} if s in codes
        )
    counts, masks = [0] * len(seeds), {}
    for gene, present in sorted(by_gene.items()):
        mask = bytearray(len(seeds) // 8)
        for i in present:
            counts[i] += 1
            mask[i // 8] |= 1 << (i % 8)
        masks[gene] = mask.hex()
    with (
        Path(destination).open("wb") as raw,
        gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as compressed,
        io.TextIOWrapper(compressed) as handle,
    ):
        json.dump(
            {
                "counts": dict(zip(seeds, counts)),
                "gene_seed_masks": masks,
                "seed_order": "ACGT_lexicographic_7mer",
                "utr3_sha256": sha256_file(utr_path),
            },
            handle,
        )


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()


def gtf_transcripts(path):
    result = {}
    with text_open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] not in ("exon", "CDS", "stop_codon"):
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', parts[8]))
            tid = attrs.get("transcript_id")
            if not tid or not attrs.get("gene_id"):
                raise ValueError("GTF feature without transcript/gene identity")
            t = result.setdefault(
                tid,
                {
                    "gene": attrs["gene_id"],
                    "name": attrs.get("gene_name", ""),
                    "strand": parts[6],
                    "exons": [],
                    "coding": [],
                },
            )
            t["exons" if parts[2] == "exon" else "coding"].append(
                (int(parts[3]) - 1, int(parts[4]))
            )
    return result


def utr3(sequence, transcript):
    """Project CDS+stop onto transcript coordinates, then take downstream exonic bases."""
    exons = sorted(set(transcript["exons"]), reverse=transcript["strand"] == "-")
    if sum(b - a for a, b in exons) != len(sequence):
        return None
    ends, offset = [], 0
    for a, b in exons:
        for c, d in transcript["coding"]:
            lo, hi = max(a, c), min(b, d)
            if lo < hi:
                ends.append(offset + (hi - a if transcript["strand"] == "+" else b - lo))
        offset += b - a
    return sequence[max(ends) :] if ends else ""


def build_reference(
    fasta,
    gtf,
    out,
    bowtie_build="bowtie-build",
    release="custom",
    hgnc=None,
    threads=4,
    index_bmax=None,
):
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    transcripts = gtf_transcripts(gtf)
    mapping, names, seen, utr_missing = {}, {}, set(), []
    target, utr_path = out / "transcripts.fa", out / "utr3.fa"
    with target.open("w") as full, utr_path.open("w") as utr, text_open(fasta) as handle:
        tid, chunks = None, []

        def emit(header, chunks):
            if header is None:
                return
            # GENCODE headers use pipes; custom headers may use bare transcript IDs.
            ident = header.split("|")[0]
            if ident in seen:
                raise ValueError("duplicate reference transcript: " + ident)
            seen.add(ident)
            tx = transcripts.get(ident)
            if not tx:
                raise ValueError("FASTA transcript missing from GTF: " + ident)
            sequence = "".join(chunks).upper().replace("U", "T")
            if set(sequence) - set("ACGTN"):
                raise ValueError("invalid reference bases")
            mapping[ident] = tx["gene"]
            names.setdefault(tx["gene"], tx["name"])
            full.write(f">{ident}\n{sequence}\n")
            part = utr3(sequence, tx)
            if part is None:
                utr_missing.append(ident)
            elif part:
                utr.write(f">{ident}\n{part}\n")

        for line in handle:
            if line.startswith(">"):
                emit(tid, chunks)
                tid, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
        emit(tid, chunks)
    if not mapping:
        raise ValueError("empty reference")
    index = out / "transcripts"
    if threads < 1 or (index_bmax is not None and index_bmax < 1):
        raise ValueError("reference build threads and index-bmax must be positive")
    command = [str(bowtie_build), "--threads", str(threads)]
    if index_bmax is not None:
        command.extend(["--noauto", "--bmax", str(index_bmax)])
    command.extend([str(target), str(index)])
    with (out / "build.log").open("w") as log:
        subprocess.run(command, stdout=log, stderr=log, check=True, timeout=7200)
    artifacts = sorted(out.glob("transcripts*.ebwt*"))
    if len(artifacts) != 6:
        raise ValueError("expected six Bowtie1 index artifacts")
    version = subprocess.run(
        [str(bowtie_build), "--version"], capture_output=True, text=True, check=True
    ).stdout.splitlines()[0]
    manifest = {
        "schema_version": 2,
        "release": release,
        "reference_fasta": target.name,
        "reference_sha256": sha256_file(target),
        "utr3_fasta": utr_path.name,
        "utr3_sha256": sha256_file(utr_path),
        "utr3_complete": not utr_missing,
        "utr3_unprojectable": utr_missing,
        "transcript_genes": mapping,
        "gene_names": names,
        "bowtie_index": index.name,
        "bowtie_index_artifacts": {p.name: sha256_file(p) for p in artifacts},
        "source_artifacts": {str(Path(p).resolve()): sha256_file(Path(p)) for p in (fasta, gtf)},
        "build_command": command,
        "build_version": version,
        "source_urls": {
            "fasta": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.transcripts.fa.gz",
            "gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz",
        }
        if release == "GENCODE_v50_GRCh38.p14"
        else {},
    }
    if hgnc:
        manifest["source_artifacts"][str(Path(hgnc).resolve())] = sha256_file(Path(hgnc))
        mapping, canonical_names, crosswalk = canonicalize_genes(mapping, names, hgnc)
        manifest.update(
            transcript_genes=mapping,
            gene_names=canonical_names,
            gene_aliases=crosswalk,
            gene_identity="pinned_HGNC_when_unambiguous; GENCODE_ID_otherwise",
        )
    seed_index = out / "seed_index.json.gz"
    build_seed_index(utr_path, mapping, seed_index)
    manifest.update(seed_index=seed_index.name, seed_index_sha256=sha256_file(seed_index))
    (out / "screen_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def add_parser(sub):
    p = sub.add_parser("build-reference", help="build pinned Bowtie1 and true 3UTR reference")
    p.add_argument("--fasta", required=True)
    p.add_argument("--gtf", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--bowtie-build", default="bowtie-build")
    p.add_argument("--release", default="custom")
    p.add_argument("--hgnc")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument(
        "--index-bmax", type=int, help="explicit Bowtie block size; large values require more RAM"
    )


def command(args):
    manifest = build_reference(
        args.fasta,
        args.gtf,
        args.out,
        args.bowtie_build,
        args.release,
        args.hgnc,
        args.threads,
        args.index_bmax,
    )
    print(json.dumps({"transcripts": len(manifest["transcript_genes"]), "output": args.out}))
    return 0
