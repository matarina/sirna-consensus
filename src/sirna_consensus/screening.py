"""Gene-aware sequence screening. A pass is a computational criterion only."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

from .fasta import sha256_file
from .types import reverse_complement


def fasta_records(path):
    name, parts = None, []
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith('>'):
                if name is not None:
                    yield name, ''.join(parts).upper().replace('U', 'T')
                name, parts = line[1:].strip().split()[0], []
            else:
                parts.append(line.strip())
    if name is not None:
        yield name, ''.join(parts).upper().replace('U', 'T')


def screen(rows: list[dict], manifest_path: Path, gene: str, executable: str,
           work: Path, timeout=1800, threads=4) -> dict:
    for row in rows:
        row['screen'] = {'status': 'unverified'}
    manifest = json.loads(manifest_path.read_text())
    mapping = manifest.get('transcript_genes', {})
    index = Path(manifest.get('bowtie_index') or '')
    reference = Path(manifest['reference_fasta'])
    if not index.is_absolute(): index = manifest_path.resolve().parent/index
    if not reference.is_absolute(): reference = manifest_path.resolve().parent/reference
    artifacts = manifest.get('bowtie_index_artifacts', {})
    suffixes = ('.1.ebwt','.2.ebwt','.3.ebwt','.4.ebwt','.rev.1.ebwt','.rev.2.ebwt')
    expected_names = {index.name + suffix for suffix in suffixes}
    if not gene or not mapping or set(artifacts) != expected_names:
        raise ValueError('screen requires target gene, transcript-gene mapping and all six pinned index files')
    if sha256_file(reference) != manifest['reference_sha256']:
        raise ValueError('reference checksum mismatch')
    for name, digest in artifacts.items():
        if Path(name).name != name or sha256_file(index.parent/name) != digest:
            raise ValueError('index checksum mismatch: ' + name)
    work.mkdir(parents=True, exist_ok=True)
    queries = work/'screen_queries.fa'
    # Query the documented target-binding region, not non-templated overhangs.
    queries.write_text(''.join(f'>q{i}\n{row["sense"]}\n' for i, row in enumerate(rows)))
    command = [executable, '-f', '-v', '2', '-a', '--norc', '--quiet', '-p', str(threads), str(index), str(queries)]
    alignments = work/'alignments.tsv'
    with alignments.open('w') as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.PIPE,
                       text=True, check=True, timeout=timeout)
    burdens = [{0: set(), 1: set(), 2: set()} for _ in rows]
    self_hit, unknown = set(), set()
    with alignments.open() as handle:
        for line in handle:
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 7 or not fields[0].startswith('q'):
                raise ValueError('malformed Bowtie alignment')
            i = int(fields[0][1:])
            if not 0 <= i < len(rows) or fields[1] != '+':
                raise ValueError('unexpected Bowtie query/strand')
            other_gene = mapping.get(fields[2])
            if not other_gene:
                unknown.add(i); continue
            n = len([x for x in fields[7].split(',') if x]) if len(fields) > 7 else 0
            if n > 2:
                raise ValueError('unexpected mismatch count')
            if other_gene == gene:
                if n == 0:
                    self_hit.add(i)
            else:
                burdens[i][n].add(other_gene)
    # Complete seed-only search across cDNA, not merely among Bowtie's <=2 hits.
    seeds = {reverse_complement(row['guide_dna'][1:8]) for row in rows}
    seed_genes = {seed: set() for seed in seeds}
    unmapped_seed = set()
    for transcript, sequence in fasta_records(reference):
        other_gene = mapping.get(transcript)
        if other_gene == gene:
            continue
        found = set(sequence[i:i+7] for i in range(len(sequence)-6)) & seeds
        for seed in found:
            if other_gene:
                seed_genes[seed].add(other_gene)
            else:
                unmapped_seed.add(seed)
    for i, row in enumerate(rows):
        seed = reverse_complement(row['guide_dna'][1:8])
        near = burdens[i][0] | burdens[i][1]
        complete = i in self_hit and i not in unknown
        row['screen'] = {'status': ('reject' if near else 'pass') if complete else 'unverified',
                         'query': row['sense'], 'query_scope': 'documented target-binding region',
                         'unintended_genes_by_mismatches': {str(n): sorted(v) for n,v in burdens[i].items()},
                         'seed_unintended_gene_count': len(seed_genes[seed]),
                         'seed_mapping_complete': seed not in unmapped_seed,
                         'intended_gene_exact_hit': i in self_hit}
    version = subprocess.run([executable, '--version'], capture_output=True, text=True, check=True, timeout=15)
    return {'command': command, 'version': version.stdout.splitlines()[0],
            'manifest_sha256': sha256_file(manifest_path), 'reference_sha256': manifest['reference_sha256'],
            'alignment_sha256': sha256_file(alignments),
            'criterion': 'no exact or one-mismatch site in an unintended gene; seed burden advisory'}
