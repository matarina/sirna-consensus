"""Frozen external scorer contract with artifact and candidate identity checks."""
from __future__ import annotations

import json
import math
import subprocess
from dataclasses import asdict
from pathlib import Path

from .fasta import sha256, sha256_file


def score_candidates(rows, transcript: str, manifest_path: Path, work: Path):
    manifest = json.loads(manifest_path.read_text())
    conv = manifest.get('input_convention')
    if conv != 'paired_antisense_19':
        raise ValueError('unsupported scorer/input convention: the frozen model requires paired_antisense_19')
    artifacts = manifest.get('artifacts', {})
    if not artifacts or not manifest.get('checkpoint_sha256') or not manifest.get('command'):
        raise ValueError('scorer manifest must pin artifacts, checkpoint and command')
    if manifest['checkpoint_sha256'] not in artifacts.values():
        raise ValueError('checkpoint must be included in checked artifacts')
    for path, digest in artifacts.items():
        artifact = Path(path) if Path(path).is_absolute() else manifest_path.parent/path
        if sha256_file(artifact) != digest:
            raise ValueError('scorer artifact checksum mismatch: ' + path)
    formats = manifest.get('supported_formats', [])
    supported = []
    for row in rows:
        core = row['core']
        # The frozen OligoFormer input is a 19-nt paired antisense core. A 21-mer
        # blunt molecule is not a supported input; it must first be canonicalised
        # to binding19 + a documented 2-nt overhang.
        if core.length == 19 and core.duplex_format in formats:
            supported.append({'key': core.key, 'core': core, 'paired_guide': core.guide[:19]})
    work.mkdir(parents=True, exist_ok=True)
    input_path, output_path = work/'scorer_input.json', work/'scorer_output.json'
    candidates = [{'key': item['key'], 'core': asdict(item['core']),
                   'paired_guide': item['paired_guide']} for item in supported]
    request = {'schema_version': 1, 'scorer_manifest_sha256': sha256_file(manifest_path),
               'transcript': transcript, 'transcript_sha256': sha256(transcript),
               'candidates': candidates}
    input_path.write_text(json.dumps(request, indent=2))
    command = [arg.replace('{input}', str(input_path.resolve())).replace('{output}', str(output_path.resolve()))
               for arg in manifest['command']]
    cache = work/'score_cache.json'
    reusable = (cache.is_file() and output_path.is_file() and
                json.loads(cache.read_text()).get('input_sha256') == sha256_file(input_path) and
                json.loads(cache.read_text()).get('output_sha256') == sha256_file(output_path))
    if supported:
        if not reusable:
            if output_path.exists():
                output_path.unlink()
            with (work/'stdout.log').open('w') as stdout, (work/'stderr.log').open('w') as stderr:
                subprocess.run(command, check=True, timeout=manifest.get('timeout_s', 1800),
                               cwd=str(manifest_path.parent), stdout=stdout, stderr=stderr)
        result = json.loads(output_path.read_text())
        if result.get('input_sha256') != sha256_file(input_path):
            raise ValueError('scorer output does not match frozen input')
        scores = result['scores']
    else:
        scores = {}
    if set(scores) != {item['key'] for item in supported}:
        raise ValueError('scorer output has missing or unexpected candidate keys')
    if any(not isinstance(v, (int, float)) or not math.isfinite(v) for v in scores.values()):
        raise ValueError('invalid scorer value')
    if supported:
        cache.write_text(json.dumps({'input_sha256':sha256_file(input_path), 'output_sha256':sha256_file(output_path)}))
    # 'replayed' is deliberately absent: it flips with cache state and would make an
    # otherwise-identical selection.json (and any freeze pinning it) drift.
    return scores, {'name': manifest['name'], 'manifest_sha256': sha256_file(manifest_path),
                    'checkpoint_sha256': manifest['checkpoint_sha256'], 'command': command,
                    'input_sha256': sha256_file(input_path), 'supported_count': len(supported),
                    'unsupported_count': len(rows)-len(supported),
                    'training_audit': manifest.get('training_audit', 'unverified'),
                    'independent_validation_passed': False}
