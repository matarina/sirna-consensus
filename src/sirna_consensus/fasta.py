"""Single-transcript FASTA input handling.

The pipeline takes exactly one transcript: the tools all slide their own windows,
so the only thing this stage must produce is a validated ``(transcript_id, sequence)``
pair that later stages can use as a coordinate axis.
"""
from __future__ import annotations

import hashlib
import re
from pathlib import Path

MIN_LENGTH = 21
MAX_LENGTH = 200_000
_AMBIGUOUS = re.compile(r"[^ACGT]")


class InputError(ValueError):
    pass


def read_transcript(path: str | Path) -> tuple[str, str]:
    """Read exactly one FASTA record. Returns ``(record_id, sequence)``.

    ``U`` is normalised to ``T``; the header's first token is the transcript id.
    """
    text = Path(path).read_text()
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines or not lines[0].startswith(">"):
        raise InputError("input must start with a '>' FASTA header")
    if sum(line.startswith(">") for line in lines) != 1:
        raise InputError("input must contain exactly one record")
    name = lines[0][1:].split()[0] if lines[0][1:].strip() else ""
    sequence = "".join(lines[1:]).upper().replace("U", "T")
    if not name:
        raise InputError("FASTA header has no transcript id")
    if not sequence:
        raise InputError("FASTA record has no sequence")
    if _AMBIGUOUS.search(sequence):
        raise InputError("sequence must be unambiguous A/C/G/T/U")
    if not MIN_LENGTH <= len(sequence) <= MAX_LENGTH:
        raise InputError(f"transcript length must be {MIN_LENGTH}..{MAX_LENGTH} nt")
    return name, sequence


def sha256(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def sha256_file(path: str | Path) -> str:
    """Streaming hash of a reference file (transcriptomes can be multi-GB)."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()
