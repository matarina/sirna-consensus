"""Adapter for SIREN (pablovargasmejia/SIREN).

Upstream: https://github.com/pablovargasmejia/SIREN (GPLv3)

What it contributes: an **off-target screen** of every 21-mer against a supplied
transcriptome, via RNAhybrid, plus a region-level construct suggestion.

SIREN writes ``other_files/sirnas.fa`` containing one record per 21-mer, named
``sirna_<start>-<end>`` (sense) and ``sirna_<start>-<end>_r`` (guide), so
coordinates come straight from the record name and never need to be searched for.
Its ``off_targets_summary.tsv`` is aggregated per off-target with the siRNA names
listed, so we invert it to get a per-21-mer off-target burden.

Because SIREN screens *every* window, its vote means "this site survived an
off-target screen", not "this site scored well". That is stated plainly in the
report rather than dressed up as a ranking.

Selector resolution (``-g``)
----------------------------
SIREN is given a transcriptome with ``-T`` and a *substring selector* with ``-g``
that must match the first whitespace token of exactly one ``-T`` header. Our
input transcript id is usually a RefSeq accession (``NM_001657.4``) while the
transcriptome we ship is Ensembl cDNA (``ENST00000395748.8``), so the input id
alone never matches and SIREN dies with *"Gene not found"*.

When a transcriptome is configured this adapter therefore **resolves a selector
that actually occurs in that transcriptome**, in this order:

1. exact transcript id,
2. transcript id without its version suffix (accepted only if unique),
3. a gene symbol / gene id taken from the transcriptome header the transcript
   maps to (matched by sequence) or from ``gene_symbol:``/``gene:`` in the input
   FASTA header — resolved back to the unique transcript id it names,
4. the longest common usable token between the input id and the transcriptome.

Membership is verified against a header index built in one streaming pass over
the transcriptome (headers plus a cheap sequence digest; sequences are never all
held in memory). The index is memoised in-process and cached on disk, so a 1.4 GB
cDNA file is parsed once per machine, not once per call.

``SIREN_GENE_SELECTOR`` overrides the whole search and is passed through verbatim.

If nothing resolves — or if SIREN still refuses the resolved selector — the run
retries against the input FASTA as a **last resort**. That "works" but the
off-target screen then only ever sees the input, so the returned
:class:`~sirna_consensus.types.ToolRun` ``detail`` says which selector was used
and whether a real transcriptome-backed screen happened, rather than degrading
silently.
"""
from __future__ import annotations

import csv
import hashlib
import os
import re
from dataclasses import dataclass
from pathlib import Path

from ..types import RawHit, ToolRun
from .base import Tool, ToolError

NAME_RE = re.compile(r"^sirna_(\d+)-(\d+)$")

#: Header tokens an Ensembl/RefSeq FASTA may use to carry the gene.
GENE_RE = re.compile(r"\b(?:gene_symbol|gene):([^\s]+)")
#: Tokens too short/too generic to be useful as a substring selector.
_MIN_TOKEN = 4

ENV_SELECTOR = "SIREN_GENE_SELECTOR"


@dataclass(frozen=True)
class Resolution:
    """Outcome of selector resolution."""

    selector: str
    source: str
    #: transcriptome record the input sequence matched, when there is one
    mapped_id: str = ""
    gene: str = ""
    #: False means "no transcriptome-backed selector": screen the input only
    matched: bool = True

    def describe(self) -> str:
        bits = [f"selector {self.selector!r} ({self.source})"]
        if self.gene:
            bits.append(f"gene {self.gene}")
        if self.mapped_id:
            bits.append(f"maps to {self.mapped_id}")
        return ", ".join(bits)


def _digest(sequence: str) -> str:
    return hashlib.sha256(sequence.upper().encode()).hexdigest()[:32]


def _strip_version(identifier: str) -> str:
    return identifier.split(".", 1)[0]


def _tokens(text: str) -> list[str]:
    return [t for t in re.split(r"[^A-Za-z0-9]+", text) if len(t) >= _MIN_TOKEN]


def _gene_tokens(header: str) -> tuple[str, ...]:
    return tuple(match for match in GENE_RE.findall(header))


class TranscriptomeIndex:
    """Header index over a cDNA FASTA: ids, gene tokens, sequence digests.

    Built in a single streaming pass. Sequences are hashed and dropped, so peak
    memory is proportional to the number of *headers*, not to the file size.
    """

    def __init__(self, path: Path) -> None:
        self.path = path
        self.ids: set[str] = set()
        self.digests: dict[str, str] = {}
        self.genes: dict[str, list[str]] = {}
        self.gene_of: dict[str, tuple[str, ...]] = {}

    # -- construction ----------------------------------------------------
    def add(self, header: str, sequence: str) -> None:
        record_id = header.split()[0]
        self.ids.add(record_id)
        self.digests.setdefault(_digest(sequence), record_id)
        genes = _gene_tokens(header)
        if genes:
            self.gene_of[record_id] = genes
            for gene in genes:
                self.genes.setdefault(gene, []).append(record_id)

    def to_rows(self):
        by_id = {record_id: digest for digest, record_id in self.digests.items()}
        for record_id in self.ids:
            genes = self.gene_of.get(record_id, ())
            yield (by_id.get(record_id, ""), record_id,
                   genes[0] if genes else "", genes[1] if len(genes) > 1 else "")

    @classmethod
    def from_rows(cls, path: Path, rows) -> "TranscriptomeIndex":
        index = cls(path)
        for digest, record_id, symbol, gene_id in rows:
            index.ids.add(record_id)
            if digest:
                index.digests[digest] = record_id
            genes = tuple(g for g in (symbol, gene_id) if g)
            if genes:
                index.gene_of[record_id] = genes
                for gene in genes:
                    index.genes.setdefault(gene, []).append(record_id)
        return index

    # -- queries ---------------------------------------------------------
    def count(self, selector: str) -> int:
        """How many transcript ids a SIREN selector would match."""
        return sum(1 for record_id in self.ids if selector in record_id)

    def unique_id(self, selector: str) -> str:
        return selector if self.count(selector) == 1 else ""


def _cache_dir() -> Path:
    base = os.environ.get("SIRNA_CONSENSUS_CACHE") or os.environ.get("XDG_CACHE_HOME")
    root = Path(base) if base else Path.home() / ".cache"
    return root / "sirna-consensus"


_INDEX_CACHE: dict[str, TranscriptomeIndex] = {}


def _iter_records(path: Path):
    header: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def load_index(path: Path) -> TranscriptomeIndex:
    """Return the header index for ``path``, using the caches when possible."""
    path = path.resolve()
    try:
        stat = path.stat()
    except OSError as exc:
        raise ToolError(f"siren: cannot stat transcriptome {path}: {exc}") from exc
    key = f"{path}:{stat.st_size}:{int(stat.st_mtime)}"
    if key in _INDEX_CACHE:
        return _INDEX_CACHE[key]

    cache_file = _cache_dir() / f"siren-index-{hashlib.sha1(key.encode()).hexdigest()}.tsv"
    index: TranscriptomeIndex | None = None
    if cache_file.is_file():
        try:
            with cache_file.open(newline="") as handle:
                rows = (line.rstrip("\n").split("\t") for line in handle)
                index = TranscriptomeIndex.from_rows(path, rows)
        except OSError:
            index = None
    if index is None:
        index = TranscriptomeIndex(path)
        for header, sequence in _iter_records(path):
            index.add(header, sequence)
        _write_cache(cache_file, index)

    _INDEX_CACHE[key] = index
    return index


def _write_cache(cache_file: Path, index: TranscriptomeIndex) -> None:
    try:
        cache_file.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache_file.with_suffix(".tmp")
        with tmp.open("w") as handle:
            for row in index.to_rows():
                handle.write("\t".join(row) + "\n")
        tmp.replace(cache_file)
    except OSError:
        pass  # a read-only home is not a reason to fail the run


def resolve_selector(fasta: Path, transcript_id: str, transcriptome: Path,
                     override: str = "") -> Resolution:
    """Pick a ``-g`` selector that occurs in ``transcriptome``.

    See the module docstring for the order. Never raises for "not found": an
    unresolved selector is returned with ``matched=False`` so the caller can take
    the documented last resort.
    """
    index = load_index(transcriptome)
    if override:
        return Resolution(override, "override",
                          mapped_id=_mapped_id(index, fasta),
                          matched=index.count(override) == 1)

    if transcript_id in index.ids:
        return Resolution(transcript_id, "exact-id")

    base = _strip_version(transcript_id)
    if base != transcript_id:
        if base in index.ids:
            return Resolution(base, "exact-id")
        if index.count(base) == 1:
            return Resolution(base, "versionless-id")

    # The transcript our input corresponds to, if any: this is what makes a
    # RefSeq id usable against an Ensembl transcriptome at all.
    mapped_id = _mapped_id(index, fasta)

    gene_candidates: list[tuple[str, str]] = []
    for gene in _gene_tokens(_input_header(fasta)):
        if gene in index.genes:
            gene_candidates.append((gene, ""))
    if mapped_id:
        for gene in index.gene_of.get(mapped_id, ()):
            gene_candidates.append((gene, mapped_id))

    for gene, prefer in gene_candidates:
        named = [rid for rid in index.genes.get(gene, ()) if rid in index.ids]
        if not named:
            continue
        if prefer and prefer in named:
            return Resolution(prefer, "gene-symbol", mapped_id=mapped_id, gene=gene)
        # A gene symbol alone is matched by SIREN only against the record *id*,
        # which for Ensembl cDNA is an ENST accession, so hand over an id even
        # when the gene names several transcripts and we cannot tell which.
        if len(named) == 1 or not prefer:
            return Resolution(named[0], "gene-symbol", mapped_id=mapped_id, gene=gene)

    # Last header-based resort: the longest token shared with the transcriptome.
    for token in sorted(set(_tokens(transcript_id)), key=len, reverse=True):
        if index.count(token) == 1:
            return Resolution(token, "common-token", mapped_id=mapped_id)

    return Resolution(transcript_id, "unresolved", matched=False)


def _input_header(fasta: Path) -> str:
    try:
        with fasta.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    return line[1:].strip()
    except OSError:
        pass
    return ""


def _input_sequence(fasta: Path) -> str:
    try:
        with fasta.open() as handle:
            return "".join(line.strip() for line in handle
                           if line and not line.startswith(">"))
    except OSError:
        return ""


def _mapped_id(index: TranscriptomeIndex, fasta: Path) -> str:
    sequence = _input_sequence(fasta)
    return index.digests.get(_digest(sequence), "") if sequence else ""


class SirenTool(Tool):
    name = "siren"
    summary = "transcriptome off-target screen of every 21-mer (RNAhybrid)"
    homepage = "https://github.com/pablovargasmejia/SIREN"
    license = "GPL-3.0"

    def _exe(self) -> str:
        return self.env.which(self.env.siren_bin or "SIREN")

    def availability(self) -> tuple[bool, str]:
        exe = self._exe()
        if not exe:
            return False, "SIREN entry point not on PATH; set SIREN_BIN"
        if not self.env.which("RNAhybrid"):
            return False, "RNAhybrid binary not found; set SIRNA_CONSENSUS_TOOL_BIN"
        return True, f"SIREN at {exe}, RNAhybrid available"

    def version(self) -> str:
        return "siren-rnai 0.1.x"

    def run(self, fasta: Path, workdir: Path, transcript_id: str) -> ToolRun:
        ok, detail = self.availability()
        if not ok:
            return self.wrap("unavailable", detail)

        transcriptome = Path(self.env.transcriptome) if self.env.transcriptome else None
        if transcriptome and not transcriptome.is_file():
            transcriptome = None

        override = self.env.siren_gene_selector
        resolution: Resolution | None = None
        if transcriptome is None:
            # No reference configured: the input is all SIREN can screen, and its
            # own id is the only selector that can match it.
            selector = override or transcript_id
            self._invoke(fasta, selector, workdir)
            note = f"input only: no transcriptome configured (selector {selector!r})"
        else:
            resolution = resolve_selector(fasta, transcript_id, transcriptome,
                                          override=override)
            if not resolution.matched:
                selector = transcript_id
                self._invoke(fasta, selector, workdir)
                note = (f"fell back to the input FASTA only: no transcriptome "
                        f"selector matched (tried {resolution.selector!r}); "
                        f"off-targets are input-scoped, not transcriptome-wide")
            else:
                selector = resolution.selector
                try:
                    self._invoke(transcriptome, selector, workdir)
                except ToolError:
                    # Matched the headers, yet SIREN still refused it: fall back
                    # to the input rather than lose the run.
                    workdir.mkdir(parents=True, exist_ok=True)
                    self._invoke(fasta, transcript_id, workdir)
                    note = (f"fell back to the input FASTA only: transcriptome "
                            f"screen {resolution.describe()} failed; no "
                            f"transcriptome-backed screen happened")
                else:
                    note = f"screened against the transcriptome ({resolution.describe()})"

        burden = self._offtarget_burden(workdir)
        hits = self._parse_sirnas(workdir, burden)
        if not hits:
            return self.wrap("empty", "SIREN produced no 21-mer records")
        screened = f"{len(burden)} siRNAs carry off-targets" if burden else "no off-targets found"
        return self.wrap("ok", f"{len(hits)} 21-mers evaluated; {note}; {screened}", hits)

    def _invoke(self, targets: Path, selector: str, workdir: Path) -> None:
        command = [self._exe(), "-T", str(targets), "-g", selector,
                   "-o", str(workdir), "-t", "8"]
        self._execute(command, workdir, self.env.child_env())

    def _offtarget_burden(self, workdir: Path) -> dict[str, int]:
        """Invert off_targets_summary.tsv into siRNA -> number of off-targets."""
        summary = workdir / "off_targets_summary.tsv"
        burden: dict[str, int] = {}
        if not summary.is_file():
            return burden
        with summary.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                names = (row.get("siRNA names") or "").split(",")
                for raw in names:
                    name = raw.strip()
                    if not name:
                        continue
                    key = name[:-2] if name.endswith("_r") else name
                    burden[key] = burden.get(key, 0) + 1
        return burden

    def _parse_sirnas(self, workdir: Path, burden: dict[str, int]) -> list[RawHit]:
        sirnas = workdir / "other_files" / "sirnas.fa"
        if not sirnas.is_file():
            sirnas = workdir / "sirnas.fa"
        if not sirnas.is_file():
            return []
        hits: list[RawHit] = []
        name: str | None = None
        chunks: list[str] = []

        def flush() -> None:
            if name is None or not chunks:
                return
            match = NAME_RE.match(name)
            if not match:
                return
            start, end = int(match.group(1)), int(match.group(2))
            sequence = "".join(chunks).upper()
            if len(sequence) != end - start + 1:
                return
            hits.append(RawHit(
                tool=self.name,
                sense=sequence,
                length=len(sequence),
                start=start,
                score=1.0 / (1.0 + burden.get(name, 0)),
                score_kind="offtarget_cleanliness",
                extra={"offtargets": burden.get(name, 0)},
            ))

        with sirnas.open() as handle:
            for line in handle:
                line = line.rstrip("\n")
                if line.startswith(">"):
                    flush()
                    raw = line[1:].strip()
                    name = raw if NAME_RE.match(raw) else None
                    chunks = []
                elif name is not None:
                    chunks.append(line.strip())
        flush()
        return hits
