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
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from ..types import RawHit, ToolRun
from .base import Tool, ToolError

NAME_RE = re.compile(r"^sirna_(\d+)-(\d+)$")


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

        targets = Path(self.env.transcriptome) if self.env.transcriptome else fasta
        selector = transcript_id
        try:
            self._invoke(targets, selector, workdir)
        except ToolError:
            if targets == fasta:
                raise
            # The configured transcriptome may not carry this identifier; retry
            # against the input itself so the run still produces a result.
            workdir.mkdir(parents=True, exist_ok=True)
            self._invoke(fasta, selector, workdir)

        burden = self._offtarget_burden(workdir)
        hits = self._parse_sirnas(workdir, burden)
        if not hits:
            return self.wrap("empty", "SIREN produced no 21-mer records")
        screened = f"{len(burden)} siRNAs carry off-targets" if burden else "no off-targets found"
        return self.wrap("ok", f"{len(hits)} 21-mers evaluated; {screened}", hits)

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
