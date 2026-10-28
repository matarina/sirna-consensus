"""Adapter for siRNAforge, run from its pinned OCI image.

Upstream: https://github.com/Austin-s-h/sirnaforge (MIT)

What it contributes: thermodynamic scoring (ViennaRNA duplex energy, end
asymmetry, melting temperature) plus optional transcriptome/miRNA off-target
counts — a family of signals the rule-based tools cannot produce.

The image is pinned by digest on purpose: an unpinned ``latest`` would silently
change the candidate set between runs.

Verbatim CLI (not guessable — it takes the FASTA positionally, writes one table
via ``-o`` and has no ``--cores``)::

    sirnaforge design <input.fa> -o /work/out.csv --top-n N --length 21

Caveats worth knowing before comparing it to the other tools:

* ``--length`` is a single value, so this adapter emits one core length per run
  (``SIRNAFORGE_LENGTH``, default 21); the others may emit several.
* the tool applies **its own** GC (30-60%) and poly-run filters first, so its hit
  list is already narrowed — a site it dropped cannot be told apart from a site
  it scored low.
* the output is named ``.tsv`` upstream but is actually comma separated, so the
  delimiter is sniffed rather than assumed.

Output decoding uses siRNAforge's own native columns: ``passenger_sequence`` is
the sense target window, ``position`` its 1-based start, ``composite_score`` the
tool score.
"""
from __future__ import annotations

import csv
import os
import re
import shutil
import subprocess
from pathlib import Path

from ..types import RawHit, ToolRun
from .base import Tool

SAFE_NAME = re.compile(r"[^A-Za-z0-9_.-]")


class SirnaforgeTool(Tool):
    name = "sirnaforge"
    summary = "thermodynamic duplex scoring and off-target counts (OCI image)"
    homepage = "https://github.com/Austin-s-h/sirnaforge"
    license = "MIT"
    role = "generator"
    acceptance = "thermodynamic candidate emitted by the pinned image"

    def availability(self) -> tuple[bool, str]:
        if not shutil.which("docker"):
            return False, "docker not on PATH"
        try:
            completed = subprocess.run(
                ["docker", "image", "inspect", self.env.sirnaforge_image],
                capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError) as exc:
            return False, f"docker unavailable: {exc}"
        if completed.returncode != 0:
            return False, (f"image {self.env.sirnaforge_image} not pulled locally "
                           f"(`docker pull` it first)")
        return True, f"image {self.env.sirnaforge_image} present"

    def version(self) -> str:
        return self.env.sirnaforge_image.rsplit(":", 1)[-1][:12]

    def run(self, fasta: Path, workdir: Path, transcript_id: str) -> ToolRun:
        ok, detail = self.availability()
        if not ok:
            return self.wrap("unavailable", detail)

        top_n = os.environ.get("SIRNAFORGE_TOP_N", "500")
        length = os.environ.get("SIRNAFORGE_LENGTH", "21")
        safe = SAFE_NAME.sub("_", transcript_id)
        staged = workdir / f"{safe}.fa"
        shutil.copyfile(fasta, staged)
        out_file = workdir / "sirnaforge_out.csv"

        # `sirnaforge design` takes the FASTA positionally, writes a single table
        # via -o, and has no --cores. The container runs as the invoking uid so
        # the output stays readable to us.
        command = [
            "docker", "run", "--rm",
            "-u", f"{os.getuid()}:{os.getgid()}",
            "-v", f"{workdir}:/work",
            self.env.sirnaforge_image,
            "sirnaforge", "design", f"/work/{staged.name}",
            "-o", f"/work/{out_file.name}",
            "--top-n", top_n,
            "--length", length,
        ]
        self._execute(command, workdir)

        rows = self._parse(out_file)
        if not rows:
            return self.wrap("empty", "no parsable candidate table produced")
        run = self.wrap("ok", f"{len(rows)} hits", rows)
        run.selection_limits = {
            "top_n": top_n,
            "length": length,
            "image": self.env.sirnaforge_image,
        }
        return run

    def _parse(self, path: Path) -> list[RawHit]:
        if not path.is_file():
            return []
        text = path.read_text()
        if not text.strip():
            return []
        delimiter = "\t" if text.splitlines()[0].count("\t") else ","
        hits: list[RawHit] = []
        reader = csv.DictReader(text.splitlines(), delimiter=delimiter)
        if "passenger_sequence" not in (reader.fieldnames or []):
            return []
        for row in reader:
            sense = (row.get("passenger_sequence") or "").strip().upper()
            position = (row.get("position") or "").strip()
            if not sense or not position:
                continue
            try:
                start = int(float(position))
                score = float(row.get("composite_score") or 0.0)
            except ValueError:
                continue
            hits.append(RawHit(
                tool=self.name,
                sense=sense,
                length=len(sense),
                start=start,
                score=score,
                score_kind="composite_score",
                extra={"off_target_count": row.get("off_target_count"),
                       "passes_filters": row.get("passes_filters")},
            ))
        return hits
