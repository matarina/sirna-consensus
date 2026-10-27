"""Adapter for siRNAforge, run from its pinned OCI image.

Upstream: https://github.com/Austin-s-h/sirnaforge (MIT)

What it contributes: thermodynamic scoring (ViennaRNA duplex energy, end
asymmetry, melting temperature) plus optional transcriptome/miRNA off-target
counts — a family of signals the rule-based tools cannot produce.

The image is pinned by digest on purpose: an unpinned `latest` would silently
change the candidate set between runs. ``sirnaforge design`` is used by default
because it skips the (expensive) transcriptome index build; set
``SIRNAFORGE_SUBCMD=workflow`` to get off-target columns as well.

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
from .base import Tool, ToolError

SAFE_NAME = re.compile(r"[^A-Za-z0-9_.-]")


class SirnaforgeTool(Tool):
    name = "sirnaforge"
    summary = "thermodynamic duplex scoring and off-target counts (OCI image)"
    homepage = "https://github.com/Austin-s-h/sirnaforge"
    license = "MIT"

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

        subcommand = os.environ.get("SIRNAFORGE_SUBCMD", "design")
        top_n = os.environ.get("SIRNAFORGE_TOP_N", "500")
        cores = os.environ.get("SIRNAFORGE_CORES", "8")
        safe = SAFE_NAME.sub("_", transcript_id)
        staged = workdir / f"{safe}.fa"
        shutil.copyfile(fasta, staged)

        command = [
            "docker", "run", "--rm",
            "-u", f"{os.getuid()}:{os.getgid()}",
            "-v", f"{workdir}:/work",
            self.env.sirnaforge_image,
            "sirnaforge", subcommand,
            "--input-fasta", f"/work/{staged.name}",
            "--output-dir", "/work",
            "--top-n", top_n,
            "--cores", cores,
        ]
        transcriptome = os.environ.get("SIRNAFORGE_TRANSCRIPTOME",
                                       self.env.transcriptome)
        if transcriptome and subcommand == "workflow":
            command += ["--transcriptome-fasta", transcriptome, "--species", "human"]
        self._execute(command, workdir)

        rows = self._parse(workdir)
        if not rows:
            return self.wrap("empty", "no parsable candidate table produced")
        return self.wrap("ok", f"{len(rows)} hits", rows)

    def _parse(self, workdir: Path) -> list[RawHit]:
        hits: list[RawHit] = []
        for path in sorted(workdir.rglob("*.csv")):
            if path.name.startswith(".") or not path.is_file():
                continue
            try:
                with path.open(newline="") as handle:
                    reader = csv.DictReader(handle)
                    fields = set(reader.fieldnames or ())
                    if "passenger_sequence" not in fields:
                        continue
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
                            extra={"source_file": path.name,
                                   "off_target_count": row.get("off_target_count")},
                        ))
            except OSError as exc:
                raise ToolError(f"sirnaforge: cannot read {path}: {exc}") from exc
        return hits
