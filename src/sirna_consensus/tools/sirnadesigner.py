"""Adapter for anazhmetdin/siRNAdesigner.

Upstream: https://github.com/anazhmetdin/siRNAdesigner (rule-based, single file)

What it contributes: a small, fully explicit published rule set (GC window, no
internal 5-mer repeats, no GC stretch, guide/passenger end preferences, A/U
composition in the seed-adjacent region) plus a hand-crafted score. Its failure
mode is the opposite of a thermodynamic tool: it is blind to real duplex energy
and to off-targets.

Output decoding (the tool's own convention, verified against its source):

* it scans a 24 nt window ``cand = transcript[i:i+24]``;
* the passenger is ``cand[2:23]`` -> target positions ``i+3 .. i+23`` (1-based);
* the guide is the reverse complement of ``cand[0:21]``;
* CSV columns are ``start, end, cand, score, guide, passenger[, extras]``.

So the sense-strand target core is the passenger, and it starts at
``reported_start + 2``.
"""
from __future__ import annotations

import csv
from pathlib import Path

from ..types import RawHit, ToolRun
from .base import Tool

CANDIDATE_SCRIPT_SUFFIXES = (
    "siRNAdesigner-main/siRNADesigner.py",
    "siRNAdesigner/siRNADesigner.py",
    "siRNADesigner.py",
)


class SirnaDesignerTool(Tool):
    name = "sirnadesigner"
    summary = "rule-based siRNA scanner with an explicit published rule set"
    homepage = "https://github.com/anazhmetdin/siRNAdesigner"
    license = "MIT (upstream); see repository"

    def _script(self) -> str:
        if self.env.sirnadesigner_script:
            return self.env.sirnadesigner_script
        roots = [Path(__file__).resolve().parents[3], Path.cwd()]
        for root in roots:
            for suffix in CANDIDATE_SCRIPT_SUFFIXES:
                candidate = root / "vendor" / suffix
                if candidate.is_file():
                    return str(candidate)
        return ""

    def availability(self) -> tuple[bool, str]:
        script = self._script()
        if not script:
            return False, ("siRNADesigner.py not found; set SIRNA_DESIGNER_SCRIPT")
        has_deps, why = self.env.interpreter_has("Bio", "numpy")
        if not has_deps:
            return False, (f"interpreter {self.env.python} lacks a dependency: {why}")
        return True, f"script found at {script}"

    def version(self) -> str:
        return "upstream main (pinned by checkout)"

    def run(self, fasta: Path, workdir: Path, transcript_id: str) -> ToolRun:
        ok, detail = self.availability()
        if not ok:
            return self.wrap("unavailable", detail)
        script = self._script()
        shim = str(Path(__file__).with_name("_shim.py"))
        out_prefix = str(workdir) + "/"
        command = [self.env.python, shim, script,
                   "-i", str(fasta), "-o", out_prefix]
        self._execute(command, workdir, self.env.child_env())

        produced = sorted(workdir.glob("_*.csv"))
        combined: list[RawHit] = []
        for path in produced:
            combined.extend(self._parse(path))
        if not combined:
            return self.wrap("empty", "no candidate survived the tool's rule set")
        return self.wrap("ok", f"{len(combined)} hits from {len(produced)} file(s)",
                         combined)

    def _parse(self, path: Path) -> list[RawHit]:
        hits: list[RawHit] = []
        with path.open(newline="") as handle:
            for row in csv.reader(handle):
                if len(row) < 6:
                    continue
                try:
                    start = int(row[0])
                    score = float(row[3])
                except ValueError:
                    continue
                passenger = row[5].strip().upper()
                if not passenger:
                    continue
                hits.append(RawHit(
                    tool=self.name,
                    sense=passenger,
                    length=len(passenger),
                    start=start + 2,
                    score=score,
                    score_kind="designer_rule_score",
                    extra={"reported_start": start, "guide": row[4]},
                ))
        return hits
