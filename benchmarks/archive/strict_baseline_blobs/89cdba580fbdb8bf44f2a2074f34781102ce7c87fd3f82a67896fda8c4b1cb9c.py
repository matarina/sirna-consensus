"""Adapter protocol for external design tools.

Each tool runs in its **own** environment and is reached over the process
boundary: FASTA in, tabular file out. Nothing is imported as a library, so a
tool's dependency pins can never leak into this package.

Adapters are responsible for their tool's quirks (strand, overhang, alphabet,
coordinate convention). Everything downstream only sees :class:`RawHit` with a
sense-strand DNA core.
"""
from __future__ import annotations

import os
import shutil
from abc import ABC, abstractmethod
from pathlib import Path

from ..types import RawHit, ToolRun


class ToolError(RuntimeError):
    pass


class Tool(ABC):
    """One external candidate generator."""

    name: str = ""
    summary: str = ""
    homepage: str = ""
    license: str = ""
    role: str = "generator"          # generator | evaluator
    acceptance: str = ""

    def __init__(self, env: "ToolEnv | None" = None) -> None:
        self.env = env or ToolEnv()

    # -- availability ----------------------------------------------------
    @abstractmethod
    def availability(self) -> tuple[bool, str]:
        """Return ``(usable, human readable reason)`` without running anything."""

    def version(self) -> str:
        return ""

    # -- execution -------------------------------------------------------
    @abstractmethod
    def run(self, fasta: Path, workdir: Path, transcript_id: str) -> list[RawHit]:
        """Run the tool and return its raw hits. Raise ToolError on failure."""

    # -- helpers ---------------------------------------------------------
    def _execute(self, command: list[str], workdir: Path,
                 extra_env: dict[str, str] | None = None) -> str:
        import subprocess
        self.last_command = list(command)
        env = dict(os.environ)
        if extra_env:
            env.update(extra_env)
        try:
            completed = subprocess.run(
                command, cwd=str(workdir), capture_output=True, text=True,
                timeout=self.env.timeout_s, env=env)
        except (OSError, subprocess.SubprocessError) as exc:
            raise ToolError(f"{self.name}: {exc}") from exc
        if completed.returncode != 0:
            tail = (completed.stderr or completed.stdout or "").strip().splitlines()
            raise ToolError(f"{self.name}: exit {completed.returncode}: "
                            + " | ".join(tail[-3:]))
        return completed.stdout

    def wrap(self, status: str, detail: str = "",
             hits: list[RawHit] | None = None) -> ToolRun:
        return ToolRun(tool=self.name, status=status, detail=detail,
                       hits=hits or [], version=self.version(),
                       command=getattr(self, "last_command", []),
                       role=self.role, acceptance=self.acceptance)


def _positive_int(raw: str | None, fallback: int) -> int:
    try:
        value = int(raw) if raw else fallback
    except ValueError:
        return fallback
    return value if value > 0 else fallback


class ToolEnv:
    """Resolved from the environment so a deployment can point at its own installs."""

    def __init__(self) -> None:
        self.python: str = os.environ.get("SIRNA_CONSENSUS_PYTHON", "") or __import__(
            "sys").executable
        self.timeout_s: int = int(os.environ.get("SIRNA_CONSENSUS_TIMEOUT", "1800"))
        self.tool_bin: str = os.environ.get("SIRNA_CONSENSUS_TOOL_BIN", "")
        self.tool_lib: str = os.environ.get("SIRNA_CONSENSUS_TOOL_LIB", "")
        self.sirnadesigner_script: str = os.environ.get("SIRNA_DESIGNER_SCRIPT", "")
        self.siren_bin: str = os.environ.get("SIREN_BIN", "")
        self.siren_threads: int = _positive_int(
            os.environ.get("SIREN_THREADS"), min(os.cpu_count() or 8, 32))
        self.sirnaforge_image: str = os.environ.get(
            "SIRNAFORGE_IMAGE", "ghcr.io/austin-s-h/sirnaforge:sha-4c05e4b")
        self.transcriptome: str = os.environ.get("SIRNA_CONSENSUS_TRANSCRIPTOME", "")
        # Explicit override for SIREN's ``-g`` selector, used verbatim.
        self.siren_gene_selector: str = os.environ.get("SIREN_GENE_SELECTOR", "")
        #: Maximum off-targets for a SIREN hit to count as clean evidence.
        self.siren_max_offtargets: int = _positive_int(
            os.environ.get("SIREN_MAX_OFFTARGETS"), 0)

    def child_env(self) -> dict[str, str]:
        """Environment additions the legacy tools need (RNAhybrid, etc.)."""
        extra: dict[str, str] = {}
        if self.tool_bin:
            extra["PATH"] = self.tool_bin + os.pathsep + os.environ.get("PATH", "")
        if self.tool_lib:
            extra["LD_LIBRARY_PATH"] = self.tool_lib + os.pathsep + os.environ.get(
                "LD_LIBRARY_PATH", "")
        return extra

    def which(self, *candidates: str) -> str:
        """Locate an executable, honouring the configured tool bin directory."""
        search = os.environ.get("PATH", "")
        if self.tool_bin:
            search = self.tool_bin + os.pathsep + search
        for candidate in candidates:
            if not candidate:
                continue
            if os.sep in candidate:
                if os.access(candidate, os.X_OK):
                    return candidate
                continue
            found = shutil.which(candidate, path=search)
            if found:
                return found
        return ""

    def interpreter_has(self, *modules: str) -> tuple[bool, str]:
        """Check the *tool* interpreter, not the one running this process."""
        import subprocess
        code = ";".join(f"import {m}" for m in modules)
        try:
            completed = subprocess.run([self.python, "-c", code],
                                       capture_output=True, text=True, timeout=60)
        except (OSError, subprocess.SubprocessError) as exc:
            return False, str(exc)
        if completed.returncode == 0:
            return True, ""
        lines = (completed.stderr or "").strip().splitlines()
        return False, lines[-1] if lines else "import failed"
