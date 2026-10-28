"""Optional smoke tests for the real external adapters.

These are marked ``external`` and run only when the caller explicitly enables
them (SIRNA_CONSENSUS_EXTERNAL_TESTS=1) because they need real tool installs,
RNAhybrid, and Docker. The hermetic test suite does not depend on them.
"""
from __future__ import annotations

import os

import pytest

from sirna_consensus.pipeline import Config, run
from sirna_consensus.tools import survey


def _external_enabled() -> bool:
    return os.environ.get("SIRNA_CONSENSUS_EXTERNAL_TESTS") == "1"


@pytest.mark.external
def test_each_registered_tool_is_available_or_explicitly_skipped():
    if not _external_enabled():
        pytest.skip("set SIRNA_CONSENSUS_EXTERNAL_TESTS=1 to run real-tool smoke tests")
    rows = survey()
    assert rows, "no registered tools"
    for _name, usable, detail in rows:
        assert isinstance(usable, bool)
        assert detail


@pytest.mark.external
def test_real_adapters_run_on_a_small_fixture(fasta_file, tmp_path):
    if not _external_enabled():
        pytest.skip("set SIRNA_CONSENSUS_EXTERNAL_TESTS=1 to run real-tool smoke tests")
    rows = survey()
    usable = [name for name, ok, _ in rows if ok]
    if not usable:
        pytest.skip("no external tool is installed")
    out = tmp_path / "external-smoke"
    result = run(Config(
        fasta=str(fasta_file),
        out_dir=str(out),
        tools=tuple(usable),
        min_votes=1,
        top_k=1,
    ))
    assert result["counts"]["hits"] >= 0
    assert (out / "run.json").is_file()
    assert (out / "top_guides.fa").is_file()
