"""Rebuild benchmarks/development_freeze.json from an explicit, present file set.

The dev lock pins generated artifacts whose bytes embed source hashes (e.g. a
selection.json records its ``implementation_hashes``), so it goes stale whenever the
package changes and needs a scripted rebuild rather than hand-editing. Pass
``--regenerate`` to first re-run the development pipeline.

Keys are stored relative to the lock directory (see freeze.py).
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]   # benchmarks/
REPO = ROOT.parent
sys.path.insert(0, str(REPO / "src"))
from sirna_consensus.fasta import sha256_file  # noqa: E402

LOCK = ROOT / "development_freeze.json"
# Explicit, ordered file set (globs expanded deterministically).
PATTERNS = [
    ROOT / "results/development/full/**/*",
    ROOT / "results/development/ablations/**/*",
    ROOT / "results/runs/huesken/*/run.json",   # native evidence the dev pipeline consumes
    ROOT / "results/runs/AREG_v2/screened/selection.json",   # screen the summarizer consumes
    ROOT / "results/development_comparison.json",
    ROOT / "results/development_results.json",
    ROOT / "results/permutation_control.json",
    ROOT / "results/selection_policy.json",
    ROOT / "results/implementation_verification.json",
    ROOT / "results/REVIEW.md",
    ROOT / "data/huesken/**/*",
    ROOT / "dataset_manifest.json", ROOT / "benchmark_config.json",
    ROOT / "README.md", ROOT / "data_sources.md",
    ROOT / "scripts/audit_exposure.py", ROOT / "scripts/audit_sources.py",
    ROOT / "scripts/verify_freezes.py", ROOT / "scripts/freeze.py",
    ROOT / "scripts/gen_development_freeze.py",
    ROOT / "scripts/run_development.py", ROOT / "scripts/run_generators.py",
    ROOT / "scripts/replay_ablations.py", ROOT / "scripts/summarize.py", ROOT / "scripts/compare.py",
    ROOT / "scripts/permutation_test.py", ROOT / "scripts/create_scorer_manifest.py",
    ROOT / "make_report_zh.py", ROOT / "report_zh.html",
    ROOT / "models/oligoformer_manifest.json",
    ROOT / "scripts/prepare_holdout.py", ROOT / "evaluate.py",
    REPO / "src/sirna_consensus/**/*.py",
    REPO / "vendor/siRNAdesigner/siRNADesigner.py",
    REPO / "README.md", REPO / "pyproject.toml", REPO / "LICENSE",
]


def files() -> list[Path]:
    import glob as globmod
    out = []
    for pattern in PATTERNS:
        text = str(pattern)
        if any(ch in text for ch in "*["):
            out += [Path(p) for p in sorted(globmod.glob(text, recursive=True)) if Path(p).is_file()]
        elif pattern.is_file():
            out.append(pattern)
    out = [f for f in out if f.resolve() != LOCK.resolve()]   # a lock cannot pin itself
    out = [f for f in out if "/work/" not in f.as_posix()]    # regenerated intermediates
    return list(dict.fromkeys(out))


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--regenerate", action="store_true",
                   help="re-run the development pipeline before re-hashing")
    p.add_argument("--data", type=Path, default=ROOT / "data/huesken")
    p.add_argument("--tool-runs", type=Path, default=ROOT / "results/runs/huesken")
    p.add_argument("--scorer-manifest", type=Path, default=ROOT / "models/oligoformer_manifest.json")
    p.add_argument("--screen", type=Path, default=ROOT / "results/runs/AREG_v2/screened/selection.json")
    a = p.parse_args()
    py = str(REPO / ".venv/bin/python")
    if a.regenerate:
        full = ROOT / "results/development/full"
        measured = a.data / "measurements.tsv"
        # Regenerate the native evidence too (its bytes are pinned, and it is consumed).
        env = {**os.environ, "SIRNA_CONSENSUS_PYTHON": str(REPO / ".scorer-venv/bin/python")}
        subprocess.run([py, str(ROOT / "scripts/run_generators.py"), "--transcripts",
                        str(a.data / "transcripts"), "--out", str(a.tool_runs), "--workers", "2"],
                       check=True, cwd=REPO, env=env)
        degraded = [(rj.parent.name, t["tool"], t["status"])
                    for rj in sorted(a.tool_runs.glob("*/run.json"))
                    for t in json.loads(rj.read_text())["tools"]
                    if t["status"] not in ("ok", "empty")]
        if degraded:
            raise SystemExit("generator not ok: " + "; ".join(f"{g}:{t}={s}" for g, t, s in degraded))
        subprocess.run([py, str(ROOT / "scripts/run_development.py"), "--data", str(a.data),
                        "--tool-runs", str(a.tool_runs), "--scorer-manifest", str(a.scorer_manifest),
                        "--out", str(full)], check=True, cwd=REPO)
        surfaced = full / "selection_policy.json"
        if not surfaced.is_file():
            raise SystemExit("no strategy won on development; selection_policy.json absent")
        (ROOT / "results/selection_policy.json").write_text(surfaced.read_text())
        # Refresh the rest of the pinned dev outputs so the lock is internally consistent.
        subprocess.run([py, str(ROOT / "scripts/compare.py"), "--candidate", str(full / "quality.json"),
                        "--baseline", str(full / "model.json"),
                        "--out", str(ROOT / "results/development_comparison.json")], check=True, cwd=REPO)
        subprocess.run([py, str(ROOT / "scripts/permutation_test.py"), "--development", str(full),
                        "--benchmark", str(measured),
                        "--out", str(ROOT / "results/permutation_control.json")], check=True, cwd=REPO)
        subprocess.run([py, str(ROOT / "scripts/replay_ablations.py"), "--development", str(full),
                        "--benchmark", str(measured),
                        "--out", str(ROOT / "results/development/ablations")], check=True, cwd=REPO)
        # summarize tolerates a missing --screen (writes screen_smoke: null), so always run it.
        subprocess.run([py, str(ROOT / "scripts/summarize.py"), "--development", str(full),
                        "--comparison", str(ROOT / "results/development_comparison.json"),
                        "--ablations", str(ROOT / "results/development/ablations"),
                        "--screen", str(a.screen), "--out", str(ROOT / "results")], check=True, cwd=REPO)
        subprocess.run([py, str(ROOT / "make_report_zh.py"), "-o", str(ROOT / "report_zh.html")],
                       check=True, cwd=REPO)
    base = LOCK.parent.resolve()
    artifacts = {os.path.relpath(f.resolve(), base): sha256_file(f) for f in files()}
    LOCK.write_text(json.dumps({"artifacts": artifacts}, indent=2) + "\n")
    print(f"{LOCK}: {len(artifacts)} entries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
