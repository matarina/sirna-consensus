"""Run the whole weight-independent final evaluation as one reproducible command.

Steps, in order: build the audited holdout, select with the development-selected
policy (candidate) and the scorer alone (baseline), freeze the lock, run the final
gate for both arms and write a paired comparison plus an honest summary.  Nothing
is re-selected after labels are read; a lock that still verifies is reused, and only
recut when its inputs actually changed.

Usage: python benchmarks/scripts/run_final.py [--top-k 10]
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
sys.path.insert(0, str(REPO / "src"))
from sirna_consensus.fasta import read_transcript, sha256  # noqa: E402

PY = str(REPO / ".venv/bin/python")
SELECT = str(REPO / ".venv/bin/sirna-consensus")
SCORER_PY = ".scorer-venv/bin/python"   # relative: a leading '/' would make REPO/p discard REPO
ARMS = {"policy": [], "model": []}


def run(cmd, env=None):
    subprocess.run(cmd, check=True, cwd=REPO, env=env)


def cut_lock(lock: Path, files) -> None:
    """Write the lock atomically so a concurrent reader never sees it absent.

    freeze.py refuses to overwrite, so build a per-process temp lock and os.replace.
    """
    tmp = lock.with_name(lock.name + f".{os.getpid()}.tmp")
    tmp.unlink(missing_ok=True)
    run([PY, str(ROOT / "scripts/freeze.py"), "--lock", str(tmp)] + [str(f) for f in files])
    os.replace(tmp, lock)


def select(out_root, fasta, candidates, scorer, policy, top_k, evidence):
    args = [SELECT, "select", "--fasta", str(fasta), "--candidates", str(candidates),
            "--out", str(out_root), "--scorer-manifest", str(scorer), "--top-k", str(top_k)]
    args += ["--policy", str(policy)] if policy else ["--strategy", "model"]
    if evidence.is_file():
        args += ["--run-evidence", str(evidence)]
    else:
        print(f"warning: no native evidence for {fasta.stem}; arm uses the scorer alone")
    run(args)


def _read_lock(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text()).get("artifacts", {})
    except (json.JSONDecodeError, AttributeError):
        return {}                       # malformed/partial lock -> treat as absent


def needs_recut(lock: Path, expected: set, summary_files) -> bool:
    """Re-cut when the file set changed OR any summary file's bytes changed.

    A key-set-only check would leave a stale REVIEW.md pinned once its content changes
    (e.g. under a different --permutations), so the digest is compared too. An absent or
    unparseable lock, or an absent summary, always recuts (no crash).
    """
    if isinstance(summary_files, Path):
        summary_files = [summary_files]
    artifacts = _read_lock(lock)
    if set(artifacts) != expected:
        return True
    for summary_file in summary_files:
        if not summary_file.is_file():
            return True
        rel = os.path.relpath(summary_file.resolve(), lock.parent.resolve())
        if artifacts.get(rel) != hashlib.sha256(summary_file.read_bytes()).hexdigest():
            return True
    return False


def gene_weighted(path):
    """Mean complete Top-5 panel precision, genes weighted equally."""
    results = json.loads(Path(path).read_text())["results"]
    per_gene = {}
    for r in results:
        if r["requested_k"] == 5 and r["complete_measured_panel"]:
            per_gene.setdefault(r["gene"], []).append(r["precision_at_k"])
    return {g: mean(v) for g, v in per_gene.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--top-k", type=int, default=10)
    p.add_argument("--permutations", type=int, default=1000,
                   help="within-unit label-permutation control per evaluation (0 disables)")
    p.add_argument("--regen-generators", action="store_true",
                   help="force regeneration of native tool runs; run id is content-derived, so byte "
                        "churn is limited to tool-version drift, but it does recut the freeze")
    a = p.parse_args()
    holdout = ROOT / "data/holdout"
    results = ROOT / "results/holdout"
    scorer = ROOT / "models/oligoformer_manifest.json"
    policy = ROOT / "results/selection_policy.json"
    cohort = ROOT / "cohort_manifest.json"
    protocol = ROOT / "protocol.json"
    measured = holdout / "measurements.tsv"
    tool_runs = ROOT / "results/runs/holdout"

    run([PY, str(ROOT / "scripts/prepare_holdout.py")])
    # Native generator evidence for the holdout (label-blind), so the dev-selected
    # policy is actually exercised rather than collapsing to the scorer alone.
    # Regenerate only when missing or forced: the run id is content-derived, but a
    # forced regeneration still recuts the freeze, so reuse keeps it stable.
    def _stale(fa: Path) -> bool:
        record = tool_runs / fa.stem / "run.json"
        if not record.is_file():
            return True
        try:
            _, sequence = read_transcript(fa)
            return json.loads(record.read_text()).get("input_sha256") != sha256(sequence)
        except (OSError, ValueError, json.JSONDecodeError):
            return True

    missing = [fa for fa in (holdout / "transcripts").glob("*.fa") if _stale(fa)]
    if a.regen_generators or missing:
        env = {**os.environ, "SIRNA_CONSENSUS_PYTHON": str(REPO / SCORER_PY)}
        run([PY, str(ROOT / "scripts/run_generators.py"), "--transcripts", str(holdout / "transcripts"),
             "--out", str(tool_runs), "--workers", "2"], env=env)
    # A degraded generator (missing interpreter/image) still writes a run.json, so
    # assert liveness: otherwise the arm silently loses evidence and looks fine.
    # "empty" is a legitimate outcome (the tool found nothing); "ok" is normal.
    degraded = []
    for fa in (holdout / "transcripts").glob("*.fa"):
        record = tool_runs / fa.stem / "run.json"
        if not record.is_file():
            degraded.append((fa.stem, "run.json", "missing"))
            continue
        degraded += [(fa.stem, t["tool"], t["status"])
                     for t in json.loads(record.read_text())["tools"]
                     if t["status"] not in ("ok", "empty")]
    if degraded:
        raise SystemExit("generator not ok: " + "; ".join(f"{g}:{tool}={st}" for g, tool, st in degraded))
    # Prune stale per-arm outputs so a transcript dropped from the cohort cannot leave a
    # selection.json that gets re-pinned while evaluating nothing. Selections are
    # deterministic, so re-selecting reproduces the same bytes.
    # Prune stale per-arm outputs so a transcript dropped from the cohort cannot leave a
    # selection.json that gets re-pinned while evaluating nothing. Selections are
    # deterministic, so re-selecting reproduces the same bytes. Only the arm dirs are
    # cleared (not REVIEW.md/final_*.json) so a valid lock still verifies and is reused.
    results.mkdir(parents=True, exist_ok=True)
    for arm in ARMS:
        shutil.rmtree(results / arm, ignore_errors=True)
    # Drop evidence for transcripts that left the cohort, or it gets re-pinned orphaned.
    live = {fa.stem for fa in (holdout / "transcripts").glob("*.fa")}
    if tool_runs.is_dir():
        for stale in tool_runs.iterdir():
            if stale.is_dir() and stale.name not in live:
                shutil.rmtree(stale, ignore_errors=True)
    for arm in ARMS:
        for fasta in sorted((holdout / "transcripts").glob("*.fa")):
            select(results / arm / fasta.stem, fasta, holdout / "candidates.tsv", scorer,
                   policy if arm == "policy" else None, a.top_k,
                   tool_runs / fasta.stem / "run.json")

    files = [measured, holdout / "candidates.tsv", cohort, ROOT / "exposure_audit.json", protocol,
             scorer, policy, ROOT / "evaluate.py",
             holdout / "dataset_manifest.json", holdout / "exclusions.json"]
    files += [ROOT / "scripts/audit_exposure.py", ROOT / "scripts/prepare_holdout.py",
              ROOT / "scripts/run_final.py", ROOT / "scripts/run_generators.py",
              ROOT / "scripts/freeze.py", ROOT / "scripts/verify_freezes.py",
              ROOT / "scripts/compare.py", ROOT / "evaluate.py"]
    files += sorted((holdout / "transcripts").glob("*.fa"))
    files += sorted(results.glob("*/[!.]*/selection.json"))
    files += sorted(tool_runs.glob("*/run.json"))
    files += sorted((REPO / "src/sirna_consensus").rglob("*.py"))   # incl. tools/ adapters
    files += [REPO / "vendor/siRNAdesigner/siRNADesigner.py"]        # decides designer hits
    files = list(dict.fromkeys(files))                              # evaluate.py was listed twice
    lock = ROOT / "final_freeze.json"
    # Never overwrite a lock that still verifies AND still covers exactly this file set;
    # a newly-added input whose bytes touch nothing pinned must not be silently unpinned.
    # REVIEW.md is written after the first freeze (its text needs the comparison), so it
    # is expected here too — otherwise covers would mismatch and recut every run.
    summary_files = [results / "REVIEW.md", results / "comparison.json"]
    covers = False
    if lock.exists():
        expected = {os.path.relpath(Path(f).resolve(), lock.parent.resolve())
                    for f in [*files, *summary_files]}
        covers = set(_read_lock(lock)) == expected
    verified = covers and subprocess.run(
        [PY, str(ROOT / "scripts/freeze.py"), "--lock", str(lock), "--verify"],
        cwd=REPO, capture_output=True).returncode == 0
    if not verified:
        cut_lock(lock, files)

    finals = {}
    for arm in ARMS:
        runs = [str(f) for f in sorted((results / arm).glob("*/selection.json"))]
        out = results / f"final_{arm}.json"
        cmd = [PY, str(ROOT / "evaluate.py"), "--phase", "final", "--freeze-lock", str(lock),
               "--cohort-manifest", str(cohort), "--protocol", str(protocol),
               "--benchmark", str(measured), "--output", str(out)]
        if a.permutations:
            cmd += ["--permutations", str(a.permutations)]
        for r in runs:
            cmd += ["--run", r]
        run(cmd)
        finals[arm] = out
    comparison = results / "comparison.json"
    run([PY, str(ROOT / "scripts/compare.py"), "--candidate", str(finals["policy"]),
         "--baseline", str(finals["model"]), "--out", str(comparison)])

    policy_genes, model_genes = gene_weighted(finals["policy"]), gene_weighted(finals["model"])
    controls = [c["p_one_sided"] for c in json.loads(finals["policy"].read_text())["controls"]
                if c.get("p_one_sided") is not None]
    shared = sorted(set(policy_genes) & set(model_genes))
    diffs = [policy_genes[g] - model_genes[g] for g in shared]
    cohort_m = json.loads(cohort.read_text())
    built = json.loads((holdout / "dataset_manifest.json").read_text())
    arm_files = sorted((results / "policy").glob("*/selection.json"))
    attached = sum(json.loads(f.read_text())["provenance"].get("evidence_attached", 0) for f in arm_files)
    candidates = sum(len(json.loads(f.read_text())["candidates"]) for f in arm_files)
    lines = [
        "# Weight-independent final evaluation", "",
        "Reproduce with `python benchmarks/scripts/run_final.py` "
        "(the lock is reused when it still verifies, else recut).", "",
        f"- Audited cohort: {len(cohort_m['holdout_studies'])} studies / {cohort_m['holdout_records']} records "
        f"/ {cohort_m['holdout_genes']} genes, criterion `{cohort_m['criterion']}`",
        f"- **Actually evaluated**: {built['records']} records / {built['genes']} genes / "
        f"{len(built['studies'])} studies ({built['exclusions']} excluded — see `data/holdout/exclusions.json`)",
        f"- native generator evidence attached to {attached}/{candidates} policy-arm candidates "
        f"({100*attached/max(1,candidates):.0f}%)",
        f"- checkpoint_training_audit: `{cohort_m['checkpoint_training_audit']}`; "
        f"independent_holdout_eligible: `{cohort_m['independent_holdout_eligible']}`",
        f"- upstream evaluation exposure: {cohort_m['reason'].split(';')[-1].strip()}", "",
        "## Development-selected policy vs frozen scorer alone", "",
        f"Complete Top-5 panels on {len(shared)} shared genes. "
        f"Paired mean gene-weighted precision difference "
        f"**{(mean(diffs)*100 if diffs else 0):+.1f} pt**.", "",
        f"Within-unit label-permutation control: {len(controls)} evaluable contrasts, "
        f"min one-sided p = {min(controls) if controls else 'n/a'}.", "",
        "| gene | policy | model | difference |", "|---|---:|---:|---:|",
    ]
    for g in shared:
        lines.append(f"| {g} | {policy_genes[g]:.2f} | {model_genes[g]:.2f} | {policy_genes[g]-model_genes[g]:+.2f} |")
    lines += ["",
              "This is a small, single-construct-per-molecule cohort; treat it as a directional "
              "holdout check, not a powered estimate. All numbers are computational.", "",
              "Native generator evidence is attached only where a generator emitted the same 19-nt "
              "binding site as an assayed guide, so coverage is partial and the `quality` strategy "
              "still leans on the scorer; a small difference can be structural, not only data.", "",
              f"Artifacts: `cohort_manifest.json`, `exposure_audit.json`, `protocol.json`, "
              f"`final_freeze.json`, `{comparison.relative_to(ROOT)}`."]
    (results / "REVIEW.md").write_text("\n".join(lines) + "\n")
    # Pin the summaries too (byte-stable across reused runs). They are written after the
    # first freeze because their text depends on the comparison, so re-cut once to add
    # them; subsequent runs already cover them and leave the lock untouched.
    expected = {os.path.relpath(Path(f).resolve(), lock.parent.resolve())
                for f in [*files, *summary_files]}
    if needs_recut(lock, expected, summary_files):
        cut_lock(lock, [*files, *summary_files])
    print(json.dumps({"shared_genes": len(shared), "mean_diff": mean(diffs) if diffs else None,
                      "review": str(summary_files[0])}))


if __name__ == "__main__":
    main()
