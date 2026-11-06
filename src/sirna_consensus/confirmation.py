"""Preregistered evaluation contracts; absence of evidence never passes a gate."""

from __future__ import annotations

import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from statistics import mean

from .fasta import sha256_file


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def load_protocol(path, registered=True):
    path = Path(path).resolve()
    protocol = json.loads(path.read_text())
    registry_path, baseline_path = path.with_name("cohorts.json"), path.with_name("baseline.json")
    for source, field in (
        (registry_path, "cohort_registry_sha256"),
        (baseline_path, "baseline_manifest_sha256"),
    ):
        if sha256_file(source) != protocol[field]:
            raise ValueError("preregistered metadata changed: " + str(source))
    registry = json.loads(registry_path.read_text())
    baseline = json.loads(baseline_path.read_text())
    if baseline["commit"] != protocol["baseline_commit"]:
        raise ValueError("baseline commit mismatch")
    registration = None
    if registered:
        package = Path(__file__).resolve().parents[2]
        relative = str(path.relative_to(package))
        registration = subprocess.check_output(
            ["git", "log", "-1", "--format=%H", "--", relative], cwd=package, text=True
        ).strip()
        if not registration:
            raise ValueError("protocol must be committed before fitting or confirmation")
        original = subprocess.check_output(
            ["git", "show", registration + ":" + relative], cwd=package
        )
        if original != path.read_bytes():
            raise ValueError("uncommitted protocol change")
    train = protocol["training"]
    for key, digest_key in (("source", "source_sha256"), ("hgnc", "hgnc_sha256")):
        if sha256_file(Path(train[key])) != train[digest_key]:
            raise ValueError("training allowlist artifact mismatch: " + key)
    if train["studies"] != ["huesken2005"]:
        raise ValueError("this attempt permits Huesken-only development")
    if sha256_file(Path(baseline["default_policy_path"])) != baseline["default_policy_sha256"]:
        raise ValueError("baseline default policy changed during this attempt")
    return protocol, registry, baseline, registration


def ensure_training_records(records, protocol):
    if not records or any(r["study"] not in protocol["training"]["studies"] for r in records):
        raise ValueError("confirmation labels or unapproved cohort in training")
    if len({r["record_id"] for r in records}) != len(records):
        raise ValueError("duplicate training record ID")


def macro_gene(rows, field):
    genes = defaultdict(list)
    for r in rows:
        if r.get(field) is not None:
            genes[r["gene_id"]].append(r[field])
    return mean(mean(v) for v in genes.values()) if genes else None


def assay_metrics(records, scores):
    """Fixed measured pool; abstention loses DCG slots, never an assay from the pool."""
    import numpy as np
    from scipy.stats import rankdata, spearmanr

    if len(records) != len(scores):
        raise ValueError("score/record count mismatch")
    if len({r["record_id"] for r in records}) != len(records):
        raise ValueError("duplicate evaluation IDs")
    units = defaultdict(list)
    for r, score in zip(records, scores):
        units[r["unit"]].append((r, score))
    results = []
    for unit, pairs in sorted(units.items()):
        rows = [r for r, _ in pairs]
        values = np.asarray([r["value"] for r in rows], dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("nonfinite experimental measurement")
        rel = (
            (rankdata(values, method="average") - 1) / (len(rows) - 1)
            if len(rows) > 1
            else np.zeros(1)
        )
        relevance = dict(zip((r["record_id"] for r in rows), rel))
        finite = [(r, float(s)) for r, s in pairs if s is not None and math.isfinite(float(s))]
        finite.sort(
            key=lambda p: (-p[1], p[0]["guide"], p[0].get("molecule_id", ""), p[0]["record_id"])
        )
        top = [r for r, _ in finite[:5]]
        active = sum(r["value"] >= r["active_threshold"] for r in top)
        ndcg = None
        if len(rows) >= 5 and np.ptp(values) > 0:
            denom = sum(
                (2**v - 1) / math.log2(i + 2) for i, v in enumerate(sorted(rel, reverse=True)[:5])
            )
            ndcg = (
                sum(
                    (2 ** relevance[r["record_id"]] - 1) / math.log2(i + 2)
                    for i, r in enumerate(top)
                )
                / denom
            )
        p5 = active / len(top) if top else None
        rho = None
        if len(finite) >= 2:
            v, s = [r["value"] for r, _ in finite], [s for _, s in finite]
            if np.ptp(v) > 0 and np.ptp(s) > 0:
                rho = float(spearmanr(v, s).statistic)
        results.append(
            dict(
                unit=unit,
                study=rows[0]["study"],
                gene_id=rows[0]["gene_id"],
                n=len(rows),
                ndcg5=ndcg,
                ndcg_eligible=bool(len(rows) >= 5 and np.ptp(values) > 0),
                precision5=p5,
                top5_active_fraction=p5,
                top3_at_least_one_active=float(
                    any(r["value"] >= r["active_threshold"] for r in top[:3])
                ),
                fill_rate=len(top) / 5,
                prediction_coverage=len(finite) / len(rows),
                abstentions=len(rows) - len(finite),
                measurement_coverage=1.0 if top else None,
                budget_yield=active / 5,
                spearman=rho,
                top_ids=[r["record_id"] for r in top],
            )
        )
    return results


def summarize(rows):
    keys = (
        "ndcg5",
        "precision5",
        "top5_active_fraction",
        "top3_at_least_one_active",
        "fill_rate",
        "prediction_coverage",
        "measurement_coverage",
        "budget_yield",
        "spearman",
    )
    return {
        **{k: macro_gene(rows, k) for k in keys},
        "assays": len(rows),
        "genes": len({r["gene_id"] for r in rows}),
        "abstentions": sum(r["abstentions"] for r in rows),
        "per_gene_assay": rows,
    }


def paired_bootstrap(candidate, baseline, groups, replicates=10000, seed=20261106):
    import numpy as np

    a, b = {r["unit"]: r for r in candidate}, {r["unit"]: r for r in baseline}
    if len(a) != len(candidate) or len(b) != len(baseline) or set(a) != set(b):
        raise ValueError("paired analysis must retain identical assay units, including abstentions")
    deltas = defaultdict(list)
    for unit in a:
        if a[unit]["ndcg5"] is not None and b[unit]["ndcg5"] is not None:
            if a[unit]["gene_id"] != b[unit]["gene_id"]:
                raise ValueError("gene identity mismatch")
            deltas[a[unit]["gene_id"]].append(a[unit]["ndcg5"] - b[unit]["ndcg5"])
    clusters = defaultdict(list)
    for gene, values in deltas.items():
        if gene not in groups:
            raise ValueError("missing independent sequence component")
        clusters[str(groups[gene])].append(mean(values))
    if not clusters:
        return dict(delta=None, ci95=[None, None], components=0, genes=0)
    values = list(clusters.values())
    rng = np.random.default_rng(seed)
    boot = [
        mean(v for i in rng.integers(len(values), size=len(values)) for v in values[int(i)])
        for _ in range(replicates)
    ]
    return dict(
        delta=mean(v for vs in values for v in vs),
        ci95=[float(x) for x in np.quantile(boot, [0.025, 0.975])],
        components=len(values),
        genes=len(deltas),
    )


def molecule_identity(row):
    """No partial/core-overlap efficacy transfer; unknown chemistry cannot match."""
    if (
        not row.get("guide")
        or not row.get("passenger")
        or row.get("chemistry") in (None, "", "unknown")
    ):
        return None
    if row.get("mapping_status") not in (
        "unique",
        "exact",
        "original_construct_exact",
        "source_accession_exact",
    ):
        return None
    return tuple(
        row.get(k, "")
        for k in ("guide", "passenger", "guide_overhang", "passenger_overhang", "chemistry")
    )


def panel_bounds(panel, measurements, k=5):
    if len({r["unit"] for r in measurements}) > 1:
        raise ValueError("never transfer activity across assay conditions")
    if len(panel) > k or len({json.dumps(x, sort_keys=True) for x in panel}) != len(panel):
        raise ValueError("panel exceeds budget or has duplicate records")
    truth = {}
    for m in measurements:
        key = molecule_identity(m)
        if key is not None:
            if key in truth:
                raise ValueError(
                    "ambiguous repeated molecule/assay; aggregate documented replicates first"
                )
            truth[key] = m["value"] >= m["active_threshold"]
    active, unknown, inactive = set(), set(), set()
    keys = set()
    for i, row in enumerate(panel):
        key = molecule_identity(row)
        if key is None:
            key = ("unresolved_identity", str(i), json.dumps(row, sort_keys=True))
        if key in keys:
            raise ValueError("duplicate molecule in panel")
        keys.add(key)
        if key in truth:
            (active if truth[key] else inactive).add(key)
        else:
            unknown.add(key)
    n = len(panel)
    first3 = panel[:3]
    known_top3 = [truth.get(molecule_identity(x)) for x in first3]
    hit3 = (
        1.0
        if any(x is True for x in known_top3)
        else (None if any(x is None for x in known_top3) else 0.0)
    )
    return dict(
        active=active,
        unknown=unknown,
        inactive=inactive,
        fill_rate=n / k,
        measurement_coverage=(n - len(unknown)) / n if n else None,
        precision5=(len(active) / n if n else None) if not unknown else None,
        top5_active_fraction=(len(active) / n if n else None) if not unknown else None,
        top3_at_least_one_active=hit3,
        yield_bounds=[len(active) / k, (len(active) + len(unknown)) / k],
    )


def paired_budget_bounds(candidate, baseline, measurements, k=5):
    a, b = panel_bounds(candidate, measurements, k), panel_bounds(baseline, measurements, k)
    # Same fully specified molecule has the same unknown outcome in the same assay.
    common = {x for x in a["unknown"] & b["unknown"] if x[0] != "unresolved_identity"}
    delta = len(a["active"]) - len(b["active"])
    return [(delta - len(b["unknown"] - common)) / k, (delta + len(a["unknown"] - common)) / k]


def decision(protocol, evidence):
    """Missing values fail; diagnostics cannot set independent eligibility."""
    gate = protocol["gate"]

    def above(x, value):
        return x is not None and math.isfinite(float(x)) and x > value

    checks = dict(
        independence=evidence.get("independence_verified") is True,
        provenance=evidence.get("provenance_verified") is True,
        minimum_studies=evidence.get("independent_studies", 0)
        >= gate["minimum_independent_studies"],
        minimum_components=evidence.get("independent_components", 0)
        >= gate["minimum_independent_components"],
        absolute_ndcg=above(evidence.get("candidate_ndcg5"), gate["ndcg_strictly_above"]),
        positive_gain=above(evidence.get("ndcg_gain"), gate["gain_strictly_above"]),
        positive_ci=above(
            evidence.get("paired_ci95_lower"), gate["paired_ci95_lower_strictly_above"]
        ),
        no_completion_loss=evidence.get("no_completion_loss") is True,
        no_coverage_loss=evidence.get("no_prediction_coverage_loss") is True,
        actual_panel_budget_gain=above(
            evidence.get("actual_budget_gain_lower"),
            gate["actual_generated_budget_gain_lower_strictly_above"],
        ),
        independent_recomputation=evidence.get("independent_recomputation_passed") is True,
    )
    passed = all(checks.values())
    return dict(
        promote=passed,
        action="eligible_for_atomic_replacement" if passed else "retain_baseline",
        independent_validation_passed=passed,
        checks=checks,
        failures=[k for k, v in checks.items() if not v],
        evidence=evidence,
        development_metrics_used_for_gate=False,
    )


def require_confirmation_ids(records, registry, protocol):
    allowed = set(protocol["confirmation_ids"])
    if not allowed or {r["record_id"] for r in records} != allowed:
        raise ValueError("no eligible preregistered confirmation IDs, or confirmation ID mismatch")
    eligible = [
        c
        for c in registry["cohorts"]
        if c.get("role") == "confirmation" and c.get("independent_eligible") is True
    ]
    if not eligible:
        raise ValueError("no audited independent confirmation cohort")
    studies = {s for c in eligible for s in c.get("studies", [])}
    if any(r["study"] not in studies for r in records):
        raise ValueError("excluded study presented as confirmation")
    for c in eligible:
        if not all(
            c.get(k) is True
            for k in (
                "baseline_training_clear",
                "candidate_training_clear",
                "historical_selection_clear",
            )
        ):
            raise ValueError("incomplete baseline/candidate/history exposure audit")
        if sha256_file(Path(c["path"])) != c["sha256"]:
            raise ValueError("confirmation source changed")
    import csv

    original = {}
    for c in eligible:
        with Path(c["path"]).open() as h:
            for row in csv.DictReader(h, delimiter="\t"):
                if row["record_id"] in original:
                    raise ValueError("duplicate authoritative confirmation ID")
                original[row["record_id"]] = row
    for row in records:
        source = original.get(row["record_id"])
        if (
            source is None
            or source["study"] != row["study"]
            or float(source["value"]) != row["value"]
        ):
            raise ValueError("confirmation values do not reproduce authoritative source")
        if source["guide_dna"].upper().replace("U", "T") != row["guide"]:
            raise ValueError("confirmation sequence differs from authoritative source")


def cross_split_audit(training, confirmation):
    """Label-free gene, forward/RC Hamming2 and shifted17 screening."""
    import numpy as np

    from .types import reverse_complement

    if any(not r.get("gene_id") or len(r.get("core", "")) != 19 for r in training + confirmation):
        raise ValueError("canonical gene and paired19 required for independence audit")
    genes = {r["gene_id"] for r in training}
    cores = np.asarray([list(r["core"]) for r in training])
    words = {
        s[i : i + 17]
        for r in training
        for s in (r["core"], reverse_complement(r["core"]))
        for i in range(3)
    }
    findings = []
    for row in confirmation:
        reasons = []
        variants = (row["core"], reverse_complement(row["core"]))
        if row["gene_id"] in genes:
            reasons.append("same_gene")
        if len(training) and any(
            np.any(np.sum(cores != np.asarray(list(s)), axis=1) <= 2) for s in variants
        ):
            reasons.append("Hamming<=2_forward_or_RC")
        if any(s[i : i + 17] in words for s in variants for i in range(3)):
            reasons.append("shared17mer")
        if reasons:
            findings.append(dict(record_id=row["record_id"], reasons=reasons))
    return dict(
        all_clear=not findings,
        findings=findings,
        training_records=len(training),
        confirmation_records=len(confirmation),
    )
