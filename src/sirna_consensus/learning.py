"""Auditable CPU rank models with nested gene/sequence-grouped validation.

ML dependencies are optional. Nothing from a held-out fold fits preprocessing,
weights, hyperparameters, or model-family selection.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import pickle
import platform
import sys
from collections import Counter, defaultdict
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "sirna_consensus"

from .fasta import sha256_file
from .rank import _quality_term
from .types import Core, reverse_complement, to_dna

SEED = 20261106
BASES = "ACGT"
FEATURE_VERSION = "paired19_sequence_v1"
PARAMS = [{"family": "ridge", "alpha": a} for a in (0.01, 0.1, 1.0, 10.0, 100.0)] + [
    {"family": "tree", "max_depth": d, "min_samples_leaf": leaf, "l2_regularization": l2}
    for d, leaf, l2 in itertools.product((2, 3), (5, 10), (0.1, 1.0))
]


def features(guides):
    import numpy as np

    rows = []
    for sequence in guides:
        g = to_dna(sequence)[:19]
        if len(g) != 19 or set(g) - set(BASES):
            raise ValueError("feature input must be a documented paired antisense 19mer")

        def gc(s):
            return sum(b in "GC" for b in s) / len(s)

        rows.append(
            [float(b == base) for b in g for base in BASES]
            + [g.count(base) / 19 for base in BASES]
            + [sum(g[i : i + 2] == a + b for i in range(18)) / 18 for a in BASES for b in BASES]
            + [
                gc(g),
                gc(g[1:8]),
                gc(g[:4]) - gc(g[-4:]),
                max(len(list(v)) for _, v in itertools.groupby(g)) / 19,
            ]
        )
    return np.asarray(rows, dtype=float)


def _write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def _hgnc(path):
    direct, aliases = {}, defaultdict(set)
    with Path(path).open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("status") != "Approved":
                continue
            symbol = row["symbol"]
            direct[symbol] = row["hgnc_id"]
            for field in ("alias_symbol", "prev_symbol"):
                for alias in row.get(field, "").split("|"):
                    if alias:
                        aliases[alias].add(symbol)
    return direct, aliases


def prepare_data(project, hgnc_path, out):
    """Use existing pinned curated artifacts; disclose unavailable upstream lineage."""
    project, out = Path(project).resolve(), Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    direct, aliases = _hgnc(hgnc_path)
    case_symbols = defaultdict(set)
    for symbol in direct:
        case_symbols[symbol.casefold()].add(symbol)
    artifacts, records, quarantine = {}, [], []
    sources = [
        (project / "sirna-consensus/benchmarks/data/huesken/measurements.tsv", "huesken2005"),
        (
            project / "independent_validation/sciabola/canonical/measurements_primary.tsv",
            "sciabola2013",
        ),
    ]
    for path, study in sources:
        artifacts[str(path)] = sha256_file(path)
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                guide = to_dna(row["guide_dna"])
                records.append(
                    {
                        "record_id": row["record_id"],
                        "study": study,
                        "gene_original": row["gene"],
                        "condition": row["assay"],
                        "endpoint": row["endpoint"],
                        "value": float(row["value"]),
                        "guide": guide,
                        "core": reverse_complement(guide[:19]),
                        "passenger": to_dna(row.get("passenger", "")),
                        "chemistry": "unknown",
                        "transcript_id": row["transcript_id"],
                        "mapping_status": row["mapping_status"],
                        "active_threshold": 0.7 if study == "huesken2005" else 70.0,
                        "source_record": row.get("source_record", row["record_id"]),
                    }
                )
    raw = (
        project
        / "independent_validation/denovo_two_benchmarks_20261105/outputs/all_source_rows.json"
    )
    artifacts[str(raw)] = sha256_file(raw)
    for row in json.loads(raw.read_text()):
        if row.get("mapping_status") != "unique":
            quarantine.append(
                {
                    "study": row["study"],
                    "record_id": row["record_id"],
                    "reason": "not_uniquely_mapped",
                }
            )
            continue
        records.append(
            {
                "record_id": row["study"] + ":" + row["record_id"],
                "study": row["study"],
                "gene_original": row["gene"],
                "condition": f"{row['study']}_primary_{row['dose_nM']}nM",
                "endpoint": "percent_knockdown",
                "value": float(row["knockdown"]),
                "guide": row["guide"],
                "core": row["core"],
                "passenger": row["passenger"],
                "chemistry": "RNA_dTdT_passenger" if row["study"] == "Wang2009" else "unknown",
                "transcript_id": row["transcript_id"],
                "start": row["mapped_positions"][0],
                "mapping_status": "unique",
                "active_threshold": 70.0,
                "source_record": row["record_id"],
            }
        )
    artifacts[str(Path(hgnc_path).resolve())] = sha256_file(Path(hgnc_path))
    by_key = {}
    for row in records:
        gene = row["gene_original"]
        row["gene_resolution"] = "approved_symbol"
        if gene not in direct:
            possible = aliases.get(gene, set())
            row["gene_resolution"] = "HGNC_previous_or_alias_symbol"
            if not possible:
                possible = case_symbols.get(gene.casefold(), set())
                row["gene_resolution"] = "HGNC_unique_case_variant"
            if len(possible) != 1:
                quarantine.append(
                    {
                        "study": row["study"],
                        "record_id": row["record_id"],
                        "reason": "gene_alias_absent_or_ambiguous",
                        "gene": gene,
                    }
                )
                continue
            gene = next(iter(possible))
        row["gene"], row["gene_id"] = gene, direct[gene]
        if len(row["core"]) != 19 or reverse_complement(row["guide"][:19]) != row["core"]:
            raise ValueError("invalid paired geometry: " + row["record_id"])
        row["chemistry_match_known"] = row["chemistry"] != "unknown"
        row["replicate_ids"] = [row["record_id"]]
        row["raw_values"] = [row["value"]]
        key = tuple(
            row.get(k, "")
            for k in ("study", "gene_id", "condition", "guide", "passenger", "chemistry")
        )
        if key in by_key:
            by_key[key]["replicate_ids"].append(row["record_id"])
            by_key[key]["raw_values"].append(row["value"])
        else:
            by_key[key] = row
    records = sorted(by_key.values(), key=lambda r: (r["study"], r["gene_id"], r["record_id"]))
    for row in records:
        row["value"] = sum(row["raw_values"]) / len(row["raw_values"])
        row["development_exposed"] = True
        row["independent_holdout"] = False
        row["unit"] = "|".join((row["study"], row["gene_id"], row["condition"]))
    from scipy.stats import rankdata

    units = defaultdict(list)
    for i, r in enumerate(records):
        units[r["unit"]].append(i)
    for indices in units.values():
        ranks = rankdata([records[i]["value"] for i in indices])
        for i, rank in zip(indices, ranks):
            records[i]["target_rank"] = (
                float((rank - 1) / (len(indices) - 1)) if len(indices) > 1 else 0.5
            )
    manifest = {
        "schema_version": 1,
        "status": "development_only",
        "sources": artifacts,
        "records": len(records),
        "raw_records": len(records)
        + sum(len(r["replicate_ids"]) - 1 for r in records)
        + len(quarantine),
        "studies": dict(Counter(r["study"] for r in records)),
        "genes": len({r["gene_id"] for r in records}),
        "quarantine": quarantine,
        "endpoint_direction": "higher_is_more_active",
        "active_thresholds": {
            "huesken2005": 0.7,
            "sciabola2013": 70.0,
            "Wang2009": 70.0,
            "Filhol2012": 70.0,
        },
        "endpoint_policy": "primary endpoint only; no cross-condition numeric pooling; technical duplicates averaged",
        "provenance_limitations": [
            "Huesken original upstream workspace is unavailable; current curated measurements are pinned.",
            "Filhol chemical identity is undocumented; Wang passenger dTdT is not all RNA.",
            "All four datasets have been inspected previously and cannot be a sealed final test.",
        ],
    }
    _write(out / "data_registry.json", manifest)
    _write(out / "development_records.json", records)
    return records, manifest


def group_records(records):
    import numpy as np

    parent = list(range(len(records)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)

    genes, mers, edges = {}, {}, []
    for i, r in enumerate(records):
        if r["gene_id"] in genes:
            union(i, genes[r["gene_id"]])
        else:
            genes[r["gene_id"]] = i
        for s in (r["core"], reverse_complement(r["core"])):
            for start in range(3):
                mer = s[start : start + 17]
                if mer in mers:
                    union(i, mers[mer])
                else:
                    mers[mer] = i
    codes = np.asarray([[BASES.index(b) for b in r["core"]] for r in records])
    rc = 3 - codes[:, ::-1]
    for i in range(len(records)):
        near = (
            np.where(
                np.minimum(
                    (codes[i + 1 :] != codes[i]).sum(axis=1), (rc[i + 1 :] != codes[i]).sum(axis=1)
                )
                <= 2
            )[0]
            + i
            + 1
        )
        for j in near:
            if records[i]["gene_id"] != records[int(j)]["gene_id"]:
                edges.append(
                    {
                        "a": records[i]["record_id"],
                        "b": records[int(j)]["record_id"],
                        "reason": "Hamming<=2",
                    }
                )
            union(i, int(j))
    groups = np.asarray([find(i) for i in range(len(records))])
    return groups, edges


def folds(indices, groups, n):
    import numpy as np

    members = defaultdict(list)
    for i in indices:
        members[int(groups[i])].append(int(i))
    n = min(n, len(members))
    if n < 2:
        raise ValueError("too few independent groups for validation")
    rng = np.random.default_rng(SEED)
    keys = list(members)
    rng.shuffle(keys)
    keys.sort(key=lambda k: -len(members[k]))
    buckets = [[] for _ in range(n)]
    for key in keys:
        bucket = min(range(n), key=lambda j: (len(buckets[j]), j))
        buckets[bucket].extend(members[key])
    all_ids = set(map(int, indices))
    return [
        (np.asarray(sorted(all_ids - set(b)), dtype=int), np.asarray(sorted(b), dtype=int))
        for b in buckets
    ]


def weights(records):
    import numpy as np

    studies = {r["study"] for r in records}
    genes = {s: {r["gene_id"] for r in records if r["study"] == s} for s in studies}
    conditions = defaultdict(set)
    counts = Counter(r["unit"] for r in records)
    for r in records:
        conditions[(r["study"], r["gene_id"])].add(r["condition"])
    w = np.asarray(
        [
            1
            / (
                len(studies)
                * len(genes[r["study"]])
                * len(conditions[(r["study"], r["gene_id"])])
                * counts[r["unit"]]
            )
            for r in records
        ]
    )
    # Mean-one weights keep ridge alpha interpretable across training fold sizes.
    return w / w.mean()


def estimator(params):
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import Ridge
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    if params["family"] == "ridge":
        return Pipeline([("scale", StandardScaler()), ("model", Ridge(alpha=params["alpha"]))])
    return HistGradientBoostingRegressor(
        **{k: v for k, v in params.items() if k != "family"},
        learning_rate=0.05,
        max_iter=200,
        early_stopping=False,
        random_state=SEED,
    )


def fit_model(params, x, y, train, records):
    model = estimator(params)
    w = weights([records[int(i)] for i in train])
    if params["family"] == "ridge":
        # Weighted preprocessing is fit solely on this training fold.
        model.fit(x[train], y[train], scale__sample_weight=w, model__sample_weight=w)
    else:
        model.fit(x[train], y[train], sample_weight=w)
    return model


def unit_metrics(records, predictions):
    import numpy as np
    from scipy.stats import spearmanr

    units = defaultdict(list)
    for i, r in enumerate(records):
        if math.isfinite(float(predictions[i])):
            units[r["unit"]].append(i)
    result = []
    for unit, ids in sorted(units.items()):
        ids = sorted(
            ids,
            key=lambda i: (-float(predictions[i]), records[i]["guide"], records[i]["record_id"]),
        )
        top, k = ids[:5], min(5, len(ids))
        ranks = np.asarray([records[i]["target_rank"] for i in ids])
        discounts = 1 / np.log2(np.arange(k) + 2)
        ideal = float(((2 ** np.sort(ranks)[::-1][:k] - 1) * discounts).sum())
        ndcg = float(((2 ** ranks[:k] - 1) * discounts).sum() / ideal) if ideal > 0 else None
        raw = np.asarray([records[i]["value"] for i in ids])
        pred = np.asarray([predictions[i] for i in ids])
        rho = (
            float(spearmanr(raw, pred).statistic)
            if len(ids) > 1 and np.ptp(raw) > 0 and np.ptp(pred) > 0
            else None
        )
        row = records[ids[0]]
        result.append(
            {
                "unit": unit,
                "study": row["study"],
                "gene": row["gene"],
                "gene_id": row["gene_id"],
                "n": len(ids),
                "ndcg5": ndcg,
                "spearman": rho,
                "precision5": sum(
                    records[i]["value"] >= records[i]["active_threshold"] for i in top
                )
                / k,
                "top_ids": [records[i]["record_id"] for i in top],
            }
        )
    return result


def macro(rows, field="ndcg5"):
    values = [r[field] for r in rows if r[field] is not None]
    return sum(values) / len(values) if values else float("-inf")


def tune(indices, x, y, groups, records):
    import numpy as np

    splits = folds(indices, groups, 3)
    trials = []
    for params in PARAMS:
        prediction = np.full(len(records), np.nan)
        for train, test in splits:
            prediction[test] = fit_model(params, x, y, train, records).predict(x[test])
        score = macro(unit_metrics(records, prediction))
        trials.append({"params": params, "ndcg5": score})
    # PARAMS order is the declared deterministic tie break.
    best = max(range(len(trials)), key=lambda i: trials[i]["ndcg5"])
    best_family = {
        family: max(
            (t for t in trials if t["params"]["family"] == family), key=lambda t: t["ndcg5"]
        )["params"]
        for family in ("ridge", "tree")
    }
    return trials[best]["params"], best_family, trials


def confidence(records, groups, a, b):
    import numpy as np

    ma, mb = (
        {r["unit"]: r for r in unit_metrics(records, a)},
        {r["unit"]: r for r in unit_metrics(records, b)},
    )
    unit_group = {r["unit"]: int(g) for r, g in zip(records, groups)}
    deltas = defaultdict(list)
    for unit, r in ma.items():
        if r["ndcg5"] is not None and mb[unit]["ndcg5"] is not None:
            deltas[unit_group[unit]].append(r["ndcg5"] - mb[unit]["ndcg5"])
    keys = sorted(deltas)
    if not keys:
        return {"delta": None, "ci95": [None, None], "groups": 0}
    rng = np.random.default_rng(SEED)
    boot = [
        np.mean([v for key in rng.choice(keys, len(keys), replace=True) for v in deltas[int(key)]])
        for _ in range(5000)
    ]
    return {
        "delta": float(np.mean([v for vs in deltas.values() for v in vs])),
        "ci95": [float(v) for v in np.quantile(boot, [0.025, 0.975])],
        "groups": len(keys),
    }


def train(project, hgnc, out):
    import numpy as np
    import scipy
    import sklearn
    from threadpoolctl import threadpool_limits

    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    records, _registry = prepare_data(project, hgnc, out)
    groups, edges = group_records(records)
    if len(set(groups)) < 3:
        raise ValueError("fewer than three independent groups; no model will be trained")
    x = features([r["guide"] for r in records])
    y = np.asarray([r["target_rank"] for r in records])
    predictions = {name: np.full(len(records), np.nan) for name in ("selected", "ridge", "tree")}
    predictions["quality"] = np.asarray(
        [
            _quality_term(
                Core(
                    r["transcript_id"],
                    1,
                    19,
                    r["core"],
                    r["guide"],
                    r["passenger"],
                    r["guide"][19:],
                    r["passenger"][19:] if r["passenger"] else "",
                    r["chemistry"],
                )
            )
            for r in records
        ]
    )
    predictions["random"] = np.random.default_rng(SEED).random(len(records))
    fold_audit = []
    with threadpool_limits(limits=4):
        for fold, (tr, te) in enumerate(folds(np.arange(len(records)), groups, 5)):
            best, by_family, trials = tune(tr, x, y, groups, records)
            audit = {
                "fold": fold,
                "train_ids": [records[int(i)]["record_id"] for i in tr],
                "test_ids": [records[int(i)]["record_id"] for i in te],
                "selected_params": best,
                "family_params": by_family,
                "train_weights": weights([records[int(i)] for i in tr]).tolist(),
                "trials": trials,
            }
            for family, params in by_family.items():
                model = fit_model(params, x, y, tr, records)
                predictions[family][te] = model.predict(x[te])
                checkpoint = out / f"fold_{fold}_{family}.pkl"
                checkpoint.write_bytes(pickle.dumps(model, protocol=4))
                audit.setdefault("checkpoints", {})[family] = {
                    "path": checkpoint.name,
                    "sha256": sha256_file(checkpoint),
                }
            predictions["selected"][te] = predictions[best["family"]][te]
            fold_audit.append(audit)
            print(
                json.dumps({"outer_fold": fold, "selected": best, "test_records": len(te)}),
                flush=True,
            )
        loso = []
        for study in sorted({r["study"] for r in records}):
            te = np.asarray([i for i, r in enumerate(records) if r["study"] == study])
            held_groups = set(groups[te])
            tr = np.asarray(
                [i for i, r in enumerate(records) if groups[i] not in held_groups], dtype=int
            )
            item = {
                "study": study,
                "test_records": len(te),
                "train_records": len(tr),
                "purged_records": len(records) - len(te) - len(tr),
                "train_ids": [records[int(i)]["record_id"] for i in tr],
                "test_ids": [records[int(i)]["record_id"] for i in te],
            }
            if len(set(groups[tr])) < 3:
                item["status"] = "insufficient_independent_training_groups"
            else:
                best, _, trials = tune(tr, x, y, groups, records)
                model = fit_model(best, x, y, tr, records)
                p = np.full(len(records), np.nan)
                p[te] = model.predict(x[te])
                item.update(
                    status="complete",
                    params=best,
                    metrics=unit_metrics(records, p),
                    quality_metrics=unit_metrics(
                        records, np.where(np.isfinite(p), predictions["quality"], np.nan)
                    ),
                    predictions={records[int(i)]["record_id"]: float(p[i]) for i in te},
                )
                path = out / f"loso_{study}.pkl"
                path.write_bytes(pickle.dumps(model, protocol=4))
                item["checkpoint"] = {"path": path.name, "sha256": sha256_file(path)}
            loso.append(item)
            print(json.dumps({"loso": study, "status": item["status"]}), flush=True)
        best, _, trials = tune(np.arange(len(records)), x, y, groups, records)
        model = fit_model(best, x, y, np.arange(len(records)), records)
    checkpoint = out / "model.pkl"
    checkpoint.write_bytes(pickle.dumps(model, protocol=4))
    improvement = confidence(records, groups, predictions["selected"], predictions["quality"])
    promoted = improvement["ci95"][0] is not None and improvement["ci95"][0] > 0
    oof = [
        {
            "record_id": r["record_id"],
            "group": int(groups[i]),
            "fold": next(a["fold"] for a in fold_audit if r["record_id"] in a["test_ids"]),
            **{name: float(values[i]) for name, values in predictions.items()},
        }
        for i, r in enumerate(records)
    ]
    _write(out / "oof_predictions.json", oof)
    _write(
        out / "split_manifest.json",
        {
            "seed": SEED,
            "group_rule": "same HGNC gene OR core forward/RC Hamming<=2 OR shared 17mer; transitive components",
            "groups": len(set(groups)),
            "records": {r["record_id"]: int(g) for r, g in zip(records, groups)},
            "cross_gene_near_edges": edges,
            "folds": fold_audit,
        },
    )
    results = {
        "schema_version": 1,
        "status": "development_only",
        "independent_validation_passed": False,
        "data_registry_sha256": sha256_file(out / "data_registry.json"),
        "primary_metric": "macro_gene_assay_NDCG@5",
        "metrics": {
            name: {
                "ndcg5": macro(m := unit_metrics(records, p)),
                "spearman": macro(m, "spearman"),
                "precision5": macro(m, "precision5"),
                "per_gene_assay": m,
            }
            for name, p in predictions.items()
        },
        "selected_vs_quality": improvement,
        "promoted": promoted,
        "final_params": best,
        "final_trials": trials,
        "loso": loso,
        "checkpoint_sha256": sha256_file(checkpoint),
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "feature_version": FEATURE_VERSION,
        "feature_count": int(x.shape[1]),
        "seed": SEED,
    }
    _write(out / "results.json", results)
    source = Path(__file__).resolve()
    manifest = {
        "name": "grouped_sequence_ranker_v1",
        "kind": "sklearn_ranker",
        "input_convention": "paired_antisense_19",
        "feature_version": FEATURE_VERSION,
        "checkpoint": "model.pkl",
        "checkpoint_sha256": sha256_file(checkpoint),
        "artifacts": {"model.pkl": sha256_file(checkpoint), str(source): sha256_file(source)},
        "command": [
            sys.executable,
            str(source),
            "--predict",
            "{input}",
            "{output}",
            str(out / "scorer_manifest.json"),
        ],
        "supported_formats": [
            "binding19-guide21-overhang2-passenger21-overhang2",
            "binding19-guide21-overhang2-passenger_unknown",
        ],
        "training_audit": "all_four_studies_exposed_development_only",
        "promoted": promoted,
        "independent_validation_passed": False,
        "target": "within_assay_gene_rank; not knockdown probability",
    }
    _write(out / "scorer_manifest.json", manifest)
    _write(
        out / "default_policy.json",
        {
            "promoted": promoted,
            "scorer_manifest": "scorer_manifest.json",
            "scorer_manifest_sha256": sha256_file(out / "scorer_manifest.json"),
            "gate": improvement,
            "default": "model" if promoted else "quality",
            "status": "development_only",
        },
    )
    verify = verify_oof(out)
    _write(out / "verification.json", verify)
    return results


def verify_oof(out):
    import numpy as np
    from threadpoolctl import threadpool_limits

    out = Path(out)
    records = json.loads((out / "development_records.json").read_text())
    splits = json.loads((out / "split_manifest.json").read_text())
    oof = {r["record_id"]: r for r in json.loads((out / "oof_predictions.json").read_text())}
    ids = {r["record_id"]: i for i, r in enumerate(records)}
    x = features([r["guide"] for r in records])
    covered = []
    errors, leaks = [], []
    with threadpool_limits(limits=4):
        for fold in splits["folds"]:
            tr, te = fold["train_ids"], fold["test_ids"]
            if {splits["records"][i] for i in tr} & {splits["records"][i] for i in te}:
                leaks.append(fold["fold"])
            covered.extend(te)
            for family, ck in fold["checkpoints"].items():
                path = out / ck["path"]
                if sha256_file(path) != ck["sha256"]:
                    raise ValueError("fold checkpoint hash mismatch")
                model = pickle.loads(path.read_bytes())
                actual = model.predict(x[[ids[i] for i in te]])
                errors.append(float(np.max(np.abs(actual - [oof[i][family] for i in te]))))
    checks = {
        "one_oof_per_record": len(covered) == len(set(covered)) == len(records)
        and set(covered) == set(ids),
        "no_group_overlap": not leaks,
        "checkpoint_prediction_max_error": max(errors),
        "checkpoints_reproduce_oof": max(errors) < 1e-12,
        "finite_predictions": all(
            math.isfinite(v)
            for r in oof.values()
            for k, v in r.items()
            if k in ("ridge", "tree", "selected")
        ),
    }
    checks["all_passed"] = (
        checks["one_oof_per_record"]
        and checks["no_group_overlap"]
        and checks["checkpoints_reproduce_oof"]
        and checks["finite_predictions"]
    )
    return checks


def predict_request(request_path, output_path, manifest_path):
    from threadpoolctl import threadpool_limits

    manifest_path = Path(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    request = json.loads(Path(request_path).read_text())
    checkpoint = manifest_path.parent / manifest["checkpoint"]
    if sha256_file(checkpoint) != manifest["checkpoint_sha256"]:
        raise ValueError("model checksum mismatch")
    model = pickle.loads(checkpoint.read_bytes())
    candidates = request["candidates"]
    with threadpool_limits(limits=4):
        p = model.predict(features([r["paired_guide"] for r in candidates])) if candidates else []
    _write(
        output_path,
        {
            "input_sha256": sha256_file(Path(request_path)),
            "scores": {r["key"]: float(v) for r, v in zip(candidates, p)},
        },
    )


def add_parser(sub):
    p = sub.add_parser("train", help="nested grouped CPU model training; development data only")
    p.add_argument("--project", required=True)
    p.add_argument("--hgnc", required=True)
    p.add_argument("--out", required=True)


def command(args):
    result = train(args.project, args.hgnc, args.out)
    print(
        json.dumps(
            {
                "output": args.out,
                "promoted": result["promoted"],
                "gain": result["selected_vs_quality"],
            }
        )
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) == 5 and sys.argv[1] == "--predict":
        # Direct script workers also need the package on sys.path in source installs.
        predict_request(sys.argv[2], sys.argv[3], sys.argv[4])
