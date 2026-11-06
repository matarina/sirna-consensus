"""Huesken-only nested model selection. Never promotes a default from development."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import pickle
import sys
from collections import defaultdict
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "sirna_consensus"

from .confirmation import assay_metrics, dump, ensure_training_records, load_protocol, summarize
from .fasta import sha256_file
from .learning import _hgnc, folds, group_records, weights
from .learning import features as original_features
from .types import reverse_complement, to_dna

FIT_AUDIT = []
FIT_CONTEXT = "test"


def features(guides, version):
    import numpy as np

    base = original_features(guides)
    if version == "sequence_v1":
        return base
    if version != "sequence_v2":
        raise ValueError("unknown preregistered feature group")
    dinucs = ["".join(x) for x in itertools.product("ACGT", repeat=2)]
    trinucs = ["".join(x) for x in itertools.product("ACGT", repeat=3)]
    extra = []
    for guide in guides:
        g = to_dna(guide)[:19]
        extra.append(
            [float(g[i : i + 2] == mer) for i in range(18) for mer in dinucs]
            + [sum(g[i : i + 3] == mer for i in range(17)) / 17 for mer in trinucs]
        )
    return np.column_stack([base, np.asarray(extra)])


def load_development(protocol):
    from scipy.stats import rankdata

    source = Path(protocol["training"]["source"])
    if sha256_file(source) != protocol["training"]["source_sha256"]:
        raise ValueError("development source checksum mismatch")
    approved, aliases = _hgnc(protocol["training"]["hgnc"])
    case_symbols = defaultdict(set)
    for symbol in approved:
        case_symbols[symbol.casefold()].add(symbol)
    records = []
    with source.open() as h:
        for row in csv.DictReader(h, delimiter="\t"):
            gene = row["gene"]
            if gene not in approved:
                options = aliases.get(gene, set()) or case_symbols.get(gene.casefold(), set())
                if len(options) != 1:
                    raise ValueError("unresolved development gene: " + gene)
                gene = next(iter(options))
            guide = to_dna(row["guide_dna"])
            if len(guide) not in (19, 21) or set(guide) - set("ACGT"):
                raise ValueError("invalid paired guide")
            threshold = protocol["active_thresholds"].get(row["endpoint"])
            if threshold is None:
                raise ValueError("unregistered endpoint")
            records.append(
                dict(
                    record_id=row["record_id"],
                    study=row["study"],
                    gene=gene,
                    gene_original=row["gene"],
                    gene_id=approved[gene],
                    condition=row["assay"],
                    unit="|".join((row["study"], approved[gene], row["assay"])),
                    endpoint=row["endpoint"],
                    value=float(row["value"]),
                    active_threshold=threshold,
                    guide=guide,
                    core=reverse_complement(guide[:19]),
                    passenger=to_dna(row.get("passenger", "")),
                    chemistry="unknown",
                    transcript_id=row["transcript_id"],
                    target_sha256=row["target_sha256"],
                    mapping_status=row["mapping_status"],
                    source_record=row["source_record"],
                )
            )
    records.sort(key=lambda r: (r["gene_id"], r["unit"], r["record_id"]))
    ensure_training_records(records, protocol)
    units = defaultdict(list)
    for r in records:
        units[r["unit"]].append(r)
    for rows in units.values():
        ranks = rankdata([r["value"] for r in rows])
        for r, rank in zip(rows, ranks):
            r["target_rank"] = float((rank - 1) / (len(rows) - 1)) if len(rows) > 1 else 0.5
    return records


def fit(params, x, records, indices, protocol):
    import numpy as np
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.preprocessing import StandardScaler

    rows = [records[int(i)] for i in indices]
    ensure_training_records(rows, protocol)
    FIT_AUDIT.append(
        dict(
            context=FIT_CONTEXT,
            params=params,
            studies=sorted({r["study"] for r in rows}),
            train_ids=[r["record_id"] for r in rows],
        )
    )
    w = weights(rows)
    values = np.asarray([r["target_rank"] for r in rows])
    matrix = x[params["feature"]][indices]
    scaler = None
    if params["family"] in ("ridge", "pairwise"):
        scaler = StandardScaler().fit(matrix, sample_weight=w)
        matrix = scaler.transform(matrix)
    if params["family"] == "ridge":
        model = Ridge(alpha=params["alpha"]).fit(matrix, values, sample_weight=w)
    elif params["family"] == "tree":
        model = HistGradientBoostingRegressor(
            max_depth=params["max_depth"],
            min_samples_leaf=params["min_samples_leaf"],
            l2_regularization=params["l2_regularization"],
            learning_rate=0.05,
            max_iter=200,
            early_stopping=False,
            random_state=protocol["seed"],
        )
        model.fit(matrix, values, sample_weight=w)
    elif params["family"] == "pairwise":
        units = defaultdict(list)
        for i, r in enumerate(rows):
            units[r["unit"]].append(i)
        left, right, units_for_pairs = [], [], []
        for unit, ids in sorted(units.items()):
            for i, j in itertools.combinations(ids, 2):
                if values[i] != values[j]:
                    lo, hi = (i, j) if values[i] < values[j] else (j, i)
                    left.extend((hi, lo))
                    right.extend((lo, hi))
                    units_for_pairs.extend((unit, unit))
        if not left:
            raise ValueError("no informative within-assay training pairs")
        pair_counts = {u: units_for_pairs.count(u) for u in units}
        genes = {r["gene_id"] for r in rows}
        conditions = defaultdict(set)
        for r in rows:
            if pair_counts[r["unit"]]:
                conditions[r["gene_id"]].add(r["unit"])
        gene_for_unit = {r["unit"]: r["gene_id"] for r in rows}
        pw = np.asarray(
            [
                1 / (len(genes) * len(conditions[gene_for_unit[u]]) * pair_counts[u])
                for u in units_for_pairs
            ]
        )
        pw /= pw.mean()
        diff = matrix[np.asarray(left)] - matrix[np.asarray(right)]
        model = LogisticRegression(C=params["C"], fit_intercept=False, max_iter=2000, tol=1e-5)
        model.fit(diff, np.tile([1, 0], len(left) // 2), sample_weight=pw)
        if int(model.n_iter_.max()) >= model.max_iter:
            raise RuntimeError("pairwise model did not converge; do not silently score")
    else:
        raise ValueError("unregistered model family")
    return dict(params=params, scaler=scaler, model=model)


def predict_matrix(bundle, matrix):
    if bundle["scaler"] is not None:
        matrix = bundle["scaler"].transform(matrix)
    model = bundle["model"]
    return (
        model.decision_function(matrix)
        if bundle["params"]["family"] == "pairwise"
        else model.predict(matrix)
    )


def predict(bundle, guides):
    return predict_matrix(bundle, features(guides, bundle["params"]["feature"]))


def choose(indices, x, records, groups, protocol):
    import numpy as np

    inner = folds(indices, groups, protocol["training"]["inner_folds"])
    trials = []
    for params in protocol["training"]["grid"]:
        prediction = np.full(len(records), np.nan)
        for train, test in inner:
            model = fit(params, x, records, train, protocol)
            prediction[test] = predict_matrix(model, x[params["feature"]][test])
        value = summarize(assay_metrics([records[int(i)] for i in indices], prediction[indices]))[
            "ndcg5"
        ]
        trials.append(dict(params=params, macro_gene_assay_NDCG5=value))
    if any(t["macro_gene_assay_NDCG5"] is None for t in trials):
        raise ValueError("undefined development selection objective")
    winner = max(range(len(trials)), key=lambda i: trials[i]["macro_gene_assay_NDCG5"])
    return trials[winner]["params"], trials


def develop(protocol_path, out):
    global FIT_CONTEXT
    import numpy as np
    import scipy
    import sklearn
    from threadpoolctl import threadpool_limits

    out = Path(out).resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("candidate output must be empty; preserve prior attempts")
    out.mkdir(parents=True, exist_ok=True)
    protocol, _, _, registration = load_protocol(protocol_path)
    records = load_development(protocol)
    groups, edges = group_records(records)
    x = {v: features([r["guide"] for r in records], v) for v in ("sequence_v1", "sequence_v2")}
    dump(out / "development_records.json", records)
    audit, prediction = [], np.full(len(records), np.nan)
    FIT_AUDIT.clear()
    with threadpool_limits(limits=protocol["training"]["threads"]):
        for n, (tr, te) in enumerate(
            folds(np.arange(len(records)), groups, protocol["training"]["outer_folds"])
        ):
            FIT_CONTEXT = f"outer_{n}"
            params, trials = choose(tr, x, records, groups, protocol)
            model = fit(params, x, records, tr, protocol)
            prediction[te] = predict_matrix(model, x[params["feature"]][te])
            checkpoint = out / f"fold_{n}.pkl"
            checkpoint.write_bytes(pickle.dumps(model, protocol=4))
            audit.append(
                dict(
                    fold=n,
                    train_ids=[records[int(i)]["record_id"] for i in tr],
                    test_ids=[records[int(i)]["record_id"] for i in te],
                    selected_params=params,
                    trials=trials,
                    train_weights=weights([records[int(i)] for i in tr]).tolist(),
                    checkpoint=checkpoint.name,
                    checkpoint_sha256=sha256_file(checkpoint),
                )
            )
            print(
                json.dumps(dict(outer_fold=n, selected=params, status="development_diagnostic")),
                flush=True,
            )
            dump(out / "progress.json", audit)
        FIT_CONTEXT = "final_development_selection"
        params, trials = choose(np.arange(len(records)), x, records, groups, protocol)
        model = fit(params, x, records, np.arange(len(records)), protocol)
    checkpoint = out / "model.pkl"
    checkpoint.write_bytes(pickle.dumps(model, protocol=4))
    oof = [
        dict(record_id=r["record_id"], score=float(v), group=int(g))
        for r, v, g in zip(records, prediction, groups)
    ]
    dump(out / "oof_predictions.json", oof)
    dump(out / "fit_audit.json", FIT_AUDIT)
    dump(
        out / "split_manifest.json",
        dict(
            folds=audit,
            groups={r["record_id"]: int(g) for r, g in zip(records, groups)},
            cross_gene_sequence_edges=edges,
        ),
    )
    result = dict(
        status="development_diagnostics_only",
        independent_validation_passed=False,
        promoted=False,
        protocol_sha256=sha256_file(Path(protocol_path)),
        registration_commit=registration,
        development_source_sha256=protocol["training"]["source_sha256"],
        records=len(records),
        genes=len({r["gene_id"] for r in records}),
        components=len(set(groups)),
        primary_metric=protocol["primary_metric"],
        metrics=summarize(assay_metrics(records, prediction)),
        selected_params=params,
        final_trials=trials,
        checkpoint_sha256=sha256_file(checkpoint),
        environment=dict(
            numpy=np.__version__, scipy=scipy.__version__, sklearn=sklearn.__version__
        ),
        loso=protocol["loso"],
        training_sources=[protocol["training"]["source"]],
        warnings=["Development scores cannot authorize a default replacement."],
    )
    dump(out / "results.json", result)
    source = Path(__file__).resolve()
    artifacts = {
        str(p): sha256_file(p)
        for p in (source, source.with_name("confirmation.py"), source.with_name("learning.py"))
    }
    artifacts["model.pkl"] = sha256_file(checkpoint)
    dump(
        out / "scorer_manifest.json",
        dict(
            name="huesken_only_challenger",
            kind="guarded_sequence_ranker",
            input_convention="paired_antisense_19",
            checkpoint_sha256=sha256_file(checkpoint),
            artifacts=artifacts,
            checkpoint="model.pkl",
            command=[
                sys.executable,
                str(source),
                "--predict",
                "{input}",
                "{output}",
                str(out / "scorer_manifest.json"),
            ],
            supported_formats=[
                "binding19-guide21-overhang2-passenger21-overhang2",
                "binding19-guide21-overhang2-passenger_unknown",
            ],
            training_audit="Huesken_only; confirmation never used for fitting or selection",
            protocol_sha256=sha256_file(Path(protocol_path)),
            promoted=False,
            independent_validation_passed=False,
        ),
    )
    print(
        json.dumps({k: result[k] for k in ("status", "records", "genes", "selected_params")}),
        flush=True,
    )
    return result


def predict_request(input_path, output_path, manifest_path):
    from threadpoolctl import threadpool_limits

    source, destination, manifest_path = Path(input_path), Path(output_path), Path(manifest_path)
    request, manifest = json.loads(source.read_text()), json.loads(manifest_path.read_text())
    if request["scorer_manifest_sha256"] != sha256_file(manifest_path):
        raise ValueError("request manifest mismatch")
    for name, sha in manifest["artifacts"].items():
        path = Path(name) if Path(name).is_absolute() else manifest_path.parent / name
        if sha256_file(path) != sha:
            raise ValueError("candidate scorer artifact mismatch")
    keys = [r["key"] for r in request["candidates"]]
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate scorer keys")
    bundle = pickle.loads((manifest_path.parent / manifest["checkpoint"]).read_bytes())
    with threadpool_limits(limits=4):
        values = predict(bundle, [r["paired_guide"] for r in request["candidates"]]) if keys else []
    dump(
        destination,
        dict(input_sha256=sha256_file(source), scores={k: float(v) for k, v in zip(keys, values)}),
    )


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--predict":
        predict_request(*sys.argv[2:])
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--protocol", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    develop(a.protocol, a.out)


if __name__ == "__main__":
    main()
