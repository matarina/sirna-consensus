"""Post-fit Huesken diagnostics only: analytic random ranking, no model selection."""

from __future__ import annotations

import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

PKG = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PKG / "src"))
from sirna_consensus.confirmation import dump  # noqa: E402
from sirna_consensus.fasta import sha256_file  # noqa: E402


def main():
    model = PKG / "models/huesken_only_20261106"
    destination = PKG / "benchmarks/results/confirmation_20261106/development_diagnostics.json"
    if destination.exists():
        raise ValueError("preserve first diagnostic")
    data = json.loads((model / "development_records.json").read_text())
    units = defaultdict(list)
    for r in data:
        if r["study"] != "huesken2005":
            raise ValueError("diagnostics allow Huesken only")
        units[r["unit"]].append(r)
    rows = []
    for unit, records in units.items():
        n = len(records)
        active = sum(r["value"] >= r["active_threshold"] for r in records)
        gains = [2 ** r["target_rank"] - 1 for r in records]
        discounts = [1 / math.log2(i + 2) for i in range(5)]
        rows.append(
            dict(
                unit=unit,
                gene_id=records[0]["gene_id"],
                random_precision5=active / n,
                random_top3_at_least_one_active=1
                - (math.comb(n - active, 3) / math.comb(n, 3) if n - active >= 3 else 0),
                random_ndcg5=mean(gains)
                * sum(discounts)
                / sum(g * d for g, d in zip(sorted(gains, reverse=True)[:5], discounts)),
            )
        )
    macro = {}
    for key in ["random_precision5", "random_top3_at_least_one_active", "random_ndcg5"]:
        genes = defaultdict(list)
        for r in rows:
            genes[r["gene_id"]].append(r[key])
        macro[key] = mean(mean(v) for v in genes.values())
    dump(
        destination,
        dict(
            kind="post_fit_development_diagnostic_only",
            used_for_selection=False,
            independent_validation=False,
            source_sha256=sha256_file(model / "development_records.json"),
            formula="uniform random permutation without replacement; analytic expectation; gene/assay macro average",
            macro=macro,
            per_gene_assay=rows,
        ),
    )
    print(json.dumps(macro))


if __name__ == "__main__":
    main()
