# Weight-independent final evaluation

Reproduce with `python benchmarks/scripts/run_final.py` (the lock is reused when it still verifies, else recut).

- Audited cohort: 6 studies / 1113 records / 12 genes, criterion `weight`
- **Actually evaluated**: 399 records / 10 genes / 5 studies (714 excluded — see `data/holdout/exclusions.json`)
- native generator evidence attached to 210/356 policy-arm candidates (59%)
- checkpoint_training_audit: `verified`; independent_holdout_eligible: `True`
- upstream evaluation exposure: all are upstream-evaluation-exposed, so this is weight-independent, not fully independent

## Development-selected policy vs frozen scorer alone

Complete Top-5 panels on 9 shared genes. Paired mean gene-weighted precision difference **+0.0 pt**.

Within-unit label-permutation control: 17 evaluable contrasts, min one-sided p = 0.000999000999000999.

| gene | policy | model | difference |
|---|---:|---:|---:|
| Cyclophilin B | 1.00 | 1.00 | +0.00 |
| DBI | 0.20 | 0.20 | +0.00 |
| Firefly luciferase | 1.00 | 1.00 | +0.00 |
| FireflyLuc | 0.80 | 0.80 | +0.00 |
| GAPDH | 1.00 | 1.00 | +0.00 |
| ICAM-1 | 0.00 | 0.00 | +0.00 |
| MMAC1 | 0.60 | 0.60 | +0.00 |
| PLK | 0.40 | 0.40 | +0.00 |
| SEAP | 0.50 | 0.50 | +0.00 |

This is a small, single-construct-per-molecule cohort; treat it as a directional holdout check, not a powered estimate. All numbers are computational.

Native generator evidence is attached only where a generator emitted the same 19-nt binding site as an assayed guide, so coverage is partial and the `quality` strategy still leans on the scorer; a small difference can be structural, not only data.

Artifacts: `cohort_manifest.json`, `exposure_audit.json`, `protocol.json`, `final_freeze.json`, `results/holdout/comparison.json`.
