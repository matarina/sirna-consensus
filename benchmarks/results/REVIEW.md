# Development benchmark results

Implementation and execution checks passed; the agreed efficacy-improvement criterion did not.

This is exposed Huesken development data: 1,941 full guides across 25 human cloned constructs. It is not independent validation or endogenous knockdown evidence. Model selection and these performance estimates use the same development cohort.

| Method | Mean Top-5 precision, complete genes | Complete genes / 25 | Mean fill rate |
|---|---:|---:|---:|
| native | 77.6% | 25 | 100.0% |
| sirnadesigner | 89.3% | 15 | 77.6% |
| sirnaforge | Unavailable | 0 | 0.0% |
| consensus | Unavailable | 0 | 0.0% |
| model | 90.4% | 25 | 100.0% |
| fusion | 89.6% | 25 | 100.0% |
| support | 90.4% | 25 | 100.0% |
| quality | 92.8% | 25 | 100.0% |

Development-selected quality fusion versus frozen model: **+2.4 percentage points**, paired 95% gene-bootstrap CI **[-0.8, +5.6] points**. The interval includes zero and the point estimate is below the +10-point target.

`native` is the union with sequence-quality fallback; it is not an individual tool. siRNAdesigner precision uses only its 15 complete genes and should not be compared directly with the 25-gene means. siRNAforge and exact consensus have no complete measured panels in this duplex-format cohort; this is missing coverage, not evidence of inefficacy.

Eight component ablations and 1,000 within-gene/assay label permutations were run. They are development diagnostics; pretrained exposure and development selection prevent interpreting their statistics as independent generalization evidence.

The real AREG sequence screen completed against the pinned Ensembl 116 cDNA reference. Its separately screened panels satisfy the computational mismatch criterion only.

See [figures](../figures/benchmark.pdf), [machine-readable results](development_results.json), [source audit](../dataset_manifest.json), and [reproduction instructions](../README.md).
