# Benchmarking `sirna-consensus`

This directory is a harness for the scientific milestone, not a claim that the
current package has been experimentally validated. No wet-lab benchmark data is
bundled here.

## What to measure

For a benchmark of experimentally measured siRNA knockdown, compare:

1. each external generator alone;
2. generator consensus alone;
3. consensus plus sequence filters;
4. consensus plus ranking.

Report:

- **recovery@K** — share of known effective sites matched by the Top-K;
- **precision@K** — share of Top-K guides that match a known effective site;
- the effect of `--min-votes`, filters, and ranking components;
- results by held-out gene or dataset to avoid leakage.

## Benchmark format

A tab-separated file needs:

```text
transcript_id	target_start	target_end	effective
NM_001657.4	343	363	1
NM_001657.4	495	515	1
```

`target_start`/`target_end` are 1-based inclusive on the transcript. `effective`
is 1/0 or a boolean-like string. To use knockdown values instead:

```bash
python benchmarks/evaluate.py \
  --benchmark benchmark.tsv \
  --run out/run.json \
  --effective-column knockdown --threshold 0.5
```

## Comparing runs

```bash
python benchmarks/evaluate.py \
  --benchmark benchmark.tsv \
  --run sirnadesigner/run.json \
  --run sirnaforge/run.json \
  --run consensus/run.json \
  --top-k 5,10 \
  --output benchmark_results.json
```

Each run should be produced from the same transcript and the same filter settings.
Because SIREN is an evaluator in v1, a SIREN-only run is not a design-consensus
comparison; use it as specificity evidence, not as a generator baseline.
