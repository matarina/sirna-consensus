# AREG example run

Real run against `NM_001657.4` (1234 nt), v1.0.0, with all three external tools
installed. The input FASTA is committed here as `NM_001657.4.fa`.

```bash
docker pull ghcr.io/austin-s-h/sirnaforge:sha-4c05e4b
# set SIREN_BIN, SIRNA_CONSENSUS_PYTHON, SIRNA_CONSENSUS_TOOL_BIN/LIB
sirna-consensus run \
  --fasta examples/AREG/NM_001657.4.fa \
  --out examples/AREG \
  --tools sirnadesigner,siren,sirnaforge \
  --min-votes 2 --top-k 10 --exclude-singletons
```

## Result

- 1,246 raw hits placed from 3 tools.
- **204 generator sites** after bounded site grouping.
- **59 sites** received two generator votes (siRNAdesigner + siRNAforge).
- **34** passed the sequence filters.
- **10** selected Top-K guides with ≥30 nt spacing.
- SIREN is an **evaluator** here, not a vote. It produced 607 input-only
  enumeration hits because no transcriptome was supplied; `specificity_status`
  is `input_only`, so this example does **not** demonstrate transcriptome-wide
  specificity.

`top_guides.fa` contains only the 10 selected guides. `all_candidates.fa`
contains all 59 candidate sites with explicit `passed=true/false` and flags.

This is a reproducibility artifact, not a wet-lab benchmark. Use
`benchmarks/evaluate.py` with a measured knockdown dataset for that comparison.

`report_zh.html` is a Chinese, paper-style write-up of this same run (abstract,
context, methods, results, figures) rendered straight from `run.json` by
`make_report_zh.py` — no recomputation, figures are inline SVG:

```bash
python examples/AREG/make_report_zh.py examples/AREG/run.json -o examples/AREG/report_zh.html
```
