# Experimental benchmark package

Primary endpoint: gene-weighted Top-5 measured precision; secondary Top-10,
coverage, fill rate, recovered hits, AUROC and average precision. Reporter activity
and endogenous inhibition are separate endpoints. Never label untested guides
ineffective. Never use ±2-nt overlap to transfer experimental activity.

The downloaded publisher-linked Figshare v1 snapshot contains 3,544 records.
`dataset_manifest.json` records source checksums and discrepancies. Its seven
studies occur in OligoFormer development/evaluation collections, so **no independent
holdout is currently established**. The Huesken cohort was also previously used by
the sibling project. All current results are development diagnostics.

## Prepare original-source development data

```bash
python benchmarks/scripts/prepare_huesken.py \
  --source ../sirnaforge-platform/data --out benchmarks/data/huesken
python benchmarks/scripts/fetch_snapshot.py
python benchmarks/scripts/audit_sources.py
```

The first command validates 1,941 human guides across 25 original cloned constructs,
preserves full 21-nt guides and explicit overhangs, retains raw normalized activity
(including values above 1), and writes label-free candidates, measurements,
transcripts, exclusions and hashes. No current full-transcript coordinate is
substituted for a construct coordinate. Unknown passenger strands remain unknown.

Measurement TSV columns:

```
record_id study gene assay transcript_id target_sha256 guide_dna duplex_format endpoint value
```

Use tab delimiters. `assay` must distinguish cell, dose, time and readout; unknown
conditions remain explicit. `endpoint` is `percent_inhibition` or
`normalized_activity`. Thresholds are 70 and 0.7 respectively, not interchangeable
units. Intermediate values stay in the denominator. Repeated same-guide/condition
observations require an explicit shared `replicate_group`; these technical
replicates are averaged before thresholding. Different conditions remain separate.

## Native tool outputs and scorer

```bash
SIRNA_CONSENSUS_PYTHON=/path/to/tools/python python benchmarks/scripts/run_generators.py \
  --transcripts benchmarks/data/huesken/transcripts --out benchmarks/results/runs/huesken
```

The OligoFormer worker environment needs PyTorch, NumPy, pandas 2.2.3, RNA-FM 0.2.2
and ptflops. Matplotlib is only needed for plots. Reuse an existing compatible
PyTorch environment; keep model dependencies outside the zero-dependency package.
Fetch RNA-FM weights from the official `cuhkaih/rnafm` release and use the audited
OligoFormer model/feature source and `best_model.pth`. No inference-time downloads.

```bash
/path/to/scorer/python benchmarks/scripts/create_scorer_manifest.py \
  --upstream ../sirnaforge-platform/data/raw/oligoformer \
  --rna-fm-checkpoint benchmarks/models/RNA-FM_pretrained.pth \
  --out benchmarks/models/oligoformer_manifest.json --device cpu
python benchmarks/scripts/run_development.py --data benchmarks/data/huesken \
  --tool-runs benchmarks/results/runs/huesken \
  --scorer-manifest benchmarks/models/oligoformer_manifest.json \
  --out benchmarks/results/development/full
```

This freezes model predictions once per construct, reuses them for strategies,
and compares individual generators, exact consensus, union with fallback, model
alone and fusion strategies. Native evidence is attached only to exact guides
with the same binding core/overhang; passenger identity not available in the
measurement table is not invented. Native generator non-emission is abstention.

Without `--tool-runs`, the script is a sequence-quality/scorer smoke benchmark,
**not** a generator or consensus comparison. Development choice is recorded in
`development_summary.json`; it is not an independent performance estimate.

The frozen model uses batch size one because upstream's LSTM layout makes changing
batch size affect inference. Eval mode disables dropout. Model output remains an
uncalibrated ranking score, not a measured knockdown percentage.

## Evaluate and compare

```bash
python benchmarks/evaluate.py --benchmark benchmarks/data/huesken/measurements.tsv \
  --run selection/selection.json --track measured --top-k 5,10 --output metrics.json
python benchmarks/scripts/compare.py --candidate candidate_metrics.json \
  --baseline baseline_metrics.json --out comparison.json
```

Measured-track ranking intersects the candidate list with the shared assay universe
**before** selecting K. Deployment evaluation uses the actual saved Top-K and
reports exact-match measurement coverage; save a budget of at least 10 to evaluate
Top-10. Zero measured matches produces unavailable precision. Missing model scores
are not converted to ineffective labels.

Comparisons require identical dataset hashes, thresholds, panel protocols,
transcripts, construct hashes and assay units. Average assays within genes, then
weight genes equally. Report completion and observed hits per five-slot budget
alongside precision. At least five complete paired genes are required for the
numeric gate; five is a minimum feasibility floor, not a power guarantee.

A numerical improvement requires >=0.10 precision gain, a paired gene-bootstrap
95% CI above zero, no completion loss, and positive budget-yield change. Independent
validation additionally requires an audited unused study/gene cohort; numerical
success alone never sets that flag. Default bootstrap replicates: 10,000.

Use `--threshold 0.6` and `0.8` only for prespecified sensitivity analyses.
`--permutations 1000` shuffles aggregated labels within gene/assay/format groups.
Random eligible emitted-pool selection is reported separately from candidate
coverage; it conditions on that method's eligibility, not an exhaustive transcript.

## Freeze and final-evaluation gate

```bash
python benchmarks/scripts/freeze.py --lock freeze.json protocol.json cohort.json \
  measurements.tsv selection/selection.json model_manifest.json source_files.py
python benchmarks/scripts/freeze.py --lock freeze.json --verify
python benchmarks/evaluate.py --phase final --freeze-lock freeze.json \
  --protocol protocol.json --cohort-manifest cohort.json \
  --benchmark measurements.tsv --run selection/selection.json --output final.json
```

Final mode verifies the lock before reading labels and requires a cohort manifest
with `independent_holdout_eligible: true` and `checkpoint_training_audit: verified`.
The supplied source audit does not satisfy this gate. Include source code,
preprocessing, policy, baseline choice, model/checkpoint manifests, splits and
label-blind predictions in the freeze. Never overwrite a frozen lock or choose a
new holdout after inspecting its results.

## Controls from frozen outputs

```bash
python benchmarks/scripts/replay_ablations.py --development benchmarks/results/development/full \
  --benchmark benchmarks/data/huesken/measurements.tsv --out benchmarks/results/development/ablations
python benchmarks/scripts/permutation_test.py --development benchmarks/results/development/full \
  --benchmark benchmarks/data/huesken/measurements.tsv --out benchmarks/results/permutation_control.json
```

These commands replay saved predictions; they do not refit the model or invent scores
for unmeasured candidates. Permutation inference is a development diagnostic and
does not remove pretrained-data exposure or account for selecting a winner on the
same development set.

## Figures and artifacts

```bash
/path/to/plot/python benchmarks/scripts/plot_results.py --result metrics.json \
  --comparison comparison.json --selection selection/selection.json \
  --measurements measurements.tsv --out figures
```

Figures export PNG, PDF and SVG. Case studies must be chosen on development-data
coverage and assay suitability, not observed advantage. Plot strata separately;
show complete-panel performance alongside coverage/completion tables. Raw source
snapshots and bulky per-run artifacts are local and ignored by Git; manifests,
protocol, scripts, compact reports and figures form the reviewable package.

`report_zh.html` is a Chinese, paper-style rendering of the frozen development
results (abstract, context, methods, results, figures). `make_report_zh.py` reads
the result JSONs verbatim — no recomputation, figures are inline SVG:

```bash
python benchmarks/make_report_zh.py -o benchmarks/report_zh.html
```
