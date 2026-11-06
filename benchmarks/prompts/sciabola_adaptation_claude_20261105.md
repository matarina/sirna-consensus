# Implement and evaluate Sciabola adaptation of the OligoFormer scorer

Work in `/data1/maxw/proj/cargo_denovo`. The user asks you to execute this plan, including implementation, appropriate checks, actual experiments, and a reviewable report. Continue through completion without asking for routine implementation choices. Do not stop at a plan or a smoke test. Make a defensible estimate of adaptation performance, whether it improves or not.

## Objective and claim

Improve the siRNA efficacy scorer on Sciabola-like assays using small adaptations of the frozen OligoFormer model. Evaluate adaptation with nested gene-held-out validation. Once Sciabola labels are used for fitting or model selection, describe Sciabola as adaptation/development data; it is no longer an independent external test. These are retrospective experiments because existing Sciabola results have already been inspected. A genuinely untouched external cohort is needed for subsequent confirmatory validation.

This task evaluates ranking of supplied measured candidates. Do not present it as validation of the full generator/consensus pipeline or of biological specificity. Do not promise improvement or tune until a numerical gate passes.

## Read these existing artifacts first

- `sirna-consensus/README.md` and `sirna-consensus/benchmarks/README.md`.
- `sirna-consensus/benchmarks/three_benchmark_report_zh.html` and `computational_defensibility_report_zh.html`.
- `sirna-consensus/benchmarks/leakage_recheck_20261105.json`.
- `sirna-consensus/benchmarks/cohort_manifest.json`, `exposure_audit.json`, `protocol.json`, and `results/benchmark_selector_policy.json`.
- `independent_validation/oligoformer/model.py`, `infer.py`, `best_model.pth`, and `independent_validation/oligoformer_manifest.json`.
- `sirna-consensus/benchmarks/scripts/oligoformer_worker.py` and the existing efficacy evaluator/aggregation code.
- `independent_validation/sciabola/prepare_canonical.py`, `rebuild_by_dose.py`, `canonical/measurements_primary.tsv`, `canonical/eval_primary.json`, canonical selections, and their saved scorer inputs/outputs.
- Original Sciabola source: DOI `10.1093/nar/gks1191`, supplementary Tables S3/S4. Local raw table: `independent_data/raw/sciabola2013_extra.csv`.

Check applicable AGENTS.md/CLAUDE.md instructions. Preserve unrelated changes, original checkpoints, source data, and historical benchmark artifacts. Write new experiments separately, for example under `independent_validation/sciabola/adaptation_20261105/`, with reusable source/scripts in an appropriate project directory. Do not overwrite existing reports or make the adapted scorer the deployment default.

## Established findings to verify and retain

1. The five-gene siRNAEfficacyDB subset contains 470 guides, all 470 of whose 19-nt antisense cores match upstream OligoFormer `Hu.csv`. It is training-exposed and cannot independently validate this model.
2. UBE2H has 70 guides with zero exact matches against the public Hu/Mix/Taka files checked. They are attributed to the Huesken source study, so provenance and exact released checkpoint lineage remain relevant. The file `snapshot_independent/measurements.tsv` actually contains 77 rows across several genes; filter explicitly to UBE2H when referring to the 70-guide cohort.
3. Sciabola has 356 primary measurements with zero exact guide-core matches against those Hu/Mix/Taka files. Zero exact overlap does not establish absence of near-duplicates, related genes, exposure in another learned component, or model-selection exposure.
4. The checked Hu/Mix/Taka hashes match the existing local audit. Some original upstream local paths recorded in old manifests are now absent. Reuse available pinned copies or obtain official files and verify hashes; do not silently use different data.
5. Existing frozen Sciabola performance is approximately +1.1 percentage points in gene-weighted Top-5 precision versus random, with a CI crossing zero. Reproduce the actual saved result and resolve discrepancies before training.
6. Primary labels are percent inhibition at 10 nM, except HPSE at 2 nM. Do not average doses. The primary table has these gene counts and >=70% active counts:
   - BIRC5: 9 / 9 active; BRAF: 5 / 1; CTNNB1: 5 / 5; EZH2: 9 / 6.
   - HIF1A: 100 / 80; HK2: 97 / 73; HPSE: 107 / 76.
   - MTOR: 10 / 5; MYC: 9 / 1; PIK3CA: 5 / 4.
   Thus 304/356 rows come from three genes. Several gene pools cannot show Top-5 enrichment over random. Using one pool per gene, the ideal-label upper bound on mean Top-5 precision improvement is about 16.59 percentage points. Recompute the bound using the exact evaluator's assay/transcript aggregation; label it an oracle diagnostic, never a predictor. Keep all prespecified genes in the primary analysis.

## Phase 1: Audit inputs, inference, and molecular identity

- Verify antisense orientation, binding19 plus documented 2-nt overhang representation, exact target placement, transcript context, boundary padding, score direction, label units, duplicate handling, and checkpoint/artifact hashes. Do not invent passenger identity or transfer activity between neighboring guides.
- Inspect the LSTM shape issue in the local model: nn.LSTM defaults to batch_first=False, while forward appears to pass batch-by-position tensors. Empirically check invariance to batch companions/order. Preserve the original checkpoint-compatible baseline and use batch-size-one frozen feature extraction. Do not silently set batch_first=True and reuse the weights as if predictions retained their original semantics. Any architectural correction is a separate retrained model and outside the initial adaptation comparison.
- Check softmax dimensions and the precise score used. Avoid applying softmax twice or silently changing the baseline's score contract. Preserve baseline reproduction and record any required corrections as separate variants with tests.
- Audit exact normalized binding cores, sequence similarity, source studies, gene aliases/homologs, and exposure in supervised learned components where evidence is available. Distinguish supervised efficacy exposure from unlabeled sequence pretraining. Record unknown checkpoint lineage honestly.

## Phase 2: Freeze the experimental protocol before fitting

Save a dated, hashed protocol/config and immutable split manifest before fitting. This is a prospective lock for the new experiment, not a claim that Sciabola was historically untouched.

- Outer evaluation: leave one Sciabola gene out, 10 folds unless sequence-related grouping requires merging folds. All transcripts/isoforms, doses, replicates, and copies of a guide belonging to an outer test gene remain excluded from fitting and validation for that fold.
- If highly similar sequences connect nominally different genes, define label-independent sequence groups and merge/purge connected groups as needed to prevent cross-fold leakage. Prespecify the similarity rule and disclose the resulting loss of coverage. Do not split adjacent tiling windows randomly.
- Inner selection: gene-grouped validation only within each outer training set, e.g. three grouped folds. Use it for scaling, feature/model selection, regularization, learning rates, loss weights, early stopping, and epoch choice. Never select using the outer test fold.
- Use a small finite search specified in the protocol, equal gene weighting, deterministic seeds, and an average of three fixed seeds for stochastic fits. Do not choose the best seed or expand the search after viewing outer-fold results.
- Any scaler, imputer, feature selector, calibrator, or fitted normalization must be trained within the appropriate training fold. Label-free extraction from an immutable frozen model may be cached for all candidates; fitted transforms may not.
- Start with primary endpoints only. Additional doses, if included at all, must be a separately prespecified secondary experiment with every guide's measurements in the same fold. Primary data confound HPSE with its 2 nM dose, so do not interpret a fitted dose effect as independently identified.
- Audit whether held-out Sciabola genes were seen by the original pretrained scorer. Distinguish held out during adaptation from genes genuinely unseen during all supervised pretraining.

## Phase 3: Implement and run the following compact comparison

Use the same frozen initial checkpoint for every applicable outer fold and keep the feature encoders frozen initially.

1. **Frozen baseline:** existing OligoFormer predictions, reproduced from canonical inputs.
2. **Small regularized regression head:** a linear/ridge head on a fixed compact representation, preferably the frozen classifier's 64-dimensional penultimate representation, with a prespecified small regularization grid. Verify the extraction point. Fit continuous percent inhibition / 100 without conflating it with Hu's normalized activity. If adding thermodynamic features, specify that feature set before fitting and validate the choice only inside inner folds.
3. **Small ranking head:** a linear or very small head on the same frozen representation, trained with within-gene, within-assay pairwise comparisons of continuous efficacy. Weight genes equally and control pair sampling so dense genes do not dominate. Pick regularization and any combined regression/ranking loss weight using only inner folds. Retain raw continuous labels; the >=70% threshold remains an evaluation definition.
4. **Limited prediction-head fine-tuning:** train a tightly restricted final portion of the classifier on frozen encoder features, using low learning rates and strong regularization with inner-fold early stopping. Keep RNA-FM and the siRNA/mRNA sequence encoders frozen. Freeze encoder dropout in eval mode. Do not accidentally train the entire large classifier as a supposedly small head. Use a common fold-specific training schedule selected without outer labels.

Run all four prespecified arms if technically viable. A documented technical incompatibility is a valid reason to mark an arm unavailable; disappointing results are not. Do not escalate to full-model fine-tuning based on outer-test outcomes. Monotonic calibration cannot improve within-pool Top-5 ordering; do not sell it as a ranking fix.

## Phase 4: Evaluate and report

- Save a single out-of-fold prediction per primary measurement per model, fold membership, validation selections, learning curves, and all fold-specific trained checkpoints/heads. Never pool in-sample fitted predictions with out-of-fold predictions.
- Primary: the existing equal-gene/equal-assay Top-5 precision and paired improvement over the frozen scorer. Also show random-pool enrichment, fill rate and measured coverage.
- Secondary: within-gene Spearman correlation, AUROC and average precision where defined, Top-10 where the pool permits, and continuous error. All-active/all-inactive genes have undefined AUROC; do not assign zero or drop them from unrelated metrics. Report pooled AUROC separately from within-gene metrics.
- Report per-gene effects and paired gene-level uncertainty. Ten genes are a small sample, and outer-fold training sets overlap; explain that a bootstrap of fixed out-of-fold predictions does not include all training/model-selection uncertainty. Report seed variability as well. Avoid significance claims unsupported by this design.
- Show the oracle Top-5 ceiling and identify nondiscriminating small/saturated pools without post-hoc removal. Do not change dose, threshold, gene inclusion, or primary metric to obtain a pass.
- Compare the fixed arms as specified. If using nested validation to choose an adaptive model, choose within each outer fold and report that selection procedure's outer predictions as a distinct arm. Do not report the best outer-scoring arm as an unbiased winner; its final superiority needs a new untouched cohort.
- Write a concise Markdown report (and an HTML report if convenient) covering dataset provenance, protocol, model changes, baseline reproduction, per-gene/aggregate results, coverage, uncertainty, leakage limitations, commands, hashes, and whether adaptation actually helped. If no model improves reliably, say so.
- Keep the 470 training-exposed guides out of any independent-validation claim. Do not repurpose already-inspected UBE2H or Sciabola results as untouched confirmatory evidence.

## Execution, checks, and handoff

Reuse existing compatible environments, including the project's `.scorer-venv` if verified. Cache frozen features to make repeated fitting cheap. CPU is adequate for a compact feature/head experiment; if GPUs or remote hosts become necessary, follow the applicable lab resource/deployment instructions, check live availability, and preserve all unrelated jobs and environments. No broad environment upgrades or system/service changes.

Add meaningful checks for split disjointness across guides/doses/genes, fold-local fitting, deterministic baseline reproduction, frozen encoder state, finite losses/gradients, and completeness of out-of-fold predictions. Ensure resumption cannot mix folds or protocol hashes. Keep commands and outputs reproducible.

At completion report what was implemented and run, the strongest supported finding without cherry-picking, exact artifact paths, unresolved limitations, and the next requirement for independent validation. If a genuine external blocker prevents completion, preserve completed work and identify it precisely. Do the work now.
