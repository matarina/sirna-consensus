"""Register an immutable, local-only confirmation attempt before candidate fitting."""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
import subprocess
from pathlib import Path

PKG = Path(__file__).resolve().parents[2]
DEST = PKG / 'benchmarks/confirmation_20261106'


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def main():
    if (DEST / 'protocol.json').exists():
        raise SystemExit('Already registered; never overwrite a preregistration.')
    baseline = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=PKG, text=True).strip()
    root = PKG.parent
    definitions = [
        ('huesken2005', PKG/'benchmarks/data/huesken/measurements.tsv', 'development',
         ['sole_development_source', 'baseline_training_exposed'], ['benchmarks/platform_v3_summary.json']),
        ('sciabola2013', root/'independent_validation/sciabola/canonical/measurements_primary.tsv', 'excluded',
         ['baseline_training_exposed', 'previous_adaptation_and_evaluation'], ['models/development_20261106/data_registry.json']),
        ('legacy_holdout', PKG/'benchmarks/data/holdout/measurements.tsv', 'excluded',
         ['previous_evaluation_exposed', 'upstream_model_selection_exposed', 'artificial_N_to_A_reference_history'],
         ['benchmarks/cohort_manifest.json', 'benchmarks/data/holdout/dataset_manifest.json']),
        ('Harborth', root/'independent_validation/harborth/measurements.tsv', 'excluded',
         ['previous_evaluation_exposed', 'check_empty_mapped_measurements'], ['benchmarks/computational_defensibility_report_zh.html']),
        ('Martinelli', root/'independent_validation/martinelli/measurements.tsv', 'excluded',
         ['previous_evaluation_exposed', 'modified_duplex_outside_primary_all_RNA_scope'], ['benchmarks/computational_defensibility_report_zh.html']),
        ('Shabalina', root/'independent_validation/shabalina/measurements.tsv', 'excluded',
         ['previous_evaluation_exposed', 'paired19_only_missing_full_duplex_identity'], ['benchmarks/computational_defensibility_report_zh.html']),
        ('siRNAEfficacyDB', root/'independent_validation/siRNAEfficacyDB_subset/measurements.tsv', 'excluded',
         ['previous_evaluation_exposed', '470_of_470_cores_in_public_Hu_training_reference'], ['benchmarks/leakage_recheck_20261105.json']),
        ('UBE2H_snapshot', root/'independent_validation/snapshot_independent/measurements.tsv', 'excluded',
         ['previous_evaluation_exposed', 'single_gene_primary_UBE2H'], ['benchmarks/leakage_recheck_20261105.json']),
    ]
    cohorts = []
    for name, path, role, reasons, evidence in definitions:
        item = dict(name=name, path=str(path.resolve()), role=role, exclusion_reasons=reasons,
                    independent_eligible=False, evidence=[])
        if path.is_file():
            with path.open() as h:
                # Census uses identifiers/metadata only; no efficacy summary is calculated.
                rows = list(csv.DictReader(h, delimiter='\t'))
            item.update(sha256=digest(path), records=len(rows),
                        studies=sorted({r['study'] for r in rows}), genes=sorted({r['gene'] for r in rows}))
        else:
            item.update(sha256=None, records=0, studies=[], genes=[])
            item['exclusion_reasons'].append('missing_upstream_artifact_no_fetch')
        for value in evidence:
            p = PKG / value
            item['evidence'].append({'path':str(p.resolve()),'sha256':digest(p) if p.is_file() else None})
        cohorts.append(item)
    oligo = root/'independent_validation/oligoformer_manifest.json'
    cohorts.append(dict(name='OligoFormer',path=str(oligo),sha256=digest(oligo),role='model_not_cohort',
                        records=0,independent_eligible=False,
                        exclusion_reasons=['model_manifest_is_not_an_independent_experimental_dataset',
                                           'upstream_checkpoint_lineage_not_fully_verified']))
    baseline_lock = json.loads((PKG/'benchmarks/platform_v3_freeze.json').read_text())
    artifacts = {}
    for key, sha in baseline_lock['artifacts'].items():
        path = (PKG/'benchmarks'/key).resolve()
        if not path.is_file() or digest(path) != sha:
            raise RuntimeError('Baseline freeze missing or mismatched: '+str(path))
        artifacts[str(path)] = sha
    policy = PKG/'models/default_policy.json'
    write(DEST/'baseline.json',dict(commit=baseline,artifacts=artifacts,default_policy_path=str(policy),
                                  default_policy_sha256=digest(policy),status='frozen_existing_development_default'))
    write(DEST/'cohorts.json',dict(schema_version=1,cohorts=cohorts,confirmation_ids=[],
        eligible_confirmation_records=0,independence='no eligible untouched cohort demonstrated',
        no_downloads=True,original_data_preserved=True,
        limitations=['Prior exposure cannot be undone by retraining on Huesken.',
                     'Unknown lineage is not evidence of absence of leakage.']))
    grid=[]
    for feature in ('sequence_v1','sequence_v2'):
        grid += [dict(feature=feature,family='ridge',alpha=a) for a in (.1,1.,10.,100.)]
        grid += [dict(feature=feature,family='tree',max_depth=d,min_samples_leaf=l,l2_regularization=r)
                 for d,l,r in itertools.product((2,3),(5,10),(.1,1.))]
        grid += [dict(feature=feature,family='pairwise',C=c) for c in (.01,.1,1.,10.)]
    train=cohorts[0]
    protocol=dict(schema_version=1,status='preregistered_before_challenger_fitting',baseline_commit=baseline,
        baseline_manifest_sha256=digest(DEST/'baseline.json'),cohort_registry_sha256=digest(DEST/'cohorts.json'),
        training=dict(studies=['huesken2005'],source=train['path'],source_sha256=train['sha256'],
                      hgnc=str(PKG/'references/gencode_v50/hgnc_complete_set.txt'),
                      hgnc_sha256=digest(PKG/'references/gencode_v50/hgnc_complete_set.txt'),
                      outer_folds=5,inner_folds=3,grid=grid,seed=20261106,threads=4,
                      weight='equal gene, equal assay within gene, equal record within assay; mean one',
                      calibration='none; relative ranking only'),
        primary_metric='macro_gene_assay_NDCG@5',
        aggregation='mean assay within canonical gene, then equal mean over genes',
        relevance='within-assay midrank (rank-1)/(n-1); gain=2**relevance-1',
        minimum_candidates=5,constant_activity='undefined_ndcg_reported_separately',
        score_ties=['guide_DNA','molecule_id','record_id'],
        active_thresholds={'normalized_activity':.7,'percent_inhibition':70.,'percent_knockdown':70.},
        bootstrap_replicates=10000,seed=20261106,
        bootstrap_unit='paired independent gene/sequence components; retain all assays of each gene',
        metrics=['macro_gene_assay_NDCG@5','precision@5','Top3_at_least_one_active','Top5_active_fraction',
                 'panel_fill_rate','prediction_coverage','measurement_coverage','abstention','budget_yield_bounds','Spearman'],
        deployment=dict(panel='diverse_top',k=5,min_spacing=30,chemistry='all_RNA',
                        match='complete guide/passenger/overhang/known chemistry and same assay; never overlap',
                        unknown='paired same-molecule unknowns cancel; remaining unknowns use worst-case bounds',
                        no_completion_loss_panels=['diverse_top','screened_diverse_top']),
        gate=dict(minimum_independent_studies=2,minimum_independent_components=5,
                  ndcg_strictly_above=.80,gain_strictly_above=0.,paired_ci95_lower_strictly_above=0.,
                  no_completion_loss=True,no_prediction_coverage_loss=True,
                  actual_generated_budget_gain_lower_strictly_above=0.,all_required=True),
        confirmation_ids=[],loso=dict(status='not_permitted_single_development_study',
                     reason='External cohorts may not be fitted even for LOSO.'),
        policy=dict(no_confirmation_training=True,no_confirmation_selection=True,no_default_change_without_all_gates=True,
                    no_new_test_after_observation=True,no_silent_download=True,baseline_retained_on_failure=True),
        preregistration_limit='Local git preregistration; no third-party timestamp or proof of undocumented history.')
    write(DEST/'protocol.json',protocol)
    print(json.dumps({'baseline_commit':baseline,'eligible_confirmation_records':0,'grid':len(grid)}))


if __name__=='__main__':main()
