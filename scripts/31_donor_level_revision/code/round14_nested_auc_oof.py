#!/usr/bin/env python3
"""
Round-14 item 2: the round-13 nested delta-AUC (+0.0279 for STATE) was
in-sample (the reduced/full models were fit and scored on the same 240
patients, no held-out data) -- round-14 flags this as unqualified in an
abstract where every other AUC is out-of-fold, and asks for the genuine
OOF version using the SAME 10-seed x 5-fold protocol as the rest of the
paper (repeated_cv.py).

Construction (unchanged from round-12/13, only the fitting protocol
changes): representation scores (PCA-200, each FM) stay FROZEN, taken from
their own already-computed 10-seed OOF ensemble (as in every prior round).
What's new: the covariate-combination models themselves --
  reduced = QC(3) + log1p(n_cells) + cohort one-hot(8) + PCA-200 score (13)
  full    = reduced + one FM's score (14)
-- are now refit INSIDE a proper 10-seed x 5-fold StratifiedGroupKFold
loop (groups=patient_id, C=0.3, same convention as repeated_cv.py's other
callers), so their predicted probabilities are genuinely out-of-fold, not
fit-and-scored-on-the-same-data. Reduced-model AUC and full-model AUC are
each summarized as repeat-mean +/- sd (repeated_cv.summarize_repeated),
and the OOF nested delta-AUC is the same-seed-paired difference
(repeated_cv.paired_delta_auc): mean, 95% CI, Wilcoxon p.
"""
import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

sys.path.insert(0, '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/code')
from repeated_cv import (DEFAULT_SEEDS, run_repeated_kfold, repeated_patient_auc,
                          summarize_repeated, paired_delta_auc)

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
QC_CSV = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
          '31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv')
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']

qc = pd.read_csv(QC_CSV)
qc = qc[qc['disease'] == 'COVID-19'].copy()
qc['patient_id'] = qc['patient_id'].astype(str)
pm = pd.read_csv(PATIENT_MAP)
pm['patient_id'] = pm['patient_id'].astype(str)
cohort_map = pm.drop_duplicates('patient_id').set_index('patient_id')['Data_NO']
qc['cohort'] = qc['patient_id'].map(cohort_map)
qc = qc.dropna(subset=['cohort'])
qc = qc.sort_values('patient_id').reset_index(drop=True)
y_obs = qc['y'].to_numpy()
assert len(qc) == 240 and int(y_obs.sum()) == 25

qc['log_n_cells'] = np.log1p(qc['n_cells'])
qc['log_counts'] = np.log1p(qc['median_counts'])
qc['log_genes'] = np.log1p(qc['median_genes'])
X_qc_only = qc[['log_counts', 'log_genes', 'median_pct_mito', 'log_n_cells']].to_numpy()
cohort_dummies = pd.get_dummies(qc['cohort'], prefix='coh', drop_first=True).astype(float)
X_with_cohort = np.hstack([X_qc_only, cohort_dummies.to_numpy()])
patient_id_arr = qc['patient_id'].to_numpy()


def load_fm_ensemble_scores(fm, patient_ids):
    seed_scores = []
    for seed in DEFAULT_SEEDS:
        if fm == 'PCA-200':
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
            df = df.merge(pm[['sample_id', 'patient_id']], on='sample_id', how='inner')
        else:
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(patient_ids)]
        seed_scores.append(df.groupby('patient_id')['p_lr'].mean())
    common = seed_scores[0].index
    for s in seed_scores[1:]:
        common = common.intersection(s.index)
    mat = np.column_stack([s.loc[common].to_numpy() for s in seed_scores])
    return pd.Series(mat.mean(axis=1), index=common)


patient_ids_set = set(qc['patient_id'])
model_vecs = {}  # FROZEN, same convention as round12/13
for fm in ['PCA-200'] + FMS:
    fm_scores = load_fm_ensemble_scores(fm, patient_ids_set)
    assert set(fm_scores.index) == patient_ids_set, f'{fm}: missing patients'
    model_vecs[fm] = qc['patient_id'].map(fm_scores).to_numpy()  # aligned to qc row order

pca_vec = model_vecs['PCA-200'].reshape(-1, 1)
X_reduced_full_feat = np.hstack([X_with_cohort, pca_vec])  # 13 features, fixed


def make_fit_predict(X_full_feat):
    def fit_predict(X_tr, y_tr, X_te):
        clf = LogisticRegression(C=0.3, solver='lbfgs', max_iter=2000)
        clf.fit(X_tr, y_tr)
        return {'p_lr': clf.predict_proba(X_te)[:, 1]}
    return fit_predict


print('Running reduced model (QC+cohort+PCA-200), 10-seed x 5-fold OOF...')
oof_reduced = run_repeated_kfold(X_reduced_full_feat, y_obs, groups=patient_id_arr,
                                  fit_predict_fn=make_fit_predict(X_reduced_full_feat))
aucs_reduced, seeds = repeated_patient_auc(oof_reduced, patient_id_arr, y_obs, score_key='p_lr')
summ_reduced = summarize_repeated(aucs_reduced)
print(f'Reduced OOF AUC: mean={summ_reduced["mean"]:.4f}, sd={summ_reduced["sd"]:.4f}')

rows = []
for fm in FMS:
    X_full_feat = np.hstack([X_reduced_full_feat, model_vecs[fm].reshape(-1, 1)])  # 14 features
    print(f'Running full model (+{fm}), 10-seed x 5-fold OOF...')
    oof_full = run_repeated_kfold(X_full_feat, y_obs, groups=patient_id_arr,
                                   fit_predict_fn=make_fit_predict(X_full_feat))
    aucs_full, seeds_full = repeated_patient_auc(oof_full, patient_id_arr, y_obs, score_key='p_lr')
    assert seeds_full == seeds
    summ_full = summarize_repeated(aucs_full)
    pair = paired_delta_auc(aucs_full, aucs_reduced)
    print(f'{fm}: full OOF AUC mean={summ_full["mean"]:.4f} sd={summ_full["sd"]:.4f}  '
          f'OOF nested dAUC={pair["mean_delta"]:+.4f}  '
          f'95%CI=[{pair["ci_lo"]:+.4f}, {pair["ci_hi"]:+.4f}]  '
          f'Wilcoxon p={pair["wilcoxon_p"]:.4f}')
    rows.append({
        'fm': fm,
        'auc_reduced_mean': summ_reduced['mean'], 'auc_reduced_sd': summ_reduced['sd'],
        'auc_full_mean': summ_full['mean'], 'auc_full_sd': summ_full['sd'],
        'oof_nested_delta_auc_mean': pair['mean_delta'],
        'oof_nested_delta_auc_ci_lo': pair['ci_lo'], 'oof_nested_delta_auc_ci_hi': pair['ci_hi'],
        'wilcoxon_p': pair['wilcoxon_p'], 'n_seeds': pair['n_seeds'],
    })

out_df = pd.DataFrame(rows)
out_path = os.path.join(OUT_DIR, 'round14_nested_delta_auc_oof.csv')
out_df.to_csv(out_path, index=False)
print(f'\nSaved {out_path}')
