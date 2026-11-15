"""Round-5 item 3: each FM + PCA-200 vs count+disease null, exact DeLong, within-COVID.
Builds the null model's own 10-seed ensemble OOF score (cell-expanded, log1p(n_cells) only,
since within-COVID is single-disease so the disease one-hot term is degenerate)."""
import numpy as np, pandas as pd, os
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
from scipy.stats import norm

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
PATIENT_MAP = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv'
QC_CSV = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results_patient_v2/qc_baseline/patient_qc_distributions.csv'
SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']

def safe_auc(y, p):
    return roc_auc_score(y, p) if len(np.unique(y)) > 1 else np.nan

qc = pd.read_csv(QC_CSV)
qc = qc[qc['disease'] == 'COVID-19'].copy().sort_values('patient_id').reset_index(drop=True)
qc['patient_id'] = qc['patient_id'].astype(str)
y = qc['y'].to_numpy()
groups = qc['patient_id'].to_numpy()
X_null = np.log1p(qc['n_cells'].to_numpy()).reshape(-1, 1)
print(f'null model cohort: n={len(y)}, fatal={int(y.sum())}')

# null model 10-seed ensemble OOF
oof_all = np.zeros((len(SEEDS), len(y)))
for si, seed in enumerate(SEEDS):
    skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
    oof = np.full(len(y), np.nan)
    for tr, te in skf.split(X_null, y, groups=groups):
        clf = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
        clf.fit(X_null[tr], y[tr])
        oof[te] = clf.predict_proba(X_null[te])[:, 1]
    oof_all[si] = oof
p_null = oof_all.mean(axis=0)
auc_null = safe_auc(y, p_null)
print(f'null model ensemble AUC = {auc_null:.4f}')

m = pd.read_csv(PATIENT_MAP)
def load_fm(fm, patient_ids):
    seed_scores = []
    for seed in SEEDS:
        if fm == 'PCA-200':
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
            pmap = m[['sample_id', 'patient_id']]
            df = df.merge(pmap, on='sample_id', how='inner')
        else:
            df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
        df['patient_id'] = df['patient_id'].astype(str)
        df = df[df['patient_id'].isin(patient_ids)]
        g = df.groupby('patient_id')['p_lr'].mean()
        seed_scores.append(g)
    common = seed_scores[0].index
    for s in seed_scores[1:]:
        common = common.intersection(s.index)
    mat = np.column_stack([s.loc[common].to_numpy() for s in seed_scores])
    return pd.Series(mat.mean(axis=1), index=common)

patient_ids = set(qc['patient_id'])
rows = []
for fm in FMS + ['PCA-200']:
    fs = load_fm(fm, patient_ids)
    idx = qc[qc['patient_id'].isin(fs.index)].index
    yb = y[idx]
    p_n = p_null[idx]
    p_f = qc.loc[idx, 'patient_id'].map(fs).to_numpy()
    pos = np.where(yb == 1)[0]; neg = np.where(yb == 0)[0]
    def placement(s):
        sx, sy = s[pos], s[neg]
        psi = (sx[:,None] > sy[None,:]).astype(float) + 0.5*(sx[:,None]==sy[None,:])
        return psi.mean(axis=1), psi.mean(axis=0), psi.mean()
    v10n, v01n, aucn = placement(p_n)
    v10f, v01f, aucf = placement(p_f)
    S10 = np.cov(np.column_stack([v10n, v10f]), rowvar=False, ddof=1)
    S01 = np.cov(np.column_stack([v01n, v01f]), rowvar=False, ddof=1)
    Cov = S10/len(pos) + S01/len(neg)
    var_diff = Cov[0,0]+Cov[1,1]-2*Cov[0,1]
    se = np.sqrt(max(var_diff,0))
    delta = aucf - aucn
    z = delta/se if se>0 else np.nan
    p = 2*(1-norm.cdf(abs(z))) if not np.isnan(z) else np.nan
    rows.append({'model': fm, 'n_patients': len(yb), 'n_fatal': int(yb.sum()),
                 'auc_null': round(aucn,4), 'auc_model': round(aucf,4),
                 'delta_vs_null': round(delta,4), 'delong_SE': round(se,5), 'delong_p': round(p,5)})
    print(f'{fm:12s} vs null: null={aucn:.4f} model={aucf:.4f} delta={delta:+.4f} p={p:.5f}')

pd.DataFrame(rows).to_csv(os.path.join(OUT_DIR, 'each_model_vs_null.csv'), index=False)
print('Saved.')
