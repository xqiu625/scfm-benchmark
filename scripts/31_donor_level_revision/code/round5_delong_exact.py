import numpy as np, pandas as pd, os

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
PATIENT_MAP = ('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/'
               '31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv')
SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']

m = pd.read_csv(PATIENT_MAP)
covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])

def build_ensemble():
    ensembled = {}
    for fm in FMS + ['PCA-200']:
        seed_scores = []
        for seed in SEEDS:
            if fm == 'PCA-200':
                df = pd.read_csv(os.path.join(OOF_DIR, f'oof_PCA-200_seed{seed}.csv'))
                pmap = m[['sample_id', 'patient_id']]
                df = df.merge(pmap, on='sample_id', how='inner')
            else:
                df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
            df = df[df['patient_id'].isin(covid_patients)]
            g = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
            seed_scores.append(g)
        common = seed_scores[0][[]]
        for s in seed_scores[1:]:
            common = common.join(s[[]], how='inner')
        y = seed_scores[0].loc[common.index, 'y']
        scores = np.mean([s.loc[common.index, 'score'].values for s in seed_scores], axis=0)
        ensembled[fm] = pd.DataFrame({'y': y.values, 'score': scores}, index=common.index)
    return ensembled

ensembled = build_ensemble()
idx = ensembled['STATE'].index
for fm in FMS + ['PCA-200']:
    ensembled[fm] = ensembled[fm].loc[idx]

y = ensembled['STATE']['y'].values
pos = np.where(y == 1)[0]
neg = np.where(y == 0)[0]
m_pos, n_neg = len(pos), len(neg)
print(f'n_patients={len(y)}, n_pos={m_pos}, n_neg={n_neg}')

models = FMS + ['PCA-200']
V10 = {}
V01 = {}
AUC = {}
for fm in models:
    s = ensembled[fm]['score'].values
    sx = s[pos]  # positive scores
    sy = s[neg]  # negative scores
    psi = (sx[:, None] > sy[None, :]).astype(float) + 0.5 * (sx[:, None] == sy[None, :])
    V10[fm] = psi.mean(axis=1)  # per positive patient
    V01[fm] = psi.mean(axis=0)  # per negative patient
    AUC[fm] = psi.mean()

V10_mat = np.column_stack([V10[fm] for fm in models])
V01_mat = np.column_stack([V01[fm] for fm in models])
S10 = np.cov(V10_mat, rowvar=False, ddof=1)
S01 = np.cov(V01_mat, rowvar=False, ddof=1)
Cov = S10 / m_pos + S01 / n_neg

print('\nAUCs:', {fm: round(AUC[fm],4) for fm in models})
print('\nPairwise DeLong (exact) vs bootstrap SD from se_check_delong.csv:')
rows = []
for i, a in enumerate(models):
    for j, b in enumerate(models):
        if i >= j: continue
        var_diff = Cov[i,i] + Cov[j,j] - 2*Cov[i,j]
        se_exact = np.sqrt(max(var_diff, 0))
        delta = AUC[a] - AUC[b]
        z = delta / se_exact if se_exact > 0 else np.nan
        from scipy.stats import norm
        p = 2*(1-norm.cdf(abs(z))) if not np.isnan(z) else np.nan
        rows.append((f'{a} - {b}', round(delta,4), round(se_exact,5), round(p,5)))
        print(f'{a:12s} - {b:12s}: delta={delta:+.4f}  DeLong_exact_SE={se_exact:.5f}  p={p:.5f}')

pd.DataFrame(rows, columns=['pair','delta_auc','delong_exact_SE','delong_exact_p']).to_csv(
    os.path.join(OOF_DIR, 'paired_inference_full', 'delong_exact.csv'), index=False)
print('\nSaved to', os.path.join(OOF_DIR, 'paired_inference_full', 'delong_exact.csv'))
