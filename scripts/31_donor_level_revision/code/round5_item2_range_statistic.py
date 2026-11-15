"""Round-5 item 2: range-statistic permutation for leaders-vs-UCE gap.
Reuses jobS2_round4_primary.py's ensemble construction + permutation loop,
adds a second statistic: max_j(AUC_j) - min_j(AUC_j) across the 5 FMs (no PCA)."""
import numpy as np, pandas as pd, os
from sklearn.metrics import roc_auc_score

OOF_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/revision_donor_level_v2'
OUT_DIR = os.path.join(OOF_DIR, 'paired_inference_full')
PATIENT_MAP = '/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/31_donor_level_revision/results/donor_meta_with_patient_id_v2.csv'
SEEDS = list(range(42, 52))
FMS = ['STATE', 'UCE', 'Geneformer', 'scPRINT-2', 'scGPT']
NPERM = 10000
RNG_SEED = 20260908

def safe_auc(y, p):
    return roc_auc_score(y, p) if len(np.unique(y)) > 1 else np.nan

m = pd.read_csv(PATIENT_MAP)
covid_patients = set(m.loc[m['disease'] == 'COVID-19', 'patient_id'])

ensembled = {}
for fm in FMS:
    seed_scores = []
    for seed in SEEDS:
        df = pd.read_csv(os.path.join(OOF_DIR, f'oof_patient_{fm}_seed{seed}.csv'))
        df = df[df['patient_id'].isin(covid_patients)]
        g = df.groupby('patient_id').agg(y=('y', 'max'), score=('p_lr', 'mean'))
        seed_scores.append(g)
    common = seed_scores[0][[]]
    for s in seed_scores[1:]:
        common = common.join(s[[]], how='inner')
    ensembled[fm] = seed_scores, common.index

common_idx = ensembled['STATE'][1]
for fm in FMS[1:]:
    common_idx = common_idx.intersection(ensembled[fm][1])

y = ensembled['STATE'][0][0].loc[common_idx, 'y'].to_numpy()
scores = {}
for fm in FMS:
    seed_scores, _ = ensembled[fm]
    mat = np.column_stack([s.loc[common_idx, 'score'].to_numpy() for s in seed_scores])
    scores[fm] = mat.mean(axis=1)

n = len(y)
print(f'n_patients={n}, n_fatal={int(y.sum())}')
assert n == 240 and int(y.sum()) == 25, 'sanity gate failed'

obs_aucs = {fm: safe_auc(y, scores[fm]) for fm in FMS}
obs_range = max(obs_aucs.values()) - min(obs_aucs.values())
print('observed AUCs:', {k: round(v,4) for k,v in obs_aucs.items()})
print(f'observed range (max-min) = {obs_range:.4f}')

rng = np.random.default_rng(RNG_SEED)
perm_ranges = np.empty(NPERM)
for i in range(NPERM):
    yp = rng.permutation(y)
    aucs = [safe_auc(yp, scores[fm]) for fm in FMS]
    perm_ranges[i] = max(aucs) - min(aucs)

p_range = (1 + np.sum(perm_ranges >= obs_range)) / (1 + NPERM)
print(f'range-statistic permutation P = {p_range:.5f} (nperm={NPERM}, rng_seed={RNG_SEED})')

pd.DataFrame([{'statistic': 'max-min AUC over 5 FMs', 'observed_range': round(obs_range,4),
               'permutation_P': round(p_range,5), 'nperm': NPERM, 'rng_seed': RNG_SEED,
               **{f'auc_{fm}': round(v,4) for fm,v in obs_aucs.items()}}]).to_csv(
    os.path.join(OUT_DIR, 'range_statistic_leaders_vs_uce.csv'), index=False)
print('Saved.')
