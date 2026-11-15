#!/usr/bin/env python3
"""
Cross-disease mortality-subspace alignment (roadmap H + N).

The manuscript claims the cross-disease transfer asymmetry has a geometric basis but
its 1D weight-cosine test was null (<0.07); it proposes a principal-angle/CCA subspace
analysis as the real test (Discussion, geometric mechanism). This runs that test.

For each FM:
  - w_covid  = LR(fatal vs survived) weight vector trained on COVID cells
  - w_sepsis = LR weight vector trained on Sepsis cells
  - cos_1d   = cosine(w_covid, w_sepsis)            -> reproduces the 1D null
              (N: for scGPT a NEGATIVE cos_1d means the mortality direction *inverts*
               between diseases -> the 0.352 anti-correlation is a learned sign-flip,
               NOT a gene-vocab artifact. scGPT vocab overlap is 13037/13066 = 99.8%.)
  - principal angles between the bootstrapped COVID vs Sepsis mortality-weight SUBSPACES
    (B donor-resamples per disease -> top-k PC subspace each) -> mean cos(principal angle)

Then correlate subspace alignment with the published forward cross-disease AUC.

Output: 05_RESULTS/option_c_primary_benchmark/subspace_alignment.csv
"""
import h5py
import numpy as np
import pandas as pd
import os
from numpy.linalg import svd
from sklearn.linear_model import LogisticRegression

ROOT = '/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c'
OUT_DIR = '/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_primary_benchmark'

EMBEDDINGS = [
    ('STATE',      'state_covid_canonical.h5ad',     'X_state'),
    ('UCE',        'uce_covid_canonical.h5ad',       'X_uce'),
    ('scPRINT-2',  'scprint2_covid_canonical.h5ad',  'X_scprint2'),
    ('Geneformer', 'geneformer_covid_canonical.h5ad', 'X_geneformer'),
    ('scGPT',      'scgpt_covid_canonical.h5ad',     'X_scgpt'),
]
# published forward (COVID->Sepsis) transfer AUC (claim C002)
FWD_AUC = {'STATE': 0.798, 'UCE': 0.807, 'scPRINT-2': 0.855, 'Geneformer': 0.773, 'scGPT': 0.352}
B = 40        # donor bootstraps per disease for the subspace
K = 5         # subspace dimension (top-k weight PCs)
KEEP_SEV = ['fatal', 'mild', 'moderate', 'recovered', 'severe']
RNG = np.random.default_rng(42)


def read_cat(g):
    if isinstance(g, h5py.Group) and 'categories' in g:
        return g['categories'][:].astype(str)[g['codes'][:]]
    return g[:].astype(str)


def fit_w(X, y):
    clf = LogisticRegression(C=1.0, solver='lbfgs', class_weight='balanced', max_iter=1000)
    clf.fit(X, y)
    w = clf.coef_.ravel()
    return w / (np.linalg.norm(w) + 1e-12)


def boot_subspace(X, y, sid, B, k, rng):
    donors = np.unique(sid)
    idx_by = {d: np.where(sid == d)[0] for d in donors}
    W = []
    for _ in range(B):
        pick = rng.choice(donors, len(donors), replace=True)
        rows = np.concatenate([idx_by[d] for d in pick])
        yy = y[rows]
        if len(np.unique(yy)) < 2:
            continue
        W.append(fit_w(X[rows], yy))
    W = np.array(W)                      # (B, d)
    # top-k right singular vectors = principal directions of the weight cloud
    _, _, vt = svd(W - W.mean(0), full_matrices=False)
    return vt[:k]                        # (k, d) orthonormal basis


def principal_angle_cosines(A, B):
    """A,B: (k,d) orthonormal row-bases. Returns cos of principal angles."""
    M = A @ B.T                          # (k,k)
    s = svd(M, compute_uv=False)
    return np.clip(s, -1, 1)


def main():
    rows = []
    for fm, fn, key in EMBEDDINGS:
        path = os.path.join(ROOT, fn)
        if not os.path.exists(path):
            print(f'SKIP {fm}', flush=True)
            continue
        with h5py.File(path, 'r') as h:
            X = h[f'obsm/{key}'][:]
            disease = read_cat(h['obs/Disease'])
            sev = read_cat(h['obs/covid_severity'])
            sid = read_cat(h['obs/Sample_ID'])
        keep = np.isin(disease, ['COVID-19', 'Sepsis']) & np.isin(sev, KEEP_SEV)
        X, disease, sev, sid = X[keep], disease[keep], sev[keep], sid[keep]
        y = (sev == 'fatal').astype(int)
        cm = disease == 'COVID-19'
        sm = disease == 'Sepsis'

        w_cov = fit_w(X[cm], y[cm])
        w_sep = fit_w(X[sm], y[sm])
        cos_1d = float(w_cov @ w_sep)

        Ac = boot_subspace(X[cm], y[cm], sid[cm], B, K, RNG)
        As = boot_subspace(X[sm], y[sm], sid[sm], B, K, RNG)
        pac = principal_angle_cosines(Ac, As)

        rows.append({
            'fm': fm, 'dim': X.shape[1],
            'cos_1d': round(cos_1d, 4),
            'subspace_mean_cos': round(float(np.mean(pac)), 4),
            'subspace_top_cos': round(float(pac.max()), 4),
            'fwd_transfer_auc': FWD_AUC.get(fm, np.nan),
        })
        print(f'{fm}: cos_1d={cos_1d:+.3f}  subspace_mean_cos={np.mean(pac):.3f}  '
              f'fwd_auc={FWD_AUC.get(fm)}', flush=True)

    df = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out = os.path.join(OUT_DIR, 'subspace_alignment.csv')
    df.to_csv(out, index=False)
    # does subspace alignment track transfer?
    if df['fwd_transfer_auc'].notna().sum() >= 3:
        r = np.corrcoef(df['subspace_mean_cos'], df['fwd_transfer_auc'])[0, 1]
        print(f'\nPearson r(subspace_mean_cos, fwd_transfer_auc) = {r:.3f}', flush=True)
    print(f'\nSAVED: {out}\n{df.to_string(index=False)}', flush=True)


if __name__ == '__main__':
    main()
