#!/usr/bin/env python3
"""Multi-classifier comparison on Option C 5-FM embeddings (mortality binary).

For each FM, evaluate 4 classifiers under identical 5-fold patient-level CV:
  1. LR        — Logistic regression (baseline; matches primary benchmark)
  2. RF        — Random forest (baseline; matches primary benchmark)
  3. DeepMLP   — 3-layer MLP with BatchNorm + Dropout, AdamW + early stopping
  4. ResidMLP  — Residual MLP with skip connections

Tests whether non-linear classifiers unlock additional performance from frozen
FM embeddings on the mortality task (legacy paper showed +6.3 pp for STATE
on the severity task with DeepMLP vs LR).

Output:
  multiclassifier_summary.csv   model × classifier → mean AUC, std, BalAcc, n_cells
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import StandardScaler

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("multiclassifier")

OBSM_KEYS = {
    "state": ["X_state", "X_emb", "X_STATE"],
    "uce": ["X_uce", "X_uce_4layer"],
    "geneformer": ["X_geneformer", "X_Geneformer"],
    "transcriptformer": ["X_transcriptformer", "X_TranscriptFormer"],
    "scprint2": ["X_scprint2", "X_scprint", "X_scPRINT", "X_scPRINT2", "scprint_emb"],
    "stack": ["X_stack", "X_STACK"],
    "scgpt": ["X_scgpt", "X_scGPT"],
}
N_FOLDS = 5
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def find_embedding(adata, model_key: str):
    """Find an FM embedding by known obsm keys. Refuses the previous unsafe
    fallback (first 2D obsm) -- that could silently grab UMAP coords or other
    auxiliary arrays. Raises a clear error if no known key is present.
    Mirrors the policy in primary_benchmark_option_c.py."""
    keys = OBSM_KEYS.get(model_key, [])
    for k in keys:
        if k in adata.obsm:
            return np.asarray(adata.obsm[k]), k
    raise RuntimeError(
        f"No known embedding key for model_key={model_key!r}. "
        f"Tried {keys}; available obsm: {list(adata.obsm.keys())}. "
        f"If a new FM is added, register its obsm key(s) in OBSM_KEYS."
    )


def build_mortality(severity: np.ndarray):
    sev = np.asarray([str(s).lower() for s in severity])
    fatal = sev == "fatal"
    survived = np.isin(sev, ["mild", "moderate", "recovered", "severe"])
    valid = fatal | survived
    y = np.zeros(len(sev), dtype=np.int8)
    y[fatal] = 1
    return y, valid


# ============================================================================
# DeepMLP
# ============================================================================
class DeepMLP(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 512, dropout: float = 0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden),
            nn.BatchNorm1d(hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden // 2),
            nn.BatchNorm1d(hidden // 2),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden // 2, 2),
        )

    def forward(self, x):
        return self.net(x)


class ResidualMLP(nn.Module):
    def __init__(self, in_dim: int, hidden: int = 512, dropout: float = 0.3):
        super().__init__()
        self.proj = nn.Linear(in_dim, hidden)
        self.bn0 = nn.BatchNorm1d(hidden)
        self.block1 = nn.Sequential(
            nn.Linear(hidden, hidden), nn.BatchNorm1d(hidden), nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.BatchNorm1d(hidden),
        )
        self.block2 = nn.Sequential(
            nn.Linear(hidden, hidden), nn.BatchNorm1d(hidden), nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden, hidden), nn.BatchNorm1d(hidden),
        )
        self.head = nn.Linear(hidden, 2)

    def forward(self, x):
        h = torch.relu(self.bn0(self.proj(x)))
        h = torch.relu(h + self.block1(h))
        h = torch.relu(h + self.block2(h))
        return self.head(h)


def _seed_torch(seed: int = 42):
    """Make the torch classifier heads deterministic run-to-run.

    Without this the DeepMLP/ResidualMLP columns of tab:multiclassifier change on
    every invocation (weight init, dropout masks and shuffling are all unseeded),
    so the published values cannot be reproduced from the committed code.
    """
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.use_deterministic_algorithms(True, warn_only=True)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    np.random.seed(seed)


def train_torch_clf(model_cls, in_dim, X_tr, y_tr, X_te,
                    epochs: int = 30, batch_size: int = 256, lr: float = 1e-3,
                    weight_decay: float = 1e-4, patience: int = 5):
    """Train a torch classifier with class-balanced loss + early stopping on
    a held-out 10% of the training fold. Returns (proba_pos, pred)."""
    _seed_torch(42)
    rng = np.random.RandomState(42)
    perm = rng.permutation(len(y_tr))
    n_val = max(int(0.1 * len(y_tr)), 100)
    val_idx = perm[:n_val]
    tr_idx = perm[n_val:]
    Xt_tr = torch.tensor(X_tr[tr_idx], dtype=torch.float32, device=DEVICE)
    yt_tr = torch.tensor(y_tr[tr_idx], dtype=torch.long, device=DEVICE)
    Xt_val = torch.tensor(X_tr[val_idx], dtype=torch.float32, device=DEVICE)
    yt_val = torch.tensor(y_tr[val_idx], dtype=torch.long, device=DEVICE)
    Xt_te = torch.tensor(X_te, dtype=torch.float32, device=DEVICE)

    # Class weights for balanced loss
    n_pos = max(int((y_tr[tr_idx] == 1).sum()), 1)
    n_neg = max(int((y_tr[tr_idx] == 0).sum()), 1)
    cw = torch.tensor([len(tr_idx) / (2 * n_neg), len(tr_idx) / (2 * n_pos)],
                      dtype=torch.float32, device=DEVICE)

    net = model_cls(in_dim).to(DEVICE)
    opt = torch.optim.AdamW(net.parameters(), lr=lr, weight_decay=weight_decay)
    loss_fn = nn.CrossEntropyLoss(weight=cw)

    best_val_loss = float("inf")
    best_state = None
    epochs_no_improve = 0

    for epoch in range(epochs):
        net.train()
        # Mini-batch
        idx = torch.randperm(len(yt_tr), device=DEVICE)
        for i in range(0, len(idx), batch_size):
            sel = idx[i:i + batch_size]
            opt.zero_grad()
            logits = net(Xt_tr[sel])
            loss = loss_fn(logits, yt_tr[sel])
            loss.backward()
            opt.step()
        # Val
        net.eval()
        with torch.no_grad():
            val_logits = net(Xt_val)
            val_loss = loss_fn(val_logits, yt_val).item()
        if val_loss < best_val_loss - 1e-4:
            best_val_loss = val_loss
            best_state = {k: v.clone() for k, v in net.state_dict().items()}
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= patience:
                break

    if best_state is not None:
        net.load_state_dict(best_state)
    net.eval()
    with torch.no_grad():
        logits = net(Xt_te)
        proba = torch.softmax(logits, dim=1)[:, 1].cpu().numpy()
        pred = logits.argmax(dim=1).cpu().numpy()
    return proba, pred


# ============================================================================
# CV runners
# ============================================================================
def cv_lr(X, y, groups):
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    aucs, baccs = [], []
    for tr, te in skf.split(X, y, groups):
        scaler = StandardScaler().fit(X[tr])
        clf = LogisticRegression(penalty="l2", C=1.0, solver="lbfgs",
                                 max_iter=2000, class_weight="balanced", n_jobs=4)
        clf.fit(scaler.transform(X[tr]), y[tr])
        p = clf.predict_proba(scaler.transform(X[te]))[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        baccs.append(balanced_accuracy_score(y[te], (p >= 0.5).astype(int)))
    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(baccs))


def cv_rf(X, y, groups):
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    aucs, baccs = [], []
    for tr, te in skf.split(X, y, groups):
        clf = RandomForestClassifier(n_estimators=100, n_jobs=8,
                                     class_weight="balanced", random_state=42)
        clf.fit(X[tr], y[tr])
        p = clf.predict_proba(X[te])[:, 1]
        aucs.append(roc_auc_score(y[te], p))
        baccs.append(balanced_accuracy_score(y[te], (p >= 0.5).astype(int)))
    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(baccs))


def cv_torch(model_cls, X, y, groups):
    skf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=42)
    aucs, baccs = [], []
    for fold, (tr, te) in enumerate(skf.split(X, y, groups)):
        scaler = StandardScaler().fit(X[tr])
        Xs_tr = scaler.transform(X[tr]).astype(np.float32)
        Xs_te = scaler.transform(X[te]).astype(np.float32)
        proba, pred = train_torch_clf(model_cls, X.shape[1], Xs_tr, y[tr], Xs_te)
        aucs.append(roc_auc_score(y[te], proba))
        baccs.append(balanced_accuracy_score(y[te], pred))
        log.info("    fold %d AUC=%.4f", fold, aucs[-1])
    return float(np.mean(aucs)), float(np.std(aucs)), float(np.mean(baccs))


def per_fm(fm_h5ad: Path):
    mk = fm_h5ad.stem.split("_")[0].lower()
    log.info("=" * 60)
    log.info("Model: %s", mk)
    log.info("=" * 60)
    a = sc.read_h5ad(fm_h5ad)
    X, key = find_embedding(a, mk)
    if "covid_severity" not in a.obs.columns or "Sample_ID" not in a.obs.columns:
        return None
    y_full, valid = build_mortality(a.obs["covid_severity"].astype(str).values)
    groups = a.obs["Sample_ID"].astype(str).values
    X_v, y_v, g_v = X[valid], y_full[valid], groups[valid]
    log.info("  obsm=%s  n=%d  dim=%d  fatal=%d", key, len(y_v), X_v.shape[1], int(y_v.sum()))
    rows = []

    log.info("  [LR]")
    auc, std, bacc = cv_lr(X_v, y_v, g_v)
    rows.append(dict(model=mk, classifier="LR", auc_mean=auc, auc_std=std, balanced_accuracy=bacc, n=len(y_v)))
    log.info("    AUC=%.4f (+/- %.4f), BalAcc=%.4f", auc, std, bacc)

    log.info("  [RF]")
    auc, std, bacc = cv_rf(X_v, y_v, g_v)
    rows.append(dict(model=mk, classifier="RF", auc_mean=auc, auc_std=std, balanced_accuracy=bacc, n=len(y_v)))
    log.info("    AUC=%.4f (+/- %.4f), BalAcc=%.4f", auc, std, bacc)

    log.info("  [DeepMLP]")
    auc, std, bacc = cv_torch(DeepMLP, X_v, y_v, g_v)
    rows.append(dict(model=mk, classifier="DeepMLP", auc_mean=auc, auc_std=std, balanced_accuracy=bacc, n=len(y_v)))
    log.info("    AUC=%.4f (+/- %.4f), BalAcc=%.4f", auc, std, bacc)

    log.info("  [ResidualMLP]")
    auc, std, bacc = cv_torch(ResidualMLP, X_v, y_v, g_v)
    rows.append(dict(model=mk, classifier="ResidualMLP", auc_mean=auc, auc_std=std, balanced_accuracy=bacc, n=len(y_v)))
    log.info("    AUC=%.4f (+/- %.4f), BalAcc=%.4f", auc, std, bacc)

    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--embedding-dir", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info("Device: %s", DEVICE)
    fm_files = sorted(args.embedding_dir.glob("*_covid_canonical.h5ad"))
    if not fm_files:
        log.error("No FMs"); sys.exit(1)

    all_rows = []
    for fm in fm_files:
        try:
            rows = per_fm(fm)
            if rows is not None:
                all_rows.extend(rows)
        except Exception as e:
            log.exception("Failed %s: %s", fm.name, e)

    df = pd.DataFrame(all_rows)
    out = args.output_dir / "multiclassifier_summary.csv"
    df.to_csv(out, index=False)
    log.info("Wrote %s", out)
    print("\n=== MULTI-CLASSIFIER SUMMARY ===")
    print(df.to_string(index=False))
    # Pivot for easy reading
    pivot = df.pivot(index="model", columns="classifier", values="auc_mean").round(3)
    print("\n=== PIVOT (AUC) ===")
    print(pivot.to_string())


if __name__ == "__main__":
    main()
