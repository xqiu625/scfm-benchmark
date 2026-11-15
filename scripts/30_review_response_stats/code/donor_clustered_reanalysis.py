#!/usr/bin/env python3
"""
Review-response statistics (roadmap P0/P1: items A, B, C, L, F).

Addresses the peer-review consensus that significance was overstated because
DeLong was computed on cell-weighted OOF predictions (~42k correlated cells from
322 donors treated as independent). Recomputes every headline comparison with the
DONOR as the unit of inference.

Inputs (in ../data/, scp'd from HPCC 05_RESULTS):
  oof_<FM>.csv               cell_id, sample_id, severity, y, p_lr, p_rf   (OOF, COVID+Sepsis subset)
  covid_sepsis_only_mortality.csv   per-fold AUCs (fold_aucs col)
  delong_pairwise.csv               original cell-level DeLong p-values
  cross_disease_fatal_vs_survived_20260504_152836.csv  (+ _20260430 for scGPT)

Outputs (../results/):
  A  donor_bootstrap_auc_ci.csv      per-FM AUC + donor-clustered 95% CI
  A  donor_bootstrap_delta_auc.csv   pairwise dAUC, donor CI, bootstrap p, vs cell-DeLong p, BH-FDR q
  B  fold_mean_ci.csv                per-FM fold-mean AUC + t-CI over 5 folds
  B  nadeau_bengio_paired.csv        corrected paired test on per-fold AUCs (top pairs)
  C  cross_disease_ci_overlap.csv    pairwise CI overlap on forward/reverse transfer
  F  severe_sensitivity.csv          eval-side AUC with vs without severe negatives
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT = os.path.join(HERE, "..", "results")
os.makedirs(OUT, exist_ok=True)

FMS = ["STATE", "scPRINT-2", "Geneformer", "scGPT", "UCE"]
B = 5000
RNG = np.random.default_rng(42)
SCORE = "p_lr"  # headline classifier


def load_oof():
    d = {}
    for fm in FMS:
        df = pd.read_csv(os.path.join(DATA, f"oof_{fm}.csv"))
        df = df.set_index("cell_id")
        d[fm] = df
    return d


def safe_auc(y, p):
    if len(np.unique(y)) < 2:
        return np.nan
    return roc_auc_score(y, p)


def donor_bootstrap_indices(sample_ids):
    """Pre-group cell row-indices by donor for cluster resampling."""
    order = np.argsort(sample_ids, kind="stable")
    sorted_ids = sample_ids[order]
    uniq, starts = np.unique(sorted_ids, return_index=True)
    groups = np.split(order, starts[1:])
    return uniq, groups  # list of arrays of row indices per donor


def cluster_bootstrap_auc(y, p, groups, B, rng):
    """Resample donors with replacement; recompute AUC each draw."""
    n = len(groups)
    aucs = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, n, n)
        idx = np.concatenate([groups[i] for i in pick])
        aucs[b] = safe_auc(y[idx], p[idx])
    return aucs


# ---------------------------------------------------------------- A: per-FM AUC CI
def analysis_A_auc_ci(oof):
    rows = []
    for fm in FMS:
        df = oof[fm]
        y = df["y"].to_numpy()
        p = df[SCORE].to_numpy()
        sid = df["sample_id"].to_numpy().astype(str)
        _, groups = donor_bootstrap_indices(sid)
        boot = cluster_bootstrap_auc(y, p, groups, B, RNG)
        boot = boot[~np.isnan(boot)]
        rows.append({
            "fm": fm, "n_cells": len(y), "n_donors": len(groups),
            "auc_point": round(safe_auc(y, p), 4),
            "auc_boot_mean": round(np.mean(boot), 4),
            "ci_lo": round(np.percentile(boot, 2.5), 4),
            "ci_hi": round(np.percentile(boot, 97.5), 4),
            "ci_width": round(np.percentile(boot, 97.5) - np.percentile(boot, 2.5), 4),
        })
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "donor_bootstrap_auc_ci.csv"), index=False)
    return out


# ------------------------------------------------- A: pairwise dAUC donor CI + FDR
def analysis_A_delta(oof, delong):
    # align pairwise on common cells; resample donors common to both
    dp = pd.read_csv(os.path.join(DATA, "delong_pairwise.csv"))
    dp_map = {(r.model_a, r.model_b): r.p_value for r in dp.itertuples()}
    rows = []
    for i, a in enumerate(FMS):
        for b_ in FMS[i + 1:]:
            da, db = oof[a], oof[b_]
            common = da.index.intersection(db.index)
            ya = da.loc[common, "y"].to_numpy()
            pa = da.loc[common, SCORE].to_numpy()
            pb = db.loc[common, SCORE].to_numpy()
            sid = da.loc[common, "sample_id"].to_numpy().astype(str)
            _, groups = donor_bootstrap_indices(sid)
            n = len(groups)
            deltas = np.empty(B)
            for k in range(B):
                pick = RNG.integers(0, n, n)
                idx = np.concatenate([groups[j] for j in pick])
                deltas[k] = safe_auc(ya[idx], pa[idx]) - safe_auc(ya[idx], pb[idx])
            deltas = deltas[~np.isnan(deltas)]
            frac_gt = np.mean(deltas > 0)
            boot_p = 2 * min(frac_gt, 1 - frac_gt)
            cell_p = dp_map.get((a, b_), dp_map.get((b_, a), np.nan))
            rows.append({
                "model_a": a, "model_b": b_, "n_align": len(common),
                "delta_auc": round(np.mean(deltas), 4),
                "ci_lo": round(np.percentile(deltas, 2.5), 4),
                "ci_hi": round(np.percentile(deltas, 97.5), 4),
                "crosses_zero": bool(np.percentile(deltas, 2.5) < 0 < np.percentile(deltas, 97.5)),
                "donor_boot_p": round(boot_p, 4),
                "cell_delong_p": cell_p,
            })
    out = pd.DataFrame(rows)
    # BH-FDR on donor bootstrap p (the honest family)
    p = out["donor_boot_p"].to_numpy()
    order = np.argsort(p)
    m = len(p)
    q = np.empty(m)
    prev = 1.0
    for rank, idx in enumerate(order[::-1]):
        i_rank = m - rank
        val = min(prev, p[idx] * m / i_rank)
        q[idx] = val
        prev = val
    out["donor_boot_q_bh"] = np.round(q, 4)
    out["sig_after_fdr_0.05"] = out["donor_boot_q_bh"] < 0.05
    out.to_csv(os.path.join(OUT, "donor_bootstrap_delta_auc.csv"), index=False)
    return out


# ---------------------------------------------------- B: fold-mean CI + Nadeau-Bengio
def parse_folds(s):
    return np.array([float(x) for x in str(s).split(";")])


def analysis_B(perfold):
    df = perfold[perfold["classifier"] == "LR"].copy()
    folds = {r.fm: parse_folds(r.fold_aucs) for r in df.itertuples()}
    # map raw fm names to our labels
    name = {"STATE": "STATE", "scPRINT-2": "scPRINT-2", "Geneformer": "Geneformer",
            "scGPT": "scGPT", "UCE": "UCE"}
    rows = []
    for raw, lbl in name.items():
        if raw not in folds:
            continue
        f = folds[raw]
        k = len(f)
        m, sd = np.mean(f), np.std(f, ddof=1)
        tcrit = stats.t.ppf(0.975, k - 1)
        half = tcrit * sd / np.sqrt(k)
        rows.append({"fm": lbl, "fold_mean": round(m, 4), "fold_sd": round(sd, 4),
                     "ci_lo": round(m - half, 4), "ci_hi": round(m + half, 4),
                     "folds": ";".join(f"{x:.3f}" for x in f)})
    ci = pd.DataFrame(rows).sort_values("fold_mean", ascending=False)
    ci.to_csv(os.path.join(OUT, "fold_mean_ci.csv"), index=False)

    # Nadeau-Bengio corrected paired t on per-fold AUC differences (k=5; n_test/n_train=1/(k-1))
    order = ci["fm"].tolist()
    nb = []
    corr = (1.0 / 5) + (1.0 / 4)  # 1/k + rho/(1-rho) with rho=1/k -> n_test/n_train=1/(k-1)
    for i, a in enumerate(order):
        for b_ in order[i + 1:]:
            ra = {v: k for k, v in name.items()}[a]
            rb = {v: k for k, v in name.items()}[b_]
            d = folds[ra] - folds[rb]
            md, vd = np.mean(d), np.var(d, ddof=1)
            t = md / np.sqrt(corr * vd) if vd > 0 else np.nan
            p = 2 * (1 - stats.t.cdf(abs(t), 4)) if np.isfinite(t) else np.nan
            nb.append({"model_a": a, "model_b": b_, "mean_fold_delta": round(md, 4),
                       "nb_t": round(t, 3) if np.isfinite(t) else np.nan,
                       "nb_p": round(p, 4) if np.isfinite(p) else np.nan,
                       "sig_0.05": bool(p < 0.05) if np.isfinite(p) else False})
    pd.DataFrame(nb).to_csv(os.path.join(OUT, "nadeau_bengio_paired.csv"), index=False)
    return ci, pd.DataFrame(nb)


# ---------------------------------------------------------- C: cross-disease CI overlap
def analysis_C():
    # canonical per claim_traceability C002: _20260504 for the 4 FMs, _20260430 for scGPT only
    main = pd.read_csv(os.path.join(DATA, "cross_disease_fatal_vs_survived_20260504_152836.csv"))
    sc = pd.read_csv(os.path.join(DATA, "cross_disease_fatal_vs_survived_20260430_193915.csv"))
    main = main[main["model"] != "scGPT"]
    sc = sc[sc["model"] == "scGPT"]
    df = pd.concat([main, sc], ignore_index=True)
    df = df[df["classifier"] == "LogReg"]
    df = df.drop_duplicates(subset=["direction", "model"], keep="first")
    rows = []
    for direction in df["direction"].unique():
        sub = df[df["direction"] == direction]
        recs = [(r.model, r.boot_auc_mean, r.boot_ci_lo, r.boot_ci_hi) for r in sub.itertuples()]
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                m1, a1, lo1, hi1 = recs[i]
                m2, a2, lo2, hi2 = recs[j]
                overlap = not (hi1 < lo2 or hi2 < lo1)
                rows.append({"direction": direction, "model_a": m1, "model_b": m2,
                             "auc_a": round(a1, 3), "auc_b": round(a2, 3),
                             "ci_a": f"[{lo1:.3f},{hi1:.3f}]", "ci_b": f"[{lo2:.3f},{hi2:.3f}]",
                             "ci_overlap": overlap})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "cross_disease_ci_overlap.csv"), index=False)
    return out


# ----------------------------------------------------------- F: severe sensitivity (eval-side)
def analysis_F(oof):
    rows = []
    for fm in FMS:
        df = oof[fm]
        y = df["y"].to_numpy()
        p = df[SCORE].to_numpy()
        sev = df["severity"].to_numpy().astype(str)
        full = safe_auc(y, p)
        # exclude severe negatives (severe is a survived/ICU category, y=0)
        keep = ~((sev == "severe") & (y == 0))
        excl = safe_auc(y[keep], p[keep])
        rows.append({"fm": fm, "auc_severe_included": round(full, 4),
                     "auc_severe_excluded": round(excl, 4),
                     "delta": round(excl - full, 4),
                     "n_severe_neg": int(((sev == "severe") & (y == 0)).sum())})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUT, "severe_sensitivity.csv"), index=False)
    return out


def main():
    oof = load_oof()
    perfold = pd.read_csv(os.path.join(DATA, "covid_sepsis_only_mortality.csv"))
    delong = pd.read_csv(os.path.join(DATA, "delong_pairwise.csv"))

    print("=== A: per-FM donor-clustered AUC CI ===")
    print(analysis_A_auc_ci(oof).to_string(index=False))
    print("\n=== A: pairwise dAUC donor CI vs cell-DeLong + BH-FDR ===")
    print(analysis_A_delta(oof, delong).to_string(index=False))
    print("\n=== B: fold-mean t-CI ===")
    ci, nb = analysis_B(perfold)
    print(ci.to_string(index=False))
    print("\n=== B: Nadeau-Bengio corrected paired test ===")
    print(nb.to_string(index=False))
    print("\n=== C: cross-disease CI overlap ===")
    print(analysis_C().to_string(index=False))
    print("\n=== F: severe-negative eval-side sensitivity ===")
    print(analysis_F(oof).to_string(index=False))
    print(f"\nSAVED -> {OUT}")


if __name__ == "__main__":
    main()
