#!/usr/bin/env python3
"""
Revision P1-A: null model — can we predict mortality from cell counts and/or
disease identity alone?

Adam 2026-08-31: "You could almost predict outcome from cell counts alone or
from identifying a disease." This trains logistic regression on ONLY
  (a) log1p(n_cells_per_donor)
  (b) log1p(n_cells_per_donor) + disease one-hot (COVID-19 / Sepsis)
under the SAME protocol as the FMs (5-fold StratifiedGroupKFold,
shuffle=True, random_state=42, Sample_ID groups, LR C=1.0 balanced),
and reports donor-level OOF AUC with a donor bootstrap CI.

Two training variants are reported:
  cell-expanded : features repeated per cell, trained on cells (mirrors the FM
                  pipeline exactly), OOF cell scores aggregated by donor mean
  donor-level   : one row per donor (the pure null; n = 322)

Inputs — one of:
  --donor-meta CSV   sample_id,disease,n_cells,y   (canonical; build on HPCC
                     with --make-donor-meta from a canonical option_c h5ad)
  --make-donor-meta H5AD --donor-meta-out CSV
                     extract donor_meta from a canonical embedding h5ad
                     (obs only, h5py; no X loaded)
  --local-oof CSV    [PRELIMINARY] build from a local oof_<FM>.csv; disease is
                     NOT available in this mode, so only model (a) runs

Outputs (../results/):
  null_model_donor_meta.csv      the donor table used (for the manuscript supplement)
  null_model_auc.csv             one row per (feature_set, variant): OOF AUC + CI

⚠ Manuscript rule (todo.md P1-A): if the (b) cell-expanded AUC is >= 0.7, the
null result becomes a HEADLINE result and the paper framing is rebuilt around
it — do not bury it in supplement.
"""
import argparse
import os

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

RNG = np.random.default_rng(42)


def safe_auc(y, p):
    return roc_auc_score(y, p) if len(np.unique(y)) > 1 else np.nan


def boot_auc_ci(y, p, B=5000, rng=RNG):
    n = len(y)
    aucs = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, n)
        aucs[b] = safe_auc(y[idx], p[idx])
    aucs = aucs[~np.isnan(aucs)]
    return np.percentile(aucs, [2.5, 97.5])


def make_donor_meta_from_h5ad(h5ad_path, out_csv):
    """obs-only extraction; mirrors stage-27 filters (disease-matched subset)."""
    import h5py

    def read_cat(g):
        if isinstance(g, h5py.Group) and "categories" in g:
            return g["categories"][:].astype(str)[g["codes"][:]]
        return g[:].astype(str)

    KEEP_DISEASE = ["COVID-19", "Sepsis"]
    KEEP_SEVERITY = ["fatal", "mild", "moderate", "recovered", "severe"]
    with h5py.File(h5ad_path, "r") as h:
        disease = read_cat(h["obs/Disease"])
        severity = read_cat(h["obs/covid_severity"])
        sample_id = read_cat(h["obs/Sample_ID"])
        extra = {}
        for col in h["obs"].keys():
            if col == "Data_NO" or "timepoint" in col.lower() or "time_point" in col.lower():
                extra[col] = read_cat(h[f"obs/{col}"])
        if extra:
            print(f"extra obs columns exported: {sorted(extra)}", flush=True)
    keep = np.isin(disease, KEEP_DISEASE) & np.isin(severity, KEEP_SEVERITY)
    df = pd.DataFrame({"sample_id": sample_id[keep],
                       "disease": disease[keep],
                       "y": (severity[keep] == "fatal").astype(int),
                       **{k: v[keep] for k, v in extra.items()}})
    meta = (df.groupby("sample_id")
              .agg(**{"disease": ("disease", "first"), "n_cells": ("y", "size"),
                      "y": ("y", "max"),
                      **{k: (k, "first") for k in extra}})
              .reset_index())
    meta.to_csv(out_csv, index=False)
    print(f"donor_meta: {len(meta)} donors, {int(meta['y'].sum())} fatal "
          f"-> {out_csv}", flush=True)
    return meta


def features(meta, with_disease):
    X = np.log1p(meta["n_cells"].to_numpy()).reshape(-1, 1)
    if with_disease:
        d = pd.get_dummies(meta["disease"], drop_first=True).to_numpy(float)
        X = np.hstack([X, d])
    return X


def run_variant(meta, with_disease, variant, seed=42):
    """OOF predictions under the stage-27 CV protocol for one seed."""
    y = meta["y"].to_numpy()
    sid = meta["sample_id"].to_numpy()
    if variant == "cell-expanded":
        # expand donor rows to per-cell rows (features constant within donor)
        counts = meta["n_cells"].to_numpy()
        Xd = features(meta, with_disease)
        X = np.repeat(Xd, counts, axis=0)
        yc = np.repeat(y, counts)
        gc = np.repeat(sid, counts)
        oof = np.full(len(yc), np.nan)
        skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        for tr, te in skf.split(X, yc, groups=gc):
            clf = LogisticRegression(C=1.0, solver="lbfgs",
                                     class_weight="balanced", max_iter=1000)
            clf.fit(X[tr], yc[tr])
            oof[te] = clf.predict_proba(X[te])[:, 1]
        # aggregate to donor by mean (same rule as P1-B)
        tmp = pd.DataFrame({"sample_id": gc, "y": yc, "p": oof})
        d = tmp.groupby("sample_id").agg(y=("y", "max"), p=("p", "mean"))
        return d["y"].to_numpy(), d["p"].to_numpy()
    else:  # donor-level
        X = features(meta, with_disease)
        oof = np.full(len(y), np.nan)
        skf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        for tr, te in skf.split(X, y, groups=sid):
            clf = LogisticRegression(C=1.0, solver="lbfgs",
                                     class_weight="balanced", max_iter=1000)
            clf.fit(X[tr], y[tr])
            oof[te] = clf.predict_proba(X[te])[:, 1]
        return y, oof


def run_variant_repeated(meta, with_disease, variant, seeds):
    """CV Protocol v2: repeated-CV, returns per-seed AUC array."""
    aucs = []
    for seed in seeds:
        y, p = run_variant(meta, with_disease, variant, seed=seed)
        aucs.append(safe_auc(y, p))
    return np.array(aucs)


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--donor-meta")
    ap.add_argument("--make-donor-meta")
    ap.add_argument("--donor-meta-out", default=os.path.join(
        here, "..", "results", "null_model_donor_meta.csv"))
    ap.add_argument("--local-oof")
    ap.add_argument("--out", default=os.path.join(here, "..", "results"))
    ap.add_argument("--n-seeds", type=int, default=1,
                    help="CV Protocol v2 (2026-09-04): repeated-CV over N seeds "
                         "(42..42+N-1), reported as repeat-mean +/- sd. Default 1 "
                         "= legacy single-split + bootstrap CI, unchanged.")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    has_disease = True
    if args.make_donor_meta:
        meta = make_donor_meta_from_h5ad(args.make_donor_meta, args.donor_meta_out)
    elif args.donor_meta:
        meta = pd.read_csv(args.donor_meta)
    elif args.local_oof:
        df = pd.read_csv(args.local_oof)
        meta = (df.groupby("sample_id")
                  .agg(n_cells=("y", "size"), y=("y", "max"))
                  .reset_index())
        has_disease = False
        print("PRELIMINARY mode: no disease column; cell-count-only model",
              flush=True)
    else:
        ap.error("one of --donor-meta / --make-donor-meta / --local-oof required")

    meta.to_csv(os.path.join(args.out, "null_model_donor_meta.csv"), index=False)

    rows = []
    fsets = [("log_cell_count", False)]
    if has_disease:
        fsets.append(("log_cell_count+disease", True))

    if args.n_seeds <= 1:
        for name, wd in fsets:
            for variant in ["cell-expanded", "donor-level"]:
                y, p = run_variant(meta, wd, variant, seed=42)
                lo, hi = boot_auc_ci(y, p)
                rows.append({"feature_set": name, "variant": variant,
                             "n_donors": len(y), "n_fatal_donors": int(y.sum()),
                             "auc_oof_donor": round(safe_auc(y, p), 4),
                             "ci_lo": round(lo, 4), "ci_hi": round(hi, 4)})
                print(f"{name:28s} {variant:14s} AUC={rows[-1]['auc_oof_donor']:.4f} "
                      f"[{lo:.4f}, {hi:.4f}]", flush=True)
        path = os.path.join(args.out, "null_model_auc.csv")
    else:
        from repeated_cv import DEFAULT_SEEDS, summarize_repeated
        seeds = DEFAULT_SEEDS[:args.n_seeds]
        for name, wd in fsets:
            for variant in ["cell-expanded", "donor-level"]:
                aucs = run_variant_repeated(meta, wd, variant, seeds)
                summ = summarize_repeated(aucs)
                rows.append({"feature_set": name, "variant": variant,
                             "n_donors": len(meta), "n_fatal_donors": int(meta["y"].sum()),
                             "n_seeds": summ["n"],
                             "auc_repeat_mean": round(summ["mean"], 4),
                             "auc_repeat_sd": round(summ["sd"], 4)})
                print(f"{name:28s} {variant:14s} "
                      f"AUC(repeat-mean+/-sd)={summ['mean']:.4f}+/-{summ['sd']:.4f} "
                      f"(n_seeds={summ['n']})", flush=True)
        path = os.path.join(args.out, "null_model_auc_repeated.csv")

    out = pd.DataFrame(rows)
    out.to_csv(path, index=False)
    print(f"SAVED -> {path}", flush=True)


if __name__ == "__main__":
    main()
