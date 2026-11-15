"""Figure S1 -- cell-count association (P1-C). Scatter + lowess + marginal density,
one small panel per FM. donor n_cells (x, log) vs predicted score (y, LR mean)."""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))
from style import MODEL_COLOR, MODEL_ORDER, check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/within_covid"
SPEC_RHO_RANGE = (0.22, 0.29)


def lowess_simple(x, y, frac=0.5):
    """Minimal lowess (no statsmodels dependency): local linear regression."""
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    n = len(xs)
    k = max(3, int(frac * n))
    yfit = np.empty(n)
    for i in range(n):
        dist = np.abs(xs - xs[i])
        idx = np.argsort(dist)[:k]
        w = 1 - (dist[idx] / (dist[idx].max() + 1e-9)) ** 3
        w = np.clip(w, 0, None) ** 3
        X = np.vstack([np.ones(k), xs[idx]]).T
        W = np.diag(w)
        try:
            beta = np.linalg.solve(X.T @ W @ X, X.T @ W @ ys[idx])
            yfit[i] = beta[0] + beta[1] * xs[i]
        except np.linalg.LinAlgError:
            yfit[i] = ys[idx].mean()
    return xs, yfit


def main():
    methods = MODEL_ORDER + ["PCA-200"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    axes = axes.flatten()
    rhos = {}
    for ax, m in zip(axes, methods):
        df = pd.read_csv(os.path.join(DATA, f"donor_predictions_{m}.csv"))
        x = np.log10(df["n_cells"].to_numpy())
        y = df["mean"].to_numpy()
        rho, p = spearmanr(df["n_cells"], y)
        rhos[m] = rho
        color = MODEL_COLOR[m]
        ax.scatter(x, y, s=14, color=color, alpha=0.5, edgecolors="none")
        xs, yfit = lowess_simple(x, y)
        ax.plot(xs, yfit, color="black", lw=1.5)
        ax.set_title(f"{m}  (Spearman ρ={rho:.3f}, p={p:.1e})", fontsize=9)
        ax.set_xlabel("log10(n_cells)")
        ax.set_ylabel("Predicted score (LR mean)")
        divider_ax = ax.inset_axes([0, 1.0, 1, 0.18], sharex=ax)
        divider_ax.hist(x, bins=20, color=color, alpha=0.5)
        divider_ax.axis("off")
    for ax in axes[len(methods):]:
        ax.axis("off")
    fig.tight_layout()
    save_fig(fig, "FigureS1_cellcount")

    fm_rhos = {m: r for m, r in rhos.items() if m in MODEL_ORDER}  # spec anchor = "every FM", excl. PCA-200
    lo_m = min(fm_rhos, key=fm_rhos.get)
    hi_m = max(fm_rhos, key=fm_rhos.get)
    check_close("FigureS1", f"Spearman rho min (across 5 FMs; {lo_m}={fm_rhos[lo_m]:.3f})",
                SPEC_RHO_RANGE[0], fm_rhos[lo_m], DATA, tol=0.02,
                note=f"P1-C anchor claims 0.22-0.29 'for every FM'; {lo_m} is actually "
                     f"{fm_rhos[lo_m]:.3f}, below the stated floor.")
    check_close("FigureS1", f"Spearman rho max (across 5 FMs; {hi_m}={fm_rhos[hi_m]:.3f})",
                SPEC_RHO_RANGE[1], fm_rhos[hi_m], DATA, tol=0.02)


if __name__ == "__main__":
    main()
