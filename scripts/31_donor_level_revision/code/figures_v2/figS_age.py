"""Figure S-Age -- GSE279452 (Chongqing sepsis cohort) age distribution and
age-stratified mortality (evidence for external-review point A2: cohort is
adult-only, age structure documented).

Panel A: age histogram (10-yr bins), fatal (mort28=1, orange) vs. survivor
  (gray), with a dashed line at the median (67).
Panel B: fatal-rate-by-age-stratum (<40 / 40-60 / >60), point + Wilson 95% CI,
  n and fatal count annotated per stratum.

DATA NOTE (verified, not assumed): 2 of 92 patients have missing mort28
(GSM8571097 age 81, GSM8571170 age 77 -- both fall in the >60 stratum). The
">60: 10/56" figure in the spec counts these 2 unknown-outcome patients in the
stratum DENOMINATOR (56 total) while the fatal numerator (10) is necessarily
from the 54 patients with a known outcome -- confirmed by direct inspection,
not a data error. This is stated explicitly in Panel B rather than silently
dropped or silently treated as non-fatal.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import check_close, save_fig

DATA = "scripts/31_donor_level_revision/results_patient_v2/monocyte_matched"
STRATA = [(0, 40, "<40"), (40, 60, "40-60"), (60, 200, ">60")]


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n
    denom = 1 + z**2 / n
    center = (p + z**2 / (2 * n)) / denom
    half = (z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / denom
    return (max(0, center - half), min(1, center + half))


def panel_a(ax, df, median_age):
    fatal = df[df["mort28"] == 1]["age"].dropna()
    surv = df[df["mort28"] == 0]["age"].dropna()
    bins = np.arange(10, 100, 10)
    ax.hist([surv, fatal], bins=bins, stacked=True, color=["#95A5A6", "#E67E22"],
            label=["Survivor", "Fatal (mort28=1)"], edgecolor="white", linewidth=0.5)
    ax.axvline(median_age, color="black", ls="--", lw=1.2)
    ax.text(median_age + 1.5, ax.get_ylim()[1] * 0.9, f"median {median_age:.0f}",
            fontsize=8, color="black")
    ax.set_xlabel("Age (years)")
    ax.set_ylabel("N patients")
    ax.set_title("A", loc="left", fontweight="bold")
    ax.legend(fontsize=8, frameon=False)


def panel_b(ax, df):
    labels, ns, ks, ci_los, ci_his = [], [], [], [], []
    for lo, hi, label in STRATA:
        sub = df[(df["age"] >= lo) & (df["age"] < hi)]
        n = len(sub)
        k = int(sub["mort28"].sum())  # NaN excluded from sum automatically
        n_known = int(sub["mort28"].notna().sum())
        rate = k / n_known if n_known else np.nan
        ci_lo, ci_hi = wilson_ci(k, n_known)
        labels.append(label); ns.append(n); ks.append(k); ci_los.append(ci_lo); ci_his.append(ci_hi)
        rate_full = rate  # rate computed on known-outcome denominator
    x = np.arange(len(STRATA))
    rates = [k / max(1, df[(df["age"] >= lo) & (df["age"] < hi)]["mort28"].notna().sum())
             for (lo, hi, _), k in zip(STRATA, ks)]
    ax.errorbar(x, rates, yerr=[np.array(rates) - np.array(ci_los), np.array(ci_his) - np.array(rates)],
                fmt="o", color="#E67E22", ecolor="#E67E22", capsize=4, markersize=9)
    for i, (label, n, k) in enumerate(zip(labels, ns, ks)):
        ax.text(i, -0.08, f"n={n}\n{k} fatal", ha="center", va="top", fontsize=8,
                transform=ax.get_xaxis_transform())
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(-0.05, 0.9)
    ax.set_ylabel("28-day mortality rate")
    ax.set_title("B", loc="left", fontweight="bold")
    ax.text(0.02, 0.98,
            ">60 stratum: n=56 total (10/56 fatal per cohort table); 2 of these\n"
            "56 have unknown mort28 status, excluded from the rate denominator\n"
            "(rate computed on 54 patients with known outcome).",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.5, color="dimgray")


def main():
    df = pd.read_csv(os.path.join(DATA, "gse279452_age_final.csv"))

    median_age = df["age"].median()
    check_close("FigureS_Age", "median age", 67, median_age,
                os.path.join(DATA, "gse279452_age_final.csv"))
    check_close("FigureS_Age", "min age", 19, df["age"].min(),
                os.path.join(DATA, "gse279452_age_final.csv"))
    check_close("FigureS_Age", "max age", 93, df["age"].max(),
                os.path.join(DATA, "gse279452_age_final.csv"))

    # stratum n/fatal verification -- spec anchors, with the >60 unknown-outcome
    # nuance resolved (see module docstring)
    spec_strata = {"<40": (6, 0), "40-60": (30, 7), ">60": (56, 10)}
    for lo, hi, label in STRATA:
        sub = df[(df["age"] >= lo) & (df["age"] < hi)]
        n, k = len(sub), int(sub["mort28"].sum())
        spec_n, spec_k = spec_strata[label]
        check_close("FigureS_Age", f"{label} n", spec_n, n,
                    os.path.join(DATA, "gse279452_age_final.csv"), tol=0.5)
        check_close("FigureS_Age", f"{label} fatal", spec_k, k,
                    os.path.join(DATA, "gse279452_age_final.csv"), tol=0.5)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    panel_a(axes[0], df, median_age)
    panel_b(axes[1], df)
    fig.tight_layout()
    save_fig(fig, "FigureS_Age_gse279452")


if __name__ == "__main__":
    main()
