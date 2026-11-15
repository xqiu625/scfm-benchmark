"""
Gate 2 diagnostic: locate Scrublet's AUTOMATIC threshold inside the observed
score distribution, and compare it against the score that would cut the
per-capture expected doublet rate.

Motivation: the first four completed cohorts detect far fewer candidate
doublets than their expected rate implies (Qiu 3x below, Wilk 210x below).
This script distinguishes the two candidate explanations --
  (a) the data really contain few doublets, versus
  (b) the automatic threshold is placed far out in the tail because the
      simulated-doublet score distribution is not cleanly bimodal --
by reporting where the automatic threshold actually falls.

Diagnostic only: writes nothing, changes no calls.
"""
import numpy as np
import pandas as pd

R = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate2_scrublet"
COHORTS = ["qiu", "cillo", "stephenson", "wilk"]
PAD = " " * 11

print("Where the AUTO threshold sits in the observed score distribution:\n")
for c in COHORTS:
    s = pd.read_csv(f"{R}/{c}/{c}_scrublet_scored.csv", low_memory=False)
    g = pd.read_csv(f"{R}/{c}/{c}_scrublet_diagnostics.csv").set_index("channel_id")
    rows = []
    for ch, sub in s.groupby("channel_id"):
        if ch not in g.index or g.loc[ch, "status"] != "ok":
            continue
        sc = sub.scrublet_score.dropna().to_numpy()
        if len(sc) < 50:
            continue
        thr = float(g.loc[ch, "threshold"])
        exp = float(g.loc[ch, "expected_rate_used"])
        rows.append({"pct_below": 100 * (sc < thr).mean(),
                     "score_at_exp": np.quantile(sc, 1 - exp),
                     "thr": thr, "exp": exp})
    d = pd.DataFrame(rows)
    ratio = d.thr.median() / d.score_at_exp.median()
    print(f"{c:11s} auto thr = percentile {d.pct_below.median():6.2f} of observed scores")
    print(f"{PAD} auto thr med = {d.thr.median():.3f}   "
          f"score cutting the expected rate = {d.score_at_exp.median():.3f}")
    print(f"{PAD} -> auto thr is {ratio:.1f}x higher than the expected-rate cut\n")
