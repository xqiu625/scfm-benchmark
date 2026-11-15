"""
Gate 2 diagnostic: determine, per cohort and per capture, whether the Scrublet
expected-rate prior was actually a PRE-QC recovered count or silently fell back
to the post-Gate-1 scored cell count.

The two Gate 2 entry points obtain the prior differently, so "prior_counts.csv
is missing" does not by itself invalidate a cohort:
  gate2_scrublet_score.py  (combined) -- reads --prior-counts, falls back when absent
  gate2_scrublet_split.py  (split)    -- counts rows in --raw-split-dir/<ch>/cells.csv,
                                         falls back when that file is absent

Both fallbacks are silent, so the only reliable evidence is in the diagnostics
themselves: a capture whose n_cells_recovered_pre_qc equals its
n_cells_scored_post_qc had no real prior (Gate 1 removed nothing for it, or the
prior was the post-QC count).

Diagnostic only: writes nothing, changes no calls.
"""
import glob
import os

import pandas as pd

R = "/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate2_scrublet"

print(f"{'cohort':<14}{'captures':>9}{'prior>scored':>14}{'prior==scored':>15}  verdict")
print("-" * 72)
for d in sorted(glob.glob(f"{R}/*")):
    f = glob.glob(f"{d}/*_scrublet_diagnostics.csv")
    if not f:
        continue
    g = pd.read_csv(f[0])
    name = os.path.basename(d).split(".")[0]
    if not {"n_cells_recovered_pre_qc", "n_cells_scored_post_qc"} <= set(g.columns):
        print(f"{name:<14}{len(g):>9}  (columns missing)")
        continue
    real = int((g.n_cells_recovered_pre_qc > g.n_cells_scored_post_qc).sum())
    same = int((g.n_cells_recovered_pre_qc == g.n_cells_scored_post_qc).sum())
    if real == len(g):
        verdict = "prior real for every capture"
    elif real == 0:
        verdict = "NO real prior -- fell back everywhere"
    else:
        verdict = f"MIXED -- {same} capture(s) without a real prior"
    print(f"{name:<14}{len(g):>9}{real:>14}{same:>15}  {verdict}")
