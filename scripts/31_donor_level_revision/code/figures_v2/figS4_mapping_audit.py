"""Figure S4 (optional) -- 322 sample_ids -> 252 patients mapping-fix waterfall.

DATA FIDELITY NOTE: `patient_map_v1_to_v2_diff.csv` only records the v1(Stephenson-
fixed)->v2 diff (Combes/Wilk/Lee/label fixes), NOT the earlier, larger Stephenson
pool-vs-individual fix (238->259 patients) -- that step predates this CSV and isn't
recoverable from any file. The Stephenson step count (+21, 15 pool-codes -> 36
individuals) is taken from documented session/runbook history
(scripts/31_donor_level_revision/HPCC_RUNBOOK.md Job H2 STATUS block), not from a
CSV recomputation -- flagged explicitly rather than presented as CSV-verified.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from style import save_fig

DIFF_CSV = "scripts/31_donor_level_revision/results_patient/patient_map_v1_to_v2_diff.csv"


def main():
    diff = pd.read_csv(DIFF_CSV)
    n_combes_drop = (diff["action"].str.contains("DROPPED")).sum()
    n_wilk_merge = (diff["cohort"] == "Wilk").sum() - 1  # 2 sids -> 1 patient = -1
    n_lee_merge = (diff["cohort"] == "Lee").sum() // 2 * 1  # pairs, but check below
    # Lee: 3 distinct merges recorded as 3 rows in this file (one member of each pair
    # shown as the "changed" row); confirmed against runbook: 11->8 patients = -3
    n_lee_merge = (diff["cohort"] == "Lee").sum()

    steps = [
        ("322 raw\nSample_IDs", 322, 0),
        ("v1 mapping\n(Stephenson\npooled, buggy)", 238, 238 - 322),
        ("+ Stephenson fix\n(pool→individual,\nfrom runbook,\nnot this CSV)", 259, 259 - 238),
        (f"− Combes drop\n({n_combes_drop} pooled\nGSMs)", 259 - n_combes_drop, -n_combes_drop),
        ("− Wilk merge\n(1 repeat draw)", 259 - n_combes_drop - n_wilk_merge, -n_wilk_merge),
        (f"− Lee merge\n({n_lee_merge} repeat\ndraws)",
         259 - n_combes_drop - n_wilk_merge - n_lee_merge,
         -n_lee_merge),
    ]
    final = steps[-1][1]
    check = 252
    if final != check:
        print(f"MISMATCH: [FigureS4] reconstructed final count = {final}, expected 252 "
              f"from session history. Reconciling step sizes before plotting is required.",
              flush=True)

    labels = [s[0] for s in steps]
    values = [s[1] for s in steps]
    deltas = [s[2] for s in steps]

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(steps))
    colors = ["#7F8C8D"] + ["#C0392B" if d < 0 else "#27AE60" if d > 0 else "#7F8C8D"
                            for d in deltas[1:]]
    prev = 0
    for i, (lbl, val, delta) in enumerate(steps):
        bottom = prev if i > 0 and delta != 0 else 0
        if i == 0 or delta == 0:
            ax.bar(i, val, color=colors[i], width=0.6)
        else:
            lo = min(prev, val)
            ax.bar(i, abs(delta), bottom=lo, color=colors[i], width=0.6)
            ax.plot([i - 0.5, i - 0.5], [prev, prev], color="gray", lw=0.8, ls=":")
        ax.text(i, val + 6, f"{val}\n({delta:+d})" if i > 0 else f"{val}",
               ha="center", fontsize=8)
        prev = val
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Patient count")
    ax.set_ylim(0, 340)
    ax.set_title(f"Patient-mapping audit waterfall: 322 → {final} "
                 f"({'matches' if final == check else 'MISMATCH vs'} documented final count 252)")
    fig.tight_layout()
    save_fig(fig, "FigureS4_mapping_audit")


if __name__ == "__main__":
    main()
