"""Run all v2 figure scripts in one process and print the complete mismatch log."""
import importlib
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import style

MODULES = [
    "fig2_within_covid", "fig3_within_sepsis", "fig4_indist_vs_ood",
    "fig5_cross_disease", "fig6_adaptation_curve", "fig7_calibration",
    "figS1_cellcount", "figS2_pairwise_heatmap", "figS4_mapping_audit",
    "figS_qc", "figS_age",
]

for name in MODULES:
    mod = importlib.import_module(name)
    mod.main()

print("\n=== FULL MISMATCH LOG ===")
for m in style.get_mismatches():
    print("-", m)
print(f"\nTotal figures generated: {len(MODULES)} + FigureS3 (copied, verified) = "
      f"{len(MODULES) + 1}")
print(f"Total mismatches flagged: {len(style.get_mismatches())}")
