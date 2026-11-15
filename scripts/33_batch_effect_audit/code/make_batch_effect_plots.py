#!/usr/bin/env python3
"""Renders the PCA/UMAP before/after-correction visual comparison a
2026-10-09 review asked for ("Produce the requested PCA/UMAP visual
comparisons with cell-ID-bound coordinates") -- deferred in earlier rounds
in favor of the coordinate CSVs alone. Reads the cell-ID-bound
`*_umap_coords.csv` files `run_batch_diagnostics.py` writes (both
`uncorrected` and, when Harmony was run, `harmony_corrected` scopes) and
renders one PNG per (lineage, embedding): a 2x2 grid of
(uncorrected, harmony_corrected) x (colored by cohort, colored by outcome).

Does not compute anything new -- pure visualization of already-validated,
already cell-ID-bound coordinates. Run on HPCC per DEVELOPMENT.md.

Usage:
    python make_batch_effect_plots.py --run-dir <diagnostic_run output dir> --out-dir <dir>
"""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger("make_batch_effect_plots")

KNOWN_LINEAGES = ["Monocyte_Myeloid", "Platelet"]  # sorted longest-first below, to anchor the
# otherwise-ambiguous split -- a prior greedy `(?P<lineage>.+)` regex split
# "Monocyte_Myeloid_X_synthFM_uncorrected_umap_coords.csv" as lineage="Monocyte_Myeloid_X",
# embedding="synthFM" (confirmed: the resulting plot title read "Platelet_X / synthFM"
# instead of "Platelet / X_synthFM"). The final OUTPUT filename was unaffected by this
# (string concatenation of the two wrong halves reassembles the same name), but the
# in-figure title text was wrong -- fixed by anchoring on the known lineage names rather
# than a greedy wildcard, since lineage values here are fixed, known strings, unlike
# embedding keys (which include arbitrary FM obsm names).
UMAP_CSV_RE = re.compile(
    r"^(?P<lineage>" + "|".join(sorted(KNOWN_LINEAGES, key=len, reverse=True)) + r")_"
    r"(?P<embedding>.+?)_(?P<scope>uncorrected|harmony_corrected)_umap_coords\.csv$"
)


def discover_pairs(run_dir: Path) -> dict:
    """Groups the coordinate CSVs in run_dir by (lineage, embedding), each
    holding whichever of {uncorrected, harmony_corrected} actually exist --
    a run without --run-harmony will only have 'uncorrected', and this is
    reported, not padded with a fabricated placeholder panel."""
    groups: dict[tuple[str, str], dict[str, Path]] = {}
    for f in sorted(run_dir.glob("*_umap_coords.csv")):
        m = UMAP_CSV_RE.match(f.name)
        if not m:
            log.warning("Skipping unrecognized file name: %s", f.name)
            continue
        key = (m["lineage"], m["embedding"])
        groups.setdefault(key, {})[m["scope"]] = f
    return groups


def plot_one(lineage: str, embedding: str, scope_paths: dict, out_dir: Path) -> Path:
    scopes = [s for s in ("uncorrected", "harmony_corrected") if s in scope_paths]
    fig, axes = plt.subplots(2, len(scopes), figsize=(6 * len(scopes), 11), squeeze=False)

    for col, scope in enumerate(scopes):
        df = pd.read_csv(scope_paths[scope])
        for row, color_col in enumerate(["cohort", "outcome"]):
            ax = axes[row][col]
            categories = sorted(df[color_col].unique())
            cmap = plt.get_cmap("tab10" if len(categories) <= 10 else "tab20")
            for i, cat in enumerate(categories):
                sub = df[df[color_col] == cat]
                ax.scatter(sub["UMAP1"], sub["UMAP2"], s=6, alpha=0.6,
                            color=cmap(i / max(len(categories) - 1, 1)), label=str(cat))
            ax.set_title(f"{scope}, colored by {color_col}\n(n={len(df)} cells, "
                          f"n_cell_ids_unique={df['cell_id'].nunique()})", fontsize=9)
            ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
            ax.legend(fontsize=6, markerscale=2, loc="best", ncol=1 if len(categories) <= 6 else 2)

    fig.suptitle(f"{lineage} / {embedding} -- batch-effect before/after comparison "
                  f"(synthetic fixture validation run)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out_path = out_dir / f"{lineage}_{embedding}_before_after.png"
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run-dir", required=True, type=Path,
                     help="A run_batch_diagnostics.py --out-dir directory to render plots from.")
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    groups = discover_pairs(args.run_dir)
    if not groups:
        raise FileNotFoundError(f"No *_umap_coords.csv files matched the expected naming "
                                  f"pattern in {args.run_dir}")

    written = []
    for (lineage, embedding), scope_paths in groups.items():
        log.info("Plotting %s / %s (scopes present: %s)", lineage, embedding,
                  sorted(scope_paths.keys()))
        out_path = plot_one(lineage, embedding, scope_paths, args.out_dir)
        written.append(out_path)
        log.info("Wrote %s", out_path)

    print(f"\n=== WROTE {len(written)} PLOT(S) ===")
    for p_ in written:
        print(p_)


if __name__ == "__main__":
    main()
