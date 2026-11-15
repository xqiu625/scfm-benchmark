#!/usr/bin/env python3
"""Build the Fig3 cross-disease snapshot from the verified C002 evidence.

The figure-data snapshot `scripts/26_figures/data/cross_disease_severe_excluded.csv`
must carry the **manuscript primary FM set** (STATE / UCE / scPRINT-2 / Geneformer / scGPT),
matching claim C002 (`docs/claim-traceability.md`). The previous snapshot was a stale cut of
the cross-disease run that kept the dropped models (TranscriptFormer, STACK) and was missing
Geneformer and scGPT.

Provenance (severe-excluded cross-disease, `tab:crossdisease`):
  - STATE / UCE / scPRINT-2 / Geneformer  <- cross_disease_fatal_vs_survived_20260504_152836.csv
  - scGPT                                  <- cross_disease_fatal_vs_survived_20260430_193915.csv
    (scGPT was dropped from the 2026-05-04 rerun; its published anti-correlated forward AUC
     0.352 lives in the 2026-04-30 file — see C002.)

This is a pure row-selection from the committed evidence — no recomputation, no fabrication.
"""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
EVID = ROOT / "scripts/06_batch_generalization/results/cross_disease_fatal_vs_survived"
OUT = ROOT / "scripts/26_figures/data/cross_disease_severe_excluded.csv"

SRC_MAIN = EVID / "cross_disease_fatal_vs_survived_20260504_152836.csv"
SRC_SCGPT = EVID / "cross_disease_fatal_vs_survived_20260430_193915.csv"

# Manuscript primary set, in display order. Each maps to (source file, model-name-in-file).
KEEP = [
    ("STATE", SRC_MAIN), ("UCE", SRC_MAIN), ("scPRINT2", SRC_MAIN),
    ("Geneformer", SRC_MAIN), ("scGPT", SRC_SCGPT),
]


def main() -> None:
    main_df = pd.read_csv(SRC_MAIN)
    scgpt_df = pd.read_csv(SRC_SCGPT)
    by_src = {SRC_MAIN: main_df, SRC_SCGPT: scgpt_df}

    parts = []
    for model, src in KEEP:
        rows = by_src[src][by_src[src]["model"] == model]
        if rows.empty:
            raise SystemExit(f"ERROR: model {model!r} not found in {src.name}")
        parts.append(rows)

    out = pd.concat(parts, ignore_index=True)
    out.to_csv(OUT, index=False)
    print(f"Wrote {OUT.relative_to(ROOT)}: {len(out)} rows, "
          f"models = {list(out['model'].unique())}")


if __name__ == "__main__":
    main()
