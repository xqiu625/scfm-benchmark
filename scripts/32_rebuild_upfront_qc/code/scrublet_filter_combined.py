"""RETIRED — superseded by the Gate 2 safety layer (gate2_run.py and friends).

This was the Phase-1 Scrublet implementation, written before the two-gate
discipline existed. `scrublet_one_channel` here is the origin of two defects
that invalidated every Gate 2 run to date:

  1. `n_recovered = int(n_cells_recovered) if n_cells_recovered is not None
     else n_cells` -- a missing pre-QC count silently became the post-QC count,
     so the expected rate was derived from what QC had already kept. Because
     `expected_doublet_rate` feeds Scrublet's score calibration, the SCORES were
     wrong, not merely the reported numbers.
  2. `adata.obs['predicted_doublet']` was read unconditionally after the
     try/except. Scanpy warns rather than raising when automatic thresholding
     fails, so a failed capture was recorded as a capture full of singlets.

Both are structurally impossible in the replacement: the rate comes only from an
authorized runtime config (gate2_runtime_config.py, fail-closed), and Scrublet's
return value is validated before any call is accepted (gate2_score_core.py),
with per-cell state carried in an explicit enum (gate2_status.py) in which
nothing but `predicted_doublet` authorizes deletion.

The body is removed rather than repaired so no old command line can reproduce
the old semantics. History is in git.
"""
import sys

sys.exit(
    "RETIRED: Phase-1 Scrublet is superseded by gate2_run.py. The rate now comes only from an "
    "authorized runtime config, and no such config exists yet -- see docs/execution-plan.md P0."
)
