"""RETIRED (combined layout) — superseded by gate2_run.py.

This script carried the two silent fallbacks that invalidated every Gate 2 run
to date:

    if a.prior_counts and os.path.exists(a.prior_counts):   # missing file -> no prior
    n_recovered = n_cells_recovered if ... is not None else n_cells   # -> post-QC count

Six of eight cohorts were scored against post-Gate-1 cell counts because of
them, and nothing raised. It also recorded a capture whose automatic threshold
failed as a capture full of singlets.

The body is removed rather than fixed, so an old command line cannot quietly
run the old semantics. Use:

    python gate2_run.py --cohort <name> --layout combined|split \
        --input <dir> --out-dir <dir> --runtime-config <csv> --base-seed <int>

Outputs produced by this script are quarantined under
05_RESULTS/rebuild_upfront_qc/gate2_scrublet/*.INVALID_* and must not be reused.
Its history is in git; see commits 52c9ba1 and earlier.
"""
import sys

sys.exit(
    "RETIRED: this Gate 2 entry point had silent prior and threshold fallbacks and has been "
    "replaced by gate2_run.py, which requires an authorized --runtime-config and an explicit "
    "--base-seed. No authorized runtime rate configuration exists yet, so Gate 2 cannot be "
    "started; see docs/data-contracts.md and docs/execution-plan.md P0."
)
