#!/bin/bash
#SBATCH --job-name=stage33_diag_run
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stage33_diag_run_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stage33_diag_run_%j.log
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --time=1:00:00
#
# NOT SUBMITTED by this script's presence -- per DEVELOPMENT.md, "No SLURM submission
# without explicit authorization, per job." Not yet deployed into a frozen,
# git-verified deployment directory either -- see
# scripts/33_batch_effect_audit/results/job2_preparation_packet_v1.md for what
# deploying and running this would require.
#
# v2 (this file) responds to review:20261011_113512 §4 item 5: the
# original step2_run_diagnostics.sh's one real defect is on the line below --
# `anndata.read_h5ad(..., backed="r")` was used ONLY to peek at the tiny
# uns['population_stage'] scalar, but a 2026-10-10 review (R1,
# review:20261010_191146) read the installed anndata source directly
# and found backed="r" still materializes every sparse `layers` entry in full
# -- so this "metadata-only" peek was silently forcing a full read of the
# whole h5ad's layers before the real diagnostic run even started. Replaced
# with a direct h5py read of the uns scalar only (the same pattern
# validate_job1_output.py's read_uns_scalar() already uses), which never
# touches X/layers. For input.h5ad's real size (933 MB as of 2026-10-10),
# this was not an OOM risk under the 16gb allocation -- it is a correctness/
# consistency fix, not an emergency one.
#
# Everything else below is unchanged from step2_run_diagnostics.sh (v1,
# reviewed 2026-10-10/2026-10-09): same DEPLOY_ROOT/EXPECTED_REVISION
# contract, same working-tree-cleanliness guard, same run_batch_diagnostics.py
# invocation and flags. This file does NOT modify the frozen deployment at
# DEPLOY_ROOT in place -- it is a new, separately reviewable wrapper version
# that a future deployment step would need to copy in, with its own
# before/after diff and checksum, exactly like any other code change to that
# deployment.

set -euo pipefail

echo "============================================================"
echo "STAGE 33 -- run_batch_diagnostics.py (wrapper v2)"
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-<not running under slurm>}"
echo "Node: $(hostname)"
echo "============================================================"

DEPLOY_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61
cd "$DEPLOY_ROOT"

PYTHON=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python
echo "Using Python: $PYTHON"
$PYTHON --version

H5AD=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/input.h5ad
OUT_DIR=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3/diagnostic_run

# v2 fix: h5py metadata-only read, never materializes X/layers (unlike
# anndata.read_h5ad(..., backed="r"), which the v1 wrapper used here).
POP_STAGE=$($PYTHON -c "
import h5py
with h5py.File('$H5AD', 'r') as f:
    node = f['uns']['population_stage']
    value = node[()]
    print(value.decode() if isinstance(value, bytes) else value)
")
echo "population_stage read from $H5AD (via h5py, metadata-only): $POP_STAGE"

EXPECTED_REVISION=aea9a6197539ddd5086e46ddd34561617ad250d7
CODE_REVISION=$(git -C "$DEPLOY_ROOT" rev-parse HEAD 2>/dev/null || echo "")
if [ -z "$CODE_REVISION" ]; then
    echo "ERROR: could not resolve a git revision for $DEPLOY_ROOT -- pass" \
         "--code-revision explicitly if deploying from a loose-file scratch copy" \
         "instead of a git checkout. Refusing to guess." >&2
    exit 1
fi
echo "code_revision: $CODE_REVISION (expected $EXPECTED_REVISION)"
if [ "$CODE_REVISION" != "$EXPECTED_REVISION" ]; then
    echo "ERROR: $DEPLOY_ROOT HEAD ($CODE_REVISION) no longer matches the" \
         "revision this wrapper was reviewed against. Re-review before running." >&2
    exit 1
fi
DIRTY=$(git -C "$DEPLOY_ROOT" status --porcelain)
if [ -n "$DIRTY" ]; then
    echo "ERROR: $DEPLOY_ROOT has uncommitted changes -- refusing to run from a" \
         "dirty deployment:" >&2
    echo "$DIRTY" >&2
    exit 1
fi

$PYTHON scripts/33_batch_effect_audit/code/run_batch_diagnostics.py \
    --h5ad "$H5AD" \
    --out-dir "$OUT_DIR" \
    --population-stage "$POP_STAGE" \
    --stage-id pre_gate2_decontx_provisional \
    --code-revision "$CODE_REVISION" \
    --expression-source decontx_counts \
    --batch-col cohort --patient-col patient_id --lineage-col lineage --outcome-col outcome \
    --run-harmony --seed 0

echo "============================================================"
echo "End time: $(date)"
echo "============================================================"
