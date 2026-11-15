#!/bin/bash
#SBATCH --job-name=stage33_prep_diag_input
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stage33_prep_diag_input_%j.log
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/stage33_prep_diag_input_%j.log
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32gb
#SBATCH --time=4:00:00
#
# NOT SUBMITTED by this script's presence -- per DEVELOPMENT.md, "No SLURM submission
# without explicit authorization, per job." This file exists so the exact command,
# interpreter, and resource request a submission would use are reviewable before
# that authorization is requested, per PROVISIONAL_JOB_SPEC_v3.md Step 1.
#
# Prepares the bounded, patient-balanced diagnostic input h5ad from the real
# (pre-Gate-2, PROVISIONAL) assembled benchmark matrices. Streams both
# 22GB/50GB counts.mtx files column-filtered; never materializes them in full.
#
# DEPLOY_ROOT/EXPECTED_REVISION below are filled in at deployment time (see
# PROVISIONAL_JOB_SPEC_v3.md's Addendum section for the exact, current
# values) -- NOT the original /bigdata/godziklab/shared/Xinru/302006/
# 302006_updategit working copy, which remains at stale commit 35131cf with
# uncommitted local modifications and is left untouched as the prior record.
# A first isolated deployment (302006_updategit_deploy_9ea4c6a, 2026-10-10)
# was superseded after a 2026-10-10 review found it had been hand-patched in
# place (two wrapper files rsync'd in for a syntax check, without a matching
# commit), leaving it dirty despite still reporting HEAD=9ea4c6a -- this is
# exactly what the working-tree-clean check below now catches. Re-deploy to
# a NEW, differently-named DEPLOY_ROOT (never overwrite or hand-patch one
# in place) if the code changes again.

set -euo pipefail

echo "============================================================"
echo "STAGE 33 -- prepare_provisional_diagnostic_input.py"
echo "Start time: $(date)"
echo "Job ID: ${SLURM_JOB_ID:-<not running under slurm>}"
echo "Node: $(hostname)"
echo "============================================================"

DEPLOY_ROOT=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_aea9a61
cd "$DEPLOY_ROOT"

EXPECTED_REVISION=aea9a6197539ddd5086e46ddd34561617ad250d7
CODE_REVISION=$(git rev-parse HEAD 2>/dev/null || echo "")
if [ -z "$CODE_REVISION" ]; then
    echo "ERROR: could not resolve a git revision for $DEPLOY_ROOT -- refusing" \
         "to run from an unverifiable deployment." >&2
    exit 1
fi
echo "deployed code_revision: $CODE_REVISION (expected $EXPECTED_REVISION)"
if [ "$CODE_REVISION" != "$EXPECTED_REVISION" ]; then
    echo "ERROR: $DEPLOY_ROOT HEAD ($CODE_REVISION) no longer matches the" \
         "revision this wrapper was reviewed against. Re-review before running." >&2
    exit 1
fi
# 2026-10-10 review: the prior guard checked HEAD only -- a clean checkout at
# the right revision that was then hand-patched (uncommitted changes) would
# still pass. Check working-tree cleanliness too, same as
# run_batch_diagnostics.py's own _resolve_code_identity() requires for
# code_revision_verified=true.
DIRTY=$(git status --porcelain)
if [ -n "$DIRTY" ]; then
    echo "ERROR: $DEPLOY_ROOT has uncommitted changes -- refusing to run from a" \
         "dirty deployment:" >&2
    echo "$DIRTY" >&2
    exit 1
fi

PYTHON=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python
echo "Using Python: $PYTHON"
$PYTHON --version

OUT_DIR=/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/batch_effect_audit/provisional_diagnostic_input_v3
mkdir -p "$OUT_DIR"

$PYTHON scripts/33_batch_effect_audit/code/prepare_provisional_diagnostic_input.py \
    --input-root /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/platelet_monocyte_benchmark_input_v2 \
    --out "$OUT_DIR/input.h5ad" \
    --n-per-donor 20 --seed 0

echo "============================================================"
echo "End time: $(date)"
echo "============================================================"
