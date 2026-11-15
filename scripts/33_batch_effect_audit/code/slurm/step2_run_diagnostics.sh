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
# without explicit authorization, per job." Depends on Step 1's output file
# existing and having passed its own acceptance checks first; requires its OWN
# separate authorization, requested only after Step 1 is confirmed (per
# PROVISIONAL_JOB_SPEC_v3.md "One job or two?").
#
# --population-stage below is passed as the LITERAL text
# prepare_provisional_diagnostic_input.py writes to uns["population_stage"] --
# copy it from that file's own output (`python -c "import anndata; print(anndata
# .read_h5ad('$H5AD').uns['population_stage'])"`) rather than retyping it, to
# avoid the exact free-text drift a 2026-10-09 review found between this
# spec's v2 and the preparer's actual output string. --stage-id is the
# machine-checkable, hard-failing version of the same check and does not have
# this risk (it is one of a small fixed vocabulary, not free text).
#
# DEPLOY_ROOT/EXPECTED_REVISION below are the same isolated, independently
# git-verified deployment Step 1 uses -- filled in at deployment time (see
# PROVISIONAL_JOB_SPEC_v3.md's Addendum section for the exact, current
# values), NOT the original /bigdata/godziklab/shared/Xinru/302006/
# 302006_updategit working copy, which remains at stale commit 35131cf with
# uncommitted local modifications and is left untouched. CODE_REVISION below
# is resolved from DEPLOY_ROOT, so --code-revision is the revision actually
# executed, not one asserted from an unrelated/stale checkout. A first
# isolated deployment (302006_updategit_deploy_9ea4c6a, 2026-10-10) was
# superseded after a 2026-10-10 review found it had been hand-patched in
# place without a matching commit, leaving it dirty despite still reporting
# HEAD=9ea4c6a -- exactly what the working-tree-clean check below now
# catches.

set -euo pipefail

echo "============================================================"
echo "STAGE 33 -- run_batch_diagnostics.py"
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

POP_STAGE=$($PYTHON -c "import anndata; print(anndata.read_h5ad('$H5AD', backed='r').uns['population_stage'])")
echo "population_stage read from $H5AD: $POP_STAGE"

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
# 2026-10-10 review: the prior guard checked HEAD only -- a clean checkout at
# the right revision that was then hand-patched (uncommitted changes) would
# still pass. Check working-tree cleanliness too, same as
# run_batch_diagnostics.py's own _resolve_code_identity() requires for
# code_revision_verified=true.
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
