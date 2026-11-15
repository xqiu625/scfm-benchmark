#!/bin/bash
#SBATCH --job-name=0006_gate2_capture_preflight_hc2
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/gate2_capture_preflight_hc2_%j.out
#SBATCH --error=logs/gate2_capture_preflight_hc2_%j.err
#
# Single-capture, score-free Gate 2 preflight CANDIDATE for GSM5102901_HC2.
# NO score, call or rate. No outcome read. Nothing deleted.
# Leaves PRODUCTION_PARAMETER_RULING = {}.
#
# NOT AUTHORISED FOR SUBMISSION. Resource ceiling copied from the reviewed
# GSM5102900_HC1 preflight (run_capture_preflight.sh), which ran the structurally
# identical score-free pipeline (same filter/normalize/HVG/manifold steps, same
# 7,666-cell-scale deposit) to completion well inside this ceiling (job 29073092).
# GSM5102901_HC2 is a single capture of comparable scale (8,345 cells, same 33,538
# genes, same preprocessing) from the same cohort loader, so the same ceiling is
# used as a starting basis, not a guarantee -- if this job is authorised and run,
# its own resource usage should be recorded and compared, not assumed to match.
#
# `logs/` must exist BEFORE sbatch, because SLURM opens the relative --output path
# at submission time; creating it inside the job is too late and the submission
# fails.
set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO"
REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"

SCRIPT="scripts/32_rebuild_upfront_qc/hc2_preflight/capture_preflight_GSM5102901_HC2.py"
WRAPPER="scripts/32_rebuild_upfront_qc/hc2_preflight/run_capture_preflight_HC2.sh"
# Every file whose content decides what runs: the script, this wrapper, and the
# modules it imports. Validating only the script would let a changed ruling or a
# changed runner reach the run unrecorded.
IMPORTED=(
  "scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py"
  "scripts/32_rebuild_upfront_qc/code/gate2_stress_design.py"
  "scripts/32_rebuild_upfront_qc/code/gate2_stress_runner.py"
  # gate2_stress_runner.py imports DIAGNOSTICS_SCHEMA_VERSION/score_distribution_summary
  # from this module at load time (confirmed by grep against the live file). Missing
  # here was review:20261012_132816 R2 -- a changed gate2_diagnostics.py would
  # have reached the run unrecorded.
  "scripts/32_rebuild_upfront_qc/code/gate2_diagnostics.py"
)

# Both the working tree AND the index must match HEAD. Checking only
# `git diff` lets a staged-but-uncommitted change through: the file on disk
# differs from HEAD while `git diff` reports clean because the change is staged.
for f in "$SCRIPT" "$WRAPPER" "${IMPORTED[@]}"; do
    git ls-files --error-unmatch "$f" >/dev/null 2>&1 || { echo "REFUSING: $f untracked" >&2; exit 2; }
    git diff --quiet -- "$f"          || { echo "REFUSING: $f has unstaged changes" >&2; exit 2; }
    git diff --cached --quiet -- "$f" || { echo "REFUSING: $f has staged changes" >&2; exit 2; }
done

COMMIT="$(git rev-parse HEAD)"
echo "repo    : $REPO"
echo "commit  : $COMMIT"
for f in "$SCRIPT" "$WRAPPER" "${IMPORTED[@]}"; do
    echo "sha256  : $(sha256sum "$f")"
done

OUT="/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_preflight/hc2_${COMMIT:0:8}_$(date +%Y%m%d_%H%M%S)"
[ -e "$OUT" ] && { echo "REFUSING: $OUT exists" >&2; exit 2; }
mkdir -p "$OUT"

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONPATH="$REPO/scripts/32_rebuild_upfront_qc/code:${PYTHONPATH:-}"
python "$SCRIPT" --out-dir "$OUT"
echo "output  : $OUT"
