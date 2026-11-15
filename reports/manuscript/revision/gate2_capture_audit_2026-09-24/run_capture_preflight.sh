#!/bin/bash
#SBATCH --job-name=0006_gate2_capture_preflight
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/gate2_capture_preflight_%j.out
#SBATCH --error=logs/gate2_capture_preflight_%j.err
#
# Single-capture, score-free Gate 2 preflight for GSM5102900_HC1.
# NO score, call or rate. No outcome read. Nothing deleted.
# Leaves PRODUCTION_PARAMETER_RULING = {}.
#
# NOT AUTHORISED FOR SUBMISSION. `logs/` must exist BEFORE sbatch, because SLURM
# opens the relative --output path at submission time; creating it inside the job
# is too late and the submission fails.
set -euo pipefail

REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO"
REPO="$(git rev-parse --show-toplevel)"
cd "$REPO"

SCRIPT="reports/manuscript/revision/gate2_capture_audit_2026-09-24/capture_preflight_GSM5102900_HC1.py"
WRAPPER="reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_preflight.sh"
# Every file whose content decides what runs: the script, this wrapper, and the
# modules it imports. Validating only the script would let a changed ruling or a
# changed runner reach the run unrecorded.
IMPORTED=(
  "scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py"
  "scripts/32_rebuild_upfront_qc/code/gate2_stress_design.py"
  "scripts/32_rebuild_upfront_qc/code/gate2_stress_runner.py"
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

OUT="/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/gate2_capture_preflight/${COMMIT:0:8}_$(date +%Y%m%d_%H%M%S)"
[ -e "$OUT" ] && { echo "REFUSING: $OUT exists" >&2; exit 2; }
mkdir -p "$OUT"

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONPATH="$REPO/scripts/32_rebuild_upfront_qc/code:${PYTHONPATH:-}"
python "$SCRIPT" --out-dir "$OUT"
echo "output  : $OUT"
