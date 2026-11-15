#!/bin/bash
#SBATCH --job-name=0006_gate2_capture_scoring
#SBATCH --partition=short
#SBATCH --time=2:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/gate2_capture_scoring_%j.out
#SBATCH --error=logs/gate2_capture_scoring_%j.err
#
# Gate 2 scoring for GSM5102900_HC1. NOT AUTHORISED AND NOT RUNNABLE.
#
# This wrapper refuses to run without a scoped pilot-permit artifact whose
# config_sha256 and input_sha256 match the audited config. No such permit
# exists. `--dry-run` performs every validation and exits WITHOUT scoring; it is
# the only mode that does anything today.
set -euo pipefail

DRY_RUN=0
[ "${1:-}" = "--dry-run" ] && DRY_RUN=1

REPO="${SLURM_SUBMIT_DIR:-$(pwd)}"
cd "$REPO"; REPO="$(git rev-parse --show-toplevel)"; cd "$REPO"

CONFIG="scripts/32_rebuild_upfront_qc/code/gate2_scoring_config_GSM5102900_HC1.py"
ENTRY="scripts/32_rebuild_upfront_qc/code/gate2_score_capture_GSM5102900_HC1.py"
WRAPPER="reports/manuscript/revision/gate2_capture_audit_2026-09-24/run_capture_scoring.sh"
IMPORTED=(
  "scripts/32_rebuild_upfront_qc/code/gate2_parameter_ruling.py"
  "scripts/32_rebuild_upfront_qc/code/gate2_stress_design.py"
  "scripts/32_rebuild_upfront_qc/code/gate2_stress_runner.py"
  # TRANSITIVE: imported by gate2_stress_runner. This list must cover exactly what
  # the config digests; the equality is asserted below rather than trusted, because
  # this list omitted this module while IMPLEMENTATION_FILES already carried it.
  "scripts/32_rebuild_upfront_qc/code/gate2_diagnostics.py"
)

for f in "$CONFIG" "$ENTRY" "$WRAPPER" "${IMPORTED[@]}"; do
    git ls-files --error-unmatch "$f" >/dev/null 2>&1 || { echo "REFUSING: $f untracked" >&2; exit 2; }
    git diff --quiet -- "$f"          || { echo "REFUSING: $f has unstaged changes" >&2; exit 2; }
    git diff --cached --quiet -- "$f" || { echo "REFUSING: $f has staged changes" >&2; exit 2; }
done

# The scoring entrypoint is newer than the preflight commit, so provenance is
# recorded by file hash rather than by a single deployment commit.
# The permit bypass was removed from the code; its presence in the environment is
# treated as an attempt to disable the authenticity check, not as a preference.
if [ -n "${GATE2_PERMIT_TEST_BYPASS:-}" ]; then
    echo "REFUSING: GATE2_PERMIT_TEST_BYPASS is set; that bypass was removed" >&2
    exit 2
fi

COMMIT="$(git rev-parse HEAD)"
echo "repo    : $REPO"
echo "commit  : $COMMIT"
for f in "$CONFIG" "$ENTRY" "$WRAPPER" "${IMPORTED[@]}"; do echo "sha256  : $(sha256sum "$f")"; done

export PATH=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin:$PATH
export PYTHONPATH="$REPO/scripts/32_rebuild_upfront_qc/code:${PYTHONPATH:-}"

# Handed to the validation so it can assert this set equals IMPLEMENTATION_FILES.
export GATE2_GUARDED_FILES="$CONFIG $ENTRY $WRAPPER ${IMPORTED[*]}"

# Validations. `set -e` means a non-zero exit from the heredoc would kill the
# shell before `rc=$?` could be read -- which made the documented dry-run success
# path unreachable, and the "validation exit: 3" reported earlier came from a
# separate direct run, not from this wrapper. errexit is suspended around exactly
# this command and restored immediately after.
set +e
python - <<'VALIDATE'
import os
import sys
import gate2_scoring_config_GSM5102900_HC1 as C
from gate2_stress_runner import Gate2RunnerError, PRODUCTION_PARAMETER_RULING

# The shell guarded a list of files; the permit binds a digest over another list.
# If they ever diverge, a file can be digested without being checked for local
# edits, or checked without being bound. Refuse rather than assume they agree.
guarded = set(os.environ.get('GATE2_GUARDED_FILES', '').split())
digested = set(C.IMPLEMENTATION_FILES)
if guarded != digested:
    print('REFUSING: the wrapper guards a different file set than the permit digests',
          file=sys.stderr)
    print(f'  guarded not digested : {sorted(guarded - digested)}', file=sys.stderr)
    print(f'  digested not guarded : {sorted(digested - guarded)}', file=sys.stderr)
    raise SystemExit(4)

cfg = C.scoring_config()
print(f'config  : {C.CONFIG_VERSION}')
print(f'          sha256 {C.config_sha256()}')
print(f'impl    : sha256 {C.implementation_sha256()}')
print(f'capture : {cfg["capture"]}  runnable={cfg["runnable"]}')
print(f'graph   : graph_preflight_verified={cfg["graph_preflight_verified"]}')
print(f'output  : {cfg["output"]["directory"]}')

C.assert_inputs_unchanged()
print('inputs  : all three deposit hashes match the audited values')

if PRODUCTION_PARAMETER_RULING:
    sys.exit('REFUSING: PRODUCTION_PARAMETER_RULING has been activated')
print('runner  : PRODUCTION_PARAMETER_RULING is still {} (not activated)')

try:
    C.assert_permitted()
except Gate2RunnerError as exc:
    print(f'permit  : ABSENT/INVALID -- {exc}')
    sys.exit(3)
print('permit  : present, tracked and matching')
VALIDATE
rc=$?
set -e

if [ $rc -eq 3 ]; then
    if [ $DRY_RUN -eq 1 ]; then
        echo
        echo "DRY RUN COMPLETE: every validation ran; scoring is blocked for want of a"
        echo "valid permit. This is the expected state today."
        exit 0
    fi
    echo "REFUSING to score without a valid scoped pilot permit." >&2
    exit 3
fi
if [ $rc -ne 0 ]; then
    echo "REFUSING: validation failed with exit $rc" >&2
    exit $rc
fi

if [ $DRY_RUN -eq 1 ]; then
    echo "DRY RUN COMPLETE: validations passed, including the permit gate."
    exit 0
fi

OUT="$(python -c 'import gate2_scoring_config_GSM5102900_HC1 as C; print(C.output_directory())')"
if [ -e "$OUT" ] && [ -n "$(ls -A "$OUT" 2>/dev/null)" ]; then
    echo "REFUSING: output reuse -- $OUT is not empty" >&2; exit 2
fi
mkdir -p "$OUT"
python scripts/32_rebuild_upfront_qc/code/gate2_score_capture_GSM5102900_HC1.py --out-dir "$OUT"
echo "output  : $OUT"
