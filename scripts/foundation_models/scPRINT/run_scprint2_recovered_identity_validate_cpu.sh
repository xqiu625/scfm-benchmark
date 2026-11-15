#!/bin/bash
# PACKET A -- NOT SUBMITTED. Prepared for review/authorization only, per
# review:20261017_223541's instruction and the owner's own
# "执行包准备好后再请你授权新作业" message (2026-10-17).
#
# Purpose: re-run the EXISTING, unmodified
# stephenson_monocyte_scprint2_pilot_validate.py against the real
# recovered_v2/ output produced this round by
# stephenson_monocyte_scprint2_pilot_recover_cell_ids.py (job 30044888's
# identity-recovered artifact), to get a formally recorded technical PASS/FAIL
# -- not just this session's own ad hoc h5py probes. CPU-only; no GPU, no
# model inference, no new checkpoint use. Never touches X/layers/obsp beyond
# whatever the unmodified validate.py script itself reads.
#
# This answers ONLY the identity question (cell ids, required QC/identity
# columns, primary embedding array shape/finiteness). It does NOT and cannot
# answer whether obsm["scprint_emb"]'s rows are numerically bound to the
# correct cell -- that is Packet B's question (row_order_verify /
# row_order_smoke_test), a different, separately gated check. Report the two
# as "identity_validation" and "producer_row_alignment" separately; do not
# merge them into one combined "PASS".
#
#SBATCH --job-name=scprint2_recovered_validate
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --no-requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_recovered_validate_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_recovered_validate_%j.err

set -euo pipefail

# ---- Fixed, reviewed inputs (no CLI args; a drifted value is a refusal, not
# a silent re-point). ----
DEPLOY_DIR=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot
EXPECTED_DEPLOY_COMMIT=3de4bca47ad4ad07d501f557d7b9065540477655
VALIDATE_PY="${DEPLOY_DIR}/scripts/34_scfm_discovery_framework/results/stephenson_monocyte_scprint2_pilot_validate.py"
EXPECTED_VALIDATE_SHA256=c134e6e5ae25d16bd0392413c042772525961b0298ff75c1a006e6a62bc65f1f

RUN_DIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/scprint2_pilot_runs/run_20261017_154430_30044888
RECOVERED_H5AD="${RUN_DIR}/recovered_v2/stephenson_monocyte_scprint2_pilot_embeddings_recovered.h5ad"
EXPECTED_RECOVERED_SHA256=dd3d4fdabc69fc8a585c0b440fe00ed2db99f3177212d53ab76482754b6407e2
EXPECTED_INPUT_H5AD="${RUN_DIR}/stephenson_monocyte_scprint2_pilot_input.h5ad"
EXPECTED_INPUT_SHA256=2bd9818773c1ad13f477979f87b645999cf59f97d6d7a75ca7649b5e150f9e1e

ENV_PY=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python3

# ---- Output: a new, versioned path. Refuses if it already exists -- never
# overwrites a prior validation record (DEVELOPMENT.md: a rerun writes to a new
# versioned directory, the old one stays as the record). ----
OUT_DIR="${RUN_DIR}/recovered_v2/formal_validation_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-local}"
if [ -e "$OUT_DIR" ]; then
    echo "ERROR: $OUT_DIR already exists -- refusing to overwrite" >&2
    exit 1
fi
mkdir -p "$OUT_DIR"
OUT_JSON="${OUT_DIR}/validation_report_recovered_v2.json"

check_sha256() {
    local path="$1" expected="$2" label="$3"
    local actual
    actual=$(sha256sum "$path" | awk '{print $1}')
    if [ "$actual" != "$expected" ]; then
        echo "ERROR: $label sha256 mismatch for $path" >&2
        echo "  expected: $expected" >&2
        echo "  actual:   $actual" >&2
        exit 1
    fi
    echo "Verified: $label sha256 matches ($actual)"
}

# ---- Deployment cleanliness/commit check, same convention as the frozen
# pilot job script -- a drifted or dirty deployment refuses, not warns. ----
if [ ! -d "$DEPLOY_DIR/.git" ]; then
    echo "ERROR: $DEPLOY_DIR is not a git checkout" >&2
    exit 1
fi
DEPLOY_DIRTY=$(git -C "$DEPLOY_DIR" status --porcelain)
if [ -n "$DEPLOY_DIRTY" ]; then
    echo "ERROR: $DEPLOY_DIR has uncommitted changes -- refusing" >&2
    echo "$DEPLOY_DIRTY" >&2
    exit 1
fi
DEPLOY_COMMIT=$(git -C "$DEPLOY_DIR" rev-parse HEAD)
if [ "$DEPLOY_COMMIT" != "$EXPECTED_DEPLOY_COMMIT" ]; then
    echo "ERROR: $DEPLOY_DIR HEAD ($DEPLOY_COMMIT) != EXPECTED_DEPLOY_COMMIT ($EXPECTED_DEPLOY_COMMIT)" >&2
    exit 1
fi
echo "Deployment: $DEPLOY_DIR @ $DEPLOY_COMMIT (clean, matches pinned EXPECTED_DEPLOY_COMMIT)"

# ---- Pin the exact 3 files this job reads, before touching any of them. ----
check_sha256 "$VALIDATE_PY" "$EXPECTED_VALIDATE_SHA256" "validate.py (unmodified since job 30044888)"
check_sha256 "$RECOVERED_H5AD" "$EXPECTED_RECOVERED_SHA256" "recovered_v2 embedded h5ad"
check_sha256 "$EXPECTED_INPUT_H5AD" "$EXPECTED_INPUT_SHA256" "frozen prepared input h5ad"

echo "=== Running the unmodified validator (no code change, no threshold change) ==="
"$ENV_PY" "$VALIDATE_PY" "$RECOVERED_H5AD" "$EXPECTED_INPUT_H5AD" "$OUT_JSON"
VALIDATE_RC=$?

echo "=== Fail-stop: validator exit code is authoritative, not inferred ==="
echo "validate.py exit code: $VALIDATE_RC"
echo "Report written to: $OUT_JSON"
echo "Reminder: a PASS here is identity_validation only, not producer_row_alignment."

exit "$VALIDATE_RC"
