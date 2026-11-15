#!/bin/bash
# DE-backend environment creation + locked install + synthetic verification.
#
# Per review:20261015_150231 F1 (rewriting the prior version this round found
# incomplete on 6 specific points -- see that review and
# DE_BACKEND_EXECUTION_PACKET_v1.md Sec.3 for the full history): this is the ACTUAL,
# concrete wrapper -- not yet run, not yet submitted to SLURM. Running it, or
# submitting it via sbatch, requires the explicit per-job authorization DEVELOPMENT.md's
# standing rules already require; writing this file is not that authorization.
#
# EXECUTION MODE, decided explicitly (the prior version's ambiguous `read -p` prompt
# is removed): this wrapper is meant to run INSIDE a real SLURM allocation (via
# `sbatch` or `salloc ... bash this_script.sh`), non-interactively end to end. The
# human review point the execution packet describes is automated here as an exact
# package-set comparison (Step 2.5) that gates Step 3 -- it exits immediately on any
# mismatch rather than pausing for a keypress a batch job has no way to supply.
#
# REQUIRED INPUTS (no defaults -- an unset/wrong value stops the script immediately,
# per `set -u` below, rather than silently proceeding against the wrong checkout):
#   REPO_CHECKOUT_DIR     -- absolute HPCC path to an already-synced git checkout of
#                            this repository (this script does not fetch/clone).
#   EXPECTED_REPO_COMMIT  -- the exact commit this specific run is pinned against;
#                            a checkout at any other commit is refused outright.
set -euo pipefail

: "${REPO_CHECKOUT_DIR:?must be set to the exact HPCC path of a synced git checkout of this repo}"
: "${EXPECTED_REPO_COMMIT:?must be set to the exact commit this run is pinned against}"
# Running inside a real SLURM allocation is required, not merely recommended --
# this is also what keeps every step off the login node without a separate,
# easy-to-forget manual check.
: "${SLURM_JOB_ID:?this wrapper must run inside a real SLURM allocation (sbatch/salloc), not interactively on a login node}"
: "${SLURM_CPUS_PER_TASK:?SLURM_CPUS_PER_TASK must be set by the allocation -- gene_ranking_providers._bounded_n_cpus() reads it}"

ENV_PREFIX="/rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1"
PROJECT_ROOT="/bigdata/godziklab/shared/Xinru/302006"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"  # measured at actual execution time, never invented
OUT_DIR="${PROJECT_ROOT}/06_LOGS/gate2_capture_scoring/claude_de_backend_v1_${TIMESTAMP}"
RESULTS_DIR="${REPO_CHECKOUT_DIR}/scripts/34_scfm_discovery_framework/results"
LOCK_FILE="${OUT_DIR}/de_backend_requirements.lock.txt"
TEST_LOCK_FILE="${OUT_DIR}/de_backend_test_tooling.lock.txt"
SYNTHETIC_TEST_NAME="test_run_pydeseq2_from_prepared_design_known_effect_and_null_ranks_correctly"

echo "=============================================================="
echo "PREFLIGHT (before any environment is created or touched)"
echo "=============================================================="

ACTUAL_COMMIT="$(git -C "${REPO_CHECKOUT_DIR}" rev-parse HEAD)"
if [[ "${ACTUAL_COMMIT}" != "${EXPECTED_REPO_COMMIT}" ]]; then
    echo "FATAL: ${REPO_CHECKOUT_DIR} is at commit ${ACTUAL_COMMIT}, expected ${EXPECTED_REPO_COMMIT}." >&2
    echo "Refusing to run the synthetic verification against an unexpected checkout/version." >&2
    exit 1
fi

if [[ -e "${ENV_PREFIX}" ]]; then
    echo "FATAL: ${ENV_PREFIX} already exists -- refusing to reuse or overwrite it." >&2
    echo "A prior attempt may be partial; investigate and remove it explicitly first," >&2
    echo "or use a different ENV_PREFIX -- this script does not decide that for you." >&2
    exit 1
fi

mkdir -p "${OUT_DIR}"
cp "${RESULTS_DIR}/de_backend_requirements.lock.txt" "${LOCK_FILE}"
cp "${RESULTS_DIR}/de_backend_test_tooling.lock.txt" "${TEST_LOCK_FILE}"
cp "${RESULTS_DIR}/de_backend_verify_native_resolution.py" "${OUT_DIR}/"
cp "${RESULTS_DIR}/de_backend_verify_synthetic_test_passed.py" "${OUT_DIR}/"
echo "Phase output directory: ${OUT_DIR}"
echo "Verified: ${REPO_CHECKOUT_DIR} is at the expected commit ${EXPECTED_REPO_COMMIT}."

# ---------------------------------------------------------------------------
# Step 1 -- environment creation.
# Running inside a real SLURM allocation (enforced by the SLURM_JOB_ID check
# above) already keeps this off the login node -- the prior OOM on the login
# node (2026-10-14) is not retried here by construction, not merely by
# instruction. `--repodata-fn current_repodata.json` further reduces the
# conda-forge metadata pulled into memory.
# ---------------------------------------------------------------------------
conda create --prefix "${ENV_PREFIX}" -c conda-forge python=3.11 pip \
    --repodata-fn current_repodata.json --yes

NEW_PY="${ENV_PREFIX}/bin/python"

# ---------------------------------------------------------------------------
# Step 2 -- native dependency check FROM the new env's own pip (a real
# resolution, not a cross-interpreter dry run).
# ---------------------------------------------------------------------------
"${NEW_PY}" -m pip install --dry-run \
    --report "${OUT_DIR}/pip_dry_run_report_native_py311.json" \
    -c "${LOCK_FILE}" \
    "pydeseq2==0.5.4" "formulaic>=1.0.2" "formulaic-contrasts>=0.2.0"

# ---------------------------------------------------------------------------
# Step 2.5 -- AUTOMATED gate, replacing a human "read -p" pause (which a
# non-interactive sbatch job has no way to answer). `-c`/constraints alone do
# NOT guarantee the native resolution matches the locked closure exactly --
# this explicit comparison does, in both directions (nothing extra, nothing
# missing, no version drift), and is what this wrapper actually stops on.
# ---------------------------------------------------------------------------
"${NEW_PY}" "${OUT_DIR}/de_backend_verify_native_resolution.py" \
    "${OUT_DIR}/pip_dry_run_report_native_py311.json" "${LOCK_FILE}"

# ---------------------------------------------------------------------------
# Step 3 -- locked installs. The first actually-irreversible-ish action in
# this sequence. TWO separate installs, from TWO separate lock files -- the
# pydeseq2/formulaic closure (de_backend_requirements.lock.txt) and the
# test-runner closure (de_backend_test_tooling.lock.txt) -- so pytest's own
# dependencies are never silently folded into pydeseq2's dependency-pin
# evidence. On ANY failure in either: this script stops immediately
# (`set -e`); pip's own stderr is captured to "${OUT_DIR}/install_stderr.log"
# -- no unpinned/"nearest compatible" substitute is installed, and Step 4
# does not run.
# ---------------------------------------------------------------------------
"${NEW_PY}" -m pip install -c "${LOCK_FILE}" \
    "pydeseq2==0.5.4" "formulaic>=1.0.2" "formulaic-contrasts>=0.2.0" \
    2> "${OUT_DIR}/install_stderr.log"
"${NEW_PY}" -m pip install -c "${TEST_LOCK_FILE}" "pytest==9.1.1" \
    2>> "${OUT_DIR}/install_stderr.log"

# ---------------------------------------------------------------------------
# Step 4 -- synthetic positive/null backend verification. Reuses the REAL
# test already written and reviewed in this repository -- never a parallel,
# duplicate fixture. Structured JUnit output is REQUIRED and checked --
# pytest's own exit code / a log's last line are NOT treated as sufficient,
# since both would look identical whether the test PASSED or was SKIPPED
# (e.g. if pydeseq2 somehow failed to import despite Step 3 succeeding).
# ---------------------------------------------------------------------------
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "${REPO_CHECKOUT_DIR}"
"${NEW_PY}" -m pytest \
    tests/test_scfm_gene_ranking_providers.py \
    -k "${SYNTHETIC_TEST_NAME}" \
    --junitxml="${OUT_DIR}/synthetic_verification_junit.xml" \
    -v 2>&1 | tee "${OUT_DIR}/synthetic_verification.log"

"${NEW_PY}" "${OUT_DIR}/de_backend_verify_synthetic_test_passed.py" \
    "${OUT_DIR}/synthetic_verification_junit.xml" "${SYNTHETIC_TEST_NAME}"

echo "Done -- the synthetic verification test genuinely passed (not merely exited 0)."
echo "See ${OUT_DIR} for every artifact from this run."
