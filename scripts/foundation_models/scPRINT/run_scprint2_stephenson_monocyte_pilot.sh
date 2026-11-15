#!/bin/bash
#SBATCH --job-name=scprint2_stephenson_monocyte_pilot
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --no-requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_stephenson_monocyte_pilot_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_stephenson_monocyte_pilot_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# v4 (2026-10-16). NOT SUBMITTED as of this version -- requires its own per-job SLURM
# authorization before `sbatch` is run, per DEVELOPMENT.md's standing rule. Fixes applied
# per review:20261016_202903 ("R1" + "R4" -- "R2"/"R3" from that review
# are in `stephenson_monocyte_scprint2_pilot_validate.py`, not this file) on top of v3
# (deployment/output-identity fixes from review:20261016_152556) and v2
# (R1-R6 from review:20261016_141242):
#
# -- R1: --time=3:00:00 exceeded short_gpu's real MaxTime -- queried directly
# (`scontrol show partition short_gpu`, 2026-10-16): `DefaultTime=02:00:00
# MaxTime=02:00:00`. Lowered to --time=2:00:00. Whether 2h/128G is actually ENOUGH for
# this larger (128,503-cell) input remains UNVERIFIED -- not measured this round, not
# assumed sufficient. If it is not enough, the run fails closed (no-requeue, no retry)
# rather than silently reducing cells or switching to an unreviewed partition.
#
# -- R4: the clean-checkout check below (carried over from v3) only checked that
# DEPLOY_DIR was clean and recorded whatever HEAD it found -- an updated-but-still-clean
# checkout would still pass, which matters if DEPLOY_DIR is fast-forwarded to a new
# commit while this job sits queued. The submitter must now pass the exact commit this
# job is pinned to via `EXPECTED_DEPLOY_COMMIT` (sbatch --export, not baked into this
# file -- baking the file's own future commit hash into itself would be self-
# referential); the job refuses to run if that variable isn't set, or doesn't match the
# deployment's actual HEAD at run time.
#
# -- Code/environment/working-directory binding, not just a correct run order --
# v2 pointed PREP_RESULTS_DIR at the LIVE `302006_updategit` checkout, which the same
# review confirmed (again, this round) is still stale at 35131cf with unrelated dirty
# changes. This script now points at DEPLOY_DIR, an INDEPENDENT clean export of a
# single fixed commit (never the dirty checkout -- that directory is never read,
# modified, or cleaned by this script or by the deployment step that produced
# DEPLOY_DIR). The script verifies DEPLOY_DIR is a clean git checkout (no uncommitted
# changes) at runtime and records its exact commit hash in the provenance file -- so a
# silent later edit to DEPLOY_DIR, not just to the dirty checkout, is also caught.
#
# -- Working directory isolation -- the installed `scprint2/tasks/cell_emb.py:183-190`
# (read directly 2026-10-16) writes intermediate per-batch prediction files to a
# RELATIVE `data/` path when no Lightning logger provides `save_dir` (true for this
# bare `embed` CLI call). Only setting the FINAL --output_filename does not isolate
# these intermediate files to this run's directory. This script now `cd`s into the
# per-run RUN_DIR before calling `scprint2 embed`, so that relative path resolves
# inside the isolated run directory, not wherever `sbatch` happened to be invoked from.
#
# -- Output identity/value validation, not just a count -- the validate step now
# receives the REAL prepared-input h5ad path (the frozen cell-id/QC manifest), not
# merely an expected integer count: `stephenson_monocyte_scprint2_pilot_validate.py`
# v3 checks exact cell-id coverage (no missing/extra/duplicate ids), the real primary
# `scprint_emb` key specifically (not just any `X_scprint*`-prefixed key), array
# dtype/shape, and QC-annotation correspondence BY ID -- see that script's own
# docstring and tests for the full list of scenarios this now rejects.
#
# -- Input/QC/checkpoint fingerprints verified BEFORE the run, not just recorded
# after -- RAW_H5AD, GATE1_FLAGS_CSV, and DECONTX_SCORES_CSV are now sha256-checked
# against the values measured directly from the real files on 2026-10-16, in addition
# to the checkpoint hash check carried over from v2. A drifted input fails closed
# before any GPU time is spent, instead of only being noted in the output summary.
#
# -- Pinned run parameters -- `--seed_everything`, `--batch_size`, `--num_workers`,
# `--how`, and `--max_len` are now passed EXPLICITLY, matching the installed
# `Embedder`/CLI class defaults read directly from source (scprint2/cli.py,
# scprint2/tasks/cell_emb.py, both read 2026-10-16) -- so a future library version
# changing its own defaults cannot silently change this run's behavior without this
# script also changing.
#
# Carried over unchanged from v2 (already fixed, not revisited this round): R1's
# `_embedded.h5ad` output-suffix handling, R2's `--is_symbol true`, R3's cell-filter
# finding, the checkpoint SHA256 runtime check, `--no-requeue`, and the fail-stop
# prep -> embed -> validate bundling into one job.

set -euo pipefail

echo "=== scPRINT-2 Stephenson monocyte TECHNICAL PILOT: prep -> embed -> validate ==="
echo "Start: $(date)"
echo "Host: $(hostname)"
echo "Job:  ${SLURM_JOB_ID:-none}"
nvidia-smi || true

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found. Run install_scprint2.sh first."
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export CONDA_DEFAULT_ENV=scprint2_env
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
echo "Python: $(which python)  ($(python --version 2>&1))"

# --- Deployment binding: an INDEPENDENT clean export of one fixed commit. Never the
# live, stale, dirty `302006_updategit` checkout -- that path is not referenced
# anywhere below. ---
# R4 (review:20261016_202903): the exact commit this job is pinned to
# must be supplied by the submitter, not read from whatever DEPLOY_DIR happens to be
# at run time -- that is what lets a fast-forward during queue wait get caught below,
# instead of only catching an uncommitted-changes drift.
: "${EXPECTED_DEPLOY_COMMIT:?EXPECTED_DEPLOY_COMMIT must be set (pass via \
'sbatch --export=ALL,EXPECTED_DEPLOY_COMMIT=<full commit hash>') -- refusing to run \
against an unpinned revision}"

DEPLOY_DIR=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot
if [ ! -d "$DEPLOY_DIR/.git" ]; then
    echo "ERROR: $DEPLOY_DIR is not a git checkout. Re-run the deployment step first."
    exit 1
fi
DEPLOY_DIRTY=$(git -C "$DEPLOY_DIR" status --porcelain)
if [ -n "$DEPLOY_DIRTY" ]; then
    echo "ERROR: $DEPLOY_DIR has uncommitted changes -- refusing to run from a drifted deployment."
    echo "$DEPLOY_DIRTY"
    exit 1
fi
DEPLOY_COMMIT=$(git -C "$DEPLOY_DIR" rev-parse HEAD)
if [ "$DEPLOY_COMMIT" != "$EXPECTED_DEPLOY_COMMIT" ]; then
    echo "ERROR: $DEPLOY_DIR HEAD ($DEPLOY_COMMIT) does not match the pinned revision this job"
    echo "was submitted for (EXPECTED_DEPLOY_COMMIT=$EXPECTED_DEPLOY_COMMIT) -- a clean-but-"
    echo "updated deployment is still a revision mismatch and must not silently run."
    exit 1
fi
echo "Deployment: $DEPLOY_DIR @ $DEPLOY_COMMIT (clean, matches pinned EXPECTED_DEPLOY_COMMIT)"
PREP_RESULTS_DIR="${DEPLOY_DIR}/scripts/34_scfm_discovery_framework/results"

RUN_ID="run_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-local}"
BASE_DIR=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/scprint2_pilot_runs
RUN_DIR="${BASE_DIR}/${RUN_ID}"
mkdir -p "$RUN_DIR"
echo "Run directory: $RUN_DIR"

# --- Input/QC fingerprints, verified BEFORE spending any GPU time -- measured
# directly from the real files on 2026-10-16 (sha256sum). ---
RAW_H5AD="/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/stephenson_monocyte_raw.h5ad"
GATE1_FLAGS_CSV="/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/gate1_qc/stephenson/qc_flags_per_cell.csv"
DECONTX_SCORES_CSV="/bigdata/godziklab/shared/Xinru/302006/05_RESULTS/rebuild_upfront_qc/stephenson_decontx_input/stephenson_decontx_contamination_scores.csv"
EXPECTED_RAW_H5AD_SHA256="11eb7450ff5c2a87158add9b906521cac43db56b71561d75308788089f0f4f4b"
EXPECTED_GATE1_FLAGS_CSV_SHA256="1034909f4c3ff3cefaba1ca5989603b09a2a64f955ca3d870a5a51b52197d4ee"
EXPECTED_DECONTX_SCORES_CSV_SHA256="6a87de6d176bddb43f3aa703d7b0af448431dcf451487bbd42578381d20ac06e"

check_sha256() {
    local path="$1" expected="$2" label="$3"
    if [ ! -f "$path" ]; then
        echo "ERROR: $label not found: $path"
        exit 1
    fi
    local actual
    actual=$(sha256sum "$path" | awk '{print $1}')
    if [ "$actual" != "$expected" ]; then
        echo "ERROR: $label sha256 mismatch -- input has drifted since this script was written."
        echo "  path:     $path"
        echo "  expected: $expected"
        echo "  actual:   $actual"
        exit 1
    fi
    echo "Verified: $label matches its recorded fingerprint ($path)"
}
check_sha256 "$RAW_H5AD" "$EXPECTED_RAW_H5AD_SHA256" "raw Stephenson monocyte h5ad"
check_sha256 "$GATE1_FLAGS_CSV" "$EXPECTED_GATE1_FLAGS_CSV_SHA256" "Gate 1 QC flags CSV"
check_sha256 "$DECONTX_SCORES_CSV" "$EXPECTED_DECONTX_SCORES_CSV_SHA256" "DecontX contamination scores CSV"

# --- Checkpoint: runtime SHA256 verification, not just a recorded value ---
EXPECTED_CKPT_SHA256="2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822"
CKPT="${HF_HOME}/scprint2_ckpts/small-v2.ckpt"
if [ ! -f "$CKPT" ]; then
    CKPT=$(find "$HF_HOME/hub" -name "small-v2.ckpt" 2>/dev/null | head -1)
fi
check_sha256 "$CKPT" "$EXPECTED_CKPT_SHA256" "scPRINT-2 small-v2.ckpt checkpoint"

PREPARED_H5AD="${RUN_DIR}/stephenson_monocyte_scprint2_pilot_input.h5ad"
PREP_SUMMARY="${PREPARED_H5AD%.h5ad}_summary.json"
EMBED_OUTPUT_STEM="${RUN_DIR}/stephenson_monocyte_scprint2_pilot_embeddings"
EMBED_OUTPUT="${EMBED_OUTPUT_STEM}_embedded.h5ad"
VALIDATION_REPORT="${RUN_DIR}/validation_report.json"
PROVENANCE_JSON="${RUN_DIR}/provenance.json"

# --- Provenance: resolved package versions, pinned run parameters, and key
# installed-source-file hashes, recorded from the ACTUAL running environment, not
# from comments or intended defaults (DEVELOPMENT.md Sec.7). ---
python -c "
import hashlib, json, sys
import importlib.metadata as m

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()

site = '${ENV_DIR}/lib/python3.11/site-packages'
provenance = {
    'deploy_dir': '${DEPLOY_DIR}',
    'deploy_commit': '${DEPLOY_COMMIT}',
    'run_id': '${RUN_ID}',
    'checkpoint_path': '${CKPT}',
    'checkpoint_sha256': '${EXPECTED_CKPT_SHA256}',
    'installed_package_versions': {
        p: m.version(p) for p in
        ['scprint2', 'scdataloader', 'anndata', 'scanpy', 'torch', 'lightning', 'numpy', 'pandas']
    },
    'installed_source_sha256': {
        'scprint2/cli.py': sha256_file(site + '/scprint2/cli.py'),
        'scdataloader/preprocess.py': sha256_file(site + '/scdataloader/preprocess.py'),
        'scprint2/tasks/cell_emb.py': sha256_file(site + '/scprint2/tasks/cell_emb.py'),
    },
    'pinned_embed_parameters': {
        'seed_everything': 42,
        'batch_size': 64,
        'num_workers': 8,
        'how': 'random expr',
        'max_len': 2000,
        'is_symbol': True,
        'species': 'NCBITaxon:9606',
    },
}
with open('${PROVENANCE_JSON}', 'w') as f:
    json.dump(provenance, f, indent=2, sort_keys=True)
print(json.dumps(provenance, indent=2, sort_keys=True))
"
echo "Provenance recorded: $PROVENANCE_JSON"

echo ""
echo "=== Step 1/3: prep (GEX-only subset + Gate1/DecontX annotation, no filtering) ==="
python -c "
import sys
sys.path.insert(0, '${PREP_RESULTS_DIR}')
import stephenson_monocyte_scprint2_pilot_prep as prep
prep.main(out_h5ad='${PREPARED_H5AD}')
"
if [ ! -f "$PREPARED_H5AD" ]; then
    echo "ERROR: prep step did not produce $PREPARED_H5AD"
    exit 2
fi
N_CELLS_IN=$(python -c "import json; print(json.load(open('${PREP_SUMMARY}'))['n_cells_in'])")
echo "Prep complete: $PREPARED_H5AD ($N_CELLS_IN cells in)"

echo ""
echo "=== Step 2/3: scprint2 embed ==="
echo "Input:  $PREPARED_H5AD"
echo "Ckpt:   $CKPT"
echo "Output: $EMBED_OUTPUT (scprint2 appends '_embedded.h5ad' to --output_filename -- verified against the installed cli.py, 2026-10-16)"

# Working directory isolation: `scprint2/tasks/cell_emb.py` writes intermediate
# per-batch prediction files to a RELATIVE `data/` path when no logger provides a
# save_dir (true here). `cd` into RUN_DIR first so that relative path, and anything
# else the library writes relative to cwd, stays inside this run's own directory.
cd "$RUN_DIR"

scprint2 embed \
    --adata "$PREPARED_H5AD" \
    --ckpt_path "$CKPT" \
    --species "NCBITaxon:9606" \
    --is_symbol true \
    --seed_everything 42 \
    --batch_size 64 \
    --num_workers 8 \
    --how "random expr" \
    --max_len 2000 \
    --output_filename "$EMBED_OUTPUT_STEM"

if [ ! -f "$EMBED_OUTPUT" ]; then
    echo "ERROR: embed step did not produce the expected output: $EMBED_OUTPUT"
    echo "Contents of run directory:"
    ls -lh "$RUN_DIR"
    exit 3
fi
echo "Embed complete:"
ls -lh "$EMBED_OUTPUT"

echo ""
echo "=== Step 3/3: validate (technical criteria only -- not a scientific judgment) ==="
echo "Checking: real cell-id coverage/duplicates against the frozen prepared-input manifest,"
echo "the primary 'scprint_emb' key, array dtype/shape, and QC-annotation correspondence by id."
VALIDATE_RC=0
python -c "
import sys
sys.path.insert(0, '${PREP_RESULTS_DIR}')
import stephenson_monocyte_scprint2_pilot_validate as val
rc = val.main('${EMBED_OUTPUT}', '${PREPARED_H5AD}', '${VALIDATION_REPORT}')
sys.exit(rc)
" || VALIDATE_RC=$?

echo ""
echo "Run directory: $RUN_DIR"
echo "Provenance: $PROVENANCE_JSON"
echo "Validation report: $VALIDATION_REPORT"
echo "Completed: $(date)"
exit $VALIDATE_RC
