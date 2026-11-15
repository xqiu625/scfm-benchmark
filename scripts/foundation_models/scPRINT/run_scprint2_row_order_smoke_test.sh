#!/bin/bash
# PACKET B -- NOT SUBMITTED. Prepared for review/authorization only.
#
# v2 of this packet (per review:20261018_000122's independent
# review of v1, commit e42086d): v1 had 3 real defects, all fixed here.
#
# B1 (fixed): v1's recovery calls passed the FULL 128,503-cell frozen input as
# the "expected manifest" for a 150-cell output -- recover_cell_ids() correctly
# refused (3 missing, 0 extra, on a 3-vs-6 equivalent probe). Fixed: both
# recovery calls now use the SAME small fixture (natural_fixture.h5ad) that
# both embed runs were built from, as their expected manifest.
#
# B2 (fixed): v1 used a pure reversal as "run B"'s order and only compared
# cross-run consistency. The independent review constructed two runs that BOTH
# apply the exact same systematic (input-order-independent) row-content bug
# and showed cross-run comparison alone returns a clean PASS
# (max_abs_diff=0.0) even though 8/9 rows were actually wrong -- a pure
# reversal is the ONE specific permutation that cancels against a
# "producer always reverses" bug. Fixed two ways: (a) "run B"'s order is now a
# fixed random derangement (seed 20261018, verified != reversal), not a pure
# reversal; (b) ONE real cell's content is duplicated at a second position in
# BOTH fixtures, and check_duplicate_row_self_consistency (an independent,
# WITHIN-one-run oracle bound to content, not position) is run on each
# recovered output in addition to the cross-run comparison. The same review
# also found two runs that were BOTH Inf at the same position passed cleanly;
# the verifier itself now refuses on any non-finite value (fixed in
# stephenson_monocyte_scprint2_pilot_row_order_verify.py, not here).
#
# B3 (fixed):
#   - This script now points at an INDEPENDENT deployment
#     (DEPLOY_DIR_B below), not the same 3de4bca checkout Packet A/the real
#     production job are pinned to -- B depends on files
#     (recover_cell_ids.py, row_order_verify.py) added AFTER 3de4bca, and
#     that production-pinned deployment must never be updated for B's sake.
#   - ATOL/RTOL are now passed EXPLICITLY to every verifier call (v1 declared
#     them but never actually passed them, silently relying on the verifier's
#     own defaults, which happened to be equal).
#   - The --genelist is now built from REAL ENSEMBL GENE IDS, not raw HGNC
#     symbols. Read directly (scdataloader/preprocess.py:268-277): with
#     --is_symbol true (this project's pinned setting), the Preprocessor
#     REPLACES adata.var.index with ensembl_gene_id (keeping the symbol only
#     as a side column) BEFORE the Collator ever runs, and Collator._setup
#     (collator.py:102-140) matches genelist against that SAME
#     ensembl-indexed load_genes() table via df.index.isin(genelist) --
#     passing symbol strings would have matched ZERO genes. This round
#     independently ran the EXACT resolution the installed Preprocessor uses
#     (genesdf.drop_duplicates("symbol").set_index("symbol")) against the
#     real frozen input's first 64 var_names: 51/64 resolve to a real
#     Ensembl ID (13 are unmapped lncRNA/clone identifiers with no HGNC
#     symbol entry in the table), and all 51 are confirmed present in the
#     table's own index. The job script re-derives and re-asserts this at
#     run time (Step 1) rather than hardcoding that count.
#   - Real batch/shard counts are now captured from each embed call's own
#     log (grepping "logging" -- model.py's own print at each shard flush --
#     and the DataLoader's own reported batch count), not just planned from
#     N_CELLS/BATCH_SIZE/SAVE_EVERY arithmetic.
#
# Deliberate, still-open deviation from frozen production parameters, stated
# plainly: uses --how "some" with a pinned --genelist instead of production's
# --how "random expr". scdataloader/collator.py's "random expr" branch draws
# from numpy's GLOBAL random state, consumed once per cell in DataLoader
# order -- changing row order changes which random genes a cell's content
# reduces to, confounding a permutation test. "some"/"all" is
# `loc = np.arange(len(expr))`, no RNG call. This test therefore answers
# "does the real DataLoader -> Collator -> _predict -> shard-write -> concat
# -> final-obsm-assignment call path preserve row identity" for this
# plumbing -- it does NOT answer whether production's actual
# --how "random expr" is itself permutation-invariant, which stays open.
#
#SBATCH --job-name=scprint2_row_order_smoke
#SBATCH --partition=short_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --no-requeue
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_row_order_smoke_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/scprint2_row_order_smoke_%j.err

set -euo pipefail

: "${EXPECTED_DEPLOY_COMMIT_B:?EXPECTED_DEPLOY_COMMIT_B must be set (pass via \
'sbatch --export=ALL,EXPECTED_DEPLOY_COMMIT_B=<full commit hash>') -- refusing to run \
against an unpinned revision}"

# Independent deployment for B -- NEVER the same checkout as the
# production-pinned 3de4bca deployment (see B3 above).
DEPLOY_DIR_B=/bigdata/godziklab/shared/Xinru/302006/302006_updategit_deploy_scprint2_pilot_rowB
if [ ! -d "$DEPLOY_DIR_B/.git" ]; then
    echo "ERROR: $DEPLOY_DIR_B is not a git checkout -- run the B-specific deploy export first" >&2
    exit 1
fi
DEPLOY_DIRTY=$(git -C "$DEPLOY_DIR_B" status --porcelain)
if [ -n "$DEPLOY_DIRTY" ]; then
    echo "ERROR: $DEPLOY_DIR_B has uncommitted changes -- refusing" >&2
    echo "$DEPLOY_DIRTY" >&2
    exit 1
fi
DEPLOY_COMMIT_B=$(git -C "$DEPLOY_DIR_B" rev-parse HEAD)
if [ "$DEPLOY_COMMIT_B" != "$EXPECTED_DEPLOY_COMMIT_B" ]; then
    echo "ERROR: $DEPLOY_DIR_B HEAD ($DEPLOY_COMMIT_B) != EXPECTED_DEPLOY_COMMIT_B ($EXPECTED_DEPLOY_COMMIT_B)" >&2
    exit 1
fi
echo "Deployment (B-only): $DEPLOY_DIR_B @ $DEPLOY_COMMIT_B (clean, matches pinned EXPECTED_DEPLOY_COMMIT_B)"
RESULTS_DIR="${DEPLOY_DIR_B}/scripts/34_scfm_discovery_framework/results"

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
if [ ! -x "$ENV_DIR/bin/python" ]; then
    echo "ERROR: $ENV_DIR/bin/python not found." >&2
    exit 1
fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
ENV_PY="$ENV_DIR/bin/python"
CLEAN_ENV_PY=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python3

HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
EXPECTED_CKPT_SHA256=2b586c144cf9a1f638b4b3e803554ebf6ce389f81c5f34179288a970addfd822
CKPT="${HF_HOME}/scprint2_ckpts/small-v2.ckpt"
ACTUAL_CKPT_SHA256=$(sha256sum "$CKPT" | awk '{print $1}')
if [ "$ACTUAL_CKPT_SHA256" != "$EXPECTED_CKPT_SHA256" ]; then
    echo "ERROR: checkpoint sha256 mismatch" >&2
    exit 1
fi
echo "Verified: checkpoint sha256 matches ($ACTUAL_CKPT_SHA256)"

FROZEN_INPUT_H5AD=/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/scprint2_pilot_runs/run_20261017_154430_30044888/stephenson_monocyte_scprint2_pilot_input.h5ad
EXPECTED_FROZEN_INPUT_SHA256=2bd9818773c1ad13f477979f87b645999cf59f97d6d7a75ca7649b5e150f9e1e
ACTUAL_FROZEN_INPUT_SHA256=$(sha256sum "$FROZEN_INPUT_H5AD" | awk '{print $1}')
if [ "$ACTUAL_FROZEN_INPUT_SHA256" != "$EXPECTED_FROZEN_INPUT_SHA256" ]; then
    echo "ERROR: frozen input sha256 mismatch" >&2
    exit 1
fi
echo "Verified: frozen prepared input sha256 matches ($ACTUAL_FROZEN_INPUT_SHA256)"

N_CELLS=150
BATCH_SIZE=16
SAVE_EVERY=50
N_SYMBOLS_TO_TRY=64
MIN_GENE_HITS_REQUIRED=20
ATOL=1e-4
RTOL=1e-3
PERMUTATION_SEED=20261018

RUN_ID="row_order_smoke_$(date +%Y%m%d_%H%M%S)_${SLURM_JOB_ID:-local}"
RUN_DIR="/bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c_monocyte/scprint2_pilot_runs/${RUN_ID}"
if [ -e "$RUN_DIR" ]; then
    echo "ERROR: $RUN_DIR already exists -- refusing to overwrite" >&2
    exit 1
fi
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

echo "=== Step 1/6: build fixtures (natural + derangement order, 1 injected content-duplicate), resolve real genelist ==="
# h5py-bounded (never anndata backed="r" on this multi-GB file -- that OOM-killed,
# exit 137, even for a 150-row slice, confirmed 2026-10-18; see that script's docstring).
"$ENV_PY" "$RESULTS_DIR/stephenson_monocyte_scprint2_pilot_row_order_build_fixtures.py" \
    "$FROZEN_INPUT_H5AD" "$N_CELLS" "$N_SYMBOLS_TO_TRY" "$MIN_GENE_HITS_REQUIRED" "$PERMUTATION_SEED" "$RUN_DIR"

GENELIST_JSON=$(cat "$RUN_DIR/genelist.json")
DUP_ORIGINAL_ID=$("$CLEAN_ENV_PY" -c "import json; print(json.load(open('$RUN_DIR/duplicate_pair.json'))['original_id'])")
DUP_DUPLICATE_ID=$("$CLEAN_ENV_PY" -c "import json; print(json.load(open('$RUN_DIR/duplicate_pair.json'))['duplicate_id'])")
echo "Pinned genelist (ensembl ids, first few): $(echo "$GENELIST_JSON" | head -c 200)..."
echo "Duplicate-content pair: $DUP_ORIGINAL_ID <-> $DUP_DUPLICATE_ID"

run_embed() {
    local order_name="$1" input_h5ad="$2" output_stem="$3" log_file="$4"
    echo "=== scprint2 embed ($order_name order) -- log: $log_file ==="
    scprint2 embed \
        --adata "$input_h5ad" \
        --ckpt_path "$CKPT" \
        --species "NCBITaxon:9606" \
        --is_symbol true \
        --seed_everything 42 \
        --batch_size "$BATCH_SIZE" \
        --num_workers 4 \
        --how "some" \
        --genelist="$GENELIST_JSON" \
        --save_every "$SAVE_EVERY" \
        --output_filename "$output_stem" \
        > "$log_file" 2>&1
    cat "$log_file"
    # Real, observed (not planned) counts for the record: model.py's own
    # `print("logging")` fires once per shard flush inside _predict(); tqdm's
    # own final progress line reports the real DataLoader batch count.
    local n_shard_flushes
    n_shard_flushes=$(grep -c '^logging$' "$log_file" || true)
    local tqdm_final_line
    tqdm_final_line=$(grep -oE '[0-9]+/[0-9]+ \[' "$log_file" | tail -1 || true)
    echo "  real shard-flush count (grep 'logging' in $log_file): $n_shard_flushes"
    echo "  real tqdm batch-progress final line: ${tqdm_final_line:-<not found>}"
}

echo "=== Step 2/6: embed, natural order ==="
run_embed natural "$RUN_DIR/fixture_natural.h5ad" "$RUN_DIR/embedded_natural" "$RUN_DIR/embed_natural.log"

echo "=== Step 3/6: embed, order_b (fixed derangement) ==="
run_embed order_b "$RUN_DIR/fixture_order_b.h5ad" "$RUN_DIR/embedded_order_b" "$RUN_DIR/embed_order_b.log"

echo "=== Step 4/6: recover real cell ids on both outputs (same fixture as expected manifest, B1 fix) ==="
"$CLEAN_ENV_PY" "$RESULTS_DIR/stephenson_monocyte_scprint2_pilot_recover_cell_ids.py" \
    "$RUN_DIR/embedded_natural_embedded.h5ad" \
    "$RUN_DIR/fixture_natural.h5ad" \
    "$RUN_DIR/embedded_natural_recovered.h5ad"
"$CLEAN_ENV_PY" "$RESULTS_DIR/stephenson_monocyte_scprint2_pilot_recover_cell_ids.py" \
    "$RUN_DIR/embedded_order_b_embedded.h5ad" \
    "$RUN_DIR/fixture_natural.h5ad" \
    "$RUN_DIR/embedded_order_b_recovered.h5ad"

echo "=== Step 5/6: cross-run consistency (necessary, NOT sufficient -- see module docstring) ==="
"$CLEAN_ENV_PY" "$RESULTS_DIR/stephenson_monocyte_scprint2_pilot_row_order_verify.py" \
    "$RUN_DIR/embedded_natural_recovered.h5ad" \
    "$RUN_DIR/embedded_order_b_recovered.h5ad" \
    "$ATOL" "$RTOL" \
    > "$RUN_DIR/cross_run_consistency_report.json"
cat "$RUN_DIR/cross_run_consistency_report.json"

echo "=== Step 6/6: independent within-run content-identity oracle (duplicate-row self-consistency) ==="
"$CLEAN_ENV_PY" -c "
import sys
sys.path.insert(0, '$RESULTS_DIR')
import json
from stephenson_monocyte_scprint2_pilot_row_order_verify import check_duplicate_row_self_consistency

pairs = [('$DUP_ORIGINAL_ID', '$DUP_DUPLICATE_ID')]
for run_name, path in [
    ('natural', '$RUN_DIR/embedded_natural_recovered.h5ad'),
    ('order_b', '$RUN_DIR/embedded_order_b_recovered.h5ad'),
]:
    result = check_duplicate_row_self_consistency(path, pairs, atol=$ATOL, rtol=$RTOL)
    result['run'] = run_name
    print(json.dumps(result, indent=2, default=str))
" | tee "$RUN_DIR/duplicate_row_self_consistency_report.json"

echo "=== DONE ==="
echo "cross_run_consistency: see $RUN_DIR/cross_run_consistency_report.json (scope: consistency between these two orderings only)"
echo "content_identity_self_consistency (x2, natural + order_b): see $RUN_DIR/duplicate_row_self_consistency_report.json"
echo "Reminder: together these are evidence about ROW-IDENTITY PLUMBING under --how=\"some\" with a fixed,"
echo "real ensembl genelist -- they do NOT validate production's --how \"random expr\", and this check is"
echo "NOT a substitute for the identity_validation check in Packet A."
