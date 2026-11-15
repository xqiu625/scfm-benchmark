#!/bin/bash
#SBATCH --job-name=0006_audit_rerun
#SBATCH --partition=batch
#SBATCH --array=0-4
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=30-00:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_audit_%A_%a.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/0006_audit_%A_%a.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

# Reruns required by the 7-module pre-submission audit (2026-08-06).
#
#   task 0  GSEA on the FULL signed ranking       -> tab:per_fm_gsea
#           gsea_per_fm_dim.py truncated the ranking to the top 200 genes BEFORE
#           gp.prerank(), which discards the whole negatively-correlated tail and
#           computes NES/FDR against a 200-gene universe. Every published NES/FDR
#           is therefore invalid. Truncation removed (--top-n 0).
#           NOTE: the script skips any output file that already exists, so this
#           writes to a FRESH directory; reusing option_c_gsea would silently
#           produce nothing.
#
#   task 1  Gene dropout on the MORTALITY endpoint -> tab:dropout
#           The published table used task 'binary', which is a SEVERITY endpoint
#           (severe+fatal vs healthy+mild+moderate+recovered). That contradicts the
#           manuscript twice: it includes the healthy controls Methods excludes, and
#           it merges severe with fatal where the mortality endpoint excludes severe
#           from y=0. New 'mortality' task matches Methods exactly.
#
#   task 2  Multi-classifier with seeded torch     -> tab:multiclassifier
#           DeepMLP/ResidualMLP had no torch.manual_seed, so those columns changed
#           on every run and the published values were not reproducible.
#
#   task 3  Fusion with nested-CV stacking         -> tab:fusion
#           The meta-learner was trained on base OOF probabilities generated using
#           the entire dataset. Emits BOTH late_stacking_lr (legacy, leaky) and
#           late_stacking_lr_nested (leakage-free) so the delta is explicit.
#
# All four are CPU-only sklearn/gseapy jobs on frozen embeddings.
# Submit:  sbatch scripts/30_review_response_stats/code/submit_audit_reruns.sh

set -euo pipefail

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found: $ENV_DIR"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1
echo "Task $SLURM_ARRAY_TASK_ID | Python: $(which python) | Start: $(date)"

REPO=/bigdata/godziklab/shared/Xinru/302006/302006_updategit
cd "$REPO"

HPCC=/bigdata/godziklab/shared/Xinru/302006
EMB=$HPCC/02_EMBEDDINGS/option_c
RES=$HPCC/05_RESULTS
C24=scripts/24_option_c_canonical_input/code

case "$SLURM_ARRAY_TASK_ID" in
  0)
    echo "=== [0a] Regenerate per-dim gene rankings WITHOUT truncation ==="
    # The truncation that invalidates GSEA lives here, not in the GSEA script:
    # biological_interpretation_option_c.py wrote only the top 200 genes by |rho|
    # per dimension, and gsea_per_fm_dim.py consumes those files. Setting --top-n 0
    # on the GSEA step alone changes nothing because its input is already truncated
    # (verified: job 27260030 task 0 still logged "n=200 genes" per dim).
    BIO=$RES/option_c_bio_interpretation_fullrank
    mkdir -p "$BIO"
    python3 -u $C24/biological_interpretation_option_c.py \
        --embedding-dir "$EMB" \
        --raw-input     "$HPCC/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad" \
        --output-dir    "$BIO"
    echo "  genes per dim now: $(tail -n +2 $(ls $BIO/top_genes_*.csv | head -1) | wc -l)"

    echo "=== [0b] GSEA on the FULL signed ranking ==="
    python -c "import gseapy" 2>/dev/null || \
        uv pip install --python "$ENV_DIR/bin/python" gseapy
    OUT=$RES/option_c_gsea_fullrank
    mkdir -p "$OUT"
    python3 -u $C24/gsea_per_fm_dim.py \
        --input-dir  "$BIO" \
        --output-dir "$OUT" \
        --libraries  KEGG_2021_Human Reactome_2022 GO_Biological_Process_2023 \
        --top-n      0 --overwrite
    ;;
  1)
    echo "=== [1] Gene dropout, in-hospital mortality endpoint ==="
    python3 -u scripts/07_gene_dropout/code/gene_dropout_robustness.py \
        --tasks mortality \
        --models STATE UCE Geneformer scPRINT2 scGPT
    ;;
  2)
    echo "=== [2] Multi-classifier, seeded torch (reproducible MLP columns) ==="
    OUT=$RES/option_c_multiclassifier_seeded
    mkdir -p "$OUT"
    python3 -u $C24/multiclassifier_option_c.py \
        --embedding-dir "$EMB" \
        --output-dir    "$OUT"
    ;;
  3)
    echo "=== [3] Fusion, nested-CV stacking (primary 5-FM set) ==="
    OUT=$RES/option_c_fusion_primary5_nested
    mkdir -p "$OUT"
    python3 -u $C24/fusion_option_c.py \
        --embedding-dir "$EMB" \
        --output-dir    "$OUT" \
        --models geneformer scgpt scprint2 state uce
    ;;
  4)
    echo "=== [4] GSEA ONLY, reusing the existing full-rank gene lists ==="
    # Same as step 0b. Split out so GSEA can be rerun without repeating the ~50 min
    # bio_interpretation recompute (which is deterministic and already done).
    python -c "import gseapy" 2>/dev/null || \
        uv pip install --python "$ENV_DIR/bin/python" gseapy
    BIO=$RES/option_c_bio_interpretation_fullrank
    [ -d "$BIO" ] || { echo "ERROR: run task 0 first, $BIO missing"; exit 1; }
    NG=$(tail -n +2 $(ls $BIO/top_genes_*.csv | head -1) | wc -l)
    echo "  genes per dim in input: $NG"
    [ "$NG" -gt 1000 ] || { echo "ERROR: input still looks truncated ($NG genes)"; exit 1; }
    OUT=$RES/option_c_gsea_fullrank
    mkdir -p "$OUT"
    python3 -u $C24/gsea_per_fm_dim.py \
        --input-dir  "$BIO" \
        --output-dir "$OUT" \
        --libraries  KEGG_2021_Human Reactome_2022 GO_Biological_Process_2023 \
        --top-n      0 --overwrite
    ;;
  *)
    echo "ERROR: unexpected SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"; exit 1
    ;;
esac

echo "Task $SLURM_ARRAY_TASK_ID done: $(date)"
