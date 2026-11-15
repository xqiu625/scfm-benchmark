#!/bin/bash
# ==============================================================================
# Option C: regenerate human COVID embeddings for 5 FMs from a SINGLE canonical
# raw AnnData. Each FM applies its OWN preprocessing internally.
#
# Models (scGPT dropped per project decision 2026-05-01):
#   STATE, UCE, Geneformer, TranscriptFormer, scPRINT
#
# Mirrors step2b_resubmit_5_fms.sh (mouse) — same pattern, COVID side.
# Output: /bigdata/godziklab/shared/Xinru/302006/02_EMBEDDINGS/option_c/<fm>_covid_canonical.h5ad
#
# Usage on HPCC:
#   bash scripts/24_option_c_canonical_input/code/step2c_human_covid_canonical_input.sh
#
# Author: Project 302006
# ==============================================================================

set -e

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
# CANONICAL RAW INPUT — verify with `python -c "import anndata; print(anndata.read_h5ad(...).shape)"`
INPUT="${INPUT:-$BASE_DIR/01_DATA/raw/human_platelet_covid_severity_raw_counts.h5ad}"
EMB_DIR="$BASE_DIR/02_EMBEDDINGS/option_c"
LOG_DIR="$BASE_DIR/06_LOGS/option_c"

mkdir -p "$EMB_DIR" "$LOG_DIR"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: canonical input not found: $INPUT"
    echo "Override with: INPUT=/path/to/raw.h5ad bash $0"
    exit 1
fi

echo "============================================================"
echo "OPTION C: 5 HUMAN COVID FM EMBEDDING JOBS (canonical raw input)"
echo "============================================================"
echo "Input:      $INPUT"
echo "Output dir: $EMB_DIR"
python -c "
import anndata as ad
a = ad.read_h5ad('$INPUT')
print(f'Input shape: {a.n_obs} cells x {a.n_vars} genes')
print(f'obs cols:    {list(a.obs.columns)[:15]}')
" 2>/dev/null || echo "(could not read input shape — make sure conda env loadable)"
echo ""

# ============================================================================
# STATE
# ============================================================================
echo "Submitting STATE..."
JOB_STATE=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=optC_STATE_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=$LOG_DIR/state_%j.log
#SBATCH --error=$LOG_DIR/state_%j.err

set -euo pipefail
echo "=== Option C: STATE on canonical COVID input ==="
echo "Start: \$(date)"

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
if [ ! -x "\$ENV_DIR/bin/python" ]; then echo "ERROR: \$ENV_DIR/bin/python not found"; exit 1; fi
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=state_env
export PYTHONNOUSERSITE=1
export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1
echo "Python: \$(which python)  (\$(python --version 2>&1))"
which state

OUTPUT=$EMB_DIR/state_covid_canonical.h5ad
if python -c "import anndata,sys; a=anndata.read_h5ad('\$OUTPUT'); sys.exit(0 if 'X_state' in a.obsm else 1)" 2>/dev/null; then
    echo "Already done. Exit."; exit 0
fi

# Use the state CLI from arc-state (state_env): state emb transform
# (one line — heredoc backslash-newlines get mangled)
state emb transform --model-folder /bigdata/godziklab/shared/Xinru/302006/03_MODELS/SE-600M --input "$INPUT" --output "\$OUTPUT" --embed-key X_state --batch-size 64
echo "Completed: \$(date)"
SLURM
)
echo "  STATE Job ID: $JOB_STATE"

# ============================================================================
# UCE
# ============================================================================
echo "Submitting UCE..."
JOB_UCE=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=optC_UCE_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=$LOG_DIR/uce_%j.log
#SBATCH --error=$LOG_DIR/uce_%j.err

set -euo pipefail
echo "=== Option C: UCE on canonical COVID input ==="
echo "Start: \$(date)"

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
if [ ! -x "\$ENV_DIR/bin/python" ]; then echo "ERROR: \$ENV_DIR/bin/python not found"; exit 1; fi
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=state_env
export PYTHONNOUSERSITE=1
echo "Python: \$(which python)  (\$(python --version 2>&1))"

OUTPUT=$EMB_DIR/uce_covid_canonical.h5ad
if python -c "import anndata,sys; a=anndata.read_h5ad('\$OUTPUT'); sys.exit(0 if 'X_uce' in a.obsm else 1)" 2>/dev/null; then
    echo "Already done. Exit."; exit 0
fi

UCE_DIR=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE
UCE_OUTDIR=$EMB_DIR/uce_run
mkdir -p "\$UCE_OUTDIR"
cd "\$UCE_DIR"
python eval_single_anndata.py --adata_path "$INPUT" --dir "\$UCE_OUTDIR" --model_loc "\$UCE_DIR/model_files/4layer_model.torch" --species human --batch_size 32

# UCE writes <input_basename>_uce_adata.h5ad under --dir. Move to canonical name.
PRODUCED=\$(find "\$UCE_OUTDIR" -name "*_uce_adata.h5ad" -newer "$INPUT" | head -1)
if [ -n "\$PRODUCED" ]; then
    mv "\$PRODUCED" "\$OUTPUT"
    rm -rf "\$UCE_OUTDIR"
    echo "Moved \$PRODUCED -> \$OUTPUT"
else
    echo "ERROR: UCE did not produce *_uce_adata.h5ad"; ls "\$UCE_OUTDIR"; exit 1
fi
echo "Completed: \$(date)"
SLURM
)
echo "  UCE Job ID: $JOB_UCE"

# ============================================================================
# Geneformer
# ============================================================================
echo "Submitting Geneformer..."
JOB_GF=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=optC_GF_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=$LOG_DIR/geneformer_%j.log
#SBATCH --error=$LOG_DIR/geneformer_%j.err

set -euo pipefail
echo "=== Option C: Geneformer on canonical COVID input ==="
echo "Start: \$(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/geneformer_py310
if [ ! -x "\$ENV_DIR/bin/python" ]; then echo "ERROR: \$ENV_DIR/bin/python not found"; exit 1; fi
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=geneformer_py310
export PYTHONNOUSERSITE=1
echo "Python: \$(which python)  (\$(python --version 2>&1))"

OUTPUT=$EMB_DIR/geneformer_covid_canonical.h5ad
if python -c "import anndata,sys; a=anndata.read_h5ad('\$OUTPUT'); sys.exit(0 if 'X_geneformer' in a.obsm else 1)" 2>/dev/null; then
    echo "Already done. Exit."; exit 0
fi

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit
python -c "
from geneformer import EmbExtractor
import scanpy as sc
import numpy as np

adata = sc.read_h5ad('$INPUT')
print(f'Loaded: {adata.shape}')

embex = EmbExtractor(
    model_type='Pretrained',
    num_classes=0,
    max_ncells=None,
    emb_layer=-1,
    emb_mode='cell',
)
embs = embex.extract_embs(
    model_directory='ctheodoris/Geneformer',
    input_data=adata,
    output_directory='/tmp/gf_optc_covid',
    output_prefix='covid_canonical',
)
adata.obsm['X_geneformer'] = np.array(embs)
print(f'Embeddings: {adata.obsm[\"X_geneformer\"].shape}')
adata.write_h5ad('$EMB_DIR/geneformer_covid_canonical.h5ad')
print(f'Saved: $EMB_DIR/geneformer_covid_canonical.h5ad')
"
echo "Completed: \$(date)"
SLURM
)
echo "  Geneformer Job ID: $JOB_GF"

# ============================================================================
# TranscriptFormer
# ============================================================================
echo "Submitting TranscriptFormer..."
JOB_TF=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=optC_TF_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=$LOG_DIR/transcriptformer_%j.log
#SBATCH --error=$LOG_DIR/transcriptformer_%j.err

set -euo pipefail
echo "=== Option C: TranscriptFormer on canonical COVID input ==="
echo "Start: \$(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/transcriptformer
if [ ! -x "\$ENV_DIR/bin/python" ]; then echo "ERROR: \$ENV_DIR/bin/python not found"; exit 1; fi
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=transcriptformer
export PYTHONNOUSERSITE=1
echo "Python: \$(which python)  (\$(python --version 2>&1))"
which transcriptformer

OUTPUT=$EMB_DIR/transcriptformer_covid_canonical.h5ad
if python -c "import anndata,sys; a=anndata.read_h5ad('\$OUTPUT'); sys.exit(0 if 'X_transcriptformer' in a.obsm else 1)" 2>/dev/null; then
    echo "Already done. Exit."; exit 0
fi

OUTPUT_DIR=$EMB_DIR/tf_temp_canonical
mkdir -p "\$OUTPUT_DIR"
# TF needs ensembl_id column in var. Use TF-ready input (canonical + ensembl_id mapping).
TF_INPUT=/bigdata/godziklab/shared/Xinru/302006/01_DATA/transcriptformer_ready/covid_transcriptformer_ready.h5ad
transcriptformer inference --checkpoint-path /bigdata/godziklab/shared/Xinru/302006/302006_updategit/checkpoints/tf_sapiens --data-file "\$TF_INPUT" --output-path "\$OUTPUT_DIR" --batch-size 16 --emb-type cell --gene-col-name ensembl_id --precision 16-mixed

python -c "
import scanpy as sc
# Use TF-prepped input (47109 cells = same n_obs as TF emb output, drops 868 cells via TF QC)
adata = sc.read_h5ad('\$TF_INPUT')
emb_adata = sc.read_h5ad('\$OUTPUT_DIR/embeddings.h5ad')
if 'embeddings' in emb_adata.obsm:
    adata.obsm['X_transcriptformer'] = emb_adata.obsm['embeddings']
elif 'X_transcriptformer' in emb_adata.obsm:
    adata.obsm['X_transcriptformer'] = emb_adata.obsm['X_transcriptformer']
print(f'Embeddings: {adata.obsm[\"X_transcriptformer\"].shape}')
adata.write_h5ad('$EMB_DIR/transcriptformer_covid_canonical.h5ad')
print(f'Saved')
"
rm -rf "\$OUTPUT_DIR"
echo "Completed: \$(date)"
SLURM
)
echo "  TranscriptFormer Job ID: $JOB_TF"

# ============================================================================
# scPRINT-2 (scprint2 CLI; uses scprint2_env + lamin ontology DB)
# ============================================================================
echo "Submitting scPRINT-2..."
JOB_SCPRINT=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=optC_scPRINT2_covid
#SBATCH --partition=short_gpu
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128gb
#SBATCH --time=2:00:00
#SBATCH --requeue
#SBATCH --output=$LOG_DIR/scprint2_%j.log
#SBATCH --error=$LOG_DIR/scprint2_%j.err

set -euo pipefail
echo "=== Option C: scPRINT-2 on canonical COVID input ==="
echo "Start: \$(date)"

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
if [ ! -x "\$ENV_DIR/bin/python" ]; then echo "ERROR: \$ENV_DIR/bin/python not found. Run install_scprint2.sh first."; exit 1; fi
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
echo "Python: \$(which python)  (\$(python --version 2>&1))"

OUTPUT=$EMB_DIR/scprint2_covid_canonical.h5ad
CKPT="\$HF_HOME/scprint2_ckpts/small-v2.ckpt"
if [ ! -f "\$CKPT" ]; then
    CKPT=\$(find "\$HF_HOME/hub" -name "small-v2.ckpt" 2>/dev/null | head -1)
fi
if [ ! -f "\$CKPT" ]; then
    echo "ERROR: small-v2.ckpt not found. Run install_scprint2.sh first."
    exit 1
fi

if python -c "import anndata,sys; a=anndata.read_h5ad('\$OUTPUT'); sys.exit(0 if any(k.startswith('X_scprint') for k in a.obsm) else 1)" 2>/dev/null; then
    echo "Already done. Exit."; exit 0
fi

scprint2 embed --adata "$INPUT" --ckpt_path "\$CKPT" --species "NCBITaxon:9606" --is_symbol true --output_filename "\$OUTPUT"
echo "Completed: \$(date)"
SLURM
)
echo "  scPRINT-2 Job ID: $JOB_SCPRINT"

echo ""
echo "============================================================"
echo "5 OPTION C JOBS SUBMITTED (scGPT excluded per 2026-05-01 decision)"
echo "============================================================"
echo "  STATE:           $JOB_STATE"
echo "  UCE:             $JOB_UCE"
echo "  Geneformer:      $JOB_GF"
echo "  TranscriptFormer: $JOB_TF"
echo "  scPRINT-2:       $JOB_SCPRINT"
echo ""
echo "Monitor: squeue -u xinruqiu --format='%.10i %.20j %.10P %.10T %.10M %R'"
echo "Outputs: $EMB_DIR/"
