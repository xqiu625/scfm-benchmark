#!/bin/bash
# ==============================================================================
# Job F Phase 3: embed the monocyte canonical object with the same 5 FMs used
# for the platelet mortality benchmark (STATE, UCE, Geneformer, scPRINT-2,
# scGPT -- matches donor_level_reanalysis.py's FMS list, NOT TranscriptFormer).
#
# Adapted from scripts/24_option_c_canonical_input/code/step2c_human_covid_canonical_input.sh
# (STATE/UCE/Geneformer/scPRINT-2 blocks) and
# scripts/foundation_models/scGPT/submit_scgpt_option_c.sh (scGPT block).
#
# Input:  02_EMBEDDINGS/option_c_monocyte/monocyte_canonical_raw.h5ad
# Output: 02_EMBEDDINGS/option_c_monocyte/<fm>_monocyte_canonical.h5ad
# ==============================================================================
set -e

BASE_DIR="/bigdata/godziklab/shared/Xinru/302006"
INPUT="$BASE_DIR/02_EMBEDDINGS/option_c_monocyte/monocyte_canonical_raw.h5ad"
EMB_DIR="$BASE_DIR/02_EMBEDDINGS/option_c_monocyte"
LOG_DIR="$BASE_DIR/06_LOGS/option_c_monocyte"
mkdir -p "$EMB_DIR" "$LOG_DIR"

if [ ! -f "$INPUT" ]; then
    echo "ERROR: canonical monocyte input not found: $INPUT"
    exit 1
fi

echo "============================================================"
echo "JOB F PHASE 3: 5 FM EMBEDDING JOBS ON MONOCYTE CANONICAL INPUT"
echo "============================================================"
echo "Input:      $INPUT"
echo "Output dir: $EMB_DIR"

# ============================================================================
# STATE
# ============================================================================
echo "Submitting STATE..."
JOB_STATE=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=jobF_STATE_mono
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
echo "=== Job F: STATE on monocyte canonical input ==="
echo "Start: \$(date)"
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=state_env
export PYTHONNOUSERSITE=1
export TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD=1
which state
OUTPUT=$EMB_DIR/state_monocyte_canonical.h5ad
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
#SBATCH --job-name=jobF_UCE_mono
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
echo "=== Job F: UCE on monocyte canonical input ==="
echo "Start: \$(date)"
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/state_env
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=state_env
export PYTHONNOUSERSITE=1

OUTPUT=$EMB_DIR/uce_monocyte_canonical.h5ad
UCE_DIR=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/uce/UCE
UCE_OUTDIR=$EMB_DIR/uce_run
mkdir -p "\$UCE_OUTDIR"
cd "\$UCE_DIR"
python eval_single_anndata.py --adata_path "$INPUT" --dir "\$UCE_OUTDIR" --model_loc "\$UCE_DIR/model_files/4layer_model.torch" --species human --batch_size 32
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
#SBATCH --job-name=jobF_GF_mono
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
echo "=== Job F: Geneformer on monocyte canonical input ==="
echo "Start: \$(date)"
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/geneformer_py310
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export CONDA_DEFAULT_ENV=geneformer_py310
export PYTHONNOUSERSITE=1

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
    output_directory='/tmp/gf_optc_monocyte',
    output_prefix='monocyte_canonical',
)
adata.obsm['X_geneformer'] = np.array(embs)
print(f'Embeddings: {adata.obsm[\"X_geneformer\"].shape}')
adata.write_h5ad('$EMB_DIR/geneformer_monocyte_canonical.h5ad')
print('Saved: $EMB_DIR/geneformer_monocyte_canonical.h5ad')
"
echo "Completed: \$(date)"
SLURM
)
echo "  Geneformer Job ID: $JOB_GF"

# ============================================================================
# scPRINT-2
# ============================================================================
echo "Submitting scPRINT-2..."
JOB_SCPRINT=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=jobF_scPRINT2_mono
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
echo "=== Job F: scPRINT-2 on monocyte canonical input ==="
echo "Start: \$(date)"
ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
export PATH="\$ENV_DIR/bin:\$PATH"
export LD_LIBRARY_PATH="\$ENV_DIR/lib:\${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="\$ENV_DIR"
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface

OUTPUT=$EMB_DIR/scprint2_monocyte_canonical.h5ad
CKPT="\$HF_HOME/scprint2_ckpts/small-v2.ckpt"
if [ ! -f "\$CKPT" ]; then
    CKPT=\$(find "\$HF_HOME/hub" -name "small-v2.ckpt" 2>/dev/null | head -1)
fi
scprint2 embed --adata "$INPUT" --ckpt_path "\$CKPT" --species "NCBITaxon:9606" --is_symbol true --output_filename "\$OUTPUT"
echo "Completed: \$(date)"
SLURM
)
echo "  scPRINT-2 Job ID: $JOB_SCPRINT"

# ============================================================================
# scGPT (singularity container, matches submit_scgpt_option_c.sh)
# ============================================================================
echo "Submitting scGPT..."
JOB_SCGPT=$(sbatch --parsable <<SLURM
#!/bin/bash
#SBATCH --job-name=jobF_scGPT_mono
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:a100:1
#SBATCH --exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=4:00:00
#SBATCH --output=$LOG_DIR/scgpt_%j.out
#SBATCH --error=$LOG_DIR/scgpt_%j.err

set -euo pipefail
echo "Host: \$(hostname) | Start: \$(date)"
BASE=/bigdata/godziklab/shared/Xinru/302006
REPO=\$BASE/302006_updategit
SIF=\$BASE/07_CONTAINERS/scgpt_0.1.7.sif
MODEL=\$BASE/models/scGPT_human/scGPT_human
SCRIPT=\$REPO/scripts/foundation_models/scGPT/embed_scgpt_option_c.py
OUTPUT=$EMB_DIR/scgpt_monocyte_canonical.h5ad

module load singularity 2>/dev/null || echo "using system singularity"
export SINGULARITYENV_PYTHONNOUSERSITE=1
singularity exec --nv --bind "\$BASE" --env PYTHONNOUSERSITE=1 "\$SIF" \
    python "\$SCRIPT" --input "$INPUT" --model "\$MODEL" --output "\$OUTPUT" --batch-size 32
echo "Done: \$(date) -> \$OUTPUT"
SLURM
)
echo "  scGPT Job ID: $JOB_SCGPT"

echo ""
echo "============================================================"
echo "5 JOB F EMBEDDING JOBS SUBMITTED"
echo "============================================================"
echo "  STATE:      $JOB_STATE"
echo "  UCE:        $JOB_UCE"
echo "  Geneformer: $JOB_GF"
echo "  scPRINT-2:  $JOB_SCPRINT"
echo "  scGPT:      $JOB_SCGPT"
echo ""
echo "Outputs: $EMB_DIR/"
