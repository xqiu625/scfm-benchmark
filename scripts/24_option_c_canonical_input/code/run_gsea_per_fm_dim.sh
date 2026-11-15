#!/bin/bash
#SBATCH --job-name=gsea_per_fm
#SBATCH --partition=batch
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gsea_per_fm_%j.out
#SBATCH --error=/bigdata/godziklab/shared/Xinru/302006/06_LOGS/gsea_per_fm_%j.err
#SBATCH --mail-user=xqiu014@ucr.edu
#SBATCH --mail-type=END,FAIL

set -euo pipefail
echo "=== Per-FM GSEA on top mortality dimensions ==="
echo "Start: $(date)"

ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
if [ ! -x "$ENV_DIR/bin/python" ]; then echo "ERROR: env not found"; exit 1; fi
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export CONDA_PREFIX="$ENV_DIR"
export PYTHONNOUSERSITE=1

# Verify gseapy installed (install if missing)
python -c "import gseapy" 2>/dev/null || \
    uv pip install --python "$ENV_DIR/bin/python" gseapy

cd /bigdata/godziklab/shared/Xinru/302006/302006_updategit

python scripts/24_option_c_canonical_input/code/gsea_per_fm_dim.py \
    --input-dir  /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_bio_interpretation \
    --output-dir /bigdata/godziklab/shared/Xinru/302006/05_RESULTS/option_c_gsea \
    --libraries  KEGG_2021_Human Reactome_2022 GO_Biological_Process_2023 \
    --top-n      0

echo "Completed: $(date)"
