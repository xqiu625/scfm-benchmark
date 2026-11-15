#!/bin/bash
# Generate all data-driven figures for the manuscript / Biohub talk.
# All scripts use hardcoded validated values from the manuscript tables —
# no SLURM job needed; runs locally on any machine with matplotlib + seaborn.
#
# Fig 1 (study overview) is a schematic — see reports/manuscript/figures/fig1/fig1_overview_spec.md
# Fig 7 (per-FM mouse UMAP grid) is a separate SLURM job (run_fig7_umap_grid.sh)
#       because it requires scanpy + UMAP computation on the mouse embeddings.

set -euo pipefail
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean
PYTHON=$ENV_DIR/bin/python
cd "$(dirname "$0")/../../.."

echo "=== Generating manuscript figures ==="
echo "Start: $(date)"

$PYTHON scripts/26_figures/code/generate_fig2_mortality.py    --output-dir reports/manuscript/figures/fig2
$PYTHON scripts/26_figures/code/generate_fig3_crossdisease.py --output-dir reports/manuscript/figures/fig3
$PYTHON scripts/26_figures/code/generate_fig4_biology.py      --output-dir reports/manuscript/figures/fig4
$PYTHON scripts/26_figures/code/generate_fig5_mouse_quant.py  --output-dir reports/manuscript/figures/fig5
$PYTHON scripts/26_figures/code/generate_fig8_lora.py         --output-dir reports/manuscript/figures/fig8
$PYTHON scripts/26_figures/code/generate_fig9_sanity.py       --output-dir reports/manuscript/figures/fig9

echo
echo "Generated:"
find reports/manuscript/figures -name "fig*.png" | sort

echo
echo "Completed: $(date)"
