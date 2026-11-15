#!/bin/bash
# Setup script for scPRINT-2 (Kalfon/Object lab, jkobject/scPRINT v2)
# Recommended: uv (~2 min). Conda also works but slower.
#
# Run on a compute node (login node will OOM during conda solve if you fall back):
#   srun -p batch --mem=32gb -c 4 --time=2:00:00 --pty bash -l
#   bash scripts/foundation_models/scPRINT/install_scprint2.sh

set -e

ENV_DIR=/bigdata/godziklab/shared/Xinru/302006/scprint2_env
LAMIN_DB_DIR=/bigdata/godziklab/shared/Xinru/302006/lamin_testdb

echo "=== Installing scPRINT-2 ==="

# 1. Try uv first (recommended by scPRINT-2 maintainers)
if command -v uv >/dev/null 2>&1; then
    echo "uv found: $(which uv)"
    uv venv "$ENV_DIR" --python 3.11
else
    echo "uv not found. Install uv first or fall back to conda."
    echo "Install uv: curl -LsSf https://astral.sh/uv/install.sh | sh"
    echo ""
    echo "Falling back to conda env creation..."
    conda create -p "$ENV_DIR" -c conda-forge --override-channels python=3.11 -y
fi

# 2. Activate via PATH (env-vars approach, robust to SLURM batch shells)
export PATH="$ENV_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$ENV_DIR/lib:${LD_LIBRARY_PATH:-}"
export PYTHONNOUSERSITE=1
export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface
export PIP_CACHE_DIR=/bigdata/godziklab/shared/Xinru/.pip-cache
mkdir -p "$HF_HOME" "$PIP_CACHE_DIR"

echo "Python: $(which python)"
python --version

# 3. Install scprint2 + supporting tools (uv venv doesn't ship pip — use uv for all installs)
if command -v uv >/dev/null 2>&1; then
    uv pip install scprint2
    # huggingface_hub + pyyaml are already pulled in as scprint2 deps; verify, install if missing
    python -c "import huggingface_hub, yaml" 2>/dev/null || uv pip install huggingface_hub pyyaml
else
    pip install --upgrade pip
    pip install scprint2 huggingface_hub pyyaml
fi

# 5. lamin init + scPRINT-2 ontology setup (one-time)
mkdir -p "$LAMIN_DB_DIR"
cd "$(dirname "$LAMIN_DB_DIR")"

if [ ! -d "$LAMIN_DB_DIR/.lamindb" ]; then
    echo "Initializing lamin..."
    lamin init --storage "$LAMIN_DB_DIR" --name testdb --modules bionty
else
    echo "lamin already initialized at $LAMIN_DB_DIR"
fi

echo "Running scprint2 easy_setup (populates ontologies — 2-10 min)..."
scprint2 easy_setup || echo "easy_setup may have warnings; continue if non-fatal"

# 6. Download scPRINT-2 checkpoint (small-v2.ckpt; v2-medium not yet released
# as of 2026-05-01 per HF API check — only small-v2.ckpt is available)
echo "Downloading small-v2.ckpt..."
hf download jkobject/scPRINT small-v2.ckpt --local-dir "$HF_HOME/scprint2_ckpts" || \
    hf download jkobject/scPRINT small-v2.ckpt

# 7. Persist env vars for future activations
mkdir -p "$ENV_DIR/etc/conda/activate.d"
{
    echo 'export PYTHONNOUSERSITE=1'
    echo 'export HF_HOME=/bigdata/godziklab/shared/Xinru/.huggingface'
    echo 'export PIP_CACHE_DIR=/bigdata/godziklab/shared/Xinru/.pip-cache'
    echo "export LAMIN_DB_DIR=$LAMIN_DB_DIR"
} > "$ENV_DIR/etc/conda/activate.d/env_vars.sh"

# 8. Sanity check
python -c "
from scprint2 import scPRINT2
print('scPRINT2 import OK')
print('class:', scPRINT2)
"

echo ""
echo "=== scPRINT-2 setup complete ==="
echo "Env:        $ENV_DIR"
echo "Lamin DB:   $LAMIN_DB_DIR"
echo "Checkpoint: $HF_HOME/scprint2_ckpts/small-v2.ckpt"
echo ""
echo "Next: bash scripts/foundation_models/scPRINT/run_scprint2_human_covid.sh"
