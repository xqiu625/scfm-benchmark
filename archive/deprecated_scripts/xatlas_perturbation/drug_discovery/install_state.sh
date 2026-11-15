#!/bin/bash
# Install STATE (Arc Institute) for perturbation prediction
# https://github.com/arc-institute/state

echo "=== Installing STATE Model ==="

# Option 1: Using uv (recommended)
if command -v uv &> /dev/null; then
    echo "Installing STATE with uv..."
    uv tool install arc-state
    echo "STATE installed via uv"
else
    echo "uv not found, using pip..."
    # Option 2: Using pip
    pip install arc-state
fi

# Verify installation
echo ""
echo "=== Verifying Installation ==="
state --version

# Download model checkpoints
echo ""
echo "=== Downloading Model Checkpoints ==="
echo "Available models:"
echo "  - SE-600M: State Embedding model (600M params)"
echo "  - ST-600M: State Transition model (600M params)"

# Create model directory
MODEL_DIR="${HOME}/.cache/state_models"
mkdir -p $MODEL_DIR

echo ""
echo "To download models, run:"
echo "  state download SE-600M --output-dir $MODEL_DIR"
echo "  state download ST-600M --output-dir $MODEL_DIR"
echo ""
echo "Installation complete!"
