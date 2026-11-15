#!/usr/bin/env python3
"""
Verify setup for STATE SE-600M True Fine-Tuning

This script checks that all required components are available and compatible
before running the full fine-tuning job.

Usage:
    python verify_true_finetune_setup.py

Author: Xinru Qiu
Date: January 2026
"""

import sys
import os
from pathlib import Path

# Patch torch.load before any imports
import torch
_original_torch_load = torch.load
def _patched_torch_load(f, *args, **kwargs):
    kwargs['weights_only'] = False
    return _original_torch_load(f, *args, **kwargs)
torch.load = _patched_torch_load

print("=" * 60)
print("STATE SE-600M True Fine-Tuning Setup Verification")
print("=" * 60)

# Check paths
BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
MODEL_DIR = BASE_DIR / "03_MODELS/SE-600M"
EMBEDDING_DIR = BASE_DIR / "02_EMBEDDINGS/state"

required_files = {
    'SE-600M Checkpoint': MODEL_DIR / "se600m_epoch16.ckpt",
    'Protein Embeddings': MODEL_DIR / "protein_embeddings.pt",
    'COVID Embeddings (raw)': EMBEDDING_DIR / "covid_state.h5ad",
}

print("\n1. Checking Required Files:")
print("-" * 40)
all_exist = True
for name, path in required_files.items():
    if path.exists():
        size = path.stat().st_size / (1024**3)  # GB
        print(f"   [OK] {name}: {path.name} ({size:.2f} GB)")
    else:
        print(f"   [MISSING] {name}: {path}")
        all_exist = False

if not all_exist:
    print("\nERROR: Some required files are missing!")
    sys.exit(1)

# Check PyTorch and CUDA
print("\n2. Checking PyTorch and CUDA:")
print("-" * 40)
print(f"   PyTorch version: {torch.__version__}")
print(f"   CUDA available: {torch.cuda.is_available()}")
if torch.cuda.is_available():
    print(f"   CUDA device: {torch.cuda.get_device_name(0)}")
    print(f"   GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")

# Check STATE imports
print("\n3. Checking STATE (arc-state) Installation:")
print("-" * 40)
try:
    from state.emb import _setup_vci_aliases
    _setup_vci_aliases()

    import state.emb as emb
    sys.modules['vci'] = emb
    sys.modules['vci.nn'] = emb.nn
    sys.modules['vci.nn.model'] = emb.nn.model

    from state.emb.nn.model import StateEmbeddingModel
    print("   [OK] StateEmbeddingModel imported successfully")
except ImportError as e:
    print(f"   [ERROR] Failed to import STATE: {e}")
    sys.exit(1)

# Load and verify protein embeddings
print("\n4. Checking Protein Embeddings Format:")
print("-" * 40)
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
protein_embeds = torch.load(str(required_files['Protein Embeddings']), map_location='cpu', weights_only=False)

if isinstance(protein_embeds, dict):
    n_genes = len(protein_embeds)
    sample_key = list(protein_embeds.keys())[0]
    sample_emb = protein_embeds[sample_key]
    print(f"   [OK] Dict format with {n_genes:,} genes")
    print(f"   Embedding dimension: {sample_emb.shape}")
    print(f"   Sample genes: {list(protein_embeds.keys())[:5]}")
elif hasattr(protein_embeds, 'shape'):
    print(f"   Tensor format: {protein_embeds.shape}")
else:
    print(f"   Unknown format: {type(protein_embeds)}")

# Load and verify COVID data
print("\n5. Checking COVID Data Format:")
print("-" * 40)
import scanpy as sc
import numpy as np

adata = sc.read_h5ad(required_files['COVID Embeddings (raw)'])
print(f"   Shape: {adata.shape}")
print(f"   obs columns: {list(adata.obs.columns)[:10]}")
print(f"   obsm keys: {list(adata.obsm.keys())}")
print(f"   var_names sample: {list(adata.var_names[:5])}")

# Check for raw counts
if hasattr(adata, 'raw') and adata.raw is not None:
    print(f"   Has .raw layer: Yes ({adata.raw.X.shape})")
else:
    print(f"   Has .raw layer: No (will use adata.X)")

# Check severity labels
severity_col = None
for col in ['covid_severity', 'severity']:
    if col in adata.obs.columns:
        severity_col = col
        break

if severity_col:
    severities = adata.obs[severity_col].value_counts()
    print(f"   Severity column: {severity_col}")
    print(f"   Label distribution:")
    for label, count in severities.items():
        print(f"      {label}: {count:,}")
else:
    print("   [WARNING] No severity column found!")

# Check gene overlap with protein embeddings
if isinstance(protein_embeds, dict):
    gene_names = [g.upper() for g in adata.var_names]
    overlap = sum(1 for g in gene_names if g in protein_embeds)
    print(f"\n6. Gene Overlap with ESM2 Embeddings:")
    print("-" * 40)
    print(f"   Genes in data: {len(gene_names):,}")
    print(f"   Genes with ESM2 embedding: {overlap:,}")
    print(f"   Coverage: {100*overlap/len(gene_names):.1f}%")

# Load SE model
print("\n7. Loading SE-600M Model:")
print("-" * 40)
try:
    checkpoint = torch.load(str(required_files['SE-600M Checkpoint']), map_location=device, weights_only=False)
    hparams = checkpoint['hyper_parameters']
    print(f"   d_model: {hparams['d_model']}")
    print(f"   nlayers: {hparams['nlayers']}")
    print(f"   nhead: {hparams['nhead']}")

    model = StateEmbeddingModel(**hparams)
    result = model.load_state_dict(checkpoint['state_dict'], strict=False)
    n_params = sum(p.numel() for p in model.parameters()) / 1e6
    print(f"   [OK] Model loaded: {n_params:.1f}M parameters")
    print(f"   Missing keys: {len(result.missing_keys)}")
    print(f"   Unexpected keys: {len(result.unexpected_keys)}")
except Exception as e:
    print(f"   [ERROR] Failed to load model: {e}")
    sys.exit(1)

# Quick forward pass test
print("\n8. Testing Forward Pass:")
print("-" * 40)
try:
    model = model.to(device)
    model.eval()

    # Create dummy input matching the expected format
    batch_size = 2
    max_genes = 100
    embed_dim = hparams.get('token_dim', 5120)
    d_model = hparams['d_model']

    # Test encoder
    dummy_gene_embeds = torch.randn(batch_size * max_genes, embed_dim).to(device)
    with torch.no_grad():
        encoded = model.encoder(dummy_gene_embeds)
    print(f"   Encoder: {embed_dim} → {encoded.shape[-1]} [OK]")

    # Test transformer
    encoded = encoded.view(batch_size, max_genes, d_model)
    encoded_t = encoded.transpose(0, 1)  # [seq, batch, dim]
    with torch.no_grad():
        transformed = model.transformer_encoder(encoded_t)
    print(f"   Transformer: {encoded_t.shape} → {transformed.shape} [OK]")

    del model, encoded, transformed, dummy_gene_embeds
    torch.cuda.empty_cache() if torch.cuda.is_available() else None

except Exception as e:
    print(f"   [ERROR] Forward pass failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# Memory estimate
print("\n9. Memory Estimate for Training:")
print("-" * 40)
n_params_total = 603e6  # 603M parameters
bytes_per_param = 4  # float32
forward_memory_multiplier = 2  # activations
backward_memory_multiplier = 4  # gradients + optimizer states

estimated_memory_gb = n_params_total * bytes_per_param * backward_memory_multiplier / 1e9
print(f"   Model parameters: 603M")
print(f"   Estimated training memory: ~{estimated_memory_gb:.1f} GB")
print(f"   Recommended GPU memory: 40GB+ (A100)")

print("\n" + "=" * 60)
print("ALL CHECKS PASSED! Ready for true fine-tuning.")
print("=" * 60)
print("\nTo run fine-tuning:")
print("  sbatch submit_true_finetune.sh")
print("\nOr interactively:")
print("  python state_true_finetune_v3.py --task binary --epochs 15")
