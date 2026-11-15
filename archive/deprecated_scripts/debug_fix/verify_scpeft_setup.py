#!/usr/bin/env python3
"""
Verify setup for scPEFT Fine-Tuning (scGPT + Geneformer PEFT)

Checks that all required components are available before running
scGPT or Geneformer PEFT fine-tuning jobs.

Usage:
    source activate /bigdata/godziklab/shared/Xinru/302006/scpeft_env
    export PYTHONPATH=/bigdata/godziklab/shared/Xinru/302006/03_MODELS/scPEFT:$PYTHONPATH
    python verify_scpeft_setup.py

Author: Xinru Qiu
Date: January 2026
"""

import sys
import os
from pathlib import Path

errors = []

print("=" * 60)
print("scPEFT Fine-Tuning Setup Verification")
print("=" * 60)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
SCPEFT_DIR = BASE_DIR / "03_MODELS/scPEFT"

# ---------------------------------------------------------------
# 1. Check scPEFT repository
# ---------------------------------------------------------------
print("\n1. Checking scPEFT Repository:")
print("-" * 40)
if SCPEFT_DIR.exists():
    print(f"   [OK] scPEFT repo: {SCPEFT_DIR}")
    # Check key subdirectories
    for subdir in ["scgpt", "tutorial_peft"]:
        p = SCPEFT_DIR / subdir
        if p.exists():
            print(f"   [OK] {subdir}/ found")
        else:
            print(f"   [WARN] {subdir}/ not found")
else:
    msg = f"scPEFT repo not found at {SCPEFT_DIR}. Run setup_scpeft.sh first."
    print(f"   [MISSING] {msg}")
    errors.append(msg)

# ---------------------------------------------------------------
# 2. Check scGPT pretrained checkpoint
# ---------------------------------------------------------------
print("\n2. Checking scGPT Pretrained Checkpoint:")
print("-" * 40)

# Check both possible locations
SCGPT_PATHS = [
    SCPEFT_DIR / "save/scGPT_human",
    BASE_DIR / "03_MODELS/scGPT/scGPT_human",
]

scgpt_model_dir = None
for p in SCGPT_PATHS:
    if p.exists():
        scgpt_model_dir = p
        break

if scgpt_model_dir:
    print(f"   [OK] scGPT model dir: {scgpt_model_dir}")
    for fname in ["best_model.pt", "args.json", "vocab.json"]:
        fpath = scgpt_model_dir / fname
        if fpath.exists():
            size_mb = fpath.stat().st_size / (1024 ** 2)
            print(f"   [OK] {fname} ({size_mb:.1f} MB)")
        else:
            msg = f"{fname} missing in {scgpt_model_dir}"
            print(f"   [MISSING] {msg}")
            errors.append(msg)
else:
    msg = "scGPT checkpoint not found. Download via setup_scpeft.sh or manually."
    print(f"   [MISSING] {msg}")
    errors.append(msg)

# ---------------------------------------------------------------
# 3. Check Geneformer pretrained model
# ---------------------------------------------------------------
print("\n3. Checking Geneformer Pretrained Model:")
print("-" * 40)

GENEFORMER_PATHS = [
    BASE_DIR / "geneformer_models/Geneformer",
    BASE_DIR / "03_MODELS/geneformer/geneformer-V2-316M-native",
]

geneformer_dir = None
for p in GENEFORMER_PATHS:
    if p.exists():
        geneformer_dir = p
        break

if geneformer_dir:
    print(f"   [OK] Geneformer model dir: {geneformer_dir}")
    for fname in ["config.json", "model.safetensors"]:
        fpath = geneformer_dir / fname
        if fpath.exists():
            size_mb = fpath.stat().st_size / (1024 ** 2)
            print(f"   [OK] {fname} ({size_mb:.1f} MB)")
        else:
            print(f"   [WARN] {fname} not found (may use alternate name)")
else:
    print("   [WARN] Geneformer pretrained model not found at expected paths")
    print("          scGPT PEFT will still work; Geneformer PEFT may need manual path config")

# Check token dictionary
TOKEN_DICT_PATHS = [
    BASE_DIR / "geneformer_models/Geneformer/geneformer/token_dictionary_gc104M.pkl",
    BASE_DIR / "03_MODELS/geneformer/geneformer-V2-316M-native/token_dictionary.pkl",
]
token_dict_found = False
for p in TOKEN_DICT_PATHS:
    if p.exists():
        print(f"   [OK] Token dictionary: {p}")
        token_dict_found = True
        break
if not token_dict_found:
    print("   [WARN] Token dictionary not found (needed for Geneformer PEFT)")

# ---------------------------------------------------------------
# 4. Check COVID data
# ---------------------------------------------------------------
print("\n4. Checking COVID Data:")
print("-" * 40)

DATA_PATHS = [
    BASE_DIR / "01_DATA/covid_pbmc_adata.h5ad",
    BASE_DIR / "02_EMBEDDINGS/state/covid_state.h5ad",
]

data_file = None
for p in DATA_PATHS:
    if p.exists():
        data_file = p
        break

if data_file:
    size_mb = data_file.stat().st_size / (1024 ** 2)
    print(f"   [OK] Data file: {data_file} ({size_mb:.1f} MB)")
else:
    msg = "COVID data not found at expected paths"
    print(f"   [MISSING] {msg}")
    errors.append(msg)

# ---------------------------------------------------------------
# 5. Check PyTorch and CUDA
# ---------------------------------------------------------------
print("\n5. Checking PyTorch and CUDA:")
print("-" * 40)

try:
    import torch
    print(f"   [OK] PyTorch: {torch.__version__}")
    print(f"   CUDA available: {torch.cuda.is_available()}")
    if torch.cuda.is_available():
        print(f"   GPU: {torch.cuda.get_device_name(0)}")
        mem_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
        print(f"   GPU memory: {mem_gb:.1f} GB")
    else:
        print("   [WARN] No GPU detected (fine-tuning will be slow on CPU)")
except ImportError:
    msg = "PyTorch not installed"
    print(f"   [ERROR] {msg}")
    errors.append(msg)

# ---------------------------------------------------------------
# 6. Check key Python packages
# ---------------------------------------------------------------
print("\n6. Checking Python Packages:")
print("-" * 40)

packages = {
    "transformers": "transformers",
    "peft": "peft",
    "accelerate": "accelerate",
    "scanpy": "scanpy",
    "sklearn": "scikit-learn",
    "loralib": "loralib",
}

for import_name, display_name in packages.items():
    try:
        mod = __import__(import_name)
        version = getattr(mod, "__version__", "?")
        print(f"   [OK] {display_name}: {version}")
    except ImportError:
        msg = f"{display_name} not installed"
        print(f"   [MISSING] {msg}")
        errors.append(msg)

# ---------------------------------------------------------------
# 7. Check scPEFT's bundled scGPT imports
# ---------------------------------------------------------------
print("\n7. Checking scPEFT scGPT Imports:")
print("-" * 40)

scgpt_checks = {
    "TransformerModel": "scgpt.model",
    "GeneVocab": "scgpt.tokenizer.gene_tokenizer",
}

for name, module_path in scgpt_checks.items():
    try:
        mod = __import__(module_path, fromlist=[name])
        cls = getattr(mod, name)
        print(f"   [OK] {name} from {module_path}")
    except Exception as e:
        print(f"   [WARN] {name}: {e}")
        print(f"          Make sure PYTHONPATH includes {SCPEFT_DIR}")

# Check PeftConfig (may or may not exist depending on scPEFT version)
try:
    from scgpt.utils import PeftConfig
    print("   [OK] PeftConfig (native scPEFT PEFT support)")
except ImportError:
    print("   [INFO] PeftConfig not found (will use manual LoRA setup)")

# ---------------------------------------------------------------
# 8. Check results/checkpoint directories
# ---------------------------------------------------------------
print("\n8. Checking Output Directories:")
print("-" * 40)

output_dirs = {
    "Results (scGPT PEFT)": BASE_DIR / "05_RESULTS/scgpt_peft",
    "Results (Geneformer PEFT)": BASE_DIR / "05_RESULTS/geneformer_peft",
    "Logs": BASE_DIR / "06_LOGS",
}

for name, path in output_dirs.items():
    if path.exists():
        print(f"   [OK] {name}: {path}")
    else:
        print(f"   [INFO] {name}: will be created at {path}")

# ---------------------------------------------------------------
# 9. Load vocabulary and test tokenization
# ---------------------------------------------------------------
print("\n9. Testing scGPT Vocabulary:")
print("-" * 40)

if scgpt_model_dir:
    vocab_file = scgpt_model_dir / "vocab.json"
    if vocab_file.exists():
        try:
            from scgpt.tokenizer.gene_tokenizer import GeneVocab
            vocab = GeneVocab.from_file(vocab_file)
            print(f"   [OK] Vocabulary loaded: {len(vocab)} genes")
            # Test a few common gene names
            test_genes = ["TP53", "BRCA1", "EGFR", "MYC", "GAPDH"]
            found = [g for g in test_genes if g in vocab]
            print(f"   Test genes found: {len(found)}/{len(test_genes)} ({', '.join(found)})")
        except Exception as e:
            print(f"   [WARN] Could not load vocabulary: {e}")
    else:
        print("   [SKIP] vocab.json not found")
else:
    print("   [SKIP] scGPT model dir not found")

# ---------------------------------------------------------------
# 10. Memory estimate
# ---------------------------------------------------------------
print("\n10. Memory Estimates for PEFT Training:")
print("-" * 40)
print("   scGPT (51M params, LoRA ~1-2% trainable):")
print("     Full model memory: ~0.5 GB")
print("     Training with LoRA: ~2-4 GB (+ data)")
print("     Recommended: 16+ GB GPU")
print("   Geneformer (30M params, LoRA ~1-2% trainable):")
print("     Full model memory: ~0.3 GB")
print("     Training with LoRA: ~1-3 GB (+ data)")
print("     Recommended: 16+ GB GPU")

# ---------------------------------------------------------------
# Summary
# ---------------------------------------------------------------
print("\n" + "=" * 60)
if errors:
    print(f"SETUP INCOMPLETE: {len(errors)} error(s) found")
    for e in errors:
        print(f"  - {e}")
    print("=" * 60)
    sys.exit(1)
else:
    print("ALL CHECKS PASSED! Ready for scPEFT fine-tuning.")
    print("=" * 60)
    print("\nTo run scGPT PEFT:")
    print("  sbatch submit_scgpt_peft.sh")
    print("  # or: python scgpt_peft_finetune.py --task binary --peft LORA --epochs 20")
    print("\nTo run Geneformer PEFT:")
    print("  sbatch submit_geneformer_peft.sh")
    print("  # or: python geneformer_peft_finetune.py --task binary --epochs 25")
    print("\nOr use the wrapper:")
    print("  python scpeft_finetune.py --backbone scgpt --task binary")
