#!/usr/bin/env python3
"""
Process HCT116 RETNLA data with scGPT - Standalone script
"""

import sys
sys.path.append('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/scGPT')

import scgpt as scg
import scanpy as sc
import numpy as np
import torch
import json
import time
import gc
import os
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

# Configuration
BASE_DIR = "/bigdata/godziklab/shared/Xinru/302006"
DATA_DIR = f"{BASE_DIR}/data/raw"
RETNLA_DIR = f"{BASE_DIR}/retnla"
MODEL_DIR = f"{BASE_DIR}/models/scGPT_human/scGPT_human"

print("Processing HCT116 RETNLA data with scGPT...")

def load_and_prepare_hct116():
    """Load and prepare HCT116 dataset"""
    
    print("Loading HCT116 dataset...")
    data_path = f"{DATA_DIR}/HCT116_retn_raw_counts.h5ad"
    
    if not os.path.exists(data_path):
        raise FileNotFoundError(f"HCT116 file not found: {data_path}")
    
    adata = sc.read_h5ad(data_path)
    adata.obs['dataset'] = 'hct116'
    adata.obs['species'] = 'human'
    adata.obs['gene_target'] = 'RETN'
    
    print(f"Loaded HCT116: {adata.shape}")
    
    # Basic preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After filtering: {adata.shape}")
    
    return adata

def process_hct116():
    """Process HCT116 dataset with scGPT"""
    
    start_time = time.time()
    
    # Load data
    adata = load_and_prepare_hct116()
    print(f"Processing {adata.n_obs} cells with {adata.n_vars} genes")
    
    # Load scGPT model components
    model_dir = Path(MODEL_DIR)
    vocab_file = model_dir / "vocab.json"
    model_file = model_dir / "best_model.pt"
    args_file = model_dir / "args.json"
    
    if not vocab_file.exists():
        raise FileNotFoundError(f"Vocabulary file not found: {vocab_file}")
    if not model_file.exists():
        raise FileNotFoundError(f"Model file not found: {model_file}")
    if not args_file.exists():
        raise FileNotFoundError(f"Args file not found: {args_file}")
    
    print(f"Loading vocabulary from {vocab_file}")
    vocab = GeneVocab.from_file(vocab_file)
    
    # Add special tokens
    special_tokens = ["<pad>", "<cls>", "<eoc>"]
    for s in special_tokens:
        if s not in vocab:
            vocab.append_token(s)
    vocab.set_default_index(vocab["<pad>"])
    
    print(f"Vocabulary size: {len(vocab):,}")
    
    # Match genes to vocabulary
    gene_names = adata.var.index.tolist()
    gene_ids = []
    valid_indices = []
    
    print("Matching genes to scGPT vocabulary...")
    for idx, gene in enumerate(gene_names):
        if gene in vocab:
            gene_ids.append(vocab[gene])
            valid_indices.append(idx)
    
    vocab_match_rate = len(gene_ids) / len(gene_names) * 100
    print(f"Genes matched to vocabulary: {len(gene_ids):,}/{len(gene_names):,} ({vocab_match_rate:.1f}%)")
    
    if len(gene_ids) < 100:
        print("WARNING: Very low vocabulary match rate!")
    
    # Filter to vocabulary genes
    adata_vocab = adata[:, valid_indices].copy()
    print(f"Using {adata_vocab.n_vars:,} vocabulary-matched genes")
    
    # Get expression data
    X = adata_vocab.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    
    print(f"Expression matrix: {X.shape}")
    print(f"Expression range: [{np.min(X):.3f}, {np.max(X):.3f}]")
    
    # Binning
    print(f"Binning {X.shape[1]:,} genes...")
    X_binned = np.zeros_like(X, dtype=int)
    
    for i in range(X.shape[1]):
        if i % 500 == 0 and i > 0:
            print(f"   Gene {i:,}/{X.shape[1]:,}")
        
        gene_expr = X[:, i]
        if np.max(gene_expr) > 0:
            non_zero = gene_expr[gene_expr > 0]
            if len(non_zero) > 10:
                # Use quantile-based binning for genes with sufficient expression
                try:
                    bins = np.quantile(non_zero, np.linspace(0, 1, min(51, len(non_zero)//2)))
                    bins = np.unique(bins)
                    if len(bins) > 1:
                        X_binned[:, i] = np.digitize(gene_expr, bins)
                    else:
                        X_binned[:, i] = (gene_expr > 0).astype(int)
                except:
                    # Fallback to simple thresholding
                    X_binned[:, i] = (gene_expr > 0).astype(int)
            else:
                # Simple thresholding for sparsely expressed genes
                X_binned[:, i] = (gene_expr > 0).astype(int)
    
    print(f"Binning complete. Range: [{np.min(X_binned)}, {np.max(X_binned)}]")
    
    # Load scGPT model
    print("Loading scGPT model...")
    with open(args_file, "r") as f:
        args = json.load(f)
    
    model = TransformerModel(
        len(vocab),
        args.get("embsize", 512),
        args.get("nheads", 8),
        args.get("d_hid", 512),
        args.get("nlayers", 12),
        vocab=vocab,
        dropout=0.0,
        pad_token="<pad>",
        pad_value=0
    )
    
    # Load model weights
    try:
        state_dict = torch.load(model_file, map_location='cpu')
        model.load_state_dict(state_dict, strict=False)
        model.eval()
        print("Model loaded successfully")
    except Exception as e:
        print(f"Error loading model: {e}")
        raise
    
    # Process cells in batches
    n_cells = X_binned.shape[0]
    embeddings = []
    batch_size = 10  # Small batches to avoid memory issues
    
    print(f"Processing {n_cells:,} cells in batches of {batch_size}...")
    
    with torch.no_grad():
        for batch_start in range(0, n_cells, batch_size):
            batch_end = min(batch_start + batch_size, n_cells)
            
            print(f"   Batch {batch_start+1:,}-{batch_end:,}/{n_cells:,}")
            
            batch_embeddings = []
            
            for i in range(batch_start, batch_end):
                try:
                    # Prepare tensors
                    gene_tensor = torch.tensor(gene_ids, dtype=torch.long).unsqueeze(0)
                    value_tensor = torch.tensor(X_binned[i:i+1], dtype=torch.float)
                    padding_mask = torch.zeros_like(value_tensor, dtype=torch.bool)
                    
                    # Generate embedding
                    output = model(gene_tensor, value_tensor, src_key_padding_mask=padding_mask)
                    
                    # Extract cell embedding
                    if hasattr(output, 'cell_emb') and output.cell_emb is not None:
                        emb = output.cell_emb.cpu().numpy()
                    elif isinstance(output, dict) and 'cell_emb' in output:
                        emb = output['cell_emb'].cpu().numpy()
                    else:
                        # Fallback: use mean of sequence embeddings
                        emb = output.mean(dim=1).cpu().numpy()
                    
                    batch_embeddings.append(emb)
                    
                    # Clean up tensors
                    del gene_tensor, value_tensor, padding_mask, output
                    
                except Exception as e:
                    print(f"   CRITICAL ERROR: Cell {i-batch_start+1} in batch failed: {e}")
                    raise RuntimeError(
                        f"Cell {i-batch_start+1} embedding generation failed. "
                        f"Cannot use synthetic data. Original error: {e}"
                    )
            
            # Add batch to main list
            embeddings.extend(batch_embeddings)
            
            # Clean up batch
            del batch_embeddings
            gc.collect()
            
            # Progress update
            elapsed = time.time() - start_time
            rate = batch_end / elapsed if elapsed > 0 else 0
            remaining = n_cells - batch_end
            if rate > 0 and remaining > 0:
                est_remaining = remaining / rate
                print(f"     {elapsed:.1f}s elapsed, {rate:.1f}/s, est {est_remaining:.1f}s remaining")
    
    # Combine embeddings
    final_embeddings = np.vstack(embeddings)
    processing_time = time.time() - start_time
    
    print(f"Generated embeddings: {final_embeddings.shape} in {processing_time:.1f}s")
    print(f"Processing speed: {n_cells/processing_time:.1f} cells/second")
    print(f"Embedding stats: mean={np.mean(final_embeddings):.3f}, std={np.std(final_embeddings):.3f}")
    
    # Create result AnnData
    result_adata = sc.AnnData(X=final_embeddings, obs=adata.obs.copy())
    result_adata.uns['scgpt_processing'] = {
        'dataset': 'hct116',
        'model': 'scGPT_human',
        'embedding_dim': final_embeddings.shape[1],
        'genes_used': len(gene_ids),
        'vocab_match_rate': vocab_match_rate,
        'processing_time_seconds': processing_time,
        'cells_per_second': n_cells / processing_time,
        'total_genes_original': adata.n_vars,
        'genes_after_filtering': len(gene_names),
        'genes_matched_vocab': len(gene_ids)
    }
    
    # Create output directory
    os.makedirs(RETNLA_DIR, exist_ok=True)
    
    # Save results
    result_path = f"{RETNLA_DIR}/hct116_scgpt_embeddings.h5ad"
    result_adata.write_h5ad(result_path)
    print(f"Saved: {result_path}")
    
    # Save embeddings as numpy array
    np.save(f"{RETNLA_DIR}/hct116_scgpt_embeddings.npy", final_embeddings)
    print(f"Saved: {RETNLA_DIR}/hct116_scgpt_embeddings.npy")
    
    # Save metadata separately
    adata.obs.to_csv(f"{RETNLA_DIR}/hct116_scgpt_metadata.csv")
    print(f"Saved: {RETNLA_DIR}/hct116_scgpt_metadata.csv")
    
    # Print summary
    print("\n" + "="*50)
    print("HCT116 scGPT Processing Summary")
    print("="*50)
    print(f"Input cells: {n_cells}")
    print(f"Original genes: {adata.n_vars}")
    print(f"Genes after filtering: {len(gene_names)}")
    print(f"Genes matched to vocab: {len(gene_ids)} ({vocab_match_rate:.1f}%)")
    print(f"Final embedding shape: {final_embeddings.shape}")
    print(f"Processing time: {processing_time:.1f} seconds")
    print(f"Processing rate: {n_cells/processing_time:.1f} cells/second")
    print(f"Results saved to: {RETNLA_DIR}/")
    print("="*50)
    
    # Cleanup
    del model, adata, embeddings, result_adata
    gc.collect()
    
    return True

# Execute processing
if __name__ == "__main__":
    try:
        success = process_hct116()
        
        if success:
            print("\nSUCCESS: HCT116 scGPT processing completed!")
            print("Next steps:")
            print("  - Process HEK293T: python process_hek293t_scgpt.py")
            print("  - Process Mouse: python process_mouse_scgpt.py")
            print("  - Combine all: python process_retnla_scgpt_individual.py --combine-only")
            
        else:
            print("\nFAILED: HCT116 scGPT processing failed")
            sys.exit(1)
            
    except Exception as e:
        print(f"\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)