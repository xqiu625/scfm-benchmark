#!/usr/bin/env python3
"""
Process HEK293T data with scGPT - Direct container execution
"""

import sys
sys.path.append('/bigdata/godziklab/shared/Xinru/302006/302006_updategit/scripts/scGPT')

import scanpy as sc
import numpy as np
import torch
import json
import time
import gc
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

from scgpt.model import TransformerModel
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.utils import set_seed

print("Processing HEK293T data with scGPT container...")

# Configuration
BASE_DIR = "/bigdata/godziklab/shared/Xinru/302006"
DATA_DIR = f"{BASE_DIR}/data/raw"
RETNLA_DIR = f"{BASE_DIR}/retnla"
MODEL_DIR = f"{BASE_DIR}/models/scGPT_human/scGPT_human"

def process_hek293t():
    """Process HEK293T dataset with scGPT"""
    
    start_time = time.time()
    
    # Load HEK293T data
    print("Loading HEK293T dataset...")
    data_path = f"{DATA_DIR}/HEK293T_retn_raw_counts.h5ad"
    adata = sc.read_h5ad(data_path)
    
    # Add metadata
    adata.obs['dataset'] = 'hek293t'
    adata.obs['species'] = 'human'
    adata.obs['gene_target'] = 'RETN'
    
    print(f"Loaded HEK293T data: {adata.shape}")
    print(f"Species: {adata.obs['species'].iloc[0]}")
    print(f"Gene target: {adata.obs['gene_target'].iloc[0]}")
    
    # Basic preprocessing
    print("Preprocessing...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"After filtering: {adata.shape}")
    
    # Load scGPT model components
    model_dir = Path(MODEL_DIR)
    vocab_file = model_dir / "vocab.json"
    model_file = model_dir / "best_model.pt"
    args_file = model_dir / "args.json"
    
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
                except Exception:
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
    
    # Process cells in smaller batches
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
            
            elapsed = time.time() - start_time
            rate = batch_end / elapsed
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
        'dataset': 'hek293t',
        'model': 'scGPT_human',
        'embedding_dim': final_embeddings.shape[1],
        'genes_used': len(gene_ids),
        'vocab_match_rate': vocab_match_rate,
        'processing_time_seconds': processing_time,
        'cells_per_second': n_cells / processing_time,
    }
    
    # Save results directly to retnla directory
    print("Saving results...")
    result_adata.write_h5ad(f"{RETNLA_DIR}/hek293t_scgpt_embeddings.h5ad")
    np.save(f"{RETNLA_DIR}/hek293t_scgpt_embeddings.npy", final_embeddings)
    
    print("Results saved:")
    print(f"  - {RETNLA_DIR}/hek293t_scgpt_embeddings.h5ad")
    print(f"  - {RETNLA_DIR}/hek293t_scgpt_embeddings.npy")
    
    # Save metadata separately
    adata.obs.to_csv(f"{RETNLA_DIR}/hek293t_scgpt_metadata.csv")
    print(f"  - {RETNLA_DIR}/hek293t_scgpt_metadata.csv")
    
    # Cleanup
    del model, adata, embeddings, result_adata
    gc.collect()
    
    return True

# Execute processing
if __name__ == "__main__":
    try:
        success = process_hek293t()
        
        if success:
            print("SUCCESS: HEK293T scGPT processing completed!")
            
        else:
            print("FAILED: HEK293T scGPT processing failed")
            sys.exit(1)
            
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)