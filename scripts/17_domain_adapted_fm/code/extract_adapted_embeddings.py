#!/usr/bin/env python3
"""
ACC 17: Extract Embeddings from Domain-Adapted Models

Extract cell embeddings from:
- Approach A: Geneformer adapted via continued pre-training -> X_Geneformer_adapted (256 dims)
- Approach B: PlateletTransformer trained from scratch -> X_platelet_tf (384 dims)

Saves embeddings into AnnData .obsm for downstream benchmarking.

Usage:
    python extract_adapted_embeddings.py --model geneformer_adapted
    python extract_adapted_embeddings.py --model platelet_tf
    python extract_adapted_embeddings.py --model both

Author: Xinru Qiu
Date: February 2026
"""

import argparse
import logging
import pickle
from pathlib import Path
import json
import numpy as np
import torch
import scanpy as sc
from scipy.sparse import issparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
EMBEDDING_DIR = BASE_DIR / "02_EMBEDDINGS"
CHECKPOINT_DIR = BASE_DIR / "05_RESULTS/checkpoints"
RESULTS_DIR = BASE_DIR / "05_RESULTS/17_domain_adapted_fm"


def extract_scgpt_adapted_embeddings():
    """Extract embeddings from domain-adapted scGPT."""
    logger.info("Extracting scGPT adapted embeddings...")

    try:
        from scgpt.model import TransformerModel
        from scgpt.tokenizer.gene_tokenizer import GeneVocab
    except ImportError:
        logger.error("scGPT not installed. Activate scpeft_env.")
        raise

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load vocabulary
    scgpt_model_dir = BASE_DIR / "03_MODELS/scPEFT/save/scGPT_human"
    if not scgpt_model_dir.exists():
        scgpt_model_dir = BASE_DIR / "03_MODELS/scGPT/scGPT_human"
    vocab = GeneVocab.from_file(scgpt_model_dir / "vocab.json")

    # Load adapted checkpoint
    ckpt_path = CHECKPOINT_DIR / "17_scgpt_adapted/scgpt_adapted_best.pt"
    if not ckpt_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")

    model = TransformerModel(
        ntoken=len(vocab), d_model=512, nhead=8, d_hid=512, nlayers=12,
        vocab=vocab, dropout=0.0, pad_token="<pad>", pad_value=-2,
        n_input_bins=51, explicit_zero_prob=False, use_fast_transformer=True,
    )

    ckpt = torch.load(ckpt_path, map_location=device)
    model.load_state_dict(ckpt['model_state_dict'], strict=False)
    model = model.to(device).eval()
    logger.info("Loaded adapted scGPT model")

    # Load COVID data
    covid_path = EMBEDDING_DIR / "state/covid_state.h5ad"
    adata = sc.read_h5ad(covid_path)
    gene_names = list(adata.var_names)

    # Extract embeddings batch by batch
    embeddings = []
    batch_size = 32
    max_seq_len = 1200

    for start in range(0, adata.n_obs, batch_size):
        end = min(start + batch_size, adata.n_obs)
        batch_gene_ids = []
        batch_values = []
        batch_lengths = []

        for i in range(start, end):
            if issparse(adata.X):
                expr = np.array(adata.X[i].todense()).flatten()
            else:
                expr = adata.X[i].flatten()

            nonzero_idx = np.where(expr > 0)[0]
            if len(nonzero_idx) == 0:
                nonzero_idx = np.arange(min(10, len(expr)))
            sorted_idx = nonzero_idx[np.argsort(-expr[nonzero_idx])][:max_seq_len]

            gene_ids = []
            values = []
            for j in sorted_idx:
                gene = gene_names[j].upper()
                if gene in vocab:
                    gene_ids.append(vocab[gene])
                    values.append(float(expr[j]))

            n = len(gene_ids)
            batch_lengths.append(n)
            pad_len = max_seq_len - n
            gene_ids.extend([0] * pad_len)
            values.extend([0.0] * pad_len)
            batch_gene_ids.append(gene_ids[:max_seq_len])
            batch_values.append(values[:max_seq_len])

        gene_ids_t = torch.tensor(batch_gene_ids, dtype=torch.long, device=device)
        values_t = torch.tensor(batch_values, dtype=torch.float32, device=device)
        padding_mask = gene_ids_t == 0

        with torch.no_grad():
            output = model(
                gene_ids_t, values_t,
                src_key_padding_mask=padding_mask,
                CLS=True, CCE=False, MVC=False, ECS=False,
            )
            # Get cell embeddings from CLS output or cell_emb
            if 'cell_emb' in output:
                emb = output['cell_emb']
            elif 'cls_output' in output:
                emb = output['cls_output']
            else:
                # Use mean of last hidden state
                emb = output.get('transformer_output', values_t)
                if emb.dim() == 3:
                    emb = emb.mean(dim=1)

            embeddings.append(emb.cpu().numpy())

        if (start // batch_size) % 50 == 0:
            logger.info(f"  Processed {end}/{adata.n_obs} cells")

    embeddings = np.concatenate(embeddings, axis=0)
    logger.info(f"scGPT adapted embeddings: {embeddings.shape}")

    # Save
    adata.obsm['X_scGPT_adapted'] = embeddings
    output_path = EMBEDDING_DIR / "scgpt/covid_scgpt_adapted.h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    logger.info(f"Saved to {output_path}")

    return embeddings


def extract_geneformer_adapted_embeddings():
    """Extract embeddings from domain-adapted Geneformer (continued pre-training)."""
    logger.info("Extracting Geneformer adapted embeddings...")

    from transformers import BertModel

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load adapted model (BertModel for hidden states, not BertForMaskedLM)
    adapted_dir = CHECKPOINT_DIR / "17_geneformer_adapted/best_model"
    if not adapted_dir.exists():
        # Try final model in results dir
        adapted_dir = RESULTS_DIR / "geneformer_adapted/model"
    if not adapted_dir.exists():
        raise FileNotFoundError(
            f"Adapted Geneformer not found. Expected at "
            f"{CHECKPOINT_DIR / '17_geneformer_adapted/best_model'} or "
            f"{RESULTS_DIR / 'geneformer_adapted/model'}"
        )

    model = BertModel.from_pretrained(adapted_dir)
    model = model.to(device).eval()
    logger.info(f"Loaded adapted Geneformer from {adapted_dir}")

    # Load token dictionary
    token_dict = None
    token_dict_candidates = [
        BASE_DIR / "03_MODELS/geneformer/geneformer-V2-316M-native/geneformer/token_dictionary_gc104M.pkl",
        BASE_DIR / "03_MODELS/geneformer/Geneformer/geneformer/token_dictionary_gc104M.pkl",
        BASE_DIR / "03_MODELS/scPEFT/geneformer_peft/geneformer/token_dictionary.pkl",
    ]
    for tdf in token_dict_candidates:
        if tdf.exists():
            with open(tdf, 'rb') as f:
                token_dict = pickle.load(f)
            logger.info(f"Loaded token dictionary: {len(token_dict)} genes")
            break
    if token_dict is None:
        raise FileNotFoundError(f"Token dictionary not found: {token_dict_candidates}")

    # Load COVID data (Geneformer-ready preferred)
    geneformer_ready = BASE_DIR / "01_DATA/geneformer_ready/covid_geneformer_ready.h5ad"
    if geneformer_ready.exists():
        adata = sc.read_h5ad(geneformer_ready)
    else:
        covid_path = EMBEDDING_DIR / "state/covid_state.h5ad"
        adata = sc.read_h5ad(covid_path)

    gene_names = list(adata.var_names)

    # Pre-compute gene token IDs
    gene_token_ids = []
    for name in gene_names:
        tid = token_dict.get(name)
        if tid is None:
            tid = token_dict.get(name.upper())
        gene_token_ids.append(tid)

    CLS, SEP, PAD = 2, 3, 0
    max_len = 2048

    # Extract embeddings batch by batch
    embeddings = []
    batch_size = 32

    for start in range(0, adata.n_obs, batch_size):
        end = min(start + batch_size, adata.n_obs)
        batch_input_ids = []
        batch_attention_masks = []

        for i in range(start, end):
            if issparse(adata.X):
                expr = np.array(adata.X[i].todense()).flatten()
            else:
                expr = adata.X[i].flatten()

            nonzero_idx = np.where(expr > 0)[0]
            if len(nonzero_idx) == 0:
                nonzero_idx = np.arange(min(10, len(expr)))
            sorted_idx = nonzero_idx[np.argsort(-expr[nonzero_idx])]

            # Tokenize
            tokens = [CLS]
            for j in sorted_idx:
                if len(tokens) >= max_len - 1:
                    break
                tid = gene_token_ids[j]
                if tid is not None:
                    tokens.append(tid)
            tokens.append(SEP)

            attn = [1] * len(tokens)
            pad_len = max_len - len(tokens)
            tokens.extend([PAD] * pad_len)
            attn.extend([0] * pad_len)

            batch_input_ids.append(tokens[:max_len])
            batch_attention_masks.append(attn[:max_len])

        input_ids_t = torch.tensor(batch_input_ids, dtype=torch.long, device=device)
        attention_mask_t = torch.tensor(batch_attention_masks, dtype=torch.long, device=device)

        with torch.no_grad():
            outputs = model(input_ids=input_ids_t, attention_mask=attention_mask_t)
            hidden_states = outputs.last_hidden_state  # (batch, seq, 256)

            # Mean pooling over non-padded tokens
            mask_expanded = attention_mask_t.unsqueeze(-1).float()  # (batch, seq, 1)
            sum_hidden = (hidden_states * mask_expanded).sum(dim=1)  # (batch, 256)
            counts = mask_expanded.sum(dim=1).clamp(min=1)  # (batch, 1)
            emb = sum_hidden / counts  # (batch, 256)

            embeddings.append(emb.cpu().numpy())

        if (start // batch_size) % 50 == 0:
            logger.info(f"  Processed {end}/{adata.n_obs} cells")

    embeddings = np.concatenate(embeddings, axis=0)
    logger.info(f"Geneformer adapted embeddings: {embeddings.shape}")

    # Save
    adata.obsm['X_Geneformer_adapted'] = embeddings
    output_path = EMBEDDING_DIR / "geneformer/covid_geneformer_adapted.h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    logger.info(f"Saved to {output_path}")

    return embeddings


def extract_platelet_tf_embeddings():
    """Extract embeddings from PlateletTransformer."""
    logger.info("Extracting PlateletTransformer embeddings...")

    # Import model
    sys.path.insert(0, str(Path(__file__).parent))
    from platelet_transformer import PlateletTransformer

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load checkpoint
    ckpt_path = CHECKPOINT_DIR / "17_platelet_transformer/platelet_tf_best.pt"
    if not ckpt_path.exists():
        raise FileNotFoundError(f"Checkpoint not found: {ckpt_path}")

    ckpt = torch.load(ckpt_path, map_location=device)
    model_config = ckpt['config']
    gene_vocab = ckpt['gene_vocab']

    model = PlateletTransformer(
        n_genes=model_config['n_genes'],
        d_model=model_config['d_model'],
        n_heads=model_config['n_heads'],
        n_layers=model_config['n_layers'],
        d_ff=model_config['d_ff'],
        dropout=0.0,
        max_seq_len=model_config['max_seq_len'],
    ).to(device)

    model.load_state_dict(ckpt['model_state_dict'])
    model.eval()
    logger.info(f"Loaded PlateletTransformer: {model_config['d_model']}-dim")

    # Load COVID data
    covid_path = EMBEDDING_DIR / "state/covid_state.h5ad"
    adata = sc.read_h5ad(covid_path)

    if issparse(adata.X):
        X = adata.X.toarray().astype(np.float32)
    else:
        X = np.array(adata.X, dtype=np.float32)

    gene_names = [g.upper() for g in adata.var_names]
    gene_to_vocab = {i: gene_vocab[gene] for i, gene in enumerate(gene_names) if gene in gene_vocab}

    # Extract embeddings
    embeddings = []
    batch_size = 64
    max_seq_len = model_config['max_seq_len']

    for start in range(0, len(X), batch_size):
        end = min(start + batch_size, len(X))
        batch_gene_ids = []
        batch_values = []
        batch_lengths = []

        for i in range(start, end):
            expr = X[i]
            nonzero = np.where(expr > 0)[0]
            if len(nonzero) == 0:
                nonzero = np.arange(min(10, len(expr)))
            sorted_idx = nonzero[np.argsort(-expr[nonzero])]

            gene_ids = []
            values = []
            for j in sorted_idx:
                if j in gene_to_vocab:
                    gene_ids.append(gene_to_vocab[j])
                    values.append(float(expr[j]))
                    if len(gene_ids) >= max_seq_len:
                        break

            n = len(gene_ids)
            if n == 0:
                gene_ids = [1]
                values = [0.0]
                n = 1

            batch_lengths.append(n)
            pad_len = max_seq_len - n
            gene_ids.extend([0] * pad_len)
            values.extend([0.0] * pad_len)
            batch_gene_ids.append(gene_ids[:max_seq_len])
            batch_values.append(values[:max_seq_len])

        gene_ids_t = torch.tensor(batch_gene_ids, dtype=torch.long, device=device)
        values_t = torch.tensor(batch_values, dtype=torch.float32, device=device)
        lengths_t = torch.tensor(batch_lengths, dtype=torch.long, device=device)

        with torch.no_grad():
            emb = model.get_cell_embeddings(gene_ids_t, values_t, lengths_t)
            embeddings.append(emb.cpu().numpy())

        if (start // batch_size) % 50 == 0:
            logger.info(f"  Processed {end}/{len(X)} cells")

    embeddings = np.concatenate(embeddings, axis=0)
    logger.info(f"PlateletTransformer embeddings: {embeddings.shape}")

    # Save
    adata.obsm['X_platelet_tf'] = embeddings
    output_path = EMBEDDING_DIR / "platelet_tf/covid_platelet_tf.h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    logger.info(f"Saved to {output_path}")

    return embeddings


def main():
    import sys as _sys
    global sys
    sys = _sys

    parser = argparse.ArgumentParser(description='Extract Adapted Model Embeddings')
    parser.add_argument('--model', default='both',
                        choices=['scgpt_adapted', 'geneformer_adapted', 'platelet_tf', 'both'])
    args = parser.parse_args()

    if args.model == 'scgpt_adapted':
        extract_scgpt_adapted_embeddings()

    if args.model in ('geneformer_adapted', 'both'):
        extract_geneformer_adapted_embeddings()

    if args.model in ('platelet_tf', 'both'):
        extract_platelet_tf_embeddings()

    logger.info("Embedding extraction complete.")


if __name__ == '__main__':
    main()
