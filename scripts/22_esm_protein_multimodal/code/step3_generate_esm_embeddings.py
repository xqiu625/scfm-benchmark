#!/usr/bin/env python3
"""
Step 3: Generate ESM-2 protein embeddings for DEG proteins.

Uses ESM-2 (650M parameter model) to generate per-protein embeddings.
Mean-pools over residues to get a single 1280-dim vector per protein.

Usage:
    python step3_generate_esm_embeddings.py [--model esm2_t33_650M_UR50D]

Requirements:
    pip install fair-esm
"""

import os
import sys
import argparse
import numpy as np
import torch
from pathlib import Path
import logging
import json

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def read_fasta(fasta_path):
    """Parse FASTA file into list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def load_esm_model(model_name="esm2_t33_650M_UR50D"):
    """Load ESM-2 model and alphabet."""
    import esm

    logger.info(f"Loading ESM model: {model_name}")
    if model_name == "esm2_t33_650M_UR50D":
        model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    elif model_name == "esm2_t36_3B_UR50D":
        model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    elif model_name == "esm2_t30_150M_UR50D":
        model, alphabet = esm.pretrained.esm2_t30_150M_UR50D()
    elif model_name == "esm2_t12_35M_UR50D":
        model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    elif model_name == "esm2_t6_8M_UR50D":
        model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
    else:
        raise ValueError(f"Unknown model: {model_name}")

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    model.eval()
    logger.info(f"  Model loaded on {device}")

    # Determine representation layer
    layer_map = {
        "esm2_t33_650M_UR50D": 33,
        "esm2_t36_3B_UR50D": 36,
        "esm2_t30_150M_UR50D": 30,
        "esm2_t12_35M_UR50D": 12,
        "esm2_t6_8M_UR50D": 6,
    }
    repr_layer = layer_map.get(model_name, 33)

    return model, alphabet, device, repr_layer


def generate_embeddings(model, alphabet, device, repr_layer, sequences, batch_size=4, max_len=1022):
    """Generate mean-pooled ESM embeddings for protein sequences."""
    batch_converter = alphabet.get_batch_converter()
    embeddings = {}

    # Process in batches
    for i in range(0, len(sequences), batch_size):
        batch_seqs = sequences[i:i + batch_size]

        # Truncate long sequences
        batch_data = []
        for header, seq in batch_seqs:
            gene = header.split('|')[1] if '|' in header else header
            truncated = seq[:max_len]
            batch_data.append((gene, truncated))

        _, _, batch_tokens = batch_converter(batch_data)
        batch_tokens = batch_tokens.to(device)

        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[repr_layer], return_contacts=False)
        representations = results["representations"][repr_layer]

        # Mean pool over sequence (exclude BOS and EOS tokens)
        for j, (gene, seq) in enumerate(batch_data):
            seq_len = min(len(seq), max_len)
            # tokens: [BOS, aa1, aa2, ..., EOS, PAD, PAD, ...]
            token_repr = representations[j, 1:seq_len + 1]  # exclude BOS, take seq_len tokens
            mean_embedding = token_repr.mean(dim=0).cpu().numpy()
            embeddings[gene] = mean_embedding

        if (i + batch_size) % 50 < batch_size:
            logger.info(f"  Processed {min(i + batch_size, len(sequences))}/{len(sequences)} proteins")

    return embeddings


def main():
    parser = argparse.ArgumentParser(description="Generate ESM-2 protein embeddings")
    parser.add_argument('--model', type=str, default='esm2_t33_650M_UR50D',
                        help='ESM model name')
    parser.add_argument('--batch-size', type=int, default=4)
    parser.add_argument('--input', type=str, default=None,
                        help='FASTA file (default: results/deg_proteins.fasta)')
    parser.add_argument('--output-dir', type=str, default=None)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    results_dir = script_dir.parent / "results" if args.output_dir is None else Path(args.output_dir)

    # Load sequences
    fasta_path = Path(args.input) if args.input else results_dir / "deg_proteins.fasta"
    sequences = read_fasta(fasta_path)
    logger.info(f"Loaded {len(sequences)} protein sequences from {fasta_path}")

    # Load model
    model, alphabet, device, repr_layer = load_esm_model(args.model)

    # Generate embeddings
    logger.info(f"Generating embeddings (layer {repr_layer})...")
    embeddings = generate_embeddings(
        model, alphabet, device, repr_layer, sequences,
        batch_size=args.batch_size
    )

    # Save
    emb_path = results_dir / "esm_protein_embeddings.npz"
    np.savez_compressed(emb_path, **embeddings)
    logger.info(f"Saved {len(embeddings)} protein embeddings to {emb_path}")

    # Save summary
    emb_dim = next(iter(embeddings.values())).shape[0]
    summary = {
        'model': args.model,
        'repr_layer': repr_layer,
        'embedding_dim': emb_dim,
        'n_proteins': len(embeddings),
        'genes': list(embeddings.keys())
    }
    summary_path = results_dir / "esm_embedding_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    logger.info(f"\nSummary:")
    logger.info(f"  Model: {args.model}")
    logger.info(f"  Embedding dim: {emb_dim}")
    logger.info(f"  Proteins embedded: {len(embeddings)}")


if __name__ == '__main__':
    main()
