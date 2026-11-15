#!/usr/bin/env python3
"""
Generate ESM-2 protein embeddings for therapeutic targets.
ESM-2 embeddings can be used for:
- Protein function prediction
- Drug-target interaction prediction
- Binding site identification

Usage:
    python generate_esm_embeddings.py \
        --fasta protein_sequences.fasta \
        --output esm_embeddings.npz

Requirements:
    pip install torch fair-esm pandas numpy
"""

import argparse
import os
import sys
from pathlib import Path
from datetime import datetime

import numpy as np
import pandas as pd
import torch


def load_esm_model(model_name: str = 'esm2_t33_650M_UR50D'):
    """
    Load ESM-2 model.

    Available models:
    - esm2_t6_8M_UR50D: 8M params, fastest
    - esm2_t12_35M_UR50D: 35M params
    - esm2_t30_150M_UR50D: 150M params
    - esm2_t33_650M_UR50D: 650M params, recommended
    - esm2_t36_3B_UR50D: 3B params, most accurate but slow
    """
    try:
        import esm
    except ImportError:
        print("ERROR: ESM not installed. Run: pip install fair-esm")
        sys.exit(1)

    print(f"Loading ESM model: {model_name}")
    model, alphabet = esm.pretrained.load_model_and_alphabet(model_name)
    model.eval()

    # Move to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    print(f"  Using device: {device}")

    return model, alphabet, device


def read_fasta(fasta_path: str) -> list:
    """Read sequences from FASTA file."""
    sequences = []
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_header:
                    sequences.append({
                        'header': current_header,
                        'sequence': ''.join(current_seq)
                    })
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)

        # Add last sequence
        if current_header:
            sequences.append({
                'header': current_header,
                'sequence': ''.join(current_seq)
            })

    return sequences


def generate_embeddings(
    sequences: list,
    model,
    alphabet,
    device,
    batch_size: int = 1
) -> dict:
    """
    Generate ESM-2 embeddings for protein sequences.

    Args:
        sequences: List of dicts with 'header' and 'sequence'
        model: ESM model
        alphabet: ESM alphabet
        device: torch device
        batch_size: Batch size for inference

    Returns:
        Dictionary mapping headers to embeddings
    """
    batch_converter = alphabet.get_batch_converter()
    embeddings = {}

    print(f"\nGenerating embeddings for {len(sequences)} proteins...")

    for i, seq_dict in enumerate(sequences):
        header = seq_dict['header']
        sequence = seq_dict['sequence']

        # Parse header to get gene name
        parts = header.split('|')
        uniprot_id = parts[0] if parts else header
        gene = parts[1] if len(parts) > 1 else uniprot_id

        print(f"  [{i+1}/{len(sequences)}] {gene} ({len(sequence)} aa)...", end=' ')

        # Truncate if too long (ESM-2 max length is 1024 for some models)
        if len(sequence) > 1022:
            print(f"truncating from {len(sequence)} to 1022...", end=' ')
            sequence = sequence[:1022]

        # Prepare batch
        batch_labels, batch_strs, batch_tokens = batch_converter([(header, sequence)])
        batch_tokens = batch_tokens.to(device)

        # Generate embeddings
        with torch.no_grad():
            results = model(batch_tokens, repr_layers=[33], return_contacts=False)

        # Extract embeddings
        # Shape: (batch, seq_len, embed_dim)
        token_embeddings = results['representations'][33]

        # Mean pooling over sequence (excluding BOS/EOS tokens)
        seq_len = len(sequence)
        mean_embedding = token_embeddings[0, 1:seq_len+1, :].mean(dim=0).cpu().numpy()

        # Also get per-residue embeddings for binding site analysis
        per_residue = token_embeddings[0, 1:seq_len+1, :].cpu().numpy()

        embeddings[gene] = {
            'uniprot_id': uniprot_id,
            'gene': gene,
            'sequence': sequence,
            'length': len(sequence),
            'mean_embedding': mean_embedding,
            'per_residue_embeddings': per_residue,
            'embedding_dim': mean_embedding.shape[0]
        }

        print(f"✓ ({mean_embedding.shape[0]}D)")

    return embeddings


def save_embeddings(embeddings: dict, output_dir: str, timestamp: str):
    """Save embeddings to files."""
    os.makedirs(output_dir, exist_ok=True)

    # Save mean embeddings as NPZ
    mean_embs = {gene: data['mean_embedding'] for gene, data in embeddings.items()}
    npz_path = os.path.join(output_dir, f'esm_embeddings_{timestamp}.npz')
    np.savez(npz_path, **mean_embs)
    print(f"\nSaved mean embeddings to {npz_path}")

    # Save per-residue embeddings (larger file)
    residue_embs = {gene: data['per_residue_embeddings'] for gene, data in embeddings.items()}
    residue_path = os.path.join(output_dir, f'esm_per_residue_{timestamp}.npz')
    np.savez_compressed(residue_path, **residue_embs)
    print(f"Saved per-residue embeddings to {residue_path}")

    # Save summary CSV
    rows = []
    for gene, data in embeddings.items():
        rows.append({
            'gene': gene,
            'uniprot_id': data['uniprot_id'],
            'length': data['length'],
            'embedding_dim': data['embedding_dim'],
        })
    summary_df = pd.DataFrame(rows)
    summary_path = os.path.join(output_dir, f'esm_summary_{timestamp}.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary to {summary_path}")

    return npz_path, residue_path, summary_path


def compute_similarity_matrix(embeddings: dict) -> pd.DataFrame:
    """Compute pairwise cosine similarity between protein embeddings."""
    from sklearn.metrics.pairwise import cosine_similarity

    genes = list(embeddings.keys())
    emb_matrix = np.array([embeddings[g]['mean_embedding'] for g in genes])

    sim_matrix = cosine_similarity(emb_matrix)
    sim_df = pd.DataFrame(sim_matrix, index=genes, columns=genes)

    return sim_df


def main():
    parser = argparse.ArgumentParser(description='Generate ESM-2 protein embeddings')
    parser.add_argument('--fasta', '-f', required=True,
                        help='Input FASTA file with protein sequences')
    parser.add_argument('--output-dir', '-o', default='results/esm_embeddings',
                        help='Output directory')
    parser.add_argument('--model', '-m', default='esm2_t33_650M_UR50D',
                        help='ESM model name')
    parser.add_argument('--batch-size', '-b', type=int, default=1,
                        help='Batch size for inference')

    args = parser.parse_args()

    print(f"\n{'='*60}")
    print("ESM-2 PROTEIN EMBEDDING GENERATION")
    print(f"{'='*60}")

    # Load model
    model, alphabet, device = load_esm_model(args.model)

    # Read sequences
    sequences = read_fasta(args.fasta)
    print(f"\nLoaded {len(sequences)} sequences from {args.fasta}")

    # Generate embeddings
    embeddings = generate_embeddings(sequences, model, alphabet, device, args.batch_size)

    # Save embeddings
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    save_embeddings(embeddings, args.output_dir, timestamp)

    # Compute similarity matrix
    print("\nComputing protein similarity matrix...")
    sim_df = compute_similarity_matrix(embeddings)
    sim_path = os.path.join(args.output_dir, f'protein_similarity_{timestamp}.csv')
    sim_df.to_csv(sim_path)
    print(f"Saved similarity matrix to {sim_path}")

    # Print similarity highlights
    print(f"\n{'='*60}")
    print("PROTEIN SIMILARITY HIGHLIGHTS")
    print(f"{'='*60}")

    genes = list(embeddings.keys())
    for i, g1 in enumerate(genes):
        for j, g2 in enumerate(genes):
            if i < j:
                sim = sim_df.loc[g1, g2]
                if sim > 0.7:  # High similarity
                    print(f"  {g1} <-> {g2}: {sim:.3f} (high similarity)")

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Proteins embedded: {len(embeddings)}")
    print(f"Embedding dimension: {list(embeddings.values())[0]['embedding_dim']}")
    print(f"Output directory: {args.output_dir}")

    print("\nDone!")


if __name__ == '__main__':
    main()
