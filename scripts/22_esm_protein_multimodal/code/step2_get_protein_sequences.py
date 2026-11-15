#!/usr/bin/env python3
"""
Step 2: Retrieve protein sequences from UniProt for DEGs.

Queries UniProt REST API for reviewed human protein entries.
Outputs FASTA file for ESM-2 embedding generation.

Usage:
    python step2_get_protein_sequences.py [--input results/deg_gene_list.txt]
"""

import os
import sys
import argparse
import requests
import time
import pandas as pd
from pathlib import Path
import logging
import json

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"


def query_uniprot(gene_name, organism_id=9606, max_retries=3):
    """Query UniProt for a human protein by gene name."""
    query = f"gene:{gene_name} AND organism_id:{organism_id} AND reviewed:true"
    params = {
        "query": query,
        "format": "json",
        "size": 1,
        "fields": "accession,gene_names,protein_name,sequence,length"
    }

    for attempt in range(max_retries):
        try:
            resp = requests.get(UNIPROT_API, params=params, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                if data.get('results'):
                    result = data['results'][0]
                    seq = result.get('sequence', {})
                    return {
                        'gene': gene_name,
                        'uniprot_id': result.get('primaryAccession', ''),
                        'protein_name': result.get('proteinDescription', {}).get(
                            'recommendedName', {}).get('fullName', {}).get('value', ''),
                        'sequence': seq.get('value', ''),
                        'length': seq.get('length', 0)
                    }
                return None
            elif resp.status_code == 429:
                time.sleep(2 ** attempt)
            else:
                logger.warning(f"  {gene_name}: HTTP {resp.status_code}")
                return None
        except requests.exceptions.RequestException as e:
            logger.warning(f"  {gene_name}: Request failed ({e}), retry {attempt+1}")
            time.sleep(1)
    return None


def write_fasta(proteins, output_path):
    """Write protein sequences to FASTA format."""
    with open(output_path, 'w') as f:
        for p in proteins:
            header = f">{p['uniprot_id']}|{p['gene']}|{p['protein_name']}"
            f.write(header + '\n')
            seq = p['sequence']
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')


def main():
    parser = argparse.ArgumentParser(description="Retrieve UniProt sequences for DEGs")
    parser.add_argument('--input', type=str, default=None,
                        help='Gene list file (one gene per line)')
    parser.add_argument('--output-dir', type=str, default=None)
    args = parser.parse_args()

    script_dir = Path(__file__).resolve().parent
    results_dir = script_dir.parent / "results" if args.output_dir is None else Path(args.output_dir)
    results_dir.mkdir(parents=True, exist_ok=True)

    # Load gene list
    if args.input is None:
        args.input = results_dir / "deg_gene_list.txt"
    genes = [g.strip() for g in open(args.input) if g.strip()]
    logger.info(f"Querying UniProt for {len(genes)} genes...")

    # Query UniProt
    proteins = []
    failed = []
    for i, gene in enumerate(genes):
        if (i + 1) % 50 == 0:
            logger.info(f"  Progress: {i+1}/{len(genes)}")
        result = query_uniprot(gene)
        if result and result['sequence']:
            proteins.append(result)
        else:
            failed.append(gene)
        time.sleep(0.2)  # Rate limiting

    logger.info(f"\nResults: {len(proteins)} found, {len(failed)} not found")

    # Save FASTA
    fasta_path = results_dir / "deg_proteins.fasta"
    write_fasta(proteins, fasta_path)
    logger.info(f"Saved FASTA: {fasta_path}")

    # Save metadata
    meta_df = pd.DataFrame(proteins)
    meta_path = results_dir / "protein_metadata.csv"
    meta_df.to_csv(meta_path, index=False)
    logger.info(f"Saved metadata: {meta_path} ({len(meta_df)} proteins)")

    # Save gene-to-uniprot mapping (needed for step 4)
    mapping = {p['gene']: p['uniprot_id'] for p in proteins}
    mapping_path = results_dir / "gene_to_uniprot.json"
    with open(mapping_path, 'w') as f:
        json.dump(mapping, f, indent=2)
    logger.info(f"Saved mapping: {mapping_path}")

    # Save failed genes
    if failed:
        failed_path = results_dir / "failed_genes.txt"
        with open(failed_path, 'w') as f:
            f.write('\n'.join(failed))
        logger.info(f"Failed genes saved to {failed_path}")

    logger.info(f"\nTop 10 proteins by sequence length:")
    for p in sorted(proteins, key=lambda x: x['length'], reverse=True)[:10]:
        logger.info(f"  {p['gene']:15s}  {p['uniprot_id']:10s}  {p['length']:5d} aa")


if __name__ == '__main__':
    main()
