#!/usr/bin/env python3
"""
Retrieve protein sequences for therapeutic target genes from UniProt.
Prepares data for ESM-2 embedding and AlphaFold structure analysis.

Usage:
    python get_protein_sequences.py \
        --genes ICMT,DUSP11,LILRB2,PIK3R2,B4GALT1,PINK1,RAMP3 \
        --output protein_sequences.fasta

Requirements:
    pip install requests biopython pandas
"""

import argparse
import requests
import json
import time
from pathlib import Path
from datetime import datetime

import pandas as pd


# Top therapeutic candidates from perturbation analysis
DEFAULT_TARGETS = [
    'ICMT',      # Isoprenylcysteine carboxyl methyltransferase
    'DUSP11',    # Dual specificity phosphatase 11
    'LILRB2',    # Leukocyte immunoglobulin-like receptor B2
    'B4GALT1',   # Beta-1,4-galactosyltransferase 1
    'PINK1',     # PTEN-induced kinase 1
    'PIK3R2',    # PI3K regulatory subunit beta
    'RAMP3',     # Receptor activity modifying protein 3
    'SLC28A1',   # Nucleoside transporter
    'ATP2B3',    # Plasma membrane Ca2+ ATPase 3
    'ESD',       # Esterase D
    'AGPAT3',    # 1-acylglycerol-3-phosphate O-acyltransferase 3
    'MED31',     # Mediator complex subunit 31
]


def query_uniprot(gene_name: str, organism: str = 'human') -> dict:
    """
    Query UniProt for protein information by gene name.

    Args:
        gene_name: Gene symbol (e.g., 'ICMT')
        organism: Organism name (default: human)

    Returns:
        Dictionary with protein info or None if not found
    """
    # UniProt REST API
    base_url = "https://rest.uniprot.org/uniprotkb/search"

    # Search for reviewed (Swiss-Prot) human proteins
    query = f'gene:{gene_name} AND organism_id:9606 AND reviewed:true'

    params = {
        'query': query,
        'format': 'json',
        'fields': 'accession,id,gene_names,protein_name,sequence,length,organism_name,function,subcellular_location',
        'size': 1
    }

    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()

        if data['results']:
            result = data['results'][0]
            return {
                'gene': gene_name,
                'uniprot_id': result.get('primaryAccession', ''),
                'entry_name': result.get('uniProtkbId', ''),
                'protein_name': result.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'sequence': result.get('sequence', {}).get('value', ''),
                'length': result.get('sequence', {}).get('length', 0),
                'organism': result.get('organism', {}).get('scientificName', ''),
                'function': extract_function(result),
                'subcellular_location': extract_location(result),
            }
        else:
            print(f"  WARNING: No UniProt entry found for {gene_name}")
            return None

    except requests.exceptions.RequestException as e:
        print(f"  ERROR querying UniProt for {gene_name}: {e}")
        return None


def extract_function(result: dict) -> str:
    """Extract function annotation from UniProt result."""
    comments = result.get('comments', [])
    for comment in comments:
        if comment.get('commentType') == 'FUNCTION':
            texts = comment.get('texts', [])
            if texts:
                return texts[0].get('value', '')
    return ''


def extract_location(result: dict) -> str:
    """Extract subcellular location from UniProt result."""
    comments = result.get('comments', [])
    for comment in comments:
        if comment.get('commentType') == 'SUBCELLULAR LOCATION':
            locations = comment.get('subcellularLocations', [])
            if locations:
                loc = locations[0].get('location', {}).get('value', '')
                return loc
    return ''


def get_alphafold_url(uniprot_id: str) -> str:
    """Get AlphaFold structure URL for a UniProt ID."""
    return f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"


def get_alphafold_pdb_url(uniprot_id: str) -> str:
    """Get AlphaFold PDB file URL."""
    return f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"


def write_fasta(proteins: list, output_path: str):
    """Write protein sequences to FASTA file."""
    with open(output_path, 'w') as f:
        for p in proteins:
            if p and p.get('sequence'):
                header = f">{p['uniprot_id']}|{p['gene']}|{p['protein_name']}"
                f.write(f"{header}\n")
                # Write sequence in 80-character lines
                seq = p['sequence']
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")

    print(f"Wrote {len([p for p in proteins if p])} sequences to {output_path}")


def create_summary_table(proteins: list, output_path: str):
    """Create summary table of protein targets."""
    rows = []
    for p in proteins:
        if p:
            rows.append({
                'Gene': p['gene'],
                'UniProt_ID': p['uniprot_id'],
                'Protein_Name': p['protein_name'],
                'Length_AA': p['length'],
                'Function': p['function'][:200] + '...' if len(p.get('function', '')) > 200 else p.get('function', ''),
                'Location': p['subcellular_location'],
                'AlphaFold_URL': get_alphafold_url(p['uniprot_id']),
                'PDB_URL': get_alphafold_pdb_url(p['uniprot_id']),
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"Wrote summary to {output_path}")

    return df


def main():
    parser = argparse.ArgumentParser(description='Retrieve protein sequences for drug targets')
    parser.add_argument('--genes', '-g', type=str,
                        help='Comma-separated list of gene symbols')
    parser.add_argument('--candidates-csv', '-c',
                        help='Path to therapeutic candidates CSV')
    parser.add_argument('--output', '-o', default='protein_sequences.fasta',
                        help='Output FASTA file')
    parser.add_argument('--output-dir', '-d', default='results/protein_data',
                        help='Output directory')
    parser.add_argument('--top-n', type=int, default=12,
                        help='Number of top candidates to retrieve')

    args = parser.parse_args()

    # Get gene list
    if args.genes:
        genes = [g.strip() for g in args.genes.split(',')]
    elif args.candidates_csv:
        df = pd.read_csv(args.candidates_csv)
        df = df[~df['perturbation'].str.contains('non-targeting', case=False, na=False)]
        df['gene'] = df['perturbation'].apply(lambda x: x.split('_')[0].split('|')[0])
        genes = df['gene'].head(args.top_n).tolist()
    else:
        genes = DEFAULT_TARGETS[:args.top_n]

    print(f"\n{'='*60}")
    print("PROTEIN SEQUENCE RETRIEVAL")
    print(f"{'='*60}")
    print(f"Genes: {', '.join(genes)}")

    # Create output directory
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Query UniProt for each gene
    proteins = []
    print("\nQuerying UniProt...")

    for gene in genes:
        print(f"  {gene}...", end=' ')
        protein = query_uniprot(gene)
        if protein:
            print(f"✓ {protein['uniprot_id']} ({protein['length']} aa)")
            proteins.append(protein)
        else:
            print("✗ Not found")
        time.sleep(0.5)  # Rate limiting

    # Write outputs
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # FASTA file
    fasta_path = Path(args.output_dir) / f'target_sequences_{timestamp}.fasta'
    write_fasta(proteins, str(fasta_path))

    # Summary table
    summary_path = Path(args.output_dir) / f'target_summary_{timestamp}.csv'
    summary_df = create_summary_table(proteins, str(summary_path))

    # Print summary
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Proteins retrieved: {len([p for p in proteins if p])}/{len(genes)}")
    print(f"\nOutput files:")
    print(f"  FASTA: {fasta_path}")
    print(f"  Summary: {summary_path}")

    print(f"\n{'='*60}")
    print("TARGET INFORMATION")
    print(f"{'='*60}")
    for p in proteins:
        if p:
            print(f"\n{p['gene']} ({p['uniprot_id']})")
            print(f"  Name: {p['protein_name']}")
            print(f"  Length: {p['length']} aa")
            print(f"  Location: {p['subcellular_location']}")
            print(f"  AlphaFold: {get_alphafold_url(p['uniprot_id'])}")

    print("\nDone!")


if __name__ == '__main__':
    main()
