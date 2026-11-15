#!/usr/bin/env python3
"""
Download AlphaFold predicted structures for therapeutic targets.
Structures are used for molecular docking and binding site analysis.

Usage:
    python download_alphafold_structures.py \
        --summary target_summary.csv \
        --output structures/

Requirements:
    pip install requests pandas biopython
"""

import argparse
import os
import requests
import json
from pathlib import Path
from datetime import datetime

import pandas as pd


ALPHAFOLD_BASE_URL = "https://alphafold.ebi.ac.uk/files"


def download_alphafold_structure(uniprot_id: str, output_dir: str, version: int = 4) -> dict:
    """
    Download AlphaFold structure files for a UniProt ID.

    Downloads:
    - PDB file (structure)
    - PAE JSON (predicted aligned error)
    - CIF file (mmCIF format)

    Args:
        uniprot_id: UniProt accession
        output_dir: Output directory
        version: AlphaFold database version (default: 4)

    Returns:
        Dictionary with download results
    """
    result = {
        'uniprot_id': uniprot_id,
        'pdb_path': None,
        'pae_path': None,
        'cif_path': None,
        'success': False,
        'error': None
    }

    os.makedirs(output_dir, exist_ok=True)

    # File URLs
    base_name = f"AF-{uniprot_id}-F1"
    files = {
        'pdb': f"{base_name}-model_v{version}.pdb",
        'pae': f"{base_name}-predicted_aligned_error_v{version}.json",
        'cif': f"{base_name}-model_v{version}.cif",
    }

    for file_type, filename in files.items():
        url = f"{ALPHAFOLD_BASE_URL}/{filename}"
        output_path = os.path.join(output_dir, filename)

        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()

            with open(output_path, 'wb') as f:
                f.write(response.content)

            result[f'{file_type}_path'] = output_path
            print(f"    ✓ {file_type.upper()}: {filename}")

        except requests.exceptions.HTTPError as e:
            if e.response.status_code == 404:
                print(f"    ✗ {file_type.upper()}: Not found")
            else:
                print(f"    ✗ {file_type.upper()}: HTTP {e.response.status_code}")
            result['error'] = str(e)

        except requests.exceptions.RequestException as e:
            print(f"    ✗ {file_type.upper()}: {e}")
            result['error'] = str(e)

    # Mark success if at least PDB was downloaded
    if result['pdb_path']:
        result['success'] = True

    return result


def extract_plddt_scores(pdb_path: str) -> dict:
    """
    Extract pLDDT confidence scores from AlphaFold PDB file.
    pLDDT is stored in the B-factor column.

    Returns:
        Dictionary with pLDDT statistics
    """
    scores = []

    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                # B-factor is columns 61-66
                try:
                    plddt = float(line[60:66].strip())
                    scores.append(plddt)
                except ValueError:
                    pass

    if scores:
        import numpy as np
        return {
            'mean_plddt': np.mean(scores),
            'min_plddt': np.min(scores),
            'max_plddt': np.max(scores),
            'high_confidence_residues': sum(1 for s in scores if s > 70),
            'very_high_confidence_residues': sum(1 for s in scores if s > 90),
            'total_residues': len(scores)
        }
    return None


def analyze_structure_quality(structures: list, output_dir: str, timestamp: str):
    """Analyze structure quality and create summary."""
    quality_data = []

    print("\nAnalyzing structure quality...")

    for struct in structures:
        if struct['success'] and struct['pdb_path']:
            plddt = extract_plddt_scores(struct['pdb_path'])
            if plddt:
                quality_data.append({
                    'uniprot_id': struct['uniprot_id'],
                    'gene': struct.get('gene', ''),
                    'mean_plddt': plddt['mean_plddt'],
                    'min_plddt': plddt['min_plddt'],
                    'max_plddt': plddt['max_plddt'],
                    'total_residues': plddt['total_residues'],
                    'high_confidence_pct': plddt['high_confidence_residues'] / plddt['total_residues'] * 100,
                    'very_high_confidence_pct': plddt['very_high_confidence_residues'] / plddt['total_residues'] * 100,
                    'pdb_path': struct['pdb_path'],
                })

    if quality_data:
        df = pd.DataFrame(quality_data)
        quality_path = os.path.join(output_dir, f'structure_quality_{timestamp}.csv')
        df.to_csv(quality_path, index=False)
        print(f"Saved quality analysis to {quality_path}")

        # Print summary
        print(f"\n{'='*60}")
        print("STRUCTURE QUALITY SUMMARY")
        print(f"{'='*60}")
        print(f"{'Gene':<12} {'pLDDT':<8} {'High Conf %':<12} {'Residues':<10}")
        print("-" * 50)
        for _, row in df.iterrows():
            print(f"{row['gene']:<12} {row['mean_plddt']:<8.1f} {row['high_confidence_pct']:<12.1f} {row['total_residues']:<10}")

        return df

    return None


def identify_binding_sites(pdb_path: str) -> list:
    """
    Identify potential binding sites based on pLDDT and structure.
    High-confidence regions are more reliable for docking.

    Note: For accurate binding site prediction, use tools like:
    - fpocket
    - SiteMap
    - DoGSiteScorer
    """
    # This is a simplified analysis
    # Real binding site prediction requires specialized tools

    binding_sites = []

    # Read pLDDT scores with residue info
    residue_scores = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line:
                try:
                    res_num = int(line[22:26].strip())
                    res_name = line[17:20].strip()
                    plddt = float(line[60:66].strip())
                    residue_scores.append({
                        'residue_num': res_num,
                        'residue_name': res_name,
                        'plddt': plddt
                    })
                except ValueError:
                    pass

    # Find high-confidence regions (potential binding sites)
    # Look for stretches of residues with pLDDT > 70
    in_region = False
    region_start = None

    for res in residue_scores:
        if res['plddt'] > 70:
            if not in_region:
                in_region = True
                region_start = res['residue_num']
        else:
            if in_region and region_start:
                region_end = res['residue_num'] - 1
                if region_end - region_start >= 5:  # At least 5 residues
                    binding_sites.append({
                        'start': region_start,
                        'end': region_end,
                        'length': region_end - region_start + 1
                    })
            in_region = False
            region_start = None

    return binding_sites


def main():
    parser = argparse.ArgumentParser(description='Download AlphaFold structures')
    parser.add_argument('--summary', '-s', required=True,
                        help='Path to target summary CSV (with UniProt IDs)')
    parser.add_argument('--output', '-o', default='results/structures',
                        help='Output directory')
    parser.add_argument('--version', '-v', type=int, default=4,
                        help='AlphaFold database version')

    args = parser.parse_args()

    print(f"\n{'='*60}")
    print("ALPHAFOLD STRUCTURE DOWNLOAD")
    print(f"{'='*60}")

    # Load target summary
    df = pd.read_csv(args.summary)
    print(f"Loaded {len(df)} targets from {args.summary}")

    # Create output directory
    os.makedirs(args.output, exist_ok=True)

    # Download structures
    structures = []
    print("\nDownloading structures...")

    for _, row in df.iterrows():
        uniprot_id = row['UniProt_ID']
        gene = row.get('Gene', uniprot_id)

        print(f"\n  {gene} ({uniprot_id}):")

        result = download_alphafold_structure(uniprot_id, args.output, args.version)
        result['gene'] = gene
        structures.append(result)

    # Analyze quality
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    quality_df = analyze_structure_quality(structures, args.output, timestamp)

    # Summary
    successful = sum(1 for s in structures if s['success'])
    print(f"\n{'='*60}")
    print("DOWNLOAD SUMMARY")
    print(f"{'='*60}")
    print(f"Structures downloaded: {successful}/{len(structures)}")
    print(f"Output directory: {args.output}")

    # List files
    print(f"\nDownloaded files:")
    for s in structures:
        if s['success']:
            print(f"  {s['gene']}: {Path(s['pdb_path']).name}")

    print("\nDone!")


if __name__ == '__main__':
    main()
