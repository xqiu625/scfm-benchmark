#!/usr/bin/env python3
"""
Search drug databases (DrugBank, ChEMBL, OpenTargets) for existing drugs
targeting therapeutic candidates identified from perturbation analysis.

Usage:
    python search_drug_databases.py \
        --genes ICMT,DUSP11,LILRB2,PIK3R2 \
        --output drug_candidates.csv

Requirements:
    pip install requests pandas chembl_webresource_client
"""

import argparse
import os
import json
import time
from pathlib import Path
from datetime import datetime

import requests
import pandas as pd


# Top therapeutic candidates
DEFAULT_TARGETS = [
    'ICMT', 'DUSP11', 'LILRB2', 'B4GALT1', 'PINK1',
    'PIK3R2', 'RAMP3', 'SLC28A1', 'ATP2B3', 'ESD'
]


def search_chembl_target(gene_name: str) -> dict:
    """
    Search ChEMBL for target information and associated compounds.

    Args:
        gene_name: Gene symbol

    Returns:
        Dictionary with target and compound information
    """
    base_url = "https://www.ebi.ac.uk/chembl/api/data"

    result = {
        'gene': gene_name,
        'chembl_target_id': None,
        'target_type': None,
        'compounds': [],
        'approved_drugs': [],
        'clinical_candidates': []
    }

    # Search for target
    try:
        search_url = f"{base_url}/target/search.json"
        params = {
            'q': gene_name,
            'limit': 10
        }
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        if data.get('targets'):
            # Find human target
            for target in data['targets']:
                if 'Homo sapiens' in str(target.get('organism', '')):
                    result['chembl_target_id'] = target.get('target_chembl_id')
                    result['target_type'] = target.get('target_type')
                    result['pref_name'] = target.get('pref_name')
                    break

        if not result['chembl_target_id']:
            return result

        # Get compounds for this target
        target_id = result['chembl_target_id']
        activity_url = f"{base_url}/activity.json"
        params = {
            'target_chembl_id': target_id,
            'limit': 100
        }
        response = requests.get(activity_url, params=params, timeout=30)
        response.raise_for_status()
        activity_data = response.json()

        compounds_seen = set()
        for activity in activity_data.get('activities', []):
            compound_id = activity.get('molecule_chembl_id')
            if compound_id and compound_id not in compounds_seen:
                compounds_seen.add(compound_id)
                result['compounds'].append({
                    'chembl_id': compound_id,
                    'name': activity.get('molecule_pref_name', ''),
                    'activity_type': activity.get('standard_type', ''),
                    'activity_value': activity.get('standard_value', ''),
                    'activity_units': activity.get('standard_units', '')
                })

        # Check for approved drugs
        for compound_id in list(compounds_seen)[:20]:
            molecule_url = f"{base_url}/molecule/{compound_id}.json"
            try:
                mol_response = requests.get(molecule_url, timeout=10)
                mol_response.raise_for_status()
                mol_data = mol_response.json()

                max_phase = mol_data.get('max_phase', 0)
                if max_phase == 4:  # Approved
                    result['approved_drugs'].append({
                        'chembl_id': compound_id,
                        'name': mol_data.get('pref_name', ''),
                        'phase': 'Approved'
                    })
                elif max_phase >= 2:  # Clinical trials
                    result['clinical_candidates'].append({
                        'chembl_id': compound_id,
                        'name': mol_data.get('pref_name', ''),
                        'phase': f'Phase {max_phase}'
                    })
            except Exception:
                pass

    except requests.exceptions.RequestException as e:
        result['error'] = str(e)

    return result


def search_opentargets(gene_name: str) -> dict:
    """
    Search OpenTargets for drug-target associations.

    Args:
        gene_name: Gene symbol

    Returns:
        Dictionary with drug associations
    """
    result = {
        'gene': gene_name,
        'ensembl_id': None,
        'drugs': [],
        'disease_associations': []
    }

    # OpenTargets GraphQL API
    api_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # First, get Ensembl ID
    query = """
    query searchTarget($queryString: String!) {
        search(queryString: $queryString, entityNames: ["target"], page: {size: 1, index: 0}) {
            hits {
                id
                name
                entity
            }
        }
    }
    """

    try:
        response = requests.post(
            api_url,
            json={'query': query, 'variables': {'queryString': gene_name}},
            timeout=30
        )
        response.raise_for_status()
        data = response.json()

        hits = data.get('data', {}).get('search', {}).get('hits', [])
        if hits:
            result['ensembl_id'] = hits[0].get('id')

        if not result['ensembl_id']:
            return result

        # Get drugs for target
        drug_query = """
        query knownDrugs($ensemblId: String!) {
            target(ensemblId: $ensemblId) {
                id
                approvedSymbol
                knownDrugs {
                    rows {
                        drug {
                            id
                            name
                            maximumClinicalTrialPhase
                            isApproved
                        }
                        mechanismOfAction
                        disease {
                            id
                            name
                        }
                    }
                }
            }
        }
        """

        response = requests.post(
            api_url,
            json={'query': drug_query, 'variables': {'ensemblId': result['ensembl_id']}},
            timeout=30
        )
        response.raise_for_status()
        drug_data = response.json()

        known_drugs = drug_data.get('data', {}).get('target', {}).get('knownDrugs', {}).get('rows', [])
        for drug_row in known_drugs:
            drug = drug_row.get('drug', {})
            result['drugs'].append({
                'drug_id': drug.get('id', ''),
                'name': drug.get('name', ''),
                'phase': drug.get('maximumClinicalTrialPhase', 0),
                'approved': drug.get('isApproved', False),
                'mechanism': drug_row.get('mechanismOfAction', ''),
                'disease': drug_row.get('disease', {}).get('name', '')
            })

    except requests.exceptions.RequestException as e:
        result['error'] = str(e)

    return result


def search_pubchem(gene_name: str) -> dict:
    """
    Search PubChem for bioactive compounds.

    Args:
        gene_name: Gene symbol

    Returns:
        Dictionary with compound information
    """
    result = {
        'gene': gene_name,
        'compounds': []
    }

    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    try:
        # Search for gene
        search_url = f"{base_url}/gene/symbol/{gene_name}/cids/JSON"
        response = requests.get(search_url, timeout=30)

        if response.status_code == 200:
            data = response.json()
            cids = data.get('InformationList', {}).get('Information', [{}])[0].get('CID', [])

            # Get compound details for top compounds
            for cid in cids[:10]:
                compound_url = f"{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName/JSON"
                comp_response = requests.get(compound_url, timeout=10)
                if comp_response.status_code == 200:
                    comp_data = comp_response.json()
                    props = comp_data.get('PropertyTable', {}).get('Properties', [{}])[0]
                    result['compounds'].append({
                        'cid': cid,
                        'name': props.get('IUPACName', ''),
                        'formula': props.get('MolecularFormula', ''),
                        'weight': props.get('MolecularWeight', '')
                    })

    except requests.exceptions.RequestException as e:
        result['error'] = str(e)

    return result


def search_clinical_trials(gene_name: str) -> dict:
    """
    Search ClinicalTrials.gov for trials targeting the gene.

    Args:
        gene_name: Gene symbol

    Returns:
        Dictionary with clinical trial information
    """
    result = {
        'gene': gene_name,
        'trials': []
    }

    base_url = "https://clinicaltrials.gov/api/v2/studies"

    try:
        params = {
            'query.term': f'{gene_name} AND (COVID-19 OR sepsis OR inflammation)',
            'pageSize': 10,
            'format': 'json'
        }
        response = requests.get(base_url, params=params, timeout=30)
        response.raise_for_status()
        data = response.json()

        for study in data.get('studies', []):
            protocol = study.get('protocolSection', {})
            id_module = protocol.get('identificationModule', {})
            status_module = protocol.get('statusModule', {})

            result['trials'].append({
                'nct_id': id_module.get('nctId', ''),
                'title': id_module.get('briefTitle', ''),
                'status': status_module.get('overallStatus', ''),
                'phase': ', '.join(protocol.get('designModule', {}).get('phases', []))
            })

    except requests.exceptions.RequestException as e:
        result['error'] = str(e)

    return result


def create_drug_summary(all_results: dict, output_dir: str, timestamp: str):
    """Create summary of drug search results."""
    rows = []

    for gene, results in all_results.items():
        chembl = results.get('chembl', {})
        opentargets = results.get('opentargets', {})
        trials = results.get('clinical_trials', {})

        # Count drugs by phase
        approved_drugs = []
        clinical_drugs = []

        # From ChEMBL
        approved_drugs.extend([d['name'] for d in chembl.get('approved_drugs', [])])
        clinical_drugs.extend([d['name'] for d in chembl.get('clinical_candidates', [])])

        # From OpenTargets
        for drug in opentargets.get('drugs', []):
            if drug.get('approved'):
                if drug['name'] not in approved_drugs:
                    approved_drugs.append(drug['name'])
            elif drug.get('phase', 0) >= 2:
                if drug['name'] not in clinical_drugs:
                    clinical_drugs.append(drug['name'])

        rows.append({
            'gene': gene,
            'chembl_target_id': chembl.get('chembl_target_id', ''),
            'n_compounds_chembl': len(chembl.get('compounds', [])),
            'n_approved_drugs': len(approved_drugs),
            'approved_drugs': '; '.join(approved_drugs[:5]),
            'n_clinical_candidates': len(clinical_drugs),
            'clinical_candidates': '; '.join(clinical_drugs[:5]),
            'n_clinical_trials': len(trials.get('trials', [])),
            'druggability': 'High' if approved_drugs else ('Medium' if clinical_drugs else 'Low')
        })

    df = pd.DataFrame(rows)
    df = df.sort_values('n_approved_drugs', ascending=False)

    summary_path = os.path.join(output_dir, f'drug_summary_{timestamp}.csv')
    df.to_csv(summary_path, index=False)
    print(f"\nSaved drug summary to {summary_path}")

    return df


def main():
    parser = argparse.ArgumentParser(description='Search drug databases for target drugs')
    parser.add_argument('--genes', '-g', type=str,
                        help='Comma-separated list of gene symbols')
    parser.add_argument('--candidates-csv', '-c',
                        help='Path to therapeutic candidates CSV')
    parser.add_argument('--output-dir', '-o', default='results/drug_search',
                        help='Output directory')
    parser.add_argument('--top-n', type=int, default=10,
                        help='Number of top candidates to search')

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
    print("DRUG DATABASE SEARCH")
    print(f"{'='*60}")
    print(f"Searching for: {', '.join(genes)}")

    os.makedirs(args.output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    all_results = {}

    for gene in genes:
        print(f"\n{'='*40}")
        print(f"Searching: {gene}")
        print(f"{'='*40}")

        results = {}

        # Search ChEMBL
        print("  ChEMBL...", end=' ')
        chembl_result = search_chembl_target(gene)
        results['chembl'] = chembl_result
        n_compounds = len(chembl_result.get('compounds', []))
        n_approved = len(chembl_result.get('approved_drugs', []))
        print(f"✓ ({n_compounds} compounds, {n_approved} approved)")
        time.sleep(0.5)

        # Search OpenTargets
        print("  OpenTargets...", end=' ')
        ot_result = search_opentargets(gene)
        results['opentargets'] = ot_result
        n_drugs = len(ot_result.get('drugs', []))
        print(f"✓ ({n_drugs} drug associations)")
        time.sleep(0.5)

        # Search Clinical Trials
        print("  ClinicalTrials.gov...", end=' ')
        ct_result = search_clinical_trials(gene)
        results['clinical_trials'] = ct_result
        n_trials = len(ct_result.get('trials', []))
        print(f"✓ ({n_trials} trials)")
        time.sleep(0.5)

        all_results[gene] = results

    # Create summary
    summary_df = create_drug_summary(all_results, args.output_dir, timestamp)

    # Save detailed results
    detailed_path = os.path.join(args.output_dir, f'drug_search_detailed_{timestamp}.json')
    with open(detailed_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Saved detailed results to {detailed_path}")

    # Print summary
    print(f"\n{'='*60}")
    print("DRUG SEARCH SUMMARY")
    print(f"{'='*60}")
    print(f"\n{'Gene':<12} {'Approved':<10} {'Clinical':<10} {'Druggability':<12}")
    print("-" * 50)
    for _, row in summary_df.iterrows():
        print(f"{row['gene']:<12} {row['n_approved_drugs']:<10} {row['n_clinical_candidates']:<10} {row['druggability']:<12}")

    # Highlight actionable targets
    print(f"\n{'='*60}")
    print("ACTIONABLE TARGETS (Existing Drugs Available)")
    print(f"{'='*60}")
    for _, row in summary_df.iterrows():
        if row['n_approved_drugs'] > 0:
            print(f"\n{row['gene']}:")
            print(f"  Approved drugs: {row['approved_drugs']}")

    print("\nDone!")


if __name__ == '__main__':
    main()
