#!/usr/bin/env python3
"""
Prepare structures for molecular docking.
Supports AutoDock Vina and DiffDock.

Usage:
    python prepare_docking.py \
        --structure protein.pdb \
        --ligands ligands.sdf \
        --output docking_inputs/

Requirements:
    pip install rdkit biopython pandas numpy

For docking:
    - AutoDock Vina: https://vina.scripps.edu/
    - DiffDock: https://github.com/gcorso/DiffDock
"""

import argparse
import os
import subprocess
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np


def prepare_receptor_vina(pdb_path: str, output_dir: str) -> str:
    """
    Prepare receptor for AutoDock Vina.
    Converts PDB to PDBQT format.

    Args:
        pdb_path: Path to input PDB file
        output_dir: Output directory

    Returns:
        Path to prepared PDBQT file
    """
    os.makedirs(output_dir, exist_ok=True)

    base_name = Path(pdb_path).stem
    pdbqt_path = os.path.join(output_dir, f"{base_name}_receptor.pdbqt")

    print(f"\nPreparing receptor: {pdb_path}")

    # Try using prepare_receptor from AutoDock Tools
    try:
        # Method 1: Using prepare_receptor4.py from MGLTools
        cmd = [
            'prepare_receptor4.py',
            '-r', pdb_path,
            '-o', pdbqt_path,
            '-A', 'hydrogens'
        ]
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"  Created: {pdbqt_path}")
        return pdbqt_path

    except (subprocess.CalledProcessError, FileNotFoundError):
        print("  WARNING: prepare_receptor4.py not found")
        print("  Creating simplified PDBQT manually...")

        # Method 2: Manual conversion (simplified)
        convert_pdb_to_pdbqt(pdb_path, pdbqt_path)
        return pdbqt_path


def convert_pdb_to_pdbqt(pdb_path: str, pdbqt_path: str):
    """
    Simple PDB to PDBQT conversion.
    Note: For production use, use proper tools like MGLTools.
    """
    # Atom type charges (simplified)
    charges = {
        'C': 0.0, 'N': -0.4, 'O': -0.4, 'S': -0.1,
        'H': 0.1, 'CA': 0.1, 'FE': 0.3, 'ZN': 0.5
    }

    with open(pdb_path, 'r') as f_in, open(pdbqt_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Get atom element
                element = line[76:78].strip() if len(line) > 76 else line[12:14].strip()[0]

                # Get charge
                charge = charges.get(element, 0.0)

                # Create PDBQT line (add charge and atom type)
                pdbqt_line = f"{line[:54]}{charge:6.3f} {element:>2}  \n"
                f_out.write(pdbqt_line)
            elif line.startswith('END'):
                f_out.write(line)

    print(f"  Created simplified PDBQT: {pdbqt_path}")


def define_binding_box(pdb_path: str, center: tuple = None, size: tuple = None) -> dict:
    """
    Define binding box for docking.

    If center/size not provided, calculates box around protein center.

    Returns:
        Dictionary with center_x, center_y, center_z, size_x, size_y, size_z
    """
    if center and size:
        return {
            'center_x': center[0],
            'center_y': center[1],
            'center_z': center[2],
            'size_x': size[0],
            'size_y': size[1],
            'size_z': size[2]
        }

    # Calculate from structure
    coords = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    pass

    if coords:
        coords = np.array(coords)
        center = coords.mean(axis=0)
        # Box size covers protein + margin
        min_coords = coords.min(axis=0)
        max_coords = coords.max(axis=0)
        box_size = max_coords - min_coords + 20  # 20 Å margin

        return {
            'center_x': float(center[0]),
            'center_y': float(center[1]),
            'center_z': float(center[2]),
            'size_x': float(box_size[0]),
            'size_y': float(box_size[1]),
            'size_z': float(box_size[2])
        }

    return None


def create_vina_config(receptor_path: str, box: dict, output_path: str, exhaustiveness: int = 8):
    """Create Vina configuration file."""
    config = f"""# AutoDock Vina configuration
receptor = {receptor_path}

center_x = {box['center_x']:.3f}
center_y = {box['center_y']:.3f}
center_z = {box['center_z']:.3f}

size_x = {box['size_x']:.1f}
size_y = {box['size_y']:.1f}
size_z = {box['size_z']:.1f}

exhaustiveness = {exhaustiveness}
num_modes = 9
energy_range = 3
"""

    with open(output_path, 'w') as f:
        f.write(config)

    print(f"  Created Vina config: {output_path}")
    return config


def prepare_diffdock_input(protein_path: str, ligand_smiles: list, output_dir: str) -> str:
    """
    Prepare input files for DiffDock.

    Args:
        protein_path: Path to protein PDB
        ligand_smiles: List of SMILES strings
        output_dir: Output directory

    Returns:
        Path to input CSV file
    """
    os.makedirs(output_dir, exist_ok=True)

    # Create input CSV for DiffDock
    rows = []
    for i, smiles in enumerate(ligand_smiles):
        rows.append({
            'complex_name': f'complex_{i}',
            'protein_path': os.path.abspath(protein_path),
            'ligand_description': smiles,
            'protein_sequence': ''  # Optional
        })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(output_dir, 'diffdock_input.csv')
    df.to_csv(csv_path, index=False)

    print(f"  Created DiffDock input: {csv_path}")
    return csv_path


def get_drugbank_compounds(gene_name: str) -> list:
    """
    Get known drug SMILES for a gene target.
    Returns list of (drug_name, smiles) tuples.
    """
    # This would query DrugBank API in production
    # For now, return example compounds for key targets

    known_drugs = {
        'PIK3R2': [
            ('Alpelisib', 'CC1=CC(=NN1C2=CC=NC=C2C3=NC=NC4=C3SC=C4)C(C)(C)C(=O)NC5=C(C=CC(=C5)F)OC'),
            ('Copanlisib', 'COc1cc2nc(N3CCOCC3)nc(Nc4ccc(cc4)S(=O)(=O)C)c2cc1OC'),
        ],
        'ICMT': [
            ('Cysmethynil', 'COC(=O)C1=CC=C(C=C1)SC2=NC(=O)C(=CN2)C'),
        ],
        'PINK1': [
            ('Kinetin', 'C1=CN=C2C(=N1)C(=O)NC(=N2)N'),
        ],
    }

    return known_drugs.get(gene_name, [])


def create_docking_script(target_name: str, output_dir: str, method: str = 'vina') -> str:
    """Create shell script to run docking."""
    script_path = os.path.join(output_dir, f'run_docking_{target_name}.sh')

    if method == 'vina':
        script = f"""#!/bin/bash
# AutoDock Vina docking script for {target_name}
# Run from the docking directory

CONFIG="{target_name}_vina_config.txt"
RECEPTOR="{target_name}_receptor.pdbqt"
LIGAND="${{1:-ligand.pdbqt}}"
OUTPUT="${{2:-docking_output}}"

echo "Running AutoDock Vina..."
echo "  Receptor: $RECEPTOR"
echo "  Ligand: $LIGAND"
echo "  Config: $CONFIG"

vina --config $CONFIG --ligand $LIGAND --out ${{OUTPUT}}.pdbqt --log ${{OUTPUT}}.log

echo "Done! Results in ${{OUTPUT}}.pdbqt"
"""
    else:  # DiffDock
        script = f"""#!/bin/bash
# DiffDock docking script for {target_name}
# Requires DiffDock installation

INPUT_CSV="diffdock_input.csv"
OUTPUT_DIR="diffdock_results"

echo "Running DiffDock..."

python -m inference \
    --protein_ligand_csv $INPUT_CSV \
    --out_dir $OUTPUT_DIR \
    --inference_steps 20 \
    --samples_per_complex 10 \
    --batch_size 10

echo "Done! Results in $OUTPUT_DIR"
"""

    with open(script_path, 'w') as f:
        f.write(script)

    os.chmod(script_path, 0o755)
    print(f"  Created docking script: {script_path}")
    return script_path


def main():
    parser = argparse.ArgumentParser(description='Prepare structures for molecular docking')
    parser.add_argument('--structures-dir', '-s', required=True,
                        help='Directory with AlphaFold PDB structures')
    parser.add_argument('--output-dir', '-o', default='results/docking',
                        help='Output directory')
    parser.add_argument('--method', '-m', choices=['vina', 'diffdock', 'both'],
                        default='both', help='Docking method to prepare')
    parser.add_argument('--exhaustiveness', type=int, default=8,
                        help='Vina exhaustiveness parameter')

    args = parser.parse_args()

    print(f"\n{'='*60}")
    print("MOLECULAR DOCKING PREPARATION")
    print(f"{'='*60}")

    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

    # Find PDB files
    pdb_files = list(Path(args.structures_dir).glob('*.pdb'))
    print(f"Found {len(pdb_files)} PDB structures")

    # Process each structure
    summary = []

    for pdb_path in pdb_files:
        pdb_name = pdb_path.stem
        # Extract gene name (format: AF-UniProtID-F1-model_v4)
        if pdb_name.startswith('AF-'):
            parts = pdb_name.split('-')
            uniprot_id = parts[1] if len(parts) > 1 else pdb_name
        else:
            uniprot_id = pdb_name

        print(f"\n{'='*40}")
        print(f"Processing: {uniprot_id}")
        print(f"{'='*40}")

        target_dir = os.path.join(args.output_dir, uniprot_id)
        os.makedirs(target_dir, exist_ok=True)

        target_info = {
            'uniprot_id': uniprot_id,
            'pdb_path': str(pdb_path),
            'docking_dir': target_dir
        }

        # Define binding box
        box = define_binding_box(str(pdb_path))
        print(f"  Binding box: center=({box['center_x']:.1f}, {box['center_y']:.1f}, {box['center_z']:.1f})")
        print(f"               size=({box['size_x']:.1f}, {box['size_y']:.1f}, {box['size_z']:.1f})")

        if args.method in ['vina', 'both']:
            # Prepare for Vina
            pdbqt_path = prepare_receptor_vina(str(pdb_path), target_dir)
            target_info['pdbqt_path'] = pdbqt_path

            config_path = os.path.join(target_dir, f'{uniprot_id}_vina_config.txt')
            create_vina_config(pdbqt_path, box, config_path, args.exhaustiveness)
            target_info['vina_config'] = config_path

            create_docking_script(uniprot_id, target_dir, 'vina')

        if args.method in ['diffdock', 'both']:
            # Get known drugs for target
            # In production, this would query DrugBank
            known_smiles = []

            if known_smiles:
                diffdock_csv = prepare_diffdock_input(str(pdb_path), known_smiles, target_dir)
                target_info['diffdock_input'] = diffdock_csv
                create_docking_script(uniprot_id, target_dir, 'diffdock')

        summary.append(target_info)

    # Save summary
    summary_df = pd.DataFrame(summary)
    summary_path = os.path.join(args.output_dir, f'docking_summary_{timestamp}.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSaved summary to {summary_path}")

    # Print instructions
    print(f"\n{'='*60}")
    print("DOCKING INSTRUCTIONS")
    print(f"{'='*60}")

    print("\n## AutoDock Vina")
    print("1. Install Vina: https://vina.scripps.edu/")
    print("2. Prepare ligand as PDBQT (use MGLTools or OpenBabel)")
    print("3. Run: ./run_docking_<target>.sh ligand.pdbqt output")

    print("\n## DiffDock")
    print("1. Install: https://github.com/gcorso/DiffDock")
    print("2. Prepare ligand SMILES in diffdock_input.csv")
    print("3. Run: python -m inference --protein_ligand_csv diffdock_input.csv")

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"Structures prepared: {len(summary)}")
    print(f"Output directory: {args.output_dir}")

    print("\nDone!")


if __name__ == '__main__':
    main()
