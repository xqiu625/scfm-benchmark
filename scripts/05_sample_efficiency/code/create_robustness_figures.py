#!/usr/bin/env python3
"""
Generate Publication-Ready Figures for Robustness Analysis

Creates figures from existing benchmark results:
1. Sample efficiency curves
2. Gene dropout robustness
3. Batch-shift generalization
4. Clinical metrics comparison
5. Summary robustness figure

Usage:
    python create_robustness_figures.py
    python create_robustness_figures.py --output-dir /path/to/figures

Author: Project 302006
Date: January 2026
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

# Model colors (consistent across figures)
MODEL_COLORS = {
    'STATE': '#E74C3C',      # Red
    'UCE': '#3498DB',        # Blue
    'Geneformer': '#2ECC71', # Green
    'TranscriptFormer': '#9B59B6',  # Purple
    'scGPT': '#F39C12',      # Orange
    'XGBoost': '#7F8C8D',    # Gray
    'PCA': '#95A5A6',        # Light gray
}


def create_sample_efficiency_figure(output_dir):
    """Create sample efficiency curves figure"""

    # Data from verified results
    data = {
        'STATE': {'fracs': [1.0, 0.5, 0.2, 0.1, 0.05], 'aucs': [0.944, 0.935, 0.910, 0.875, 0.821]},
        'UCE': {'fracs': [1.0, 0.5, 0.2, 0.1, 0.05], 'aucs': [0.910, 0.905, 0.890, 0.870, 0.847]},
        'Geneformer': {'fracs': [1.0, 0.5, 0.2, 0.1, 0.05], 'aucs': [0.845, 0.830, 0.800, 0.760, 0.718]},
        'scGPT': {'fracs': [1.0, 0.5, 0.2, 0.1, 0.05], 'aucs': [0.735, 0.725, 0.710, 0.695, 0.706]},
    }

    fig, ax = plt.subplots(figsize=(10, 6))

    for model, vals in data.items():
        color = MODEL_COLORS.get(model, '#333333')
        ax.plot(vals['fracs'], vals['aucs'], 'o-', label=model, color=color,
                linewidth=2.5, markersize=8)

    ax.set_xlabel('Training Data Fraction', fontsize=12)
    ax.set_ylabel('AUC-ROC', fontsize=12)
    ax.set_title('Sample Efficiency: Performance vs Training Data Size', fontsize=14, fontweight='bold')
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0.65, 1.0)
    ax.legend(loc='lower right', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add annotations
    ax.axhline(y=0.897, color='gray', linestyle='--', alpha=0.5, label='XGBoost baseline')
    ax.text(0.02, 0.90, 'XGBoost baseline', fontsize=9, color='gray')

    plt.tight_layout()

    # Save
    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig_sample_efficiency.{fmt}', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: fig_sample_efficiency.png/pdf")


def create_gene_dropout_figure(output_dir):
    """Create gene dropout robustness figure"""

    # Data from verified results
    dropout_levels = [0, 10, 30, 50, 70]
    data = {
        'UCE': [0.876, 0.874, 0.870, 0.867, 0.861],
        'STATE': [0.893, 0.890, 0.882, 0.873, 0.862],
    }

    fig, ax = plt.subplots(figsize=(10, 6))

    for model, aucs in data.items():
        color = MODEL_COLORS.get(model, '#333333')
        ax.plot(dropout_levels, aucs, 'o-', label=model, color=color,
                linewidth=2.5, markersize=8)

    ax.set_xlabel('Gene Dropout Rate (%)', fontsize=12)
    ax.set_ylabel('AUC-ROC', fontsize=12)
    ax.set_title('Gene Dropout Robustness: Performance vs Missing Genes', fontsize=14, fontweight='bold')
    ax.set_xlim(-5, 75)
    ax.set_ylim(0.80, 0.95)
    ax.legend(loc='lower left', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Add drop annotations
    ax.annotate('UCE: -1.7%', xy=(70, 0.861), xytext=(55, 0.84),
                fontsize=10, color=MODEL_COLORS['UCE'],
                arrowprops=dict(arrowstyle='->', color=MODEL_COLORS['UCE']))
    ax.annotate('STATE: -3.4%', xy=(70, 0.862), xytext=(55, 0.89),
                fontsize=10, color=MODEL_COLORS['STATE'],
                arrowprops=dict(arrowstyle='->', color=MODEL_COLORS['STATE']))

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig_gene_dropout.{fmt}', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: fig_gene_dropout.png/pdf")


def create_batch_shift_figure(output_dir):
    """Create batch-shift generalization figure"""

    models = ['STATE', 'UCE', 'Geneformer', 'scGPT']
    aucs = [0.789, 0.750, 0.730, 0.524]
    colors = [MODEL_COLORS[m] for m in models]

    fig, ax = plt.subplots(figsize=(8, 6))

    bars = ax.bar(models, aucs, color=colors, edgecolor='black', linewidth=1.5)

    ax.set_ylabel('AUC-ROC', fontsize=12)
    ax.set_title('Batch-Shift Generalization\n(Train: COVID → Test: Sepsis)', fontsize=14, fontweight='bold')
    ax.set_ylim(0, 1.0)

    # Add value labels
    for bar, auc in zip(bars, aucs):
        height = bar.get_height()
        ax.annotate(f'{auc:.3f}',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    # Add random baseline
    ax.axhline(y=0.5, color='red', linestyle='--', alpha=0.7, linewidth=2)
    ax.text(3.5, 0.52, 'Random', fontsize=10, color='red')

    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig_batch_generalization.{fmt}', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: fig_batch_generalization.png/pdf")


def create_clinical_metrics_figure(output_dir):
    """Create clinical metrics comparison figure"""

    models = ['STATE', 'UCE', 'TranscriptFormer', 'Geneformer', 'scGPT']
    metrics = {
        'AUC-ROC': [0.894, 0.876, 0.838, 0.824, 0.776],
        'AUC-PR': [0.934, 0.923, 0.895, 0.888, 0.836],
        'Sens@90%Spec': [0.717, 0.674, 0.585, 0.564, 0.408],
        "Cohen's κ": [0.611, 0.568, 0.501, 0.475, 0.402],
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, (metric_name, values) in zip(axes.flat, metrics.items()):
        colors = [MODEL_COLORS[m] for m in models]
        bars = ax.bar(models, values, color=colors, edgecolor='black', linewidth=1)

        ax.set_ylabel(metric_name, fontsize=11)
        ax.set_title(metric_name, fontsize=12, fontweight='bold')

        if 'Sens' in metric_name or 'κ' in metric_name:
            ax.set_ylim(0, 1.0)
        else:
            ax.set_ylim(0.7, 1.0)

        # Add value labels
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax.annotate(f'{val:.3f}',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 2),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=9)

        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')

    fig.suptitle('Clinical Utility Metrics (Binary Classification)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig_clinical_metrics.{fmt}', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: fig_clinical_metrics.png/pdf")


def create_robustness_summary_figure(output_dir):
    """Create comprehensive robustness summary figure"""

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Sample Efficiency (drop at 5%)
    models = ['UCE', 'scGPT', 'STATE', 'Geneformer']
    drops = [-7, -4, -13, -15]
    colors = [MODEL_COLORS[m] for m in models]

    ax = axes[0]
    bars = ax.barh(models, drops, color=colors, edgecolor='black')
    ax.set_xlabel('% Drop at 5% Training Data', fontsize=11)
    ax.set_title('A. Sample Efficiency', fontsize=12, fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.set_xlim(-20, 5)
    for bar, drop in zip(bars, drops):
        ax.annotate(f'{drop}%', xy=(drop - 1, bar.get_y() + bar.get_height()/2),
                    ha='right', va='center', fontsize=10, fontweight='bold')

    # Panel B: Gene Dropout (drop at 70%)
    models = ['UCE', 'STATE']
    drops = [-1.7, -3.4]
    colors = [MODEL_COLORS[m] for m in models]

    ax = axes[1]
    bars = ax.barh(models, drops, color=colors, edgecolor='black')
    ax.set_xlabel('% Drop at 70% Gene Dropout', fontsize=11)
    ax.set_title('B. Gene Dropout Robustness', fontsize=12, fontweight='bold')
    ax.axvline(x=0, color='black', linewidth=0.5)
    ax.set_xlim(-5, 1)
    for bar, drop in zip(bars, drops):
        ax.annotate(f'{drop}%', xy=(drop - 0.3, bar.get_y() + bar.get_height()/2),
                    ha='right', va='center', fontsize=10, fontweight='bold')

    # Panel C: Batch Generalization
    models = ['STATE', 'UCE', 'Geneformer', 'scGPT']
    aucs = [0.789, 0.750, 0.730, 0.524]
    colors = [MODEL_COLORS[m] for m in models]

    ax = axes[2]
    bars = ax.barh(models, aucs, color=colors, edgecolor='black')
    ax.set_xlabel('AUC (COVID → Sepsis)', fontsize=11)
    ax.set_title('C. Batch-Shift Generalization', fontsize=12, fontweight='bold')
    ax.axvline(x=0.5, color='red', linestyle='--', linewidth=1.5, label='Random')
    ax.set_xlim(0.4, 0.85)
    for bar, auc in zip(bars, aucs):
        ax.annotate(f'{auc:.3f}', xy=(auc + 0.01, bar.get_y() + bar.get_height()/2),
                    ha='left', va='center', fontsize=10, fontweight='bold')

    fig.suptitle('Foundation Model Robustness Analysis', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()

    for fmt in ['png', 'pdf']:
        fig.savefig(output_dir / f'fig_robustness_summary.{fmt}', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: fig_robustness_summary.png/pdf")


def main():
    parser = argparse.ArgumentParser(description='Generate Robustness Figures')
    parser.add_argument('--output-dir', type=str,
                        default="/bigdata/godziklab/shared/Xinru/302006/04_FIGURES",
                        help='Output directory for figures')
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("="*60)
    print("GENERATING ROBUSTNESS FIGURES")
    print("="*60)
    print(f"Output directory: {output_dir}")
    print()

    # Generate all figures
    create_sample_efficiency_figure(output_dir)
    create_gene_dropout_figure(output_dir)
    create_batch_shift_figure(output_dir)
    create_clinical_metrics_figure(output_dir)
    create_robustness_summary_figure(output_dir)

    print()
    print("="*60)
    print("All figures generated successfully!")
    print("="*60)

    return 0


if __name__ == "__main__":
    sys.exit(main())
