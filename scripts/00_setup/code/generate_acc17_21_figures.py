#!/usr/bin/env python3
"""
Generate publication figures for ACC 17-21 extension modules.
Uses local result data (CSVs/MDs) — no GPU required.

Figures generated:
1. Domain adaptation benchmark (ACC 17) — bar chart
2. Multimodal fusion strategies (ACC 18) — bar chart
3. Model optimization (ACC 19) — multi-panel
4. Drug target precision@K (ACC 20) — line plot
5. Conditional generation metrics (ACC 21) — summary panel
6. Differential diagnosis per-disease AUC (ACC 15) — heatmap

Usage:
    python scripts/generate_acc17_21_figures.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Output directory
FIGURE_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)

# Publication style
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

COLORS = {
    'Geneformer_adapted': '#2ecc71',
    'STATE': '#e74c3c',
    'UCE': '#3498db',
    'TranscriptFormer': '#f39c12',
    'Geneformer': '#9b59b6',
    'scGPT': '#95a5a6',
    'scGPT_PCA256': '#7f8c8d',
    'highlight': '#e74c3c',
    'baseline': '#bdc3c7',
}


def fig_domain_adaptation():
    """ACC 17: Domain-adapted model benchmark — grouped bar chart."""
    models = ['Geneformer\n_adapted', 'STATE', 'UCE', 'Transcript\nFormer', 'Geneformer\n(base)', 'scGPT']
    binary = [0.879, 0.866, 0.859, 0.784, 0.773, 0.773]
    three_class = [0.891, 0.881, 0.874, 0.799, 0.781, 0.765]
    six_class = [0.884, 0.874, 0.871, 0.796, 0.771, 0.806]
    colors = ['#2ecc71', '#e74c3c', '#3498db', '#f39c12', '#9b59b6', '#95a5a6']

    x = np.arange(len(models))
    width = 0.25

    fig, ax = plt.subplots(figsize=(10, 5))
    bars1 = ax.bar(x - width, binary, width, label='Binary', color=colors, alpha=0.9, edgecolor='white')
    bars2 = ax.bar(x, three_class, width, label='3-Class', color=colors, alpha=0.7, edgecolor='white')
    bars3 = ax.bar(x + width, six_class, width, label='6-Class', color=colors, alpha=0.5, edgecolor='white')

    ax.set_ylabel('AUC-ROC')
    ax.set_title('Domain-Adapted Foundation Model Benchmark (PCA-256 + LogReg, 5-fold CV)')
    ax.set_xticks(x)
    ax.set_xticklabels(models)
    ax.legend(loc='lower right')
    ax.set_ylim(0.72, 0.92)
    ax.axhline(y=0.773, color='#9b59b6', linestyle='--', alpha=0.3, label='Base Geneformer')

    # Annotate improvement
    ax.annotate('+10.6%', xy=(0, 0.879), xytext=(0.3, 0.895),
                fontsize=10, fontweight='bold', color='#2ecc71',
                arrowprops=dict(arrowstyle='->', color='#2ecc71'))

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_domain_adaptation.png'))
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_domain_adaptation.pdf'))
    plt.close()
    print('  Generated: fig_domain_adaptation.png/pdf')


def fig_multimodal_fusion():
    """ACC 18: Fusion strategy comparison — horizontal bar chart."""
    strategies = ['Attention fusion', 'Late stacking', 'Late average', 'Late rank',
                  'Early concat', 'Expr + embedding', 'STATE baseline', 'Early PCA-256']
    aucs = [0.937, 0.915, 0.910, 0.909, 0.899, 0.895, 0.895, 0.883]
    colors = ['#2ecc71' if s == 'Attention fusion' else
              '#e74c3c' if 'STATE' in s else '#3498db' for s in strategies]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.barh(range(len(strategies)), aucs, color=colors, edgecolor='white', height=0.6)
    ax.set_yticks(range(len(strategies)))
    ax.set_yticklabels(strategies)
    ax.set_xlabel('AUC-ROC (Binary Severity)')
    ax.set_title('Multimodal Fusion Strategy Comparison')
    ax.set_xlim(0.87, 0.945)
    ax.invert_yaxis()

    # Annotate best
    ax.annotate('+4.2% vs STATE', xy=(0.937, 0), xytext=(0.94, 1.5),
                fontsize=10, fontweight='bold', color='#2ecc71',
                arrowprops=dict(arrowstyle='->', color='#2ecc71'))

    # Add value labels
    for i, (bar, v) in enumerate(zip(bars, aucs)):
        ax.text(v - 0.003, i, f'{v:.3f}', ha='right', va='center',
                fontweight='bold', color='white', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_multimodal_fusion.png'))
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_multimodal_fusion.pdf'))
    plt.close()
    print('  Generated: fig_multimodal_fusion.png/pdf')


def fig_optimization():
    """ACC 19: Model optimization — multi-panel figure."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    # Panel A: Quantization
    ax = axes[0]
    precisions = ['FP32', 'FP16', 'INT8\ndynamic', 'INT8\nstatic']
    aucs = [0.8938, 0.8938, 0.8934, 0.8899]
    speedups = [1.0, 2.0, 2.9, 3.3]
    colors_q = ['#95a5a6', '#3498db', '#2ecc71', '#f39c12']

    ax2 = ax.twinx()
    bars = ax.bar(range(len(precisions)), aucs, color=colors_q, alpha=0.8, edgecolor='white')
    ax2.plot(range(len(precisions)), speedups, 'ko-', markersize=8, linewidth=2)
    ax.set_xticks(range(len(precisions)))
    ax.set_xticklabels(precisions)
    ax.set_ylabel('AUC-ROC')
    ax2.set_ylabel('Speedup (x)')
    ax.set_title('A. Quantization')
    ax.set_ylim(0.885, 0.900)

    # Panel B: FSDP
    ax = axes[1]
    configs = ['1-GPU\nbaseline', '2-GPU\nFSDP']
    fsdp_auc = [0.930, 0.926]
    memory = [4.2, 3.7]
    colors_f = ['#3498db', '#2ecc71']

    ax2 = ax.twinx()
    bars = ax.bar(range(len(configs)), fsdp_auc, color=colors_f, alpha=0.8, width=0.4, edgecolor='white')
    ax2.bar([x + 0.4 for x in range(len(configs))], memory, color=colors_f, alpha=0.4, width=0.4, edgecolor='white')
    ax.set_xticks([x + 0.2 for x in range(len(configs))])
    ax.set_xticklabels(configs)
    ax.set_ylabel('AUC-ROC')
    ax2.set_ylabel('Peak Memory (GB)')
    ax.set_title('B. FSDP Distributed Training')
    ax.set_ylim(0.91, 0.94)
    ax2.set_ylim(3.0, 5.0)

    # Annotate memory savings
    ax.annotate('-12% memory', xy=(1.4, 0.926), xytext=(1.2, 0.935),
                fontsize=9, fontweight='bold', color='#2ecc71')

    # Panel C: torch.compile
    ax = axes[2]
    models_c = ['UCE', 'Gene-\nformer', 'Transcript\nFormer', 'STATE', 'scGPT']
    speedups_c = [1.13, 1.13, 1.05, 1.02, 0.83]
    colors_c = ['#3498db', '#9b59b6', '#f39c12', '#e74c3c', '#95a5a6']

    bars = ax.bar(range(len(models_c)), speedups_c, color=colors_c, edgecolor='white')
    ax.axhline(y=1.0, color='black', linestyle='--', alpha=0.3)
    ax.set_xticks(range(len(models_c)))
    ax.set_xticklabels(models_c)
    ax.set_ylabel('Speedup (x)')
    ax.set_title('C. torch.compile (max-autotune)')
    ax.set_ylim(0.7, 1.25)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_optimization.png'))
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_optimization.pdf'))
    plt.close()
    print('  Generated: fig_optimization.png/pdf')


def fig_drug_targets():
    """ACC 20: Precision@K curves."""
    k_values = [5, 10, 25, 50, 100, 200, 500, 1000]

    # STATE precision@K
    state_precision = [0.0, 0.0, 0.0, 0.0, 0.0, 0.010, 0.006, 0.004]
    # Random baseline: 9/2561 ≈ 0.00351
    random_precision = [9/2561] * len(k_values)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Panel A: Precision@K
    ax1.plot(k_values, state_precision, 'o-', color='#e74c3c', linewidth=2,
             markersize=6, label='STATE')
    ax1.plot(k_values, random_precision, '--', color='#95a5a6', linewidth=1.5,
             label='Random baseline')
    ax1.set_xlabel('K (top candidates)')
    ax1.set_ylabel('Precision@K')
    ax1.set_title('A. Precision@K for Known Drug Targets')
    ax1.set_xscale('log')
    ax1.legend()
    ax1.set_xlim(3, 1500)

    # Panel B: Top-10 validation
    models_v = ['UCE', 'STATE']
    plausible = [5, 2]
    novel = [5, 8]

    x = np.arange(len(models_v))
    width = 0.35
    ax2.bar(x - width/2, plausible, width, label='Plausible (ChEMBL)', color='#2ecc71')
    ax2.bar(x + width/2, novel, width, label='Novel', color='#3498db')
    ax2.set_xticks(x)
    ax2.set_xticklabels(models_v)
    ax2.set_ylabel('Count (of top 10)')
    ax2.set_title('B. Top-10 Candidate ChEMBL Validation')
    ax2.legend()
    ax2.set_ylim(0, 10)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_drug_targets.png'))
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_drug_targets.pdf'))
    plt.close()
    print('  Generated: fig_drug_targets.png/pdf')


def fig_conditional_generation():
    """ACC 21: Conditional generation summary metrics."""
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    # Panel A: Cosine similarity gauge
    ax = axes[0]
    val = 0.784
    colors_seg = ['#e74c3c', '#f39c12', '#2ecc71']
    ax.barh(0, val, color='#2ecc71', height=0.4, edgecolor='white')
    ax.barh(0, 1.0 - val, left=val, color='#ecf0f1', height=0.4, edgecolor='white')
    ax.set_xlim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel('Cosine Similarity')
    ax.set_title('A. Recovery Direction Alignment')
    ax.text(val/2, 0, f'{val:.3f}', ha='center', va='center', fontweight='bold',
            fontsize=14, color='white')

    # Panel B: Classifier agreement
    ax = axes[1]
    classified = 64.0
    not_classified = 36.0
    ax.pie([classified, not_classified], labels=['Recovered', 'Other'],
           colors=['#2ecc71', '#e74c3c'], autopct='%1.1f%%',
           startangle=90, textprops={'fontsize': 11})
    ax.set_title('B. Independent Classifier\nAgreement')

    # Panel C: Key metrics summary
    ax = axes[2]
    metrics = ['Severe cells\nencoded', 'Real recovered\ncells', 'Cosine\nsimilarity',
               'Classified as\nrecovered', 'Mean\nP(recovered)']
    values = ['19,805', '3,695', '0.784', '64.0%', '0.592']
    ax.axis('off')
    ax.set_title('C. Generation Quality Summary')
    for i, (m, v) in enumerate(zip(metrics, values)):
        y = 0.85 - i * 0.18
        ax.text(0.1, y, m, fontsize=10, va='center', ha='left')
        ax.text(0.85, y, v, fontsize=11, va='center', ha='right', fontweight='bold',
                color='#2ecc71' if i >= 2 else '#34495e')

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_conditional_generation.png'))
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_conditional_generation.pdf'))
    plt.close()
    print('  Generated: fig_conditional_generation.png/pdf')


def fig_differential_diagnosis():
    """ACC 15: Per-disease AUC heatmap."""
    models = ['STATE', 'UCE', 'TranscriptFormer', 'scGPT (PCA-256)', 'Geneformer', 'scGPT (raw)']
    diseases = ['COVID-19', 'Healthy', 'Pneumonia', 'SLE', 'Sepsis']
    data = np.array([
        [0.958, 0.961, 0.973, 0.999, 0.998],
        [0.941, 0.948, 0.963, 0.999, 0.995],
        [0.923, 0.919, 0.947, 0.998, 0.989],
        [0.909, 0.889, 0.980, 0.999, 0.994],
        [0.892, 0.915, 0.908, 0.996, 0.984],
        [0.827, 0.715, 0.954, 0.993, 0.976],
    ])

    fig, ax = plt.subplots(figsize=(8, 5))
    im = ax.imshow(data, cmap='RdYlGn', aspect='auto', vmin=0.7, vmax=1.0)
    ax.set_xticks(range(len(diseases)))
    ax.set_xticklabels(diseases, rotation=45, ha='right')
    ax.set_yticks(range(len(models)))
    ax.set_yticklabels(models)
    ax.set_title('Per-Disease One-vs-Rest AUC (5-Class Differential Diagnosis)')

    # Annotate each cell
    for i in range(len(models)):
        for j in range(len(diseases)):
            text_color = 'white' if data[i, j] < 0.85 else 'black'
            ax.text(j, i, f'{data[i, j]:.3f}', ha='center', va='center',
                    fontsize=9, color=text_color, fontweight='bold')

    plt.colorbar(im, ax=ax, label='AUC-ROC', shrink=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_differential_diagnosis.png'))
    plt.savefig(os.path.join(FIGURE_DIR, 'fig_differential_diagnosis.pdf'))
    plt.close()
    print('  Generated: fig_differential_diagnosis.png/pdf')


if __name__ == '__main__':
    print('Generating ACC 15 & 17-21 publication figures...')
    print(f'Output directory: {FIGURE_DIR}')
    print()

    fig_differential_diagnosis()
    fig_domain_adaptation()
    fig_multimodal_fusion()
    fig_optimization()
    fig_drug_targets()
    fig_conditional_generation()

    print()
    print('All figures generated successfully.')
    print(f'Total: 6 figure sets (PNG + PDF) in {FIGURE_DIR}/')
