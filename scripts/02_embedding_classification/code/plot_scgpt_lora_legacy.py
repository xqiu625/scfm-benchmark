#!/usr/bin/env python3
"""
Legacy scPEFT/LoRA-classifier benchmark figure (Jan 2026).

Reproduces the parameter-efficient fine-tuning result referenced in the CV bullet:
  "Performed parameter-efficient fine-tuning (LoRA) of scGPT via scPEFT."

Source: scripts/02_embedding_classification/results/benchmark_all_models_LoRA_Classifier.csv
        (LoRA-style adapters on frozen embeddings, classifier head only;
         scPEFT-inspired protocol; 5-fold stratified CV; focal loss)

Two-panel figure:
  Panel A: Zero-shot AUC vs LoRA-Classifier AUC per FM (binary task)
           Grouped bars; Delta-AUC annotated above each pair.
  Panel B: Delta-AUC across the three classification tasks (binary, 3-class, 6-class)
           Bars colored red for negative (scGPT family), green for positive (other FMs).
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DATA = '/Users/xinruqiu/Desktop/Files/302006/302006_updategit/scripts/02_embedding_classification/results/benchmark_all_models_LoRA_Classifier.csv'
OUT_DIR = '/Users/xinruqiu/Desktop/Files/302006/Platelet-FM-Benchmark/figures'
os.makedirs(OUT_DIR, exist_ok=True)

MODEL_ORDER = ['STATE', 'UCE', 'TranscriptFormer', 'Geneformer', 'scGPT', 'scGPT_BP']
TASK_ORDER  = ['binary', '3class', '6class']
TASK_LABEL  = {'binary': 'Binary', '3class': '3-class', '6class': '6-class'}

# Subtle FM color palette (Tableau-ish)
COLORS = {
    'STATE':            '#1F77B4',
    'UCE':              '#FF7F0E',
    'TranscriptFormer': '#2CA02C',
    'Geneformer':       '#9467BD',
    'scGPT':            '#D62728',
    'scGPT_BP':         '#8C564B',
}


def panel_a(ax, df):
    """Zero-shot vs LoRA grouped bars for BINARY task."""
    d = df[df['Task'] == 'binary'].set_index('Model').reindex(MODEL_ORDER)
    x = np.arange(len(MODEL_ORDER))
    width = 0.38

    zs = d['Zero_Shot_AUC'].values
    lr = d['LoRA_Classifier_AUC'].values
    delta = lr - zs

    bars_zs = ax.bar(x - width / 2, zs, width, label='Zero-shot (LR on frozen embedding)',
                     color='#BBBBBB', edgecolor='black', linewidth=0.5)
    bars_lr = ax.bar(x + width / 2, lr, width, label='LoRA-Classifier (scPEFT-style)',
                     color=[COLORS[m] for m in MODEL_ORDER],
                     edgecolor='black', linewidth=0.5, alpha=0.9)

    # Delta annotation above each pair
    for xi, (z, l, dv) in enumerate(zip(zs, lr, delta)):
        ytop = max(z, l) + 0.02
        color = '#2CA02C' if dv > 0 else '#D62728'
        sym = '+' if dv > 0 else ''
        ax.text(xi, ytop, f'{sym}{dv*100:.1f}%', ha='center', va='bottom',
                fontsize=8.5, color=color, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels(MODEL_ORDER, rotation=20, ha='right', fontsize=9)
    ax.set_ylim(0.40, 1.00)
    ax.set_ylabel('Binary classification AUC')
    ax.set_title('A. scPEFT-style LoRA-Classifier vs zero-shot (binary task)',
                 loc='left', fontsize=11)
    ax.legend(loc='lower right', fontsize=8, frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.axhline(0.5, color='#888888', linestyle=':', linewidth=0.7, alpha=0.5)


def panel_b(ax, df):
    """Delta AUC across the 3 tasks, per FM."""
    d = df.copy()
    d['delta'] = d['LoRA_Classifier_AUC'] - d['Zero_Shot_AUC']

    x = np.arange(len(MODEL_ORDER))
    width = 0.26
    for i, task in enumerate(TASK_ORDER):
        vals = (d[d['Task'] == task].set_index('Model').reindex(MODEL_ORDER)['delta'].values)
        offset = (i - 1) * width
        colors = ['#D62728' if v < 0 else '#2CA02C' for v in vals]
        alpha = [0.45, 0.7, 1.0][i]
        ax.bar(x + offset, vals * 100, width, color=colors, alpha=alpha,
               edgecolor='black', linewidth=0.5, label=TASK_LABEL[task])
        for xi, v in zip(x + offset, vals):
            sign = '+' if v > 0 else ''
            # Pull the value label closer to the bar tip and away from neighbors
            pad = 0.6
            ytop = v * 100 + (pad if v >= 0 else -pad)
            ax.text(xi, ytop, f'{sign}{v*100:.1f}', ha='center',
                    va='bottom' if v >= 0 else 'top', fontsize=7.5, color='#333333')

    ax.axhline(0, color='black', linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(MODEL_ORDER, rotation=20, ha='right', fontsize=9)
    # Expanded y-range so the -24.3 label and the failure-mode annotation both fit cleanly
    ax.set_ylim(-32, 12)
    ax.set_ylabel('$\\Delta$ AUC vs zero-shot (%)')
    ax.set_title('B. $\\Delta$ AUC across tasks (alpha-shaded by task)',
                 loc='left', fontsize=11)
    # Legend in upper-right corner where there is empty whitespace
    ax.legend(loc='upper right', fontsize=8, frameon=False, ncol=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Annotation: scGPT family failure mode — placed in lower-left dead zone
    # well below any bar, with arrow tip on the scGPT_BP 6-class bar
    ax.annotate('scGPT family: LoRA-classifier degrades performance\n'
                '(generative-pretrained embeddings do not compose with adapter heads)',
                xy=(5 + width, -22), xytext=(0.1, -30),
                fontsize=8, color='#7A0000', ha='left',
                arrowprops=dict(arrowstyle='->', color='#7A0000', lw=0.8))


def main():
    df = pd.read_csv(DATA)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(14.5, 5.8),
                                   gridspec_kw={'width_ratios': [1.0, 1.25], 'wspace': 0.24})
    panel_a(axA, df)
    panel_b(axB, df)
    fig.suptitle('Parameter-efficient fine-tuning of scGPT via scPEFT '
                 '(legacy severity benchmark, Jan 2026)',
                 fontsize=12, y=1.02)

    pdf = os.path.join(OUT_DIR, 'fig_scgpt_lora_scpeft.pdf')
    png = os.path.join(OUT_DIR, 'fig_scgpt_lora_scpeft.png')
    fig.savefig(pdf, bbox_inches='tight')
    fig.savefig(png, bbox_inches='tight', dpi=200)
    print(f'Wrote {pdf}')
    print(f'Wrote {png}')


if __name__ == '__main__':
    main()
