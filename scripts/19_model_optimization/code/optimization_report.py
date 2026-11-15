#!/usr/bin/env python3
"""
ACC 19: Generate Optimization Report

Consolidate results from quantization, FSDP, and torch.compile benchmarks
into publication-ready figures and tables.

Usage:
    python optimization_report.py

Author: Xinru Qiu
Date: February 2026
"""

import logging
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
OPT_DIR = BASE_DIR / "05_RESULTS/19_model_optimization"
FIGURES_DIR = OPT_DIR / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def plot_quantization_results():
    """Plot precision vs accuracy vs latency vs memory."""
    quant_file = OPT_DIR / "quantization/quantization_results.csv"
    if not quant_file.exists():
        logger.warning("Quantization results not found, skipping plot")
        return

    df = pd.read_csv(quant_file)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1. Accuracy vs Precision
    ax = axes[0]
    colors = {'fp32': '#2196F3', 'fp16': '#4CAF50', 'int8_dynamic': '#FF9800', 'int8_static': '#F44336'}
    for _, row in df.iterrows():
        ax.bar(row['precision'], row['auc'], color=colors.get(row['precision'], '#999'),
              edgecolor='black', linewidth=0.5)
        ax.text(row['precision'], row['auc'] + 0.002, f"{row['auc']:.4f}",
               ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('AUC')
    ax.set_title('Classification Accuracy')
    ax.set_ylim(0.9, 1.0)

    # 2. Speedup vs Precision
    ax = axes[1]
    for _, row in df.iterrows():
        ax.bar(row['precision'], row['speedup_vs_fp32'], color=colors.get(row['precision'], '#999'),
              edgecolor='black', linewidth=0.5)
        ax.text(row['precision'], row['speedup_vs_fp32'] + 0.05,
               f"{row['speedup_vs_fp32']:.1f}x", ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('Speedup vs FP32')
    ax.set_title('Inference Speedup')
    ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)

    # 3. Memory vs Precision
    ax = axes[2]
    for _, row in df.iterrows():
        ax.bar(row['precision'], row['model_size_mb'], color=colors.get(row['precision'], '#999'),
              edgecolor='black', linewidth=0.5)
        ax.text(row['precision'], row['model_size_mb'] + 10,
               f"{row['model_size_mb']:.0f}MB", ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('Model Size (MB)')
    ax.set_title('Memory Footprint')

    plt.suptitle('SE-600M Quantization Benchmark', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'quantization_benchmark.png', dpi=150, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'quantization_benchmark.pdf', bbox_inches='tight')
    plt.close()
    logger.info("Saved quantization figure")


def plot_fsdp_results():
    """Plot FSDP scaling curve."""
    fsdp_files = list((OPT_DIR / "fsdp").glob("fsdp_results_*.json"))
    if not fsdp_files:
        logger.warning("FSDP results not found, skipping plot")
        return

    results = []
    for f in fsdp_files:
        with open(f) as fp:
            results.append(json.load(fp))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Extract data
    n_gpus = [r['n_gpus'] for r in results]
    throughputs = [r['final_throughput'] for r in results]
    memories = [r['final_peak_memory_gb'] for r in results]
    aucs = [r['best_auc'] for r in results]

    # 1. Throughput scaling
    ax = axes[0]
    ax.plot(n_gpus, throughputs, 'o-', color='#2196F3', linewidth=2, markersize=8)
    ax.set_xlabel('Number of GPUs')
    ax.set_ylabel('Throughput (cells/sec)')
    ax.set_title('FSDP Throughput Scaling')

    # 2. Memory per GPU
    ax = axes[1]
    ax.plot(n_gpus, memories, 'o-', color='#4CAF50', linewidth=2, markersize=8)
    ax.set_xlabel('Number of GPUs')
    ax.set_ylabel('Peak Memory (GB)')
    ax.set_title('Memory per GPU')

    # 3. AUC (should be stable)
    ax = axes[2]
    ax.plot(n_gpus, aucs, 'o-', color='#FF9800', linewidth=2, markersize=8)
    ax.set_xlabel('Number of GPUs')
    ax.set_ylabel('Validation AUC')
    ax.set_title('Model Quality')
    ax.set_ylim(0.9, 1.0)

    plt.suptitle('FSDP Distributed Training Scaling', fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'fsdp_scaling.png', dpi=150, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'fsdp_scaling.pdf', bbox_inches='tight')
    plt.close()
    logger.info("Saved FSDP scaling figure")


def plot_compile_results():
    """Plot torch.compile speedup chart."""
    compile_file = OPT_DIR / "torch_compile/torch_compile_results.csv"
    if not compile_file.exists():
        logger.warning("torch.compile results not found, skipping plot")
        return

    df = pd.read_csv(compile_file)
    if 'error' in df.columns:
        df = df[df['error'].isna() | (df['error'] == '')]

    models = df['model'].unique()
    modes = [m for m in df['mode'].unique() if m != 'eager']

    fig, ax = plt.subplots(figsize=(12, 6))

    x = np.arange(len(models))
    width = 0.35 / max(1, len(modes))

    # Eager baseline bars
    eager = df[df['mode'] == 'eager'].set_index('model')
    ax.bar(x - 0.2, eager.loc[models, 'mean_latency_ms'].values,
          width=0.35, label='Eager', color='#BDBDBD', edgecolor='black', linewidth=0.5)

    # Compiled bars
    colors = ['#2196F3', '#4CAF50', '#FF9800']
    for i, mode in enumerate(modes):
        compiled = df[df['mode'] == mode].set_index('model')
        if len(compiled) == 0:
            continue
        vals = []
        for m in models:
            if m in compiled.index:
                vals.append(compiled.loc[m, 'mean_latency_ms'])
            else:
                vals.append(0)

        offset = 0.2 + i * width
        ax.bar(x + offset, vals, width=width, label=f'torch.compile({mode})',
              color=colors[i % len(colors)], edgecolor='black', linewidth=0.5)

    ax.set_xlabel('Foundation Model')
    ax.set_ylabel('Latency (ms/batch)')
    ax.set_title('torch.compile Inference Speedup')
    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=30, ha='right')
    ax.legend()

    plt.tight_layout()
    plt.savefig(FIGURES_DIR / 'torch_compile_speedup.png', dpi=150, bbox_inches='tight')
    plt.savefig(FIGURES_DIR / 'torch_compile_speedup.pdf', bbox_inches='tight')
    plt.close()
    logger.info("Saved torch.compile figure")


def generate_report():
    """Generate consolidated optimization report."""
    logger.info("=" * 60)
    logger.info("ACC 19: Optimization Report")
    logger.info("=" * 60)

    plot_quantization_results()
    plot_fsdp_results()
    plot_compile_results()

    logger.info(f"\nAll figures saved to {FIGURES_DIR}")


if __name__ == '__main__':
    generate_report()
