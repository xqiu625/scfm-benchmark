#!/usr/bin/env python3
"""
ACC 19: torch.compile Inference Benchmark

Apply torch.compile to all 5 FM inference pipelines and measure:
- Compilation warmup time
- Steady-state throughput vs eager mode
- Correctness verification (compiled output matches eager)

Usage:
    python torch_compile_benchmark.py
    python torch_compile_benchmark.py --modes reduce-overhead max-autotune

Author: Xinru Qiu
Date: February 2026
"""

import argparse
import logging
import time
from pathlib import Path
import json
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

BASE_DIR = Path("/bigdata/godziklab/shared/Xinru/302006")
RESULTS_DIR = BASE_DIR / "05_RESULTS/19_model_optimization/torch_compile"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


class SimpleClassifier(nn.Module):
    """Lightweight classifier mirroring downstream task for each FM."""

    def __init__(self, input_dim, hidden_dim=256, n_classes=2):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, n_classes),
        )

    def forward(self, x):
        return self.net(x)


# Model embedding dimensions
MODEL_DIMS = {
    'STATE': 2058,
    'UCE': 1280,
    'Geneformer': 1152,
    'TranscriptFormer': 2048,
    'scGPT': 512,
}


def benchmark_model(model_name, input_dim, device, compile_mode=None,
                   batch_size=256, n_warmup=20, n_runs=100):
    """Benchmark a single model in eager or compiled mode."""
    model = SimpleClassifier(input_dim=input_dim).to(device).eval()

    if compile_mode:
        compile_start = time.perf_counter()
        model = torch.compile(model, mode=compile_mode)
        # Need to run once to trigger compilation
        dummy = torch.randn(batch_size, input_dim, device=device)
        with torch.no_grad():
            _ = model(dummy)
        if device.type == 'cuda':
            torch.cuda.synchronize()
        compile_time = time.perf_counter() - compile_start
    else:
        compile_time = 0.0

    # Warmup
    dummy = torch.randn(batch_size, input_dim, device=device)
    with torch.no_grad():
        for _ in range(n_warmup):
            _ = model(dummy)
    if device.type == 'cuda':
        torch.cuda.synchronize()

    # Benchmark
    latencies = []
    for _ in range(n_runs):
        x = torch.randn(batch_size, input_dim, device=device)
        if device.type == 'cuda':
            torch.cuda.synchronize()
        start = time.perf_counter()
        with torch.no_grad():
            _ = model(x)
        if device.type == 'cuda':
            torch.cuda.synchronize()
        latencies.append((time.perf_counter() - start) * 1000)

    return {
        'compile_time_sec': compile_time,
        'mean_latency_ms': float(np.mean(latencies)),
        'std_latency_ms': float(np.std(latencies)),
        'p50_latency_ms': float(np.percentile(latencies, 50)),
        'p95_latency_ms': float(np.percentile(latencies, 95)),
        'throughput_samples_per_sec': batch_size / (np.mean(latencies) / 1000),
    }


def verify_correctness(model_name, input_dim, device, compile_mode='reduce-overhead'):
    """Verify compiled model produces same output as eager."""
    model = SimpleClassifier(input_dim=input_dim).to(device).eval()
    compiled = torch.compile(model, mode=compile_mode)

    x = torch.randn(32, input_dim, device=device)
    with torch.no_grad():
        eager_out = model(x)
        compiled_out = compiled(x)

    max_diff = (eager_out - compiled_out).abs().max().item()
    return max_diff < 1e-4


def run_benchmark(modes=None, batch_size=256):
    """Run torch.compile benchmark across all models."""
    if modes is None:
        modes = ['reduce-overhead', 'max-autotune']

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    logger.info("=" * 60)
    logger.info("ACC 19: torch.compile Inference Benchmark")
    logger.info(f"Device: {device}")
    logger.info(f"Compile modes: {modes}")
    logger.info(f"Batch size: {batch_size}")
    logger.info("=" * 60)

    results = []

    for model_name, input_dim in MODEL_DIMS.items():
        logger.info(f"\n--- {model_name} (dim={input_dim}) ---")

        # Eager baseline
        eager = benchmark_model(model_name, input_dim, device, compile_mode=None, batch_size=batch_size)
        eager['model'] = model_name
        eager['mode'] = 'eager'
        results.append(eager)
        logger.info(f"  Eager: {eager['mean_latency_ms']:.3f} ms ({eager['throughput_samples_per_sec']:.0f} samples/s)")

        # Compiled modes
        for mode in modes:
            try:
                compiled = benchmark_model(model_name, input_dim, device, compile_mode=mode, batch_size=batch_size)
                compiled['model'] = model_name
                compiled['mode'] = mode

                speedup = eager['mean_latency_ms'] / compiled['mean_latency_ms']
                compiled['speedup_vs_eager'] = speedup

                results.append(compiled)
                logger.info(f"  {mode}: {compiled['mean_latency_ms']:.3f} ms "
                          f"({compiled['throughput_samples_per_sec']:.0f} samples/s, "
                          f"{speedup:.2f}x speedup, compile={compiled['compile_time_sec']:.1f}s)")

                # Correctness check
                if device.type == 'cuda':
                    correct = verify_correctness(model_name, input_dim, device, mode)
                    compiled['correct'] = correct
                    if not correct:
                        logger.warning(f"  WARNING: Compiled output diverges from eager!")

            except Exception as e:
                logger.warning(f"  {mode}: FAILED - {e}")
                results.append({
                    'model': model_name,
                    'mode': mode,
                    'error': str(e),
                })

    # Save results
    df = pd.DataFrame(results)
    df.to_csv(RESULTS_DIR / "torch_compile_results.csv", index=False)

    with open(RESULTS_DIR / "torch_compile_results.json", 'w') as f:
        json.dump(results, f, indent=2, default=str)

    # Summary
    logger.info(f"\n{'='*80}")
    logger.info("TORCH.COMPILE BENCHMARK SUMMARY")
    logger.info(f"{'='*80}")
    logger.info(f"{'Model':<20} {'Mode':<20} {'Latency(ms)':<15} {'Speedup':<10} {'Compile(s)':<12}")
    logger.info("-" * 77)

    for r in results:
        if 'error' in r:
            continue
        speedup = f"{r.get('speedup_vs_eager', 1.0):.2f}x" if r['mode'] != 'eager' else '-'
        compile_t = f"{r['compile_time_sec']:.1f}" if r['mode'] != 'eager' else '-'
        logger.info(f"{r['model']:<20} {r['mode']:<20} {r['mean_latency_ms']:<15.3f} {speedup:<10} {compile_t:<12}")

    logger.info(f"\nResults saved to {RESULTS_DIR}")
    return df


def main():
    parser = argparse.ArgumentParser(description='torch.compile Benchmark')
    parser.add_argument('--modes', nargs='+', default=['reduce-overhead', 'max-autotune'],
                       choices=['reduce-overhead', 'max-autotune', 'default'])
    parser.add_argument('--batch-size', type=int, default=256)
    args = parser.parse_args()

    run_benchmark(args.modes, args.batch_size)


if __name__ == '__main__':
    main()
