#!/usr/bin/env python3
"""
scPEFT Fine-Tuning Wrapper

This is a wrapper script that calls the backbone-specific PEFT fine-tuning scripts.

For actual implementation, see:
    - scgpt_peft_finetune.py: scGPT with PEFT adapters (LoRA, PREFIX, ENCODER, TOKEN)
    - geneformer_peft_finetune.py: Geneformer with LoRA adapters

Based on scPEFT methodology (Nature Machine Intelligence, 2025):
    Paper: https://doi.org/10.1038/s42256-025-01170-z
    Code: https://github.com/coffee19850519/scPEFT

Usage:
    python scpeft_finetune.py --backbone scgpt --task binary
    python scpeft_finetune.py --backbone geneformer --task binary

Or use backbone-specific scripts directly:
    python scgpt_peft_finetune.py --task binary --peft LORA
    python geneformer_peft_finetune.py --task binary --lora-r 16

Author: Xinru Qiu
Date: January 2026
"""

import subprocess
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(
        description='scPEFT Fine-Tuning Wrapper',
        epilog='For more options, use the backbone-specific scripts directly.'
    )
    parser.add_argument('--backbone', required=True, choices=['scgpt', 'geneformer'],
                        help='Backbone model')
    parser.add_argument('--task', default='binary', choices=['binary', '3class', '6class'])
    parser.add_argument('--epochs', type=int, default=20)
    parser.add_argument('--all-tasks', action='store_true')

    args = parser.parse_args()

    if args.backbone == 'scgpt':
        script = 'scgpt_peft_finetune.py'
        cmd = [sys.executable, script, '--task', args.task, '--epochs', str(args.epochs)]
        if args.all_tasks:
            cmd.append('--all-tasks')
    else:
        script = 'geneformer_peft_finetune.py'
        cmd = [sys.executable, script, '--task', args.task, '--epochs', str(args.epochs)]
        if args.all_tasks:
            cmd.append('--all-tasks')

    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd)


if __name__ == '__main__':
    main()
