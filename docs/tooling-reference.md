# Tooling Reference

Non-analysis code: foundation-model embedding generation, fine-tuning,
experiment tracking, and HPCC utilities.

Analysis code lives in `../scripts/` — see `../scripts/README.md` for
the per-ACC registry and status markers.

## Path resolution — `scfm_benchmark.config`

`src/scfm_benchmark/config.py` is the single source of truth for the base path.
Resolve paths through it instead of hardcoding `/bigdata/...`:

```python
from scfm_benchmark.config import paths, settings

paths.project_root          # env-resolved root (hpcc | local | cloud)
paths.data_root             # …/data
paths.resolve("results", "primary_benchmark.csv")
settings["random_seed"]     # merged configs/defaults.yaml + <env>.yaml
```

Environment is auto-detected (`/bigdata` present → `hpcc`, else `local`) or forced
with `SCFM_ENV=hpcc|local|cloud`; `PROJECT_ROOT` / `SCFM_DATA_ROOT` / `SCFM_LOGS_DIR`
override individual roots. Inspect with `python -m scfm_benchmark.config`. The
numbered stage scripts still carry hardcoded absolutes as a fallback and are migrated
to `paths` opportunistically (see `../docs/decisions.md`, 2026-06-06).

## Structure

```
scripts/
├── agent/              # Optional research orchestration prototype
├── finetuning/         # Parameter-efficient fine-tuning (Geneformer LoRA, scGPT PEFT, STATE LoRA)
├── foundation_models/  # Per-FM embedding-generation SLURM submitters (one folder per FM)
└── tracking/           # Weights & Biases integration (wandb_tracker + backfill)
```

## Subfolder details

| Path | Purpose |
|---|---|
| `agent/scfm_benchmark_agent.py` | Optional research orchestration prototype; separate from the package CLI in `src/scfm_benchmark/cli.py`. |
| `finetuning/geneformer_peft_finetune.py` | Geneformer LoRA fine-tune driver (the May 6 chain that populated `tab:domainadapt`). |
| `finetuning/state_lora_finetune.py` | STATE LoRA fine-tuning (auxiliary; not in v1 manuscript primary). |
| `finetuning/scgpt_peft_finetune.py` | scGPT PEFT (scGPT was dropped from the v1 active set on 2026-05-01; retained for reference). |
| `foundation_models/<FM>/` | Per-FM SLURM submitter for embedding generation against the canonical Option C input. One folder per FM: `Geneformer/`, `scGPT/`, `scPRINT/`, `STACK/`. |
| `tracking/wandb_tracker.py` | Logs benchmark results to Weights & Biases for cross-run comparison. |
| `tracking/backfill_experiments.py` | One-time backfill of pre-existing benchmark CSVs into W&B. |

## Top-level utility files

| File | Purpose |
|---|---|
| `download_results_from_hpcc.sh` | Helper for pulling result files from HPCC to local. Largely superseded by the per-ACC `rsync` workflow used in the 2026-05-24 reorganization (`../docs/submission-checklist.md`). |
| `generate_acc17_21_figures.py` | Figure-generation helper for ACCs 17–22 (domain adaptation, multimodal fusion, model optimization, drug-target validation, conditional generation, ESM multimodal). |
| `inspect_cell_indices.py` | Debugging helper to inspect cell-index alignment across FM embeddings. |
