# Deprecated Scripts Archive

Archived on 2026-04-10 during code review cleanup.

## duplicates/
Exact copies that existed in both `analyses/` and `scripts/`. The `analyses/` versions are kept as canonical.

| Archived File | Kept Version |
|---------------|-------------|
| `generate_geneformer_embeddings.py` | `analyses/16_embedding_generation_5models/code/` |
| `generate_transcriptformer_embeddings.py` | `analyses/16_embedding_generation_5models/code/` |
| `generate_perturbation_embeddings.py` | `analyses/16_embedding_generation_5models/code/generate_state_embeddings.py` |
| `benchmark_single_model_02.py` | `analyses/03_zeroshot_benchmark/code/benchmark_single_model.py` (has scGPT_BP) |

## old_versions/
Superseded by optimized `_fast` or `_extended` replacements.

| Archived File | Replacement | Speedup |
|---------------|------------|---------|
| `statistical_significance.py` | `statistical_significance_fast.py` | 1,440x |
| `embedding_ablations.py` | `embedding_ablations_fast.py` | 60x |
| `run_stats_{state,uce,scgpt,geneformer}.sh` | `run_stats_*_fast.sh` | N/A |
| `clinical_metrics.py` | `clinical_metrics_extended.py` | N/A |

## debug_fix/
One-time diagnostic, fix, setup, and verification scripts. Purpose already served.

- TranscriptFormer gene ID fixes (`fix_*.py`, `diagnose_*.py`, `debug_*.py`, `final_fix_entrez_ids.sh`)
- Environment verification (`verify_scpeft_setup.py`, `verify_true_finetune_setup.py`)
- Data cleaning (`clean_scgpt_embeddings.py`)
- Metadata inspection (`inspect_metadata.py`)
- Thin wrapper (`scpeft_finetune.py`)

## archived_from_filename/
Files previously marked with `(archive)` in their filenames, now properly relocated and renamed.

- UCE cell-type extraction scripts (`extract_hct116_*.py`, `extract_hek293t_*.py`)
- scGPT batch processing (`scgpt_all_platelets.py`, `scgpt_singularity_batch.sh`)
