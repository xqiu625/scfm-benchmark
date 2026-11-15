# Environment Manifests

Conda / venv environment definitions for the 302006 platelet-OOD scFM benchmark.

## Core analysis environment

- `scfm_benchmark.yml` — the downstream analysis/reporting environment
  (embedding classification, evaluation, figures, `scfm-bench` CLI). This is a
  copy of the repository-root `environment.yml`; a Dockerfile for the same core
  environment is provided at the repository root.

Every result, table, and figure in the manuscript is reproducible with this
core environment plus the pre-computed embeddings. The per-model environments
below are needed only to regenerate the embeddings themselves.

## Per-foundation-model HPCC environments

| Manifest | Environment | Type | Python | Used for |
|---|---|---|---|---|
| `state_env.yml` | `state_env` | conda (prefix) | 3.11 | STATE and UCE embedding inference |
| `scprint2_env-requirements.txt` | `scprint2_env` | uv venv | 3.11 | scPRINT-2 embedding inference (`scprint2==1.0.3`) |
| `scpeft_env.yml` | `scpeft_env` | conda (prefix) | 3.9 | scGPT inference; Geneformer zero-shot extraction (V2-316M) and LoRA fine-tuning (V1-10M) |
| `scgpt_finetune.yml` | `scgpt_finetune` | conda (named) | 3.10 | legacy scGPT fine-tuning arm |

These environments are large GPU stacks maintained on the UCR HPCC (see
`configs/hpcc.yaml` and `docs/tooling-reference.md`). They are **not** uniformly
conda environments and are **not** all registered by name, so a single
`conda env export --name <env>` loop does not work. Regenerate them with:

```bash
# from the repo root on the HPCC
S=/bigdata/godziklab/shared/Xinru/302006

# conda environments referenced by prefix
conda env export --prefix $S/state_env  --no-builds > environments/state_env.yml
conda env export --prefix $S/scpeft_env --no-builds > environments/scpeft_env.yml

# conda environment registered by name
conda env export --prefix /rhome/xinruqiu/bigdata/.conda/envs/scgpt_finetune \
    --no-builds > environments/scgpt_finetune.yml

# scprint2_env is a uv venv, not conda, and ships no pip module
uv pip freeze --python $S/scprint2_env/bin/python \
    > environments/scprint2_env-requirements.txt
```

Exported 2026-08-05.
