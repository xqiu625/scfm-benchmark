# Runbook — HPCC (UCR)

Extracted 2026-09-21 from `docs/memory-bank/techcontext.md` (archived under
`archive/superseded_2026-09-21/`). Only entries verified in use as of that date are carried over;
the archived file additionally contains retired models, old repo paths and superseded commands.

---

## Access

```bash
ssh -x xinruqiu@cluster.hpcc.ucr.edu     # works non-interactively, no Duo prompt
```

Hostname is `cluster.hpcc.ucr.edu` — not `skylark`, not `bluejay`.

`/tmp` is **not** shared between login and compute nodes. Stage anything a job needs under
`/bigdata`, never `/tmp`.

---

## Paths

| What | Where |
|---|---|
| Project root | `/bigdata/godziklab/shared/Xinru/302006` |
| Repo mirror | `…/302006/302006_updategit` |
| Raw data | `…/302006/01_DATA/raw` |
| Results | `…/302006/05_RESULTS` |
| Logs | `…/302006/06_LOGS` |
| Shared R library | `/bigdata/godziklab/shared/Xinru/R` |

The repo mirror is a *mirror*. Code written there and not synced back to the local checkout is
invisible to git — this has already caused one full reproducibility failure (24 untracked
scripts, fixed in `f0c1070`). Write locally, `scp` out; never the reverse.

---

## Conda environments

| Env | Use |
|---|---|
| `302006_env_clean` | scanpy / scrublet / sklearn — QC, Gate 1, Gate 2 |
| `scgpt_finetune` | Geneformer |
| `scpeft_env` | scGPT |
| `state_env` | STATE, UCE |

Activate without `conda activate` inside batch scripts:

```bash
ENV_DIR=/rhome/xinruqiu/bigdata/.conda/envs/<env>
export PATH="$ENV_DIR/bin:$PATH"
export PYTHONNOUSERSITE=1
```

`PYTHONNOUSERSITE=1` is required; without it user site-packages leak in and shadow the env.

---

## R

```bash
module load gcc/16.2.0        # 12.2.0 is NOT enough — celda needs CXXABI_1.3.15
module load R/4.5.2
export R_LIBS_USER=/bigdata/godziklab/shared/Xinru/R
```

Do not modify the shared `cellrank_env`.

---

## SLURM

| Partition | Walltime cap | Use |
|---|---|---|
| `short` | 2h | default for CPU work |
| `short_gpu` | 2h | default for GPU work |
| `gpu` + `--qos=gpu` | longer | only when a GPU job genuinely cannot be chunked |
| `highmem` | — | minimum 100G memory request |

**Default to `short`/`short_gpu` with checkpointed array jobs.** Long unchunked jobs on
`batch`/`highmem` have repeatedly been lost to walltime — a Chongqing DecontX run died at 8h
mid-write and had to be redone at 20h.

A checkpoint must bind to its configuration, not just its channel name, or a resumed run will
silently reuse results computed under different settings. See `runbooks/gate1-gate2.md`.

GPU exclusions: `--exclude=gpu01,gpu02,gpu03,gpu04,gpu13,gpu14`.

---

## Gotchas found the hard way

- **R / Matrix 32-bit ceiling.** Large sparse matrices overflow; chunk or stream instead.
- **Double-gzipped RDS** (Wilk): `zcat file | zcat | ...` before `readRDS`.
- **`np.char.upper` needs true `str_` dtype**, not object — `np.asarray(x, dtype=np.str_)`.
- **NaN outcome silently casts to int.** Filter with `event.notna()` explicitly. This is an
  endpoint-availability filter and must stay separate from outcome-blind QC.
- **Write checkpoint outputs unconditionally.** A `genes.csv` omitted from a success branch
  produced a checkpoint that resumed into a missing-file error.
- **Write to a temp name then `os.rename`** so an interrupted write never leaves a half-file
  that a resume accepts as complete.
