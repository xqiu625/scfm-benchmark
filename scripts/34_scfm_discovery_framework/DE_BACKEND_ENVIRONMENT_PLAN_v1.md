# DE statistical backend: environment plan v1 (revised 2026-10-14; not installed this round)

Prepared per `review:20261014_151747` step 3, revised per `review:20261014_203915`
F3: v1 pointed its only install-command draft at the EXISTING shared `302006_env_clean` and quoted
`lfc_shrink(coeff=None, adapt=True)` as though `coeff=None` were a usable default, neither of which
survives inspection. This revision gives a concrete NEW isolated environment path (never
`302006_env_clean` as an install target), a REAL (executed, read-only) dependency-resolution
dry-run for two candidate Python/`pydeseq2` pairings, and a corrected backend-call draft using
`pydeseq2==0.5.2`'s/`0.5.4`'s actual source, not its stale docstring text. This remains a **plan
only** — nothing below was installed, upgraded, or run against real data this round.
`ExpressionDERankingProvider` (`gene_ranking_providers.py`) still has no real backend after this
document; it still always returns `status in ('failed', 'not_run')`, never a fabricated
`'complete'`. Facts and recommendations are kept separate per the convention this framework
already uses in `MODEL_SELECTION_v1.md`.

## What this round actually did to produce the facts below

- Read the live `302006_env_clean` conda environment directly (`pip list` over SSH, on
  `cluster.hpcc.ucr.edu`, 2026-10-14) — not from a requirements file or a prior claim.
- Queried PyPI's own JSON metadata API directly for `pydeseq2` (`https://pypi.org/pypi/pydeseq2/json`
  and per-version `.../pydeseq2/<version>/json`, fetched 2026-10-14) for `requires_python` and
  `requires_dist` — not from a cached search snippet.
- Checked R's actual availability on the cluster via an interactive login shell's module system
  (`bash -l -c 'module avail r'`), **not** from a non-interactive SSH shell's `which Rscript`
  returning nothing — per this review's own explicit caution ("不能从一个非交互 shell 找不到
  Rscript 推断整个集群没有 R") and `DEVELOPMENT.md` rule 3 (absence of a guarantee is not a guarantee
  of the opposite). The non-interactive check genuinely found no `Rscript` on `PATH` and no
  `module` command at all; the login-shell check found a real Lmod module tree with R present.
  This revision did not re-run that check (per the review's own note that it did not re-verify
  the module list either) — it is restated from the prior round's own verified command/output,
  cited exactly as before, not re-claimed as freshly re-checked.
- Fetched `pydeseq2`'s actual source at tag `v0.5.2` directly from GitHub (`raw.githubusercontent.com
  /owkin/PyDESeq2/v0.5.2/pydeseq2/{ds.py,dds.py}`, fetched 2026-10-14) to read `DeseqStats.
  lfc_shrink`'s and `DeseqDataSet.__init__`'s REAL signatures, rather than trusting the
  ReadTheDocs-rendered docstring text (which, per below, is stale/misleading on one parameter).
- Ran a REAL (executed, read-only) `pip install --dry-run --report <file>.json` dependency
  resolution from `302006_env_clean`'s own `pip` (25.1.1), targeting `--python-version 3.11
  --implementation cp --abi cp311 --platform manylinux2014_x86_64 --only-binary=:all:` for
  `pydeseq2==0.5.4` and the equivalent `--python-version 3.10 ... --abi cp310` for
  `pydeseq2==0.5.2`, both on `cluster.hpcc.ucr.edu`, 2026-10-14 (see "Dependency-resolution
  dry-run" below for the important caveat on what this does and does not prove).
- Attempted a REAL `conda create --dry-run -c conda-forge python=3.11 pip` check (also read-only —
  `--dry-run` creates nothing) to verify a Python 3.11 conda-forge environment is resolvable on
  this cluster; **this process was killed (OOM) on the login node** before completing metadata
  resolution, consistent with this project's existing caution against running memory-heavy
  operations on a login node. This is reported as attempted-and-inconclusive, not as evidence
  either way — see "New isolated environment" below.

## Fact: R is available on this cluster via modules; not yet bridged into this project's env

Confirmed 2026-10-14, login shell, `module avail r`:

```
R/4.1.2   R/4.1.3   R/4.2.0   R/4.2.2   R/4.3.0   R/4.4.0   R/4.4.2   R/4.5.0
R/4.5.2 (default) <L>   R/4.6.0
```

So an R/DESeq2/edgeR/limma path is genuinely possible here (`module load R/4.5.2` or similar),
contrary to what a naive non-interactive-shell check would have wrongly concluded. This is a real
option, not a dead end, and it does NOT require `rpy2` specifically: per
`review:20261014_203915` F3, bridging to R can equally be done via a separate `Rscript`
process invoked as a subprocess, with the design matrix/counts/contrast passed through plain
files (CSV/TSV/RDS written by Python, read by the R script, results written back to a file Python
reads) — never importing `rpy2`'s in-process R<->Python object conversion at all. Either bridge
still requires a SEPARATE R package install (`DESeq2`/`edgeR`/`apeglm` via Bioconductor, not yet
attempted, not yet verified for this R version) and crosses a second-language failure surface
(permissions on a shared module tree, R session/subprocess lifecycle). None of that has been
attempted or verified this round; the Python-native PyDESeq2 path below avoids this boundary
entirely and remains the primary recommendation.

## Fact: this project's own environment is 100% Python; the DE-stream interface is already Python

Every file this B increment touches (`pseudobulk_providers.py`, `gene_ranking_providers.py`,
`de_design_providers.py`, `artifact_schema.py`) is plain Python with no R dependency anywhere.
`302006_env_clean` already has (confirmed via live `pip list`, 2026-10-14) — shown here only as
the baseline that will remain UNCHANGED (the DE backend installs into the new isolated
`302006_de_backend_v1` below, never here):

| package | installed version |
|---|---|
| python | 3.10.18 |
| numpy | 2.2.6 |
| pandas | 2.3.1 |
| scipy | 1.15.3 |
| scikit-learn | 1.7.1 |
| matplotlib | 3.10.3 |
| anndata | 0.11.4 |
| scanpy | 1.11.3 |
| statsmodels | 0.14.5 |
| pydeseq2 | **not installed** |
| rpy2 | **not installed** |
| formulaic / formulaic-contrasts | **not installed** |

## New isolated environment (NEVER `302006_env_clean` as an install target)

Per `review:20261014_203915` F3: v1's only install-command draft pointed at the SHARED
`302006_env_clean`, which this project uses across multiple stages — installing a DE backend
there risks every other already-pinned package in that env, and conflates "the version PyDESeq2
happens to still support" with "the version this project is stuck with". A NEW, ISOLATED conda
environment removes that constraint entirely: it can have its OWN Python version, chosen for the
DE backend alone, never the shared env's.

- **Proposed absolute path**: `/rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1` — a
  sibling of `302006_env_clean` under the same `.conda/envs` directory (so it resolves via the
  same `conda env list`/`envs_dirs` convention), but a DISTINCT environment `302006_env_clean`
  is never touched by any command in this document.
- **Draft creation command (NOT RUN this round — environment creation, even empty, is an
  authorized-install-equivalent action per `DEVELOPMENT.md`)**:
  ```bash
  conda create -n 302006_de_backend_v1 -c conda-forge python=3.11 pip
  ```
- An actual (read-only) `conda create --dry-run -c conda-forge python=3.11 pip` check was
  attempted on `cluster.hpcc.ucr.edu`'s login node 2026-10-14 to verify this resolves at all; it
  was **killed (OOM) before completing conda-forge metadata resolution** — this project's
  existing caution against memory-heavy operations on a login node (`DEVELOPMENT.md`'s "do not run
  the real 5 GB input on the login node" note) generalises to conda's own metadata solve, which
  is a genuinely memory-hungry operation. **This is reported as attempted-and-inconclusive**:
  it neither confirms nor rules out that a Python 3.11 conda-forge environment resolves cleanly
  here. The real check (if authorized) should run from an interactive compute allocation, not
  the login node, or via `mamba`/`conda`'s `--repodata-fn current_repodata.json` to reduce the
  metadata it has to load.

## Recommendation: `pydeseq2==0.5.4` on Python 3.11 in the NEW isolated env (primary); `pydeseq2==0.5.2` on Python 3.10 (documented fallback)

**Why the isolated env changes the version choice from v1.** v1 recommended pinning to `0.5.2`
specifically because `302006_env_clean` runs Python 3.10.18 and `pydeseq2>=0.5.3` requires
`>=3.11`. In a NEW, ISOLATED environment that constraint does not exist — there is no reason to
deliberately install an older, no-longer-actively-supported `pydeseq2` release merely because a
DIFFERENT, shared environment happens to run an older Python. The isolated env should run Python
**3.11** and install the current PyPI release, **`pydeseq2==0.5.4`**.

| version | `requires_python` | upload date (PyPI) | `anndata`/`numpy`/`pandas`/`scikit-learn`/`scipy`/`matplotlib` lower bounds |
|---|---|---|---|
| 0.5.2 (fallback) | `>=3.10.0` | 2025-06-26 | `>=0.8.0` / `>=1.23.0` / `>=1.4.0` / `>=1.1.0` / `>=1.11.0` / `>=3.6.2` |
| 0.5.4 (primary) | `>=3.11` | 2026-01-23 | `>=0.11.0` / `>=2.0.0` / `>=2.2.0` / `>=1.4.0` / `>=1.12.0` / `>=3.9.0` |

(`0.5.4`'s metadata fetched directly from `https://pypi.org/pypi/pydeseq2/0.5.4/json`, 2026-10-14;
`0.5.2`'s from the equivalent per-version endpoint, same date, repeated from v1.) Both versions
declare `formulaic>=1.0.2` and `formulaic-contrasts>=0.2.0` as dependencies (introduced at
`0.5.0`); neither is installed anywhere on this cluster today.

**Fallback condition, stated explicitly**: if the isolated env's Python 3.11 provisioning is
later found NOT to resolve cleanly (the conda dry-run above was inconclusive, not negative), the
fallback is Python 3.10 + `pydeseq2==0.5.2` in the SAME new isolated env (never in
`302006_env_clean`) — not a reason to abandon the isolated-env approach itself.

**Why PyDESeq2 over the R path at all, given R is confirmed available.** PyDESeq2 implements
DESeq2's own negative-binomial GLM + dispersion-shrinkage model natively in Python (it is a
from-scratch reimplementation of DESeq2's statistical method, not a wrapper around it), so it
reaches the SAME estimator family the module docstring already commits to ("donor/animal-aware
pseudobulk aggregation + a shrunken log-fold-change-style estimator") without crossing a second
language boundary and without needing a SEPARATE Bioconductor install on the shared module R
tree (permissions/versioning unverified). The R path (via `Rscript` + a file interface, or via
`rpy2` — neither is mandatory, see above) remains a DOCUMENTED, available alternative if a future
round needs apeglm's own reference R implementation specifically, or a method PyDESeq2 does not
cover — but it is not this plan's primary recommendation.

## Dependency-resolution dry-run (REAL, executed, read-only — NOT an install)

**Revised per `review:20261014_221155` R4**: the first dry-run (below, kept for its own
transitive-chain evidence) reused `302006_env_clean`'s ALREADY-INSTALLED packages as "satisfied,"
understating a fresh environment's footprint. Per the review's own pointer to `pip install
--dry-run`'s official docs (`pip.pypa.io/en/stable/cli/pip_install`), adding **`--ignore-installed`**
(`-I`) forces pip to resolve the FULL dependency closure from scratch, ignoring what happens to
already be installed — still read-only, still no actual install, confirmed via `pip install
--help`'s own text ("Can be used in combination with `--ignore-installed` to 'resolve' the
requirements") before running it.

```bash
# Full closure (REAL, executed, 2026-10-14):
/rhome/xinruqiu/bigdata/.conda/envs/302006_env_clean/bin/python -m pip install --dry-run \
    --ignore-installed --report pip_dry_run_report_py311_full.json \
    --python-version 3.11 --implementation cp --abi cp311 \
    --platform manylinux2014_x86_64 --only-binary=:all: \
    "pydeseq2==0.5.4" "formulaic>=1.0.2" "formulaic-contrasts>=0.2.0"
```

**Result**: exit 0, no resolver conflict. **40 packages** in the full closure (`pip_dry_run_
report_py311_full.json`'s own `install` list), including every direct dependency PyPI's metadata
already named PLUS their own transitive chain this round newly surfaced: `anndata==0.12.19` (a
newer release than the `>=0.11.0` floor), `numpy==2.2.6`, `pandas==2.3.2`, `scipy==1.16.3`,
`scikit-learn==1.7.2`, `matplotlib==3.11.2`, `h5py==3.14.0`, plus `anndata`'s own real transitive
set this round newly surfaced: `zarr`, `numcodecs`, `donfig`, `array-api-compat`,
`legacy-api-wrap`, `scverse-misc`, `session-info2`, `google-crc32c`, `PyYAML`. This is now the
AUTHORITATIVE dependency-pin evidence for `302006_de_backend_v1`, superseding the earlier
7-package list below.

**Residual caveat, stated plainly (not resolved, per `DEVELOPMENT.md` rule 3)**: the review's own
pointer flags that a cross-interpreter dry-run still carries a boundary — `--python-version
3.11`/`--abi cp311`/`--platform` steer which WHEEL FILES pip considers, but the report's own
`environment` block still records pip's OWN running interpreter (`python_full_version: 3.10.18`,
confirmed by reading the report), because pip evaluates PEP 508 environment markers (e.g. a
dependency declared only `; python_version<"3.11"`) against the interpreter actually running pip,
not the target. No marker-conditional dependency was found in this closure that visibly differs
between 3.10 and 3.11, but this was not exhaustively checked against every transitive package's
own metadata — the authoritative confirmation is still a resolution run FROM an actual Python
3.11 interpreter once `302006_de_backend_v1` exists, not simulated from `302006_env_clean`'s 3.10
`pip`, however complete the `--ignore-installed` closure looks.

**Earlier, incomplete dry-run (kept for its own evidence value)**: run without `--ignore-
installed`, from the SAME `302006_env_clean` `pip` (25.1.1): `Would install formulaic-1.2.2
formulaic-contrasts-1.0.0 interface_meta-2.0.1 narwhals-2.27.1 pydeseq2-0.5.4 session_info-1.0.1
stdlib-list-0.12.0` (and the identical 7-package list for `pydeseq2==0.5.2`/Python 3.10) — this
treated `302006_env_clean`'s own numpy/pandas/scipy/etc. as "already satisfied," which is why it
undercounted; its value is that `narwhals`/`interface_meta`/`session_info`/`stdlib-list` were not
independently visible from either run alone until compared.

Raw JSON reports (`pip_dry_run_report_py311_full.json`, 40-package closure; `pip_dry_run_report_
py311.json`/`pip_dry_run_report_py310_pydeseq052.json`, the earlier incomplete runs) are kept
under the HPCC scratch directory this round used, not committed to this repository
(machine-generated, large, and reproducible from the commands above).

## Draft install command for the new isolated environment (NOT RUN this round)

Per `review:20261014_221155` R4: uses an explicit `--prefix <absolute path>` rather than
only `-n <name>`, since `-n` alone resolves through conda's own `envs_dirs` search order and does
not itself GUARANTEE the absolute location — `--prefix` states it explicitly and unambiguously.

```bash
# Draft only -- environment CREATION and any install into it each require explicit
# authorization before execution, per DEVELOPMENT.md's standing rule that a code review is not
# itself an execution/environment-mutation authorization. The dry-run above is read-only and
# needed no such authorization; this command does. Per the review's own correction, this
# read-only-vs-authorized distinction applies to EVERY check in this document, not only
# installs -- DEVELOPMENT.md's "no SLURM submission without explicit authorization" is about SLURM
# jobs specifically, not a blanket requirement to re-ask before every metadata read.
conda create --prefix /rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1 \
    -c conda-forge python=3.11 pip
conda run --prefix /rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1 python -m pip install \
    "pydeseq2==0.5.4" "formulaic>=1.0.2" "formulaic-contrasts>=0.2.0"
# Fallback, same new env (its own --prefix), only if Python 3.11 provisioning is later found not
# to resolve cleanly -- this has NOT been retried on the login node after the earlier conda
# solver OOM below; it would need an interactive compute allocation or a repodata-size reduction
# (`--repodata-fn current_repodata.json`), per this round's own explicit instruction not to
# repeat an already-OOM'd solve on the login node:
#   conda create --prefix /rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1 \
#       -c conda-forge python=3.10 pip
#   conda run --prefix /rhome/xinruqiu/bigdata/.conda/envs/302006_de_backend_v1 python -m pip \
#       install "pydeseq2==0.5.2" "formulaic>=1.0.2" "formulaic-contrasts>=0.2.0"
```

`302006_env_clean` is never named in either branch above — this document's own v1 pointed its
install draft at that shared environment, which is exactly the gap this revision closes.

## Backend call draft: DataFrame counts + the SAME frozen design DataFrame + a numeric contrast

Per `review:20261014_203915` F3, RE-VERIFIED per `review:20261014_221155` R4
against **`v0.5.4`'s own source directly** (not carried over from `v0.5.2`, since the plan now
recommends `0.5.4` as primary): fetched `dds.py`/`ds.py` at BOTH tags from
`raw.githubusercontent.com/owkin/PyDESeq2`, 2026-10-14. `lfc_shrink`'s signature and
`DeseqDataSet`'s `counts`/`metadata`/`design` parameters are IDENTICAL between `v0.5.2` and
`v0.5.4` (confirmed by diffing the relevant sections directly) — the facts below hold for the
actually-recommended version, not merely inferred from the fallback's.

- `DeseqDataSet(counts=<pandas.DataFrame, samples x genes>, metadata=<pandas.DataFrame>,
  design=<pandas.DataFrame, samples x design columns>)` — `dds.py` accepts `design` as EITHER a
  formula string OR a `pandas.DataFrame` directly (`self.obsm["design_matrix"] = self.design`
  when a DataFrame is given, `dds.py` lines ~295-301); passing `de_design_providers.
  PreparedDEDesign.design_matrix`/`.column_names` directly as that DataFrame means PyDESeq2 never
  re-derives the design from a formula string, and the column order/encoding this framework
  already froze (and hash-binds, per F2's fix) is exactly what gets fit — never silently
  re-parsed or re-ordered by a separate formula engine.
- `DeseqStats(dds, contrast=<numpy.ndarray, length == len(design.columns)>)` — `ds.py` accepts a
  numeric contrast VECTOR aligned to the design matrix's own columns (not only a 3-element
  formula-style string contrast); `PreparedDEDesign.contrast_vector` is exactly this shape
  already.
- **`lfc_shrink` correction (the actual bug in v1)**: `ds.py`'s real signature is
  `def lfc_shrink(self, coeff: str, adapt: bool = True)` — `coeff` has **no usable default**
  despite its docstring's own text ("default: ``None``"); calling `lfc_shrink()` with no `coeff`
  argument at all is a `TypeError` (missing a required positional argument), and calling it with
  `coeff=None` explicitly raises `KeyError` from `if coeff not in self.LFC.columns:` (`None` is
  never a column name). **The plan must pass an EXPLICIT, verified target coefficient column
  name** — e.g. `f"condition[{comparison.comparison_level}]"`, read back from
  `PreparedDEDesign.contrast_name` (already recorded, already hash/recompute-bound per F2) — NOT
  rely on any default. Once a real backend exists, EVERY `'complete'` result must record both the
  unshrunken Wald LFC/p-value AND the shrunken LFC from `lfc_shrink(coeff=<verified name>)`, with
  the shrunken value never silently overwriting or replacing the raw one in the persisted record
  — the two are different quantities with different appropriate uses (shrunken LFC for
  ranking/reporting magnitude, raw p-value for significance), and conflating them would repeat
  the general pattern `DEVELOPMENT.md` rule 3 warns about (one quantity silently standing in for a
  different one). Gene filtering (e.g. a minimum total/mean count floor before fitting) and
  missing-value handling must be decided and recorded explicitly when this is actually
  implemented — NOT left as an unstated default; this plan does not pick those values now.

## Positive/null synthetic verification plan (not yet built)

Before any real backend call is trusted, this framework's existing convention (every other
provider in this round was proved against a hand-computable synthetic fixture, never only
against real data) must be extended to the actual DE statistic:

- **Positive control**: a synthetic pseudobulk counts matrix built with a KNOWN, injected
  log-fold-change for a known gene subset (e.g. multiply one group's counts for genes
  `{G0, G1}` by a fixed factor before adding Poisson/negative-binomial noise) and the REST of the
  genes held at equal expected counts between groups. The real backend, run through
  `de_design_providers.PreparedDEDesign`'s own design matrix/contrast (never a hand-rolled
  alternative), must recover the injected genes as the most significant/largest-magnitude hits,
  in the correct sign (per the design's own `contrast_direction`, "comparison minus reference").
- **Null control**: a synthetic matrix with NO injected difference between groups (same
  generating distribution for every gene, every group) run through the SAME real
  `PreparedDEDesign` path. The review's own instruction is explicit here and must not be
  relaxed: **do not require a random null's estimate to be exactly 0** (sampling noise means a
  null gene's estimated LFC will scatter around 0, not sit on it) — the actual check is
  calibration-shaped (e.g. the fraction of genes called significant at a fixed p-value threshold
  stays near that threshold's nominal rate across many resampled null draws), not a single
  point-estimate equality assertion.
- **Paired-null caution**: the review's instruction is also explicit that a PAIRED design's null
  must preserve pairing structure — label-shuffling that ignores which samples came from the SAME
  donor is NOT a valid null for a paired/donor-fixed-effect design (it would inject spurious
  between-donor variance into the "null", making the calibration check meaningless). Any null
  construction for `"paired_or_repeated_within_donor"` designs must shuffle group labels WITHIN
  each donor's own pair, never across donors.
- Both controls run through the real `DEDesignPreparationProvider`/future backend call, never a
  standalone statistics script that bypasses this framework's own design-matrix/contrast object —
  per `DEVELOPMENT.md` rule 5 ("test the real producer-consumer interface").

None of the above has been built or run this round; this section is the plan for a FUTURE round,
explicitly gated on (a) an owner decision to actually create `302006_de_backend_v1` and install
`pydeseq2==0.5.4` (or the `0.5.2`/Python 3.10 fallback) into it, each a separate explicit step,
and (b) this document being read and not objected to.

## What this document is NOT

- Not an installation and not an environment creation. Nothing in `302006_env_clean` changed
  this round and `302006_de_backend_v1` does not exist yet; `pydeseq2`/`rpy2`/`formulaic`/
  `formulaic-contrasts` remain absent from `302006_env_clean`, confirmed by the same live
  `pip list` query this document's facts came from.
- Not a claim that R is unavailable anywhere on the cluster — the opposite was verified directly.
- **Not a claim that per-SLURM-job authorization extends to every read-only check in this
  document.** Per `review:20261014_203915` F3's correction: `DEVELOPMENT.md`'s "no SLURM
  submission without explicit authorization" rule is about SLURM jobs specifically; it does not
  mean every version-metadata read, dependency-resolution dry-run, or draft command here needed
  a separate authorization round before being RUN — the `pip install --dry-run`/`conda create
  --dry-run` checks in this document were executed directly (nothing they do is reversible-risk
  the way an actual install or SLURM job is), and only environment CREATION and actual installs
  (explicitly marked "NOT RUN this round" above) are gated on a separate ask.
- Not a decision about which comparison, donor/visit structure, or masking score the real
  backend will eventually be run against — those remain the owner's decisions, per this
  framework's standing human-in-the-loop convention, and are entirely unaffected by which DE
  statistical package this document recommends.
