# Single-cell Foundation Model Robustness Benchmark — runtime container.
#
# Builds a CUDA-enabled environment that supports both the downstream
# benchmark scripts (logistic regression, calibration, fusion, etc.) and
# optional in-container fine-tuning experiments. Foundation-model
# embedding-generation jobs themselves run on HPCC GPUs in their own
# per-FM conda envs (state_env, scprint2_env, scpeft_env, stack_env,
# transcriptformer_env, saturn_env) and are not built into this image —
# see scripts/foundation_models/ for those.
#
# Build:   docker build -t scfm-benchmark:latest .
# Run:     docker run --rm -v $(pwd)/data:/app/data scfm-benchmark:latest pytest tests/
# CLI:     docker run --rm scfm-benchmark:latest scfm-bench --help

FROM pytorch/pytorch:2.1.0-cuda11.8-cudnn8-runtime

LABEL maintainer="Xinru Qiu <xqiu014@ucr.edu>"
LABEL description="Single-cell foundation model benchmark — analysis runtime"
LABEL repository="https://github.com/xqiu625/scfm-benchmark"

WORKDIR /app

# System dependencies (minimal — git for VCS-installed deps, build tools for any sdist)
RUN apt-get update && apt-get install -y --no-install-recommends \
        git \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy package metadata + source first so the editable install layer can be
# cached when only scripts/ changes. src/scfm_benchmark/ is the
# importable package referenced by pyproject.toml [tool.hatch.build.targets.wheel].
COPY pyproject.toml README.md ./
COPY src/ src/

# Install package + ML + tracking + dev extras directly from pyproject (no requirements.txt needed).
# `dev` brings pytest/ruff/black/mypy so tests run in the container.
RUN pip install --upgrade pip && \
    pip install --no-cache-dir -e ".[ml,tracking,dev]"

# Copy non-package source (scripts/ has ~67 MB of results CSVs; the .dockerignore
# excludes scripts/*/figures/ and scripts/*/results/ so only code/ is brought in).
COPY scripts/ scripts/
COPY tests/ tests/

# Strict import + CLI sanity check at build time — fails the build if either is broken.
RUN python -c "import scfm_benchmark; print('scfm_benchmark', scfm_benchmark.__version__)" && \
    scfm-bench --version

# Default: print CLI help. Overridable via `docker run ... <cmd>`.
CMD ["scfm-bench", "--help"]
