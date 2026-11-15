.PHONY: help test test-cov test-fast lint format typecheck \
        docker-build docker-run docker-test \
        benchmark benchmark-stats benchmark-fusion \
        option-c option-c-mortality option-c-calibration option-c-activation \
        clean install-dev install-full

help:  ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-25s\033[0m %s\n", $$1, $$2}'

# ─── Testing ───────────────────────────────────────────────────────────────────

test:  ## Run unit tests
	pytest tests/ -v --tb=short

test-cov:  ## Run tests with coverage report
	pytest tests/ --cov=scripts --cov=src --cov-report=term-missing --cov-report=html

test-fast:  ## Run fast tests only (skip slow/integration markers)
	pytest tests/ -v -m "not slow and not integration"

# ─── Code Quality ──────────────────────────────────────────────────────────────

lint:  ## Run linter (ruff)
	ruff check scripts/ src/ tests/

format:  ## Format code (black + ruff --fix)
	black scripts/ src/ tests/
	ruff check --fix scripts/ src/ tests/

typecheck:  ## Run type checking (mypy)
	mypy scripts/ src/ --ignore-missing-imports

# ─── Docker ────────────────────────────────────────────────────────────────────

docker-build:  ## Build Docker image
	docker build -t scfm-benchmark:latest .

docker-run:  ## Run image with default CMD (prints help)
	docker run --rm scfm-benchmark:latest

docker-test:  ## Run unit tests inside the container
	docker run --rm scfm-benchmark:latest pytest tests/ -v -m "not slow and not integration"

# ─── Legacy Benchmarks (pre-Option-C, kept for reproducibility) ────────────────

benchmark:  ## Run legacy COVID severity benchmark (cell-level CV; superseded by Option C)
	python scripts/02_embedding_classification/code/covid_severity_benchmark.py

benchmark-stats:  ## Run statistical significance analysis (DeLong + bootstrap)
	python scripts/08_statistical_significance/code/statistical_significance_fast.py

benchmark-fusion:  ## Run multimodal fusion (all 7 strategies)
	python scripts/18_multimodal_fusion/code/multimodal_fusion.py --strategy all

# ─── Option C (Active 2026-05-01 pipeline) ─────────────────────────────────────
# These targets generate the SLURM submission commands; they don't sbatch
# directly because that would run on the wrong machine. Pipe to `sh` on HPCC.

option-c:  ## Print sbatch commands for Option C 5-FM re-run on canonical raw input
	@echo "# On HPCC:"
	@echo "bash scripts/24_option_c_canonical_input/code/step2c_human_covid_canonical_input.sh"

option-c-mortality:  ## Print sbatch for cross-disease mortality + Geneformer fine-tune
	@echo "# On HPCC:"
	@echo "sbatch scripts/06_batch_generalization/code/run_cross_disease_fatal.sh"
	@echo "sbatch scripts/finetuning/submit_geneformer_finetune_mortality.sh"

option-c-calibration:  ## Print sbatch for per-cohort calibration + per-cohort AUC
	@echo "# On HPCC (after Option C embeddings land):"
	@echo "sbatch scripts/24_option_c_canonical_input/code/run_calibration_per_cohort.sh"

option-c-activation:  ## Print sbatch for platelet activation marker analysis
	@echo "# On HPCC (after Option C embeddings land):"
	@echo "sbatch scripts/24_option_c_canonical_input/code/run_activation_markers.sh"

# ─── Utilities ─────────────────────────────────────────────────────────────────

clean:  ## Remove build artifacts and caches
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name htmlcov -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name .DS_Store -delete 2>/dev/null || true
	rm -rf .mypy_cache/ .ruff_cache/ dist/ build/ *.egg-info/

install-dev:  ## Install development dependencies
	pip install -e ".[dev]"

install-full:  ## Install all dependencies (ml + tracking + dev + esm + scprint2)
	pip install -e ".[dev,ml,tracking,esm,scprint2]"
