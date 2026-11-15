# ACC 20: Drug Target Validation

## Goal
Validate FM-nominated therapeutic targets against known drug databases, demonstrating drug discovery relevance.

## Pipeline
1. Load perturbation recovery rankings from ACC 11 (STATE + UCE)
2. Build ground truth: curated anti-inflammatory/anti-COVID/anti-thrombosis targets from ChEMBL, OpenTargets
3. Compute precision@k (k=5,10,20,50,100) and recall@k
4. Fisher's exact test for enrichment of known targets in top-k
5. Classify top-10 candidates: "validated" / "plausible" / "novel"

## Execution
```bash
# First: run interactively to cache API results
python analyses/20_drug_target_validation/code/validate_therapeutic_targets.py

# Then: submit batch job
sbatch analyses/20_drug_target_validation/code/submit_20_validation.sh
```

## Expected Results
- Top-10 includes 3-5 known targets
- Precision@10 significantly above random baseline (Fisher p < 0.05)
- STATE rankings show higher enrichment than UCE

## Env
- `302006_env_clean` | batch partition, CPU-only, 32GB, 4h
