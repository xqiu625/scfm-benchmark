# 27_covid_sepsis_only_rerun

Re-run of the primary in-hospital mortality binary benchmark on the **COVID + Sepsis subset only**.

## Why

The original Option C 6-FM mortality benchmark was trained on a 5-disease-class pool (COVID-19 + Flu + Lung_infect + Non_covid + Sepsis; HC and SLE excluded). After audit of the dataset composition (`tab:disease_composition` in main.tex), we restrict the analysis to the two clinically-paired inflammatory diseases (COVID-19 and Sepsis) for a cleaner cross-disease story.

## Filter

- `Disease in {'COVID-19', 'Sepsis'}` (drops 5,882 cells / 91 donors: Flu, HC, Lung_infect, Non_covid, SLE)
- `covid_severity in {'fatal', 'mild', 'moderate', 'recovered', 'severe'}` (drops nothing further within COVID/Sepsis since neither has 'healthy' or 'unknown')

After filter: **42,095 cells / 322 donors / 8,850 fatal / 33,245 survived**.
(Per-FM cell counts vary by ~13–868 depending on each FM's native preprocessing.)

## Pipeline

5-fold patient-level `StratifiedGroupKFold` on `Sample_ID`; LR (C=1.0, lbfgs, balanced) + RF (n=100, balanced); seed 42; n_jobs=-1 for RF.

## How to run

```bash
# On HPCC bluejay/skylark, after git pull:
sbatch /bigdata/godziklab/shared/Xinru/302006/scripts/27_covid_sepsis_only_rerun/code/rerun_covid_sepsis.sh
```

## Outputs

- `05_RESULTS/option_c_primary_benchmark/covid_sepsis_only_mortality.csv` — main result table (per-FM × LR/RF AUC, BalAcc, fold AUCs)
- `scripts/27_covid_sepsis_only_rerun/logs/rerun_<JOBID>.{out,err}` — SLURM logs

## Expected runtime

~30 min on epyc with 16 cpus + 64G mem (LR fits are <10s each; RF is the bottleneck at ~1–2 min × 5 folds × 6 FMs).
