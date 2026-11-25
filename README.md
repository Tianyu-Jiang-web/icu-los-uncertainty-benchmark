# Uncertainty Quantification for ICU Length of Stay — modular R pipeline

Code for my MSc Health Data Science dissertation, *Uncertainty Quantification for ICU Length of Stay: A Comparative Benchmark on MIMIC-IV*. It benchmarks ten probabilistic models on 38,181 first ICU stays and evaluates them on accuracy, calibration, subgroup heterogeneity, fairness and robustness to missing inputs.

This repository is a **modular refactor** of the original research scripts. Behaviour, seeds, train/test split, hyper-parameter grids and quantile grids are unchanged; the code was reorganised into small, single-purpose modules with a config file and three entry points.

(中文说明: [README.zh-CN.md](README.zh-CN.md))

**Live demos** (GitHub Pages, served from `docs/`):
[Results dashboard](https://tianyu-jiang-web.github.io/icu-los-uncertainty-benchmark/result/) ·
[CRPS explorer](https://tianyu-jiang-web.github.io/icu-los-uncertainty-benchmark/crps-explorer/)

## Models compared

DRF · Quantile Regression Forest · Random Forest (ranger) · XGBoost (one pinball-loss model per quantile) · Ensemble DRF (stacking) · Hybrid NN + DRF · Deep Ensemble · MC Dropout · Bayesian NN (Bayes by Backprop) · Distributional DNN. Primary metric: CRPS on a 39-level quantile grid (0.025–0.975); also MAE, RMSE, 90% coverage/width, PIT, coverage calibration, subgroup CRPS, fairness and missingness sensitivity.

Headline result (test set, 80:20 temporal split): the Hybrid NN–DRF model has the lowest CRPS (≈1.05), followed by MC Dropout and QRF.

## Project layout

```
config.R                    paths, quantile grids, switches (env-var overridable)
run_train.R                 pipeline 1: tuning + training + 19-level evaluation
run_recovery.R              pipeline 2: retrain torch models with saved hyper-parameters
run_predict_visualize.R     pipeline 3: 39-level prediction, tables and all figures
R/
  pipeline.R                module loader, stage runner, device selection
  utils_metrics.R           pure functions: pinball loss, CRPS, weighted quantiles
  data_prep.R               feature construction + temporal split (shared)
  models/                   one file per model (train + predict stages)
  analysis/                 evaluation, heterogeneity, fairness, calibration, sensitivity
  plots/                    publication-style (Lancet) figures
  recovery.R, outputs.R, predict_setup.R, table1.R, eda_plots.R, ...
```

Each `stage_*` function is one block of the original analysis; pipelines are ordered lists of stages, so any subset can be run, e.g.

```bash
DRF_SKIP="stage_eval39_sensitivity,stage_plot39_sensitivity" Rscript run_predict_visualize.R
```

## How the refactor was checked

- All 73 stage functions were compared with the original scripts at the parsed-expression level: **0 differences** (apart from documented path/device substitutions).
- The shared data-prep module reproduces the original train/test matrices exactly (`identical()` on all objects).
- Dependency check along each pipeline: no undefined objects.
- Re-running the prediction/evaluation stages from saved models reproduced the original per-model CRPS: tree models and Ensemble DRF identical, ENS/DDNN to 1e-8; Hybrid/MC-Dropout/BNN within 3e-3 (stochastic inference, CPU vs. MPS).

## Running

Developed and tested with R 4.5; packages: dplyr, readr, lubridate, tidyr, ggplot2, patchwork, scales, ggridges, caret, drf, quantregForest, ranger, xgboost, torch, scoringutils, coro, quantreg, ggforce, gt, gtsummary.

```bash
Rscript run_train.R               # hours: tuning + training
Rscript run_recovery.R            # retrain torch models, writes outputs/final/torch_models
Rscript run_predict_visualize.R   # needs all_model_results.rds (+ torch_models/)
```

Paths are set in `config.R` (`DRF_DATA_PATH`, `DRF_OUTPUT_DIR`, `DRF_MODELS_DIR`, `DRF_FORCE_CPU=1`).

## Data

**No patient data or trained models are included.** MIMIC-IV is a credentialed-access dataset (PhysioNet); obtain access at https://physionet.org/content/mimiciv/ and build the analysis table as described in the dissertation. The pipeline expects it at `DRF_DATA_PATH`.

Method for DRF: Ćevid et al., *Distributional Random Forests*, JMLR 2022; Ensemble DRF follows Papacharalampous et al. License: MIT.
