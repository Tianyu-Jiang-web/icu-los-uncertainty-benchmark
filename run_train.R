# ============================================================
# PIPELINE 1 - Full training (modular version of scripts.R)
#   EDA -> split -> HP tuning + training of 10 models -> 19-level
#   evaluation, heterogeneity, fairness, calibration, sensitivity
#   -> all_model_results.rds
# WARNING: tuning + training is computationally heavy (hours).
# Usage (from the project root):  Rscript run_train.R
# ============================================================
source("config.R")
source("R/pipeline.R")
load_modules()

dir.create(cfg$final_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- cfg$final_dir                      # variable name used by the original code

train_stages <- c(
  "stage_packages_train", "stage_prepare_features", "stage_eda_plots", "stage_split_data",
  "stage_train_setup_cv", "stage_train_define_crps_helpers",
  "stage_train_drf", "stage_train_qrf", "stage_train_ranger", "stage_train_xgboost",
  "stage_train_torch_common", "stage_train_gamma_helpers",
  "stage_train_ens", "stage_train_mcd", "stage_train_ddnn",
  "stage_train_edrf", "stage_train_bnn", "stage_train_hybrid",
  "stage_train_collect_hyperparameters",
  "stage_eval19_summary", "stage_eval19_build_eval_df",
  "stage_eval19_heterogeneity", "stage_eval19_fairness",
  "stage_eval19_calibration", "stage_eval19_sensitivity",
  "stage_train_save_outputs", "stage_train_package_report"
)
run_stages(train_stages)
