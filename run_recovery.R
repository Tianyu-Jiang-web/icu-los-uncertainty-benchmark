# ============================================================
# PIPELINE 2 - Torch recovery (modular version of script_recovery.R)
#   Reads all_model_results.rds, retrains the 5 torch-based models with
#   the saved best hyper-parameters (no tuning) and writes them with
#   torch_save() to <output>/final/torch_models/.
# Usage (from the project root):  Rscript run_recovery.R
# ============================================================
source("config.R")
source("R/pipeline.R")
load_modules()

dir.create(cfg$final_dir, recursive = TRUE, showWarnings = FALSE)
out_dir  <- cfg$final_dir                                       # torch_models/ is written here
rds_path <- file.path(cfg$models_dir, "all_model_results.rds")  # rds is read from here

recovery_stages <- c(
  "stage_packages_recovery", "stage_prepare_features", "stage_split_data",
  "stage_recovery_helpers", "stage_recovery_load_rds", "stage_recovery_torch_infra",
  "stage_recovery_trees", "stage_recovery_ens", "stage_recovery_mcd", "stage_recovery_ddnn",
  "stage_recovery_edrf", "stage_recovery_bnn", "stage_recovery_hybrid",
  "stage_recovery_assemble", "stage_recovery_check", "stage_recovery_save_torch"
)
run_stages(recovery_stages)
