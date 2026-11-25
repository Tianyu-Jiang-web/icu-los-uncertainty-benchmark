# ============================================================
# PIPELINE 3 - Predict & visualise on the 39-level quantile grid
#   (modular version of 39quantiles.R, + crps_guidance_lancet.R)
#   Loads saved models (rds + torch_models/), re-predicts the test set,
#   then produces Table 1, evaluation summary, heterogeneity, fairness,
#   calibration and sensitivity analyses with Lancet-style figures.
# Usage (from the project root):  Rscript run_predict_visualize.R
# ============================================================
source("config.R")
source("R/pipeline.R")
load_modules()

# ---- paths used by the original code ----
base_fig_dir <- cfg$output_dir
input_dir    <- cfg$models_dir
out_dir      <- cfg$pred_dir
rds_path     <- file.path(cfg$models_dir, "all_model_results.rds")
torch_dir    <- cfg$torch_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

predict_stages <- c(
  "stage_packages_predict", "stage_prepare_features", "stage_split_data", "stage_table1",
  "stage_predict39_qgrid", "stage_predict39_define_helpers",
  "stage_predict39_load_rds", "stage_predict39_load_torch", "stage_predict39_setup",
  "stage_predict39_drf", "stage_predict39_qrf", "stage_predict39_ranger",
  "stage_predict39_xgboost", "stage_predict39_ens", "stage_predict39_mcd",
  "stage_predict39_ddnn", "stage_predict39_edrf", "stage_predict39_bnn",
  "stage_predict39_hybrid", "stage_predict39_assemble",
  "stage_eval39_summary",
  "stage_eval39_heterogeneity", "stage_plot39_heterogeneity_heatmap",
  "stage_plot39_heterogeneity_stability",
  "stage_eval39_fairness", "stage_plot39_fairness_disparity",
  "stage_eval39_calibration_pit", "stage_eval39_calibration_coverage",
  "stage_eval39_predictive_density", "stage_eval39_calibration_conditional",
  "stage_eval39_sensitivity", "stage_plot39_sensitivity",
  "stage_predict39_save_summary"
)
if (isTRUE(cfg$run_crps_guide)) predict_stages <- c(predict_stages, "stage_crps_guide")
run_stages(predict_stages)
