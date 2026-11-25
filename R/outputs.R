# ============================================================
# Saving outputs and reporting.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Collect the selected hyper-parameters of all models.
# Source: scripts.R:1897-1908
stage_train_collect_hyperparameters <- function(env = .GlobalEnv) {
  local({
    best_hyperparameters <- list(
      DRF          = list(num.trees = num_trees_drf,    min.node.size = min_node_size_drf),
      QRF          = list(ntree = ntree_qrf,            nodesize = nodesize_qrf),
      Ranger       = list(num.trees = num_trees_ranger, min.node.size = min_node_size_ranger),
      XGBoost      = list(max_depth = max_depth_xgb,    eta = eta_xgb, nrounds = nrounds_xgb),
      ENS_Gamma    = best_ens_hp,
      MCD_Gamma    = best_mcd_hp,
      DDNN_Gamma   = best_ddnn_hp,
      Ensemble_DRF = list(K = K_drf, num.trees = num_trees_drf,
                          min.node.size = min_node_size_drf, combiner = "linear_QR"),
      BNN_Gamma    = best_bnn_hp
    )
  }, envir = env)
}

# ---- Write CSVs and all_model_results.rds.
# Source: scripts.R:2874-2924
stage_train_save_outputs <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # SAVE OUTPUTS
    # ============================================================
    write.csv(summary_all, file.path(out_dir, "summary_all_models.csv"), row.names=FALSE)

    # ---- Heterogeneity subgroup outputs ----
    if (!is.null(out_vasopressin_crps))
      write.csv(out_vasopressin_crps, file.path(out_dir, "heterogeneity_vasopressin_crps.csv"),  row.names=FALSE)
    if (!is.null(out_vaso_crps))
      write.csv(out_vaso_crps,        file.path(out_dir, "heterogeneity_vaso_crps.csv"),         row.names=FALSE)
    if (!is.null(out_vent_crps))
      write.csv(out_vent_crps,        file.path(out_dir, "heterogeneity_vent_crps.csv"),         row.names=FALSE)
    if (!is.null(out_sepsis_crps))
      write.csv(out_sepsis_crps,      file.path(out_dir, "heterogeneity_sepsis_crps.csv"),       row.names=FALSE)
    if (!is.null(out_stroke_crps))
      write.csv(out_stroke_crps,      file.path(out_dir, "heterogeneity_stroke_crps.csv"),       row.names=FALSE)
    if (!is.null(out_unit_crps))
      write.csv(out_unit_crps,        file.path(out_dir, "heterogeneity_unit_crps.csv"),         row.names=FALSE)
    if (!is.null(out_age_crps))
      write.csv(out_age_crps,         file.path(out_dir, "heterogeneity_age_crps.csv"),          row.names=FALSE)

    # ---- Fairness subgroup outputs ----
    if (!is.null(out_gender_crps))
      write.csv(out_gender_crps,      file.path(out_dir, "fairness_gender_crps.csv"),            row.names=FALSE)
    if (!is.null(out_race_crps))
      write.csv(out_race_crps,        file.path(out_dir, "fairness_race_crps.csv"),              row.names=FALSE)
    if (!is.null(out_insurance_crps))
      write.csv(out_insurance_crps,   file.path(out_dir, "fairness_insurance_crps.csv"),         row.names=FALSE)

    saveRDS(list(
      summary          = summary_all,
      crps             = crps_vec,
      best_hyperparams = best_hyperparameters,
      # model objects
      drf_fit          = drf_fit,
      qrf_fit          = qrf_fit,
      rf_q_fit         = rf_q_fit,
      xgb_quantile_models = xgb_quantile_models,
      ens_members      = ens_members,
      mc_dropout_net   = mc_dropout_net,
      ddnn_runs        = ddnn_runs,
      edrf_members     = drf_base_full,
      bnn_net          = bnn_net,
      hybrid_drf       = hybrid_drf,
      nn_feat_extractor = nn_feat_extractor
    ), file.path(out_dir, "all_model_results.rds"))

    cat("\n✓ All results saved to", out_dir, "\n")
    cat("\n========================================\n")
    cat("✓ ✓ ✓  ALL 10 MODELS COMPLETE  ✓ ✓ ✓\n")
    cat("========================================\n\n")
  }, envir = env)
}

# ---- Package / session version report.
# Source: scripts.R:2927-2972
stage_train_package_report <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # PART: Package version report (req. 5)
    # ============================================================
    cat("\n========================================\n")
    cat("PACKAGE VERSION REPORT\n")
    cat("========================================\n\n")

    # Collect all packages used in this script
    pkg_list <- c(
      # Data wrangling / EDA
      "dplyr", "readr", "lubridate", "tidyr",
      # Visualisation
      "ggplot2", "patchwork", "scales", "ggridges",
      # ML / modelling
      "caret",
      "drf",
      "quantregForest",
      "ranger",
      "xgboost",
      "torch",
      "scoringutils",
      # DL utilities
      "coro"
    )

    # R itself
    cat(sprintf("R version:  %s\n\n", R.version.string))

    # Per-package version
    pkg_versions <- lapply(pkg_list, function(pkg) {
      tryCatch({
        ver <- as.character(packageVersion(pkg))
        data.frame(Package = pkg, Version = ver, stringsAsFactors = FALSE)
      }, error = function(e) {
        data.frame(Package = pkg, Version = "NOT INSTALLED", stringsAsFactors = FALSE)
      })
    })

    pkg_df <- do.call(rbind, pkg_versions)
    rownames(pkg_df) <- NULL
    print(pkg_df, row.names = FALSE)

    cat("\n")
    # Also print as sessionInfo() for full reproducibility record
    cat("--- Full sessionInfo() ---\n")
    print(sessionInfo())
  }, envir = env)
}

# ---- Write summary_all_models.csv and print output inventory.
# Source: 39quantiles.R:3027-3048
stage_predict39_save_summary <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # SAVE ALL OUTPUTS SUMMARY
    # ============================================================
    write.csv(summary_all, file.path(out_dir, "summary_all_models.csv"), row.names = FALSE)

    cat("\n========================================\n")
    cat("✓ ALL OUTPUTS SAVED TO:", out_dir, "\n")
    cat("  q_grid: seq(0.025, 0.975, by=0.025) — 39 levels\n\n")
    cat("  Calibration plots:\n")
    cat("    pit_histogram.png\n")
    cat("    transformed_pit_histogram.png\n")
    cat("    coverage_calibration_plot.png\n")
    cat("    conditional_calibration_uncertainty_quintile.png\n")
    cat("    conditional_calibration_decile_table.csv\n\n")
    cat("  Heterogeneity CSVs (7): vasopressin, vaso, vent,\n")
    cat("    sepsis, stroke, unit, age\n\n")
    cat("  Fairness CSVs (3): gender, race, insurance\n\n")
    cat("  Sensitivity plots + tables:\n")
    cat("    sensitivity_missing_crps.png\n")
    cat("    sensitivity_missing_crps_pct.png\n")
    cat("    sensitivity_missing_results.csv\n")
    cat("========================================\n\n")
  }, envir = env)
}
