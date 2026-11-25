# ============================================================
# Overall evaluation (coverage, width, CRPS, MAE, RMSE).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Assemble all_q / all_qgrid, coverage, width, CRPS, MAE, RMSE (19 levels).
# Source: scripts.R:1914-2011
stage_eval19_summary <- function(env = .GlobalEnv) {
  local({
    cat("\n========================================\n")
    cat("UNIFIED EVALUATION — ALL MODELS\n")
    cat("========================================\n\n")

    eval_interval <- function(q, y) {
      covered <- (y >= q$q05) & (y <= q$q95)
      list(coverage = mean(covered, na.rm = TRUE),
           width    = mean(q$q95 - q$q05, na.rm = TRUE))
    }

    # CRPS (primary metric): per-observation, then averaged (req. 7)
    crps_score_full <- function(qmat, y, na.rm = TRUE) {
      crps_from_quantiles(y, as.matrix(qmat), q_grid)
    }



    # Named list: model label -> prediction data.frame (q05, q50, q95)
    all_q <- list(
      DRF          = drf_q,
      QRF          = qrf_q,
      RF           = rf_q,
      XGB          = xgb_q,
      ENS_Gamma    = ens_q,
      MCD_Gamma    = mcd_q,
      DDNN_Gamma   = ddnn_q,
      Ensemble_DRF = edrf_q,
      BNN_Gamma    = bnn_q,
      Hybrid_NN_DRF = hybrid_q
    )
    all_qgrid <- list(
      DRF          = drf_qgrid,
      QRF          = qrf_qgrid,
      RF           = rf_qgrid,
      XGB          = xgb_qgrid,
      ENS_Gamma    = ens_qgrid,
      MCD_Gamma    = mcd_qgrid,
      DDNN_Gamma   = ddnn_qgrid,
      Ensemble_DRF = edrf_qgrid,
      BNN_Gamma    = bnn_qgrid,
      Hybrid_NN_DRF = hybrid_qgrid
    )
    model_names <- names(all_q)

    # MODEL_COLOR palette (req. 4: remove density distribution plot but keep colors for other plots)
    MODEL_COLOR <- c(
      Actual        = "gray20",
      DRF           = "#D65C5C", QRF         = "#5A9BD4", RF            = "#66BB66",
      XGB           = "#E8A838", ENS_Gamma   = "#9B59B6", MCD_Gamma     = "#E67E22",
      DDNN_Gamma    = "#1ABC9C", Ensemble_DRF = "#C0392B", BNN_Gamma    = "#2471A3",
      Hybrid_NN_DRF = "#6E2F1A"
    )

    # Coverage & width
    cov_width <- do.call(rbind, lapply(model_names, function(m) {
      r <- eval_interval(all_q[[m]], y_te)
      data.frame(Model = m, Coverage = r$coverage, Width = r$width)
    }))

    cat("Coverage & Width:\n")
    cov_width[, c("Coverage", "Width")] <- round(cov_width[, c("Coverage", "Width")], 4)
    print(cov_width)

    # ============================================================
    # PRIMARY METRIC: CRPS (per-observation then averaged, req. 7)
    # ============================================================
    crps_vec <- sapply(model_names, function(m)
      crps_score_full(all_qgrid[[m]], y_te))
    cat("\nCRPS (lower = better; per-obs then averaged):\n"); print(round(crps_vec, 4))


    # MAE & RMSE (based on median prediction q50)
    mae_vec <- sapply(model_names, function(m) {
      q50 <- all_q[[m]]$q50
      mean(abs(y_te - q50), na.rm = TRUE)
    })
    rmse_vec <- sapply(model_names, function(m) {
      q50 <- all_q[[m]]$q50
      sqrt(mean((y_te - q50)^2, na.rm = TRUE))
    })
    cat("\nMAE (lower = better):\n");  print(round(mae_vec,  4))
    cat("\nRMSE (lower = better):\n"); print(round(rmse_vec, 4))

    # Summary table (ranked by CRPS — primary metric)
    summary_all <- data.frame(
      Model    = model_names,
      Coverage = cov_width$Coverage,
      Width    = cov_width$Width,
      MAE      = mae_vec,
      RMSE     = rmse_vec,
      CRPS     = crps_vec
      )
    summary_all$CRPS_Rank <- rank(summary_all$CRPS, na.last = "keep")
    cat("\n=== SUMMARY TABLE (ranked by CRPS) ===\n")
    tmp <- summary_all[order(summary_all$CRPS_Rank), ]
    numeric_cols_tmp <- sapply(tmp, is.numeric)
    tmp[, numeric_cols_tmp] <- round(tmp[, numeric_cols_tmp], 4)
    print(tmp)
  }, envir = env)
}

# ---- eval_df (wide) and subgroup CRPS helper (19 levels).
# Source: scripts.R:2014-2047
stage_eval19_build_eval_df <- function(env = .GlobalEnv) {
  local({
    # ---- eval_df: wide format for subgroup analysis ----
    eval_df <- do.call(data.frame, lapply(names(all_q), function(m) {
      mn <- tolower(gsub("[^a-zA-Z0-9]","_",m))
      q  <- all_q[[m]]
      setNames(q, paste0(mn, c("_q05","_q50","_q95")))
    }))
    eval_df$y <- y_te

    eval_df$los_group  <- cut(y_te, c(-Inf,3,7,Inf), labels=c("Short(≤3d)","Medium(3-7d)","Long(>7d)"))

    # ---- CRPS helper for subgroup (19 quantiles, per-obs averaged) ----
    make_idx_list <- function(group_var) {
      grp <- eval_df[[group_var]]
      lvls <- levels(grp)
      idx_list <- lapply(lvls, function(lv) which(grp == lv))
      names(idx_list) <- lvls
      list(idx_list = idx_list, labels = lvls)
    }

    sg_models <- model_names

    # Subgroup CRPS table (replaces previous WIS-based subgroup tables, req. 4)
    subgroup_crps_table <- function(df, group_var) {
      res      <- make_idx_list(group_var)
      crps_mat <- sapply(sg_models, function(m) {
        qmat <- as.matrix(all_qgrid[[m]])
        sapply(res$idx_list, function(idx) {
          crps_from_quantiles(y_te[idx], qmat[idx, , drop = FALSE], q_grid)
        })
      })
      df_out <- as.data.frame(crps_mat)
      cbind(Group = res$labels, df_out)
    }
  }, envir = env)
}

# ---- Coverage, width, MAE, RMSE, CRPS summary + eval_df + subgroup helper (39 levels).
# Source: 39quantiles.R:802-853
stage_eval39_summary <- function(env = .GlobalEnv) {
  local({
    cat("\n========================================\n")
    cat("EVALUATION\n")
    cat("========================================\n\n")

    eval_interval <- function(q, y) {
      covered <- (y >= q$q05) & (y <= q$q95)
      list(coverage = mean(covered, na.rm = TRUE),
           width    = mean(q$q95 - q$q05, na.rm = TRUE))
    }

    cov_width <- do.call(rbind, lapply(model_names, function(m) {
      r <- eval_interval(all_q[[m]], y_te)
      data.frame(Model = m, Coverage = r$coverage, Width = r$width)
    }))

    crps_vec <- sapply(model_names, function(m) crps_score_full(all_qgrid[[m]], y_te))
    mae_vec  <- sapply(model_names, function(m) mean(abs(y_te - all_q[[m]]$q50), na.rm = TRUE))
    rmse_vec <- sapply(model_names, function(m) sqrt(mean((y_te - all_q[[m]]$q50)^2, na.rm = TRUE)))

    summary_all <- data.frame(
      Model    = model_names, Coverage = cov_width$Coverage,
      Width    = cov_width$Width, MAE = mae_vec,
      RMSE     = rmse_vec, CRPS = crps_vec
    )
    summary_all$CRPS_Rank <- rank(summary_all$CRPS, na.last = "keep")
    cat("\n=== SUMMARY TABLE (ranked by CRPS) ===\n")
    tmp <- summary_all[order(summary_all$CRPS_Rank), ]
    tmp[, sapply(tmp, is.numeric)] <- round(tmp[, sapply(tmp, is.numeric)], 4)
    print(tmp)

    write.csv(summary_all, file.path(out_dir, "summary_all_models.csv"), row.names = FALSE)
    cat("✓ Saved: summary_all_models.csv\n")

    # eval_df: wide format used throughout subgroup analyses
    eval_df <- do.call(data.frame, lapply(names(all_q), function(m) {
      mn <- tolower(gsub("[^a-zA-Z0-9]", "_", m))
      setNames(all_q[[m]], paste0(mn, c("_q05", "_q50", "_q95")))
    }))
    eval_df$y <- y_te

    # Subgroup CRPS helper (mirrors subgroup_crps_table in scripts.R)
    subgroup_crps_table <- function(df_sub, group_var) {
      grp    <- df_sub[[group_var]]
      lvls   <- levels(grp)
      idx_ls <- lapply(lvls, function(lv) which(grp == lv)); names(idx_ls) <- lvls
      crps_mat <- sapply(model_names, function(m) {
        qmat <- as.matrix(all_qgrid[[m]])
        sapply(idx_ls, function(idx) crps_from_quantiles(y_te[idx], qmat[idx, , drop = FALSE], q_grid))
      })
      cbind(Group = lvls, as.data.frame(crps_mat))
    }
  }, envir = env)
}
