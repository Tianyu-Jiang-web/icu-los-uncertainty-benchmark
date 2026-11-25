# ============================================================
# Sensitivity to introduced test-time missingness.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Robustness to introduced missingness - 19 levels.
# Source: scripts.R:2577-2872
stage_eval19_sensitivity <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # SENSITIVITY ANALYSIS: Robustness to Introduced Missingness
    # Metric changed to CRPS (req. 1 & 2)
    # ============================================================
    cat("\n========================================\n")
    cat("SENSITIVITY ANALYSIS: Robustness to Missing Data (CRPS-based)\n")
    cat("========================================\n\n")

    num_cols_te <- names(Xte_tree)[sapply(Xte_tree, is.numeric)]
    num_cols_te <- grep("_miss$", num_cols_te, value = TRUE, invert = TRUE)
    cat(sprintf("Numeric variables for missingness injection: %d\n", length(num_cols_te)))

    train_medians <- sapply(num_cols_te, function(v) {
      if (v %in% names(Xtr_tree)) median(Xtr_tree[[v]], na.rm = TRUE) else NA_real_
    })

    miss_rates <- c(0.10, 0.20, 0.30)
    n_seeds    <- 3
    sens_seeds <- c(42L, 123L, 999L)

    inject_missing <- function(Xte_orig, miss_rate, cols, seed) {
      set.seed(seed)
      Xte_new <- Xte_orig
      n_rows   <- nrow(Xte_new)
      for (v in cols) {
        if (!v %in% names(Xte_new)) next
        n_mask   <- round(miss_rate * n_rows)
        if (n_mask == 0) next
        idx_mask <- sample(n_rows, n_mask, replace = FALSE)
        Xte_new[[v]][idx_mask] <- NA_real_
      }
      Xte_new
    }

    impute_with_train_median <- function(Xte_df, medians) {
      Xte_imp <- Xte_df
      for (v in names(medians)) {
        if (!v %in% names(Xte_imp)) next
        if (!is.numeric(Xte_imp[[v]])) next
        na_idx <- is.na(Xte_imp[[v]])
        if (any(na_idx) && !is.na(medians[v]))
          Xte_imp[[v]][na_idx] <- medians[v]
      }
      Xte_imp
    }

    prepare_nn_features <- function(Xte_imp_all) {
      mm <- model.matrix(~ . - 1, data = Xte_imp_all)
      for (cn in setdiff(common_nn, colnames(mm))) {
        mm <- cbind(mm, matrix(0, nrow(mm), 1, dimnames = list(NULL, cn)))
      }
      mm <- mm[, common_nn, drop = FALSE]
      mm[is.na(mm)] <- 0
      sc <- scale(mm, center = x_mean, scale = x_sd)
      sc[sc >  10] <-  10
      sc[sc < -10] <- -10
      sc
    }

    # CRPS helper for sensitivity loop (19 quantiles, per-obs averaged)
    crps19 <- function(pm) {
      pm <- pmax(t(apply(pm, 1, sort)), 1)
      crps_from_quantiles(y_te, pm, q_grid)
    }

    # ---- Baseline (0% missing) from already-computed results ----
    baseline_crps <- sapply(names(all_qgrid), function(m) {
      crps_score_full(all_qgrid[[m]], y_te)
    })

    sens_results <- list()
    sens_results[["0.00"]] <- data.frame(
      model      = names(baseline_crps),
      miss_rate  = 0.00,
      crps_mean  = as.numeric(baseline_crps),
      crps_sd    = 0.0
    )

    # ---- Loop over miss rates ----
    for (mr in miss_rates) {
      cat(sprintf("\n--- Miss rate = %.0f%% ---\n", mr * 100))

      crps_by_seed <- lapply(sens_seeds, function(sd) {
        Xte_tree_m <- inject_missing(Xte_tree, mr, num_cols_te, sd)

        num_cols_all <- intersect(num_cols_te, names(Xte_all))
        Xte_all_m    <- inject_missing(Xte_all, mr, num_cols_all, sd)
        Xte_all_imp  <- impute_with_train_median(Xte_all_m, train_medians)

        Xte_nn_sc <- prepare_nn_features(Xte_all_imp)

        # DRF: native missing handling
        drf_pred_m  <- predict(drf_fit, Xte_tree_m)
        drf_q_m     <- wquant_matrix(as.numeric(drf_pred_m$y), drf_pred_m$weights, q_grid)
        crps_drf    <- crps19(drf_q_m)

        # QRF: imputed
        qrf_q_m  <- as.matrix(predict(qrf_fit, Xte_all_imp, what = q_grid))
        crps_qrf <- crps19(qrf_q_m)

        # Ranger: imputed
        rf_q_m   <- as.matrix(predict(rf_q_fit, Xte_all_imp,
                                      type = "quantiles", quantiles = q_grid)$predictions)
        crps_rf  <- crps19(rf_q_m)

        # XGBoost: imputed + one-hot, predict from each quantile model
        Xte_mm_m <- model.matrix(~ . - 1, data = Xte_all_imp)
        for (cn in setdiff(common_cols, colnames(Xte_mm_m)))
          Xte_mm_m <- cbind(Xte_mm_m, matrix(0, nrow(Xte_mm_m), 1, dimnames=list(NULL,cn)))
        Xte_mm_m <- Xte_mm_m[, common_cols, drop = FALSE]; Xte_mm_m[is.na(Xte_mm_m)] <- 0
        dtest_m  <- xgb.DMatrix(Xte_mm_m)
        xgb_qgrid_m <- matrix(NA_real_, nrow(Xte_mm_m), length(q_grid))
        for (ki in seq_along(q_grid)) {
          xgb_qgrid_m[, ki] <- predict(xgb_quantile_models[[ki]], dtest_m)
        }
        xgb_qgrid_m <- pmax(t(apply(xgb_qgrid_m, 1, sort)), min(y_tr))
        crps_xgb    <- crps19(xgb_qgrid_m)

        # ENS: Gamma, imputed + NN features (parameter-space averaging)
        ens_params_m <- lapply(ens_members, function(net) get_gamma_params(net, Xte_nn_sc))
        ens_shape_m  <- Reduce(`+`, lapply(ens_params_m, `[[`, "shape")) / M_ens
        ens_rate_m   <- Reduce(`+`, lapply(ens_params_m, `[[`, "rate"))  / M_ens
        ens_qgrid_m  <- pmax(gamma_qgrid_fn(ens_shape_m, ens_rate_m, q_grid), 1)
        crps_ens     <- crps19(ens_qgrid_m)

        # MCD: Gamma, imputed + NN features
        mcd_q_m     <- mcd_gamma_quantiles(mc_dropout_net, Xte_nn_sc, M_pass = 50)
        crps_mcd    <- crps19(mcd_q_m)

        # DDNN: Gamma, imputed + NN features
        ddnn_params_m <- lapply(ddnn_runs, function(r) get_gamma_params(r$net, Xte_nn_sc))
        ddnn_shape_m  <- Reduce(`+`, lapply(ddnn_params_m, `[[`, "shape")) / N_runs_ddnn
        ddnn_rate_m   <- Reduce(`+`, lapply(ddnn_params_m, `[[`, "rate"))  / N_runs_ddnn
        ddnn_q_m      <- pmax(gamma_qgrid_fn(ddnn_shape_m, ddnn_rate_m, q_grid), 1)
        crps_ddnn     <- crps19(ddnn_q_m)

        # Ensemble DRF
        te_qpreds_m <- lapply(drf_base_full, function(fit) {
          pr <- predict(fit, Xte_tree_m)
          wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
        })
        edrf_mat_m <- matrix(NA_real_, nrow(Xte_tree_m), length(q_grid))
        for (qi in seq_along(q_grid)) {
          X_comb  <- do.call(cbind, lapply(te_qpreds_m, function(qm) qm[, qi]))
          df_te_c <- as.data.frame(X_comb)
          colnames(df_te_c) <- paste0("bl", seq_len(K_drf))
          if (!is.null(qr_combiners[[qi]])) {
            edrf_mat_m[, qi] <- pmax(predict(qr_combiners[[qi]], newdata = df_te_c), 1)
          } else {
            edrf_mat_m[, qi] <- pmax(rowMeans(X_comb), 1)
          }
        }
        for (i in seq_len(nrow(edrf_mat_m))) edrf_mat_m[i, ] <- sort(edrf_mat_m[i, ])
        crps_edrf <- crps19(edrf_mat_m)

        # BNN: Gamma, imputed + NN features
        bnn_q_m   <- bnn_gamma_quantiles(bnn_net, Xte_nn_sc, M_pass = 50)
        crps_bnn  <- crps19(bnn_q_m)

        # Hybrid NN+DRF
        nn_feat_extractor$eval()
        X_te_t_m <- torch_tensor(Xte_nn_sc, dtype = torch_float(), device = device)
        with_no_grad({ neural_te_m <- as.matrix(nn_feat_extractor(X_te_t_m)) })
        colnames(neural_te_m) <- paste0("nn_feat_", seq_len(ncol(neural_te_m)))
        hyb_pred_m    <- predict(hybrid_drf, newdata = as.data.frame(neural_te_m))
        hyb_q_m       <- wquant_matrix(as.numeric(hyb_pred_m$y), hyb_pred_m$weights, q_grid)
        hyb_q_m       <- pmax(hyb_q_m, min(y_tr))
        crps_hyb      <- crps19(hyb_q_m)

        c(DRF          = crps_drf,
          QRF          = crps_qrf,
          RF           = crps_rf,
          XGB          = crps_xgb,
          ENS_Gamma    = crps_ens,
          MCD_Gamma    = crps_mcd,
          DDNN_Gamma   = crps_ddnn,
          Ensemble_DRF = crps_edrf,
          BNN_Gamma    = crps_bnn,
          Hybrid_NN_DRF = crps_hyb)
      })

      crps_mat <- do.call(rbind, crps_by_seed)
      cat(sprintf("  CRPS means:\n"))
      print(round(colMeans(crps_mat), 4))

      sens_results[[sprintf("%.2f", mr)]] <- data.frame(
        model     = colnames(crps_mat),
        miss_rate = mr,
        crps_mean = colMeans(crps_mat),
        crps_sd   = apply(crps_mat, 2, sd)
      )
    }

    # ---- Combine into long data frame ----
    sens_df <- bind_rows(sens_results)
    rownames(sens_df) <- NULL

    cat("\n=== Sensitivity Analysis — CRPS by Miss Rate ===\n")
    tmp_sens <- tidyr::pivot_wider(sens_df[,c("model","miss_rate","crps_mean")],
                                   names_from = miss_rate, values_from = crps_mean)
    num_cols_s <- sapply(tmp_sens, is.numeric)
    tmp_sens[, num_cols_s] <- round(tmp_sens[, num_cols_s], 4)
    print(tmp_sens)

    # ---- CRPS degradation relative to baseline ----
    baseline_df_s <- sens_df %>% filter(miss_rate == 0) %>%
      select(model, crps_baseline = crps_mean)
    sens_df <- sens_df %>%
      left_join(baseline_df_s, by = "model") %>%
      mutate(crps_delta      = crps_mean - crps_baseline,
             crps_pct_change = 100 * crps_delta / crps_baseline)

    cat("\n=== CRPS % change from baseline (lower slope = more robust) ===\n")
    tmp_pct <- tidyr::pivot_wider(
      sens_df %>% filter(miss_rate > 0) %>% select(model, miss_rate, crps_pct_change),
      names_from = miss_rate, values_from = crps_pct_change)
    num_cols_p <- sapply(tmp_pct, is.numeric)
    tmp_pct[, num_cols_p] <- round(tmp_pct[, num_cols_p], 2)
    print(tmp_pct)

    # ---- Sensitivity plot: CRPS vs miss rate ----
    drf_family <- c("DRF", "Ensemble_DRF", "Hybrid_NN_DRF")

    sens_plot_df <- sens_df %>%
      mutate(
        model      = factor(model, levels = model_names),
        is_drf_fam = model %in% drf_family,
        line_size  = ifelse(model == "DRF", 1.5, 0.8),
        alpha_val  = ifelse(is_drf_fam, 1.0, 0.6)
      )

    p_sens <- ggplot(sens_plot_df,
                     aes(x = miss_rate * 100, y = crps_mean,
                         colour = model, group = model,
                         linewidth = I(line_size), alpha = I(alpha_val))) +
      geom_ribbon(aes(ymin = crps_mean - crps_sd, ymax = crps_mean + crps_sd,
                      fill = model), alpha = 0.08, colour = NA) +
      geom_line() +
      geom_point(size = 2.2) +
      scale_colour_manual(values = MODEL_COLOR) +
      scale_fill_manual(values   = MODEL_COLOR) +
      scale_x_continuous(breaks = c(0, 10, 20, 30),
                         labels = c("0%\n(baseline)", "10%", "20%", "30%")) +
      labs(
        title    = "Sensitivity Analysis: CRPS Degradation under Introduced Missingness",
        subtitle = paste0(
          "Missingness injected uniformly at random across all numeric variables\n",
          "DRF: native missing-value handling  |  All other models: median imputation\n",
          "Shaded band = ±1 SD across 3 random seeds  |  Flatter slope = more robust"),
        x      = "Additional missingness injected into test set (% of values per numeric variable)",
        y      = "CRPS (lower = better)",
        colour = "Model", fill = "Model"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, colour = "grey40"),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )

    print(p_sens)
    ggsave(file.path(out_dir, "sensitivity_missing_crps.png"),
           plot = p_sens, width = 12, height = 7, dpi = 300)
    cat("✓ Saved: sensitivity_missing_crps.png\n")

    # ---- CRPS % change plot ----
    p_sens_delta <- ggplot(
      sens_df %>% filter(miss_rate > 0) %>%
        mutate(model = factor(model, levels = model_names)),
      aes(x = miss_rate * 100, y = crps_pct_change, colour = model, group = model)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(aes(linewidth = I(ifelse(model == "DRF", 1.5, 0.8)))) +
      geom_point(size = 2.2) +
      scale_colour_manual(values = MODEL_COLOR) +
      scale_x_continuous(breaks = c(10, 20, 30), labels = c("10%","20%","30%")) +
      labs(
        title    = "CRPS Degradation (% change from baseline) under Introduced Missingness",
        subtitle = "Dashed line = no change from baseline  |  Lower / flatter = more robust",
        x        = "Additional missingness injected (%)",
        y        = "CRPS % change from baseline",
        colour   = "Model"
      ) +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank())

    print(p_sens_delta)
    ggsave(file.path(out_dir, "sensitivity_missing_crps_pct.png"),
           plot = p_sens_delta, width = 12, height = 7, dpi = 300)
    cat("✓ Saved: sensitivity_missing_crps_pct.png\n")

    # Save sensitivity table
    write.csv(sens_df, file.path(out_dir, "sensitivity_missing_results.csv"), row.names = FALSE)
    cat("✓ Saved: sensitivity_missing_results.csv\n\n")
  }, envir = env)
}

# ---- Missingness-injection sensitivity analysis (39 levels).
# Source: 39quantiles.R:2697-2874
stage_eval39_sensitivity <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # STEP 14: SENSITIVITY ANALYSIS — robustness to missingness
    # ============================================================
    cat("\n========================================\n")
    cat("SENSITIVITY ANALYSIS: Robustness to Missing Data (CRPS-based)\n")
    cat("========================================\n\n")

    num_cols_te <- names(Xte_tree)[sapply(Xte_tree, is.numeric)]
    num_cols_te <- grep("_miss$", num_cols_te, value = TRUE, invert = TRUE)
    cat(sprintf("Numeric variables for missingness injection: %d\n", length(num_cols_te)))

    train_medians <- sapply(num_cols_te, function(v) {
      if (v %in% names(Xtr_tree)) median(Xtr_tree[[v]], na.rm = TRUE) else NA_real_
    })

    miss_rates <- c(0.10, 0.20, 0.30)
    sens_seeds <- c(42L, 123L, 999L)

    inject_missing <- function(Xte_orig, miss_rate, cols, seed) {
      set.seed(seed); Xte_new <- Xte_orig; n_rows <- nrow(Xte_new)
      for (v in cols) {
        if (!v %in% names(Xte_new)) next
        n_mask <- round(miss_rate * n_rows)
        if (n_mask == 0) next
        Xte_new[[v]][sample(n_rows, n_mask, replace = FALSE)] <- NA_real_
      }
      Xte_new
    }

    impute_with_train_median <- function(Xte_df, medians) {
      Xte_imp <- Xte_df
      for (v in names(medians)) {
        if (!v %in% names(Xte_imp) || !is.numeric(Xte_imp[[v]])) next
        na_idx <- is.na(Xte_imp[[v]])
        if (any(na_idx) && !is.na(medians[v])) Xte_imp[[v]][na_idx] <- medians[v]
      }
      Xte_imp
    }

    # Baseline (0% missing) from already-computed all_qgrid
    baseline_crps <- sapply(names(all_qgrid), function(m) crps_score_full(all_qgrid[[m]], y_te))
    sens_results  <- list()
    sens_results[["0.00"]] <- data.frame(
      model     = names(baseline_crps), miss_rate = 0.00,
      crps_mean = as.numeric(baseline_crps), crps_sd = 0.0
    )

    for (mr in miss_rates) {
      cat(sprintf("\n--- Miss rate = %.0f%% ---\n", mr * 100))

      crps_by_seed <- lapply(sens_seeds, function(sd) {
        Xte_tree_m <- inject_missing(Xte_tree, mr, num_cols_te, sd)
        num_cols_all <- intersect(num_cols_te, names(Xte_all))
        Xte_all_m    <- inject_missing(Xte_all, mr, num_cols_all, sd)
        Xte_all_imp  <- impute_with_train_median(Xte_all_m, train_medians)
        Xte_nn_sc    <- prepare_nn_features(Xte_all_imp)
        X_te_t_m     <- torch_tensor(Xte_nn_sc, dtype = torch_float(), device = device)

        # DRF: native missing handling
        drf_pred_m <- predict(drf_fit, Xte_tree_m)
        drf_q_m    <- wquant_matrix(as.numeric(drf_pred_m$y), drf_pred_m$weights, q_grid)
        crps_drf   <- crps_qmat_safe(drf_q_m, y_te)

        # QRF: imputed
        qrf_q_m  <- as.matrix(predict(qrf_fit, Xte_all_imp, what = q_grid))
        crps_qrf <- crps_qmat_safe(qrf_q_m, y_te)

        # Ranger: imputed
        rf_q_m  <- as.matrix(predict(rf_q_fit, Xte_all_imp,
                                     type = "quantiles", quantiles = q_grid)$predictions)
        crps_rf <- crps_qmat_safe(rf_q_m, y_te)

        # XGBoost: imputed + one-hot + interpolate
        Xte_mm_m <- model.matrix(~ . - 1, data = Xte_all_imp)
        for (cn in setdiff(common_cols, colnames(Xte_mm_m)))
          Xte_mm_m <- cbind(Xte_mm_m, matrix(0, nrow(Xte_mm_m), 1, dimnames = list(NULL, cn)))
        Xte_mm_m <- Xte_mm_m[, common_cols, drop = FALSE]; Xte_mm_m[is.na(Xte_mm_m)] <- 0
        dtest_m  <- xgb.DMatrix(Xte_mm_m)
        xgb_orig_m <- matrix(NA_real_, nrow(Xte_mm_m), length(q_grid_orig))
        for (ki in seq_along(q_grid_orig))
          xgb_orig_m[, ki] <- predict(xgb_quantile_models[[ki]], dtest_m)
        xgb_orig_m <- pmax(t(apply(xgb_orig_m, 1, sort)), min(y_tr))
        xgb_q_m    <- t(apply(xgb_orig_m, 1, function(row)
          approx(x = q_grid_orig, y = row, xout = q_grid, method = "linear", rule = 2)$y))
        crps_xgb   <- crps_qmat_safe(xgb_q_m, y_te)

        # ENS Gamma
        ens_params_m <- lapply(ens_members, function(net) get_gamma_params(net, Xte_nn_sc))
        ens_shape_m  <- Reduce(`+`, lapply(ens_params_m, `[[`, "shape")) / M_ens
        ens_rate_m   <- Reduce(`+`, lapply(ens_params_m, `[[`, "rate"))  / M_ens
        ens_q_m      <- pmax(gamma_qgrid_fn(ens_shape_m, ens_rate_m, q_grid), 1)
        crps_ens     <- crps_qmat_safe(ens_q_m, y_te)

        # MCD Gamma
        mc_dropout_net$train()
        acc_m <- matrix(0, nrow(Xte_nn_sc), 2)
        for (i in seq_len(50)) with_no_grad({ acc_m <- acc_m + as.matrix(mc_dropout_net(X_te_t_m)$cpu()) })
        mcd_q_m  <- pmax(gamma_qgrid_fn(exp(acc_m[, 1] / 50), exp(acc_m[, 2] / 50), q_grid), 1)
        crps_mcd <- crps_qmat_safe(mcd_q_m, y_te)

        # DDNN Gamma
        ddnn_params_m <- lapply(ddnn_nets, function(net) get_gamma_params_ddnn(net, Xte_nn_sc))
        ddnn_shape_m  <- Reduce(`+`, lapply(ddnn_params_m, `[[`, "shape")) / N_runs_ddnn
        ddnn_rate_m   <- Reduce(`+`, lapply(ddnn_params_m, `[[`, "rate"))  / N_runs_ddnn
        ddnn_q_m      <- pmax(gamma_qgrid_fn(ddnn_shape_m, ddnn_rate_m, q_grid), 1)
        crps_ddnn     <- crps_qmat_safe(ddnn_q_m, y_te)

        # Ensemble DRF: native missing handling
        te_qpreds_m <- lapply(drf_base_full, function(fit) {
          pr <- predict(fit, Xte_tree_m)
          wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
        })
        edrf_mat_m <- matrix(NA_real_, nrow(Xte_tree_m), n_q)
        for (qi in seq_len(n_q)) {
          X_comb  <- do.call(cbind, lapply(te_qpreds_m, function(qm) qm[, qi]))
          df_te_c <- as.data.frame(X_comb)
          colnames(df_te_c) <- paste0("bl", seq_len(K_drf))
          if (!is.null(qr_combiners[[qi]])) {
            edrf_mat_m[, qi] <- pmax(predict(qr_combiners[[qi]], newdata = df_te_c), 1)
          } else {
            edrf_mat_m[, qi] <- pmax(rowMeans(X_comb), 1)
          }
        }
        for (i in seq_len(nrow(edrf_mat_m))) edrf_mat_m[i, ] <- sort(edrf_mat_m[i, ])
        crps_edrf <- crps_qmat_safe(edrf_mat_m, y_te)

        # BNN Gamma
        bnn_net$train()
        acc_bm <- matrix(0, nrow(Xte_nn_sc), 2)
        for (i in seq_len(50)) with_no_grad({ acc_bm <- acc_bm + as.matrix(bnn_net(X_te_t_m)$cpu()) })
        bnn_q_m  <- pmax(gamma_qgrid_fn(exp(acc_bm[, 1] / 50), exp(acc_bm[, 2] / 50), q_grid), 1)
        crps_bnn <- crps_qmat_safe(bnn_q_m, y_te)

        # Hybrid NN+DRF
        nn_feat_extractor$eval()
        with_no_grad({ neural_te_m <- as.matrix(nn_feat_extractor(X_te_t_m)$cpu()) })
        colnames(neural_te_m) <- paste0("nn_feat_", seq_len(ncol(neural_te_m)))
        hyb_pred_m <- predict(hybrid_drf, newdata = as.data.frame(neural_te_m))
        hyb_q_m    <- pmax(wquant_matrix(as.numeric(hyb_pred_m$y),
                                         hyb_pred_m$weights, q_grid), min(y_tr))
        crps_hyb   <- crps_qmat_safe(hyb_q_m, y_te)

        c(DRF = crps_drf, QRF = crps_qrf, RF = crps_rf, XGB = crps_xgb,
          ENS_Gamma = crps_ens, MCD_Gamma = crps_mcd, DDNN_Gamma = crps_ddnn,
          Ensemble_DRF = crps_edrf, BNN_Gamma = crps_bnn, Hybrid_NN_DRF = crps_hyb)
      })

      crps_mat <- do.call(rbind, crps_by_seed)
      cat(sprintf("  CRPS means:\n")); print(round(colMeans(crps_mat), 4))

      sens_results[[sprintf("%.2f", mr)]] <- data.frame(
        model     = colnames(crps_mat), miss_rate = mr,
        crps_mean = colMeans(crps_mat), crps_sd   = apply(crps_mat, 2, sd)
      )
    }

    sens_df <- bind_rows(sens_results); rownames(sens_df) <- NULL

    cat("\n=== Sensitivity Analysis — CRPS by Miss Rate ===\n")
    tmp_sens <- tidyr::pivot_wider(sens_df[, c("model", "miss_rate", "crps_mean")],
                                   names_from = miss_rate, values_from = crps_mean)
    tmp_sens[, sapply(tmp_sens, is.numeric)] <- round(tmp_sens[, sapply(tmp_sens, is.numeric)], 4)
    print(tmp_sens)

    baseline_df_s <- sens_df %>% filter(miss_rate == 0) %>%
      select(model, crps_baseline = crps_mean)
    sens_df <- sens_df %>%
      left_join(baseline_df_s, by = "model") %>%
      mutate(crps_delta      = crps_mean - crps_baseline,
             crps_pct_change = 100 * crps_delta / crps_baseline)

    cat("\n=== CRPS % change from baseline (lower slope = more robust) ===\n")
    tmp_pct <- tidyr::pivot_wider(
      sens_df %>% filter(miss_rate > 0) %>% select(model, miss_rate, crps_pct_change),
      names_from = miss_rate, values_from = crps_pct_change)
    tmp_pct[, sapply(tmp_pct, is.numeric)] <- round(tmp_pct[, sapply(tmp_pct, is.numeric)], 2)
    print(tmp_pct)
  }, envir = env)
}
