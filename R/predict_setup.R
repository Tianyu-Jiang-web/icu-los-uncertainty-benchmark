# ============================================================
# Set-up of the 39-quantile predict pipeline: load saved models, predict, assemble.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Define the 39-level q_grid and quantile indices.
# Source: 39quantiles.R:452-459
stage_predict39_qgrid <- function(env = .GlobalEnv) {
  local({
    # STEP 4: Denser q_grid (39-point journal standard)
    # ============================================================
    q_grid  <- cfg$q_grid_predict   # 39 levels
    n_q     <- length(q_grid)
    q05_idx <- which.min(abs(q_grid - 0.05))
    q50_idx <- which.min(abs(q_grid - 0.50))
    q95_idx <- which.min(abs(q_grid - 0.95))
    cat(sprintf("✓ q_grid: %d levels [%.3f, %.3f]\n", n_q, q_grid[1], q_grid[n_q]))
  }, envir = env)
}

# ---- crps_score_full / make_qgrid_df / crps_qmat_safe bound to the 39-level grid.
# Source: 39quantiles.R:482-483; 39quantiles.R:507-517
stage_predict39_define_helpers <- function(env = .GlobalEnv) {
  local({
    crps_from_qmat  <- function(y, qmat) crps_from_quantiles(y, as.matrix(qmat), q_grid)
    crps_score_full <- function(qmat, y) crps_from_quantiles(y, as.matrix(qmat), q_grid)

    make_qgrid_df <- function(mat) {
      df <- as.data.frame(mat)
      colnames(df) <- paste0("q", sprintf("%03d", round(1000 * q_grid)))
      df
    }

    # Enforce monotonicity + clip to 1, then compute CRPS (used in sensitivity loop)
    crps_qmat_safe <- function(pm, y_obs) {
      pm <- pmax(t(apply(pm, 1, sort)), 1)
      crps_from_quantiles(y_obs, pm, q_grid)
    }
  }, envir = env)
}

# ---- Load non-torch fitted models and hyper-parameters from all_model_results.rds.
# Source: 39quantiles.R:532-542
stage_predict39_load_rds <- function(env = .GlobalEnv) {
  local({
    cat(sprintf("\nLoading RDS: %s\n", rds_path))
    results <- readRDS(rds_path)

    drf_fit              <- results$drf_fit
    qrf_fit              <- results$qrf_fit
    rf_q_fit             <- results$rf_q_fit
    xgb_quantile_models  <- results$xgb_quantile_models
    drf_base_full        <- results$edrf_members
    best_hyperparameters <- results$best_hyperparams

    cat("✓ Non-torch models loaded\n")
  }, envir = env)
}

# ---- Neural feature matrices, helper functions and torch_load() of the saved .pt models.
# Source: 39quantiles.R:548-607
stage_predict39_load_torch <- function(env = .GlobalEnv) {
  local({
    device <- torch_device(select_torch_device())
    cat("Device:", device$type, "\n")

    # One-hot encode + standardise (train-led)
    Xtr_ens_mm <- model.matrix(~ . - 1, data = Xtr_all)
    Xte_ens_mm <- model.matrix(~ . - 1, data = Xte_all)
    common_nn  <- union(colnames(Xtr_ens_mm), colnames(Xte_ens_mm))
    Xtr_ens_mm <- Xtr_ens_mm[, common_nn, drop = FALSE]; Xtr_ens_mm[is.na(Xtr_ens_mm)] <- 0
    Xte_ens_mm <- Xte_ens_mm[, common_nn, drop = FALSE]; Xte_ens_mm[is.na(Xte_ens_mm)] <- 0
    x_mean     <- colMeans(Xtr_ens_mm)
    x_sd       <- apply(Xtr_ens_mm, 2, sd); x_sd[x_sd == 0] <- 1
    Xtr_ens_sc <- scale(Xtr_ens_mm, center = x_mean, scale = x_sd)
    Xte_ens_sc <- scale(Xte_ens_mm, center = x_mean, scale = x_sd)
    Xtr_ens_sc[Xtr_ens_sc >  10] <-  10; Xtr_ens_sc[Xtr_ens_sc < -10] <- -10
    Xte_ens_sc[Xte_ens_sc >  10] <-  10; Xte_ens_sc[Xte_ens_sc < -10] <- -10
    n_input <- ncol(Xtr_ens_sc)

    # Gamma helpers
    gamma_qgrid_fn <- function(shape, rate, qs)
      t(mapply(function(sh, ra) qgamma(qs, shape = sh, rate = ra), shape, rate))

    get_gamma_params <- function(net, X_sc) {
      X_sc <- as.matrix(X_sc); net$eval()
      Xv   <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      with_no_grad({ p <- as.matrix(net(Xv)$cpu()) })
      list(shape = exp(p[, 1]), rate = exp(p[, 2]))
    }

    get_gamma_params_ddnn <- function(net, X_sc) {
      X_sc <- as.matrix(X_sc); net$eval()
      Xv   <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      with_no_grad({ p <- as.matrix(net(Xv)$cpu()) })
      log_mu  <- p[, 1]; log_var <- pmin(pmax(p[, 2], -3), 4)
      mu <- exp(log_mu); v <- exp(log_var)
      list(shape = mu^2 / v, rate = mu / v)
    }

    # Align any perturbed test set to the training feature space + standardise
    prepare_nn_features <- function(Xte_imp_all) {
      mm <- model.matrix(~ . - 1, data = Xte_imp_all)
      for (cn in setdiff(common_nn, colnames(mm)))
        mm <- cbind(mm, matrix(0, nrow(mm), 1, dimnames = list(NULL, cn)))
      mm <- mm[, common_nn, drop = FALSE]; mm[is.na(mm)] <- 0
      sc <- scale(mm, center = x_mean, scale = x_sd)
      sc[sc >  10] <-  10; sc[sc < -10] <- -10
      sc
    }

    # Load torch models (script_recovery.R saved these via torch_save)
    ens_members <- lapply(seq_len(5), function(m)
      torch_load(file.path(torch_dir, sprintf("ens_member_%d.pt", m)), device = device))
    mc_dropout_net <- torch_load(file.path(torch_dir, "mcd_net.pt"), device = device)
    ddnn_nets      <- lapply(seq_len(4), function(run)
      torch_load(file.path(torch_dir, sprintf("ddnn_run_%d.pt", run)), device = device))
    bnn_net           <- torch_load(file.path(torch_dir, "bnn_net.pt"), device = device)
    nn_feat_extractor <- torch_load(file.path(torch_dir, "hybrid_feat_extractor.pt"),
                                    device = device)
    M_ens       <- 5
    N_runs_ddnn <- 4
    cat("✓ All torch models loaded\n")
  }, envir = env)
}

# ---- Test tensor and XGBoost design matrices / interpolation helper.
# Source: 39quantiles.R:613-633
stage_predict39_setup <- function(env = .GlobalEnv) {
  local({
    cat("\n--- Predicting all models ---\n")
    X_te_t <- torch_tensor(as.matrix(Xte_ens_sc), dtype = torch_float(), device = device)

    # XGBoost — original 19 models interpolated to 39 levels
    q_grid_orig <- seq(0.05, 0.95, by = 0.05)   # original training grid
    Xtr_all_mm  <- model.matrix(~ . - 1, data = Xtr_all)
    Xte_all_mm  <- model.matrix(~ . - 1, data = Xte_all)
    common_cols <- union(colnames(Xtr_all_mm), colnames(Xte_all_mm))
    Xtr_all_mm  <- Xtr_all_mm[, common_cols, drop = FALSE]; Xtr_all_mm[is.na(Xtr_all_mm)] <- 0
    Xte_all_mm  <- Xte_all_mm[, common_cols, drop = FALSE]; Xte_all_mm[is.na(Xte_all_mm)] <- 0
    dtest_full  <- xgb.DMatrix(Xte_all_mm)

    xgb_interp <- function(dmat, models_19) {
      orig_mat <- matrix(NA_real_, nrow(Xte_all_mm), length(q_grid_orig))
      for (ki in seq_along(q_grid_orig))
        orig_mat[, ki] <- predict(models_19[[ki]], dmat)
      orig_mat <- pmax(t(apply(orig_mat, 1, sort)), min(y_tr, na.rm = TRUE))
      new_mat  <- t(apply(orig_mat, 1, function(row)
        approx(x = q_grid_orig, y = row, xout = q_grid, method = "linear", rule = 2)$y))
      t(apply(new_mat, 1, sort))
    }
  }, envir = env)
}

# ---- Assemble all_q, all_qgrid, model_names, MODEL_COLOR.
# Source: 39quantiles.R:759-797
stage_predict39_assemble <- function(env = .GlobalEnv) {
  local({
    make_q3 <- function(qmat)
      data.frame(q05 = qmat[, q05_idx], q50 = qmat[, q50_idx], q95 = qmat[, q95_idx])

    all_q <- list(
      DRF           = make_q3(drf_qmat),
      QRF           = make_q3(qrf_qmat),
      RF            = make_q3(rf_qmat),
      XGB           = make_q3(xgb_qmat),
      ENS_Gamma     = make_q3(ens_qmat),
      MCD_Gamma     = make_q3(mcd_qmat),
      DDNN_Gamma    = make_q3(ddnn_qmat),
      Ensemble_DRF  = make_q3(edrf_qmat),
      BNN_Gamma     = make_q3(bnn_qmat),
      Hybrid_NN_DRF = make_q3(hybrid_qmat)
    )
    all_qgrid <- list(
      DRF           = make_qgrid_df(drf_qmat),
      QRF           = make_qgrid_df(qrf_qmat),
      RF            = make_qgrid_df(rf_qmat),
      XGB           = make_qgrid_df(xgb_qmat),
      ENS_Gamma     = make_qgrid_df(ens_qmat),
      MCD_Gamma     = make_qgrid_df(mcd_qmat),
      DDNN_Gamma    = make_qgrid_df(ddnn_qmat),
      Ensemble_DRF  = make_qgrid_df(edrf_qmat),
      BNN_Gamma     = make_qgrid_df(bnn_qmat),
      Hybrid_NN_DRF = make_qgrid_df(hybrid_qmat)
    )
    model_names <- names(all_q)

    MODEL_COLOR <- c(
      Actual        = "gray20",
      DRF           = "#D65C5C", QRF          = "#5A9BD4", RF            = "#66BB66",
      XGB           = "#E8A838", ENS_Gamma    = "#9B59B6", MCD_Gamma     = "#E67E22",
      DDNN_Gamma    = "#1ABC9C", Ensemble_DRF = "#C0392B", BNN_Gamma     = "#2471A3",
      Hybrid_NN_DRF = "#6E2F1A"
    )

    cat("\n✓ all_q and all_qgrid assembled\n")
  }, envir = env)
}
