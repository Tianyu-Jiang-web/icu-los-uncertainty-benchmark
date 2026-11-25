# ============================================================
# XGBoost direct quantile regression.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- XGBoost: one pinball-loss model per quantile level (tuning + final fit).
# Source: scripts.R:845-961
stage_train_xgboost <- function(env = .GlobalEnv) {
  local({
    library(xgboost)
    set.seed(2026)

    stopifnot(exists("Xtr_all"), exists("Xte_all"), exists("y_tr"))
    stopifnot(nrow(Xtr_all) == length(y_tr))

    # One-hot encoding for XGBoost
    Xtr_all_mm_hp <- model.matrix(~ . - 1, data = Xtr_all)
    Xte_all_mm    <- model.matrix(~ . - 1, data = Xte_all)

    # Align columns (req. 3: same approach as before for feature alignment)

    common_cols <- union(colnames(Xtr_all_mm_hp), colnames(Xte_all_mm))
    Xtr_all_mm_hp <- Xtr_all_mm_hp[, common_cols, drop = FALSE]
    Xte_all_mm    <- Xte_all_mm[,    common_cols, drop = FALSE]

    Xtr_all_mm_hp[is.na(Xtr_all_mm_hp)] <- 0
    Xte_all_mm[is.na(Xte_all_mm)]       <- 0

    cat("\n========================================\n")
    cat("XGBoost: Tuning HPs using CRPS + quantile loss (req. 1, 2, 3)\n")
    cat("========================================\n")

    # ---- HP tuning: use median regression (alpha=0.5) for speed, then evaluate CRPS ----
    val_idx_xgb   <- cv_folds_hp[[1]]
    train_idx_xgb <- setdiff(seq_len(nrow(Xtr_all_mm_hp)), val_idx_xgb)

    dtrain_cv <- xgb.DMatrix(Xtr_all_mm_hp[train_idx_xgb, , drop = FALSE],
                             label = y_tr[train_idx_xgb])
    dval_cv   <- xgb.DMatrix(Xtr_all_mm_hp[val_idx_xgb,   , drop = FALSE],
                             label = y_tr[val_idx_xgb])

    xgb_hp_grid <- expand.grid(
      max_depth = c(4, 6),
      eta       = c(0.05, 0.1),
      nrounds   = c(300, 500),
      stringsAsFactors = FALSE
    )

    xgb_cv_crps <- sapply(seq_len(nrow(xgb_hp_grid)), function(i) {
      hp <- xgb_hp_grid[i, ]
      cat(sprintf("  [%d/%d] depth=%d eta=%.2f nrounds=%d ... ",
                  i, nrow(xgb_hp_grid), hp$max_depth, hp$eta, hp$nrounds))
      tryCatch({
        # Train K quantile models on CV train split, predict on val, compute CRPS
        qpreds_val <- matrix(NA_real_, length(val_idx_xgb), length(q_grid))
        for (ki in seq_along(q_grid)) {
          tau_k <- q_grid[ki]
          params_k <- list(
            objective  = "reg:quantileerror",
            quantile_alpha = tau_k,
            max_depth  = hp$max_depth,
            eta        = hp$eta,
            verbosity  = 0
          )
          fit_k <- xgb.train(params = params_k, data = dtrain_cv,
                             nrounds = hp$nrounds, verbose = 0)
          qpreds_val[, ki] <- predict(fit_k, dval_cv)
        }
        # Enforce monotonicity across quantiles
        qpreds_val <- t(apply(qpreds_val, 1, sort))
        cr <- crps_from_qmat(y_tr[val_idx_xgb], qpreds_val)
        cat(sprintf("CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat(sprintf("ERROR: %s\n", conditionMessage(e))); Inf })
    })

    best_xgb_i  <- which.min(xgb_cv_crps)
    best_xgb_hp <- as.list(xgb_hp_grid[best_xgb_i, ])
    cat(sprintf("✓ Best XGBoost HP: depth=%d eta=%.2f nrounds=%d (CRPS=%.4f)\n",
                best_xgb_hp$max_depth, best_xgb_hp$eta, best_xgb_hp$nrounds,
                xgb_cv_crps[best_xgb_i]))

    max_depth_xgb <- best_xgb_hp$max_depth
    eta_xgb       <- best_xgb_hp$eta
    nrounds_xgb   <- best_xgb_hp$nrounds

    # ---- Final fit: train one model per quantile on full training set (req. 3) ----
    cat(sprintf("\nXGBoost: Training %d quantile models on full train set...\n", length(q_grid)))

    dtrain_full <- xgb.DMatrix(Xtr_all_mm_hp, label = y_tr)
    dtest_full  <- xgb.DMatrix(Xte_all_mm)

    xgb_quantile_models <- vector("list", length(q_grid))
    xgb_qgrid_mat       <- matrix(NA_real_, nrow(Xte_all_mm), length(q_grid))

    for (ki in seq_along(q_grid)) {
      tau_k  <- q_grid[ki]
      cat(sprintf("  Fitting quantile tau=%.2f (%d/%d)...\n", tau_k, ki, length(q_grid)))
      params_k <- list(
        objective      = "reg:quantileerror",
        quantile_alpha = tau_k,
        max_depth      = max_depth_xgb,
        eta            = eta_xgb,
        verbosity      = 0
      )
      fit_k <- xgb.train(params = params_k, data = dtrain_full,
                         nrounds = nrounds_xgb, verbose = 0)
      xgb_quantile_models[[ki]] <- fit_k
      xgb_qgrid_mat[, ki]       <- predict(fit_k, dtest_full)
    }

    # Enforce monotonicity and clip to training range lower bound
    xgb_qgrid_mat <- t(apply(xgb_qgrid_mat, 1, sort))
    xgb_qgrid_mat <- pmax(xgb_qgrid_mat, min(y_tr, na.rm = TRUE))

    q05_idx_xgb <- which.min(abs(q_grid - 0.05))
    q50_idx_xgb <- which.min(abs(q_grid - 0.50))
    q95_idx_xgb <- which.min(abs(q_grid - 0.95))

    xgb_q <- data.frame(
      q05 = xgb_qgrid_mat[, q05_idx_xgb],
      q50 = xgb_qgrid_mat[, q50_idx_xgb],
      q95 = xgb_qgrid_mat[, q95_idx_xgb]
    )
    xgb_qgrid <- as.data.frame(xgb_qgrid_mat)
    colnames(xgb_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    cat("✓ XGBoost (direct quantile regression) complete.\n")
  }, envir = env)
}

# ---- XGBoost: predict with the 19 saved models and linearly interpolate to 39 levels.
# Source: 39quantiles.R:649-651
stage_predict39_xgboost <- function(env = .GlobalEnv) {
  local({
    cat("  XGBoost (interpolating 19 -> 39 levels)...")
    xgb_qmat <- xgb_interp(dtest_full, xgb_quantile_models)
    cat(" ✓\n")
  }, envir = env)
}
