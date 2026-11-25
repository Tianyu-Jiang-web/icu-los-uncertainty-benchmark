# ============================================================
# Distributional Random Forest.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- DRF: CRPS-based tuning + final fit + test quantiles.
# Source: scripts.R:593-595; scripts.R:622-691
stage_train_drf <- function(env = .GlobalEnv) {
  local({
    library(drf)

    set.seed(2026)

    # ============================================================
    # HP TUNING: DRF — CRPS-based (req. 1 & 2)
    # Grid: num.trees x min.node.size = 2 x 4 = 8 configs, 1 CV fold
    # ============================================================
    cat("\n========================================\n")
    cat("DRF: Tuning HPs using CRPS (re-tuned, req. 1 & 2)\n")
    cat("========================================\n")

    drf_hp_grid <- expand.grid(
      num_trees     = c(500, 1000),
      min_node_size = c(5, 10, 20, 50),
      stringsAsFactors = FALSE
    )

    val_idx_drf   <- cv_folds_hp[[1]]
    train_idx_drf <- setdiff(seq_len(nrow(Xtr_tree)), val_idx_drf)

    drf_cv_crps <- sapply(seq_len(nrow(drf_hp_grid)), function(i) {
      hp <- drf_hp_grid[i, ]
      cat(sprintf("  [%d/%d] trees=%d node=%d ... ", i, nrow(drf_hp_grid), hp$num_trees, hp$min_node_size))
      tryCatch({
        fit <- drf(
          X = Xtr_tree[train_idx_drf, , drop = FALSE],
          Y = y_tr[train_idx_drf],
          num.trees = hp$num_trees,
          min.node.size = hp$min_node_size,
          mtry = floor(sqrt(ncol(Xtr_tree))),
          sample.fraction = 0.5
        )
        pr  <- predict(fit, Xtr_tree[val_idx_drf, , drop = FALSE])
        qm  <- wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
        cr  <- crps_from_qmat(y_tr[val_idx_drf], qm)
        cat(sprintf("CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat("ERROR\n"); Inf })
    })

    best_drf_i  <- which.min(drf_cv_crps)
    best_drf_hp <- as.list(drf_hp_grid[best_drf_i, ])
    cat(sprintf("✓ Best DRF HP: trees=%d node=%d (CRPS=%.4f)\n",
                best_drf_hp$num_trees, best_drf_hp$min_node_size, drf_cv_crps[best_drf_i]))

    num_trees_drf     <- best_drf_hp$num_trees
    min_node_size_drf <- best_drf_hp$min_node_size

    # --- Final Fit ---
    set.seed(2026)
    drf_fit <- drf(
      X = Xtr_tree,
      Y = y_tr,
      num.trees = num_trees_drf,
      min.node.size = min_node_size_drf,
      mtry = floor(sqrt(ncol(Xtr_tree))),
      sample.fraction = 0.5
    )

    drf_pred <- predict(drf_fit, Xte_tree)
    W <- drf_pred$weights
    y_train <- as.numeric(drf_pred$y)

    # compute quantiles
    drf_q_mat <- wquant_matrix(y_train, W, c(0.05, 0.5, 0.95))
    drf_q <- as.data.frame(drf_q_mat)
    colnames(drf_q) <- c("q05", "q50", "q95")

    # full quantile grid
    drf_qgrid_mat <- wquant_matrix(y_train, W, q_grid)
    drf_qgrid <- as.data.frame(drf_qgrid_mat)
    colnames(drf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))

    cat("\n✓ DRF complete.\n")
  }, envir = env)
}

# ---- DRF: re-predict test set on the 39-level grid from the saved fit.
# Source: 39quantiles.R:635-638
stage_predict39_drf <- function(env = .GlobalEnv) {
  local({
    cat("  DRF...")
    drf_pred <- predict(drf_fit, Xte_tree)
    drf_qmat <- wquant_matrix(as.numeric(drf_pred$y), drf_pred$weights, q_grid)
    cat(" ✓\n")
  }, envir = env)
}
