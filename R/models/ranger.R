# ============================================================
# Random Forest quantile regression (ranger).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Ranger quantile RF: CRPS-based tuning + final fit + test quantiles.
# Source: scripts.R:763-827
stage_train_ranger <- function(env = .GlobalEnv) {
  local({
    library(ranger)

    # ============================================================
    # HP TUNING: Ranger — CRPS-based (req. 1 & 2)
    # Grid: num.trees x min.node.size = 2 x 3 = 6 configs
    # ============================================================
    cat("\n========================================\n")
    cat("Ranger RF: Tuning HPs using CRPS (re-tuned, req. 1 & 2)\n")
    cat("========================================\n")

    ranger_hp_grid <- expand.grid(
      num_trees     = c(500, 1000),
      min_node_size = c(5, 10, 20),
      stringsAsFactors = FALSE
    )

    val_idx_rng   <- cv_folds_hp[[1]]
    train_idx_rng <- setdiff(seq_len(nrow(Xtr_tree)), val_idx_rng)

    ranger_cv_crps <- sapply(seq_len(nrow(ranger_hp_grid)), function(i) {
      hp <- ranger_hp_grid[i, ]
      cat(sprintf("  [%d/%d] trees=%d node=%d ... ", i, nrow(ranger_hp_grid), hp$num_trees, hp$min_node_size))
      tryCatch({
        fit <- ranger(
          x = Xtr_tree[train_idx_rng, , drop = FALSE],
          y = y_tr[train_idx_rng],
          num.trees = hp$num_trees,
          min.node.size = hp$min_node_size,
          quantreg = TRUE,
          keep.inbag = TRUE
        )
        qm <- as.matrix(predict(fit, Xtr_tree[val_idx_rng, , drop = FALSE],
                                type = "quantiles", quantiles = q_grid)$predictions)
        cr <- crps_from_qmat(y_tr[val_idx_rng], qm)
        cat(sprintf("CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat("ERROR\n"); Inf })
    })

    best_rng_i  <- which.min(ranger_cv_crps)
    best_rng_hp <- as.list(ranger_hp_grid[best_rng_i, ])
    cat(sprintf("✓ Best Ranger HP: trees=%d node=%d (CRPS=%.4f)\n",
                best_rng_hp$num_trees, best_rng_hp$min_node_size, ranger_cv_crps[best_rng_i]))

    num_trees_ranger     <- best_rng_hp$num_trees
    min_node_size_ranger <- best_rng_hp$min_node_size

    rf_q_fit <- ranger(
      x = Xtr_tree,
      y = y_tr,
      num.trees = num_trees_ranger,
      min.node.size = min_node_size_ranger,
      quantreg = TRUE,
      keep.inbag = TRUE
    )

    rf_q <- as.data.frame(
      predict(rf_q_fit, Xte_tree, type = "quantiles", quantiles = qs)$predictions
    )
    colnames(rf_q) <- c("q05","q50","q95")

    rf_qgrid <- as.data.frame(
      predict(rf_q_fit, Xte_tree, type = "quantiles", quantiles = q_grid)$predictions
    )
    colnames(rf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    cat("✓ Ranger RF complete.\n")
  }, envir = env)
}

# ---- Ranger RF: 39-level prediction from the saved fit.
# Source: 39quantiles.R:644-647
stage_predict39_ranger <- function(env = .GlobalEnv) {
  local({
    cat("  Ranger RF...")
    rf_qmat <- as.matrix(predict(rf_q_fit, Xte_tree, type = "quantiles",
                                 quantiles = q_grid)$predictions)
    cat(" ✓\n")
  }, envir = env)
}
