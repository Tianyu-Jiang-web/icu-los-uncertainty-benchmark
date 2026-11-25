# ============================================================
# Quantile Regression Forest.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- QRF: CRPS-based tuning + final fit + test quantiles.
# Source: scripts.R:697-757
stage_train_qrf <- function(env = .GlobalEnv) {
  local({
    library(quantregForest)

    # ============================================================
    # HP TUNING: QRF — CRPS-based (req. 1 & 2)
    # Grid: ntree x nodesize = 2 x 3 = 6 configs
    # ============================================================
    cat("\n========================================\n")
    cat("QRF: Tuning HPs using CRPS (re-tuned, req. 1 & 2)\n")
    cat("========================================\n")

    qrf_hp_grid <- expand.grid(
      ntree    = c(500, 1000),
      nodesize = c(10, 20, 50),
      stringsAsFactors = FALSE
    )

    val_idx_qrf   <- cv_folds_hp[[1]]
    train_idx_qrf <- setdiff(seq_len(nrow(Xtr_all)), val_idx_qrf)

    qrf_cv_crps <- sapply(seq_len(nrow(qrf_hp_grid)), function(i) {
      hp <- qrf_hp_grid[i, ]
      cat(sprintf("  [%d/%d] ntree=%d nodesize=%d ... ", i, nrow(qrf_hp_grid), hp$ntree, hp$nodesize))
      tryCatch({
        fit <- quantregForest(
          x = Xtr_all[train_idx_qrf, , drop = FALSE],
          y = y_tr[train_idx_qrf],
          ntree = hp$ntree,
          nodesize = hp$nodesize
        )
        qm <- as.matrix(predict(fit, Xtr_all[val_idx_qrf, , drop = FALSE], what = q_grid))
        cr <- crps_from_qmat(y_tr[val_idx_qrf], qm)
        cat(sprintf("CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat("ERROR\n"); Inf })
    })

    best_qrf_i  <- which.min(qrf_cv_crps)
    best_qrf_hp <- as.list(qrf_hp_grid[best_qrf_i, ])
    cat(sprintf("✓ Best QRF HP: ntree=%d nodesize=%d (CRPS=%.4f)\n",
                best_qrf_hp$ntree, best_qrf_hp$nodesize, qrf_cv_crps[best_qrf_i]))

    ntree_qrf    <- best_qrf_hp$ntree
    nodesize_qrf <- best_qrf_hp$nodesize

    qs <- c(0.05, 0.5, 0.95)

    # fit qrf model
    qrf_fit <- quantregForest(
      x = Xtr_all,
      y = y_tr,
      ntree = ntree_qrf,
      nodesize = nodesize_qrf
    )

    qrf_q <- as.data.frame(
      predict(qrf_fit, Xte_all, what = qs)
    )
    colnames(qrf_q) <- c("q05","q50","q95")

    qrf_qgrid <- as.data.frame(predict(qrf_fit, Xte_all, what = q_grid))
    colnames(qrf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    cat("✓ QRF complete.\n")
  }, envir = env)
}

# ---- QRF: 39-level prediction from the saved fit.
# Source: 39quantiles.R:640-642
stage_predict39_qrf <- function(env = .GlobalEnv) {
  local({
    cat("  QRF...")
    qrf_qmat <- as.matrix(predict(qrf_fit, Xte_all, what = q_grid))
    cat(" ✓\n")
  }, envir = env)
}
