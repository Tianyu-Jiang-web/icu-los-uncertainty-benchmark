# ============================================================
# Ensemble DRF (stacked generalisation).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Ensemble DRF: stacked DRF base learners + quantile-regression combiners.
# Source: scripts.R:1483-1593
stage_train_edrf <- function(env = .GlobalEnv) {
  local({
    cat("========================================\n")
    cat("MODEL 9: Ensemble DRF (Stacking)\n")
    cat("========================================\n")

    # Stacking design follows Papacharalampous et al. (2025):
    # Split training set into Set1 (train base learners) and Set2 (train combiner).
    # This avoids predicting on the full training set, keeping weight matrices
    # at ~15k × 15k — well within memory limits.
    library(quantreg)

    set.seed(2026 + 900)
    n_tr     <- nrow(Xtr_tree)
    idx_set1 <- sample(n_tr, floor(n_tr / 2), replace = FALSE)
    idx_set2 <- setdiff(seq_len(n_tr), idx_set1)
    cat(sprintf("Stacking split: Set1=%d, Set2=%d\n", length(idx_set1), length(idx_set2)))

    X_set1 <- Xtr_tree[idx_set1, , drop = FALSE];  y_set1 <- y_tr[idx_set1]
    X_set2 <- Xtr_tree[idx_set2, , drop = FALSE];  y_set2 <- y_tr[idx_set2]

    K_drf <- 2
    cat(sprintf("Training %d DRF base learners on Set1...\n", K_drf))

    drf_base_s1 <- lapply(seq_len(K_drf), function(k) {
      mtry_k <- if (k == 1) floor(sqrt(ncol(Xtr_tree))) else floor(ncol(Xtr_tree) / 3)
      cat(sprintf("  Base learner %d: mtry=%d\n", k, mtry_k))
      drf(
        X               = X_set1,
        Y               = y_set1,
        num.trees       = num_trees_drf,
        min.node.size   = min_node_size_drf,
        mtry            = mtry_k,
        sample.fraction = 0.5,
        seed            = 2026 + k * 7
      )
    })

    cat("Predicting quantiles for Set2 (combiner training data)...\n")
    set2_qpreds <- lapply(drf_base_s1, function(fit) {
      pr <- predict(fit, X_set2)
      wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
    })

    cat("Fitting QR combiners for all quantile levels (Set2)...\n")
    qr_combiners <- vector("list", length(q_grid))

    for (qi in seq_along(q_grid)) {
      tau_i  <- q_grid[qi]
      X_comb <- do.call(cbind, lapply(set2_qpreds, function(qm) qm[, qi]))
      df_comb <- as.data.frame(X_comb)
      colnames(df_comb) <- paste0("bl", seq_len(K_drf))
      df_comb$y <- y_set2
      fmla <- as.formula(paste("y ~", paste(paste0("bl", seq_len(K_drf)), collapse = " + ")))
      qr_combiners[[qi]] <- tryCatch(
        rq(fmla, tau = tau_i, data = df_comb, method = "br"),
        error = function(e) NULL
      )
    }
    cat(sprintf("  Fitted QR combiners for %d quantile levels.\n",
                sum(!sapply(qr_combiners, is.null))))

    cat("Retraining DRF base learners on full training set...\n")
    drf_base_full <- lapply(seq_len(K_drf), function(k) {
      mtry_k <- if (k == 1) floor(sqrt(ncol(Xtr_tree))) else floor(ncol(Xtr_tree) / 3)
      cat(sprintf("  Base learner %d: mtry=%d (full train n=%d)\n",
                  k, mtry_k, nrow(Xtr_tree)))
      drf(
        X               = Xtr_tree,
        Y               = y_tr,
        num.trees       = num_trees_drf,
        min.node.size   = min_node_size_drf,
        mtry            = mtry_k,
        sample.fraction = 0.5,
        seed            = 2026 + k * 7 + 100
      )
    })

    cat("Predicting test quantiles from updated base learners...\n")
    te_qpreds <- lapply(drf_base_full, function(fit) {
      pr <- predict(fit, Xte_tree)
      wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
    })

    cat("Applying stacking combiner to test predictions...\n")
    edrf_qgrid_mat <- matrix(NA_real_, nrow(Xte_tree), length(q_grid))

    for (qi in seq_along(q_grid)) {
      X_te_comb <- do.call(cbind, lapply(te_qpreds, function(qm) qm[, qi]))
      df_te_comb <- as.data.frame(X_te_comb)
      colnames(df_te_comb) <- paste0("bl", seq_len(K_drf))

      if (!is.null(qr_combiners[[qi]])) {
        preds <- predict(qr_combiners[[qi]], newdata = df_te_comb)
        edrf_qgrid_mat[, qi] <- pmax(preds, 1)
      } else {
        edrf_qgrid_mat[, qi] <- rowMeans(X_te_comb)
      }
    }

    for (i in seq_len(nrow(edrf_qgrid_mat))) {
      edrf_qgrid_mat[i, ] <- sort(edrf_qgrid_mat[i, ])
    }

    edrf_q <- data.frame(
      q05 = edrf_qgrid_mat[, q05_idx],
      q50 = edrf_qgrid_mat[, q50_idx],
      q95 = edrf_qgrid_mat[, q95_idx]
    )
    edrf_qgrid <- as.data.frame(edrf_qgrid_mat)
    colnames(edrf_qgrid) <- paste0("q", sprintf("%02d", round(100 * q_grid)))

    cat("Ensemble DRF (stacking) complete.\n\n")
  }, envir = env)
}

# ---- Ensemble DRF: rebuild QR combiners and predict on the 39-level grid.
# Source: 39quantiles.R:677-722
stage_predict39_edrf <- function(env = .GlobalEnv) {
  local({
    cat("  Ensemble DRF (rebuilding QR combiners)...")
    K_drf <- length(drf_base_full)
    set.seed(2026 + 900)
    n_tr     <- nrow(Xtr_tree)
    idx_set1 <- sample(n_tr, floor(n_tr / 2), replace = FALSE)
    idx_set2 <- setdiff(seq_len(n_tr), idx_set1)
    X_set2   <- Xtr_tree[idx_set2, , drop = FALSE]; y_set2 <- y_tr[idx_set2]
    te_qpreds <- lapply(drf_base_full, function(fit) {
      pr <- predict(fit, Xte_tree)
      wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
    })
    drf_base_s1 <- lapply(seq_len(K_drf), function(k) {
      mtry_k <- if (k == 1) floor(sqrt(ncol(Xtr_tree))) else floor(ncol(Xtr_tree) / 3)
      drf(X = Xtr_tree[idx_set1, , drop = FALSE], Y = y_tr[idx_set1],
          num.trees     = best_hyperparameters$Ensemble_DRF$num.trees,
          min.node.size = best_hyperparameters$Ensemble_DRF$min.node.size,
          mtry = mtry_k, sample.fraction = 0.5, seed = 2026 + k * 7)
    })
    set2_qpreds <- lapply(drf_base_s1, function(fit) {
      pr <- predict(fit, X_set2)
      wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
    })
    qr_combiners <- vector("list", n_q)
    for (qi in seq_len(n_q)) {
      tau_i   <- q_grid[qi]
      X_comb  <- do.call(cbind, lapply(set2_qpreds, function(qm) qm[, qi]))
      df_comb <- as.data.frame(X_comb)
      colnames(df_comb) <- paste0("bl", seq_len(K_drf))
      df_comb$y <- y_set2
      fmla <- as.formula(paste("y ~", paste(paste0("bl", seq_len(K_drf)), collapse = " + ")))
      qr_combiners[[qi]] <- tryCatch(rq(fmla, tau = tau_i, data = df_comb, method = "br"),
                                     error = function(e) NULL)
    }
    edrf_qmat <- matrix(NA_real_, nrow(Xte_tree), n_q)
    for (qi in seq_len(n_q)) {
      X_te_comb  <- do.call(cbind, lapply(te_qpreds, function(qm) qm[, qi]))
      df_te_comb <- as.data.frame(X_te_comb)
      colnames(df_te_comb) <- paste0("bl", seq_len(K_drf))
      if (!is.null(qr_combiners[[qi]])) {
        edrf_qmat[, qi] <- pmax(predict(qr_combiners[[qi]], newdata = df_te_comb), 1)
      } else {
        edrf_qmat[, qi] <- rowMeans(X_te_comb)
      }
    }
    for (i in seq_len(nrow(edrf_qmat))) edrf_qmat[i, ] <- sort(edrf_qmat[i, ])
    cat(" ✓\n")
  }, envir = env)
}
