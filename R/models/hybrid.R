# ============================================================
# Hybrid neural network + DRF.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Hybrid NN+DRF: QRNN feature extractor then DRF on neural features.
# Source: scripts.R:1795-1891
stage_train_hybrid <- function(env = .GlobalEnv) {
  local({
    cat("========================================\n")
    cat("MODEL 11: Hybrid NN + DRF\n")
    cat("========================================\n")

    cat("Stage 1: Training neural feature extractor (QRNN, 150 epochs)...\n")
    set.seed(2026)

    y_mean_ens  <- mean(y_tr)
    y_sd_ens    <- sd(y_tr)
    y_tr_sc     <- (y_tr - y_mean_ens) / y_sd_ens

    pinball_loss <- function(pred, y_true, tau) {
      err  <- y_true - pred
      loss <- torch_where(err >= 0, tau * err, (tau - 1) * err)
      torch_mean(loss)
    }

    nn_feat_extractor <- nn_sequential(
      nn_linear(n_input, 128), nn_relu(),
      nn_linear(128, 32),      nn_relu()
    )
    nn_out_head <- nn_linear(32, n_quantiles)

    nn_feat_extractor$to(device = device)
    nn_out_head$to(device = device)

    opt_hyb <- optim_adam(
      c(nn_feat_extractor$parameters, nn_out_head$parameters),
      lr = 5e-5
    )

    dl_hyb <- make_dl(Xtr_ens_sc, y_tr_sc, batch_size = 512)

    nn_feat_extractor$train()
    nn_out_head$train()

    for (epoch in seq_len(150)) {
      coro::loop(for (b in dl_hyb) {
        x_batch <- b[[1]]
        y_batch <- b[[2]]
        opt_hyb$zero_grad()
        feats <- nn_feat_extractor(x_batch)
        preds <- nn_out_head(feats)
        loss  <- pinball_loss(preds, y_batch, tau_t)
        loss$backward()
        opt_hyb$step()
      })
      if (epoch %% 50 == 0) {
        cat(sprintf("    Epoch %d, Loss: %.4f\n", epoch, as.numeric(loss$item())))
      }
    }
    cat("  Feature extractor training complete.\n\n")

    cat("Stage 2: Extracting 32-dim neural features...\n")
    nn_feat_extractor$eval()

    X_tr_t <- torch_tensor(as.matrix(Xtr_ens_sc), dtype = torch_float(), device = device)
    X_te_t <- torch_tensor(as.matrix(Xte_ens_sc), dtype = torch_float(), device = device)

    with_no_grad({
      neural_tr <- as.matrix(nn_feat_extractor(X_tr_t)$cpu())
      neural_te <- as.matrix(nn_feat_extractor(X_te_t)$cpu())
    })

    colnames(neural_tr) <- colnames(neural_te) <- paste0("nn_feat_", seq_len(ncol(neural_tr)))
    df_neural_tr <- as.data.frame(neural_tr)
    df_neural_te <- as.data.frame(neural_te)

    cat(sprintf("  Features: %d train x %d dims,  %d test x %d dims\n",
                nrow(df_neural_tr), ncol(df_neural_tr), nrow(df_neural_te), ncol(df_neural_te)))

    cat("Stage 3: Training DRF on neural features...\n")
    hybrid_drf <- drf(
      X               = df_neural_tr,
      Y               = y_tr,
      num.trees       = num_trees_drf,
      min.node.size   = min_node_size_drf,
      mtry            = floor(sqrt(ncol(df_neural_tr))),
      sample.fraction = 0.5,
      seed            = 2026
    )

    cat("Stage 4: Computing Hybrid predictions (wquant_matrix)...\n")
    hybrid_pred      <- predict(hybrid_drf, newdata = df_neural_te)
    hybrid_qgrid_mat <- wquant_matrix(as.numeric(hybrid_pred$y), hybrid_pred$weights, q_grid)
    lower_hyb        <- min(y_tr)
    hybrid_qgrid_mat <- pmax(hybrid_qgrid_mat, lower_hyb)

    hybrid_q <- data.frame(
      q05 = hybrid_qgrid_mat[, q05_idx],
      q50 = hybrid_qgrid_mat[, q50_idx],
      q95 = hybrid_qgrid_mat[, q95_idx]
    )

    hybrid_qgrid <- as.data.frame(hybrid_qgrid_mat)
    colnames(hybrid_qgrid) <- paste0("q", sprintf("%02d", round(100 * q_grid)))
    cat("Hybrid NN+DRF complete.\n\n")
  }, envir = env)
}

# ---- Hybrid NN+DRF: extract features with loaded extractor, refit DRF, predict 39 levels.
# Source: 39quantiles.R:733-753
stage_predict39_hybrid <- function(env = .GlobalEnv) {
  local({
    cat("  Hybrid NN+DRF...")
    nn_feat_extractor$eval()
    X_tr_t <- torch_tensor(as.matrix(Xtr_ens_sc), dtype = torch_float(), device = device)
    with_no_grad({
      neural_tr <- as.matrix(nn_feat_extractor(X_tr_t)$cpu())
      neural_te <- as.matrix(nn_feat_extractor(X_te_t)$cpu())
    })
    colnames(neural_tr) <- colnames(neural_te) <- paste0("nn_feat_", seq_len(ncol(neural_tr)))
    df_neural_tr <- as.data.frame(neural_tr)
    df_neural_te <- as.data.frame(neural_te)
    hybrid_drf <- drf(
      X = df_neural_tr, Y = y_tr,
      num.trees     = best_hyperparameters$DRF$num.trees,
      min.node.size = best_hyperparameters$DRF$min.node.size,
      mtry = floor(sqrt(ncol(df_neural_tr))), sample.fraction = 0.5, seed = 2026
    )
    hybrid_pred <- predict(hybrid_drf, newdata = df_neural_te)
    hybrid_qmat <- pmax(wquant_matrix(as.numeric(hybrid_pred$y),
                                      hybrid_pred$weights, q_grid),
                        min(y_tr, na.rm = TRUE))
    cat(" ✓\n")
  }, envir = env)
}
