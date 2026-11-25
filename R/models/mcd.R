# ============================================================
# Monte Carlo Dropout (Gamma).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- MC Dropout (Gamma): tuning, final net, MC inference.
# Source: scripts.R:1183-1275
stage_train_mcd <- function(env = .GlobalEnv) {
  local({
    cat("========================================\n")
    cat("MODEL 7: MC Dropout (Gamma parametric, req. 6)\n")
    cat("========================================\n")



    # MCD network architecture: outputs 2 Gamma params
    make_mcd_gamma_net <- function(n_in, h1 = 128, h2 = 64, p_drop = 0.1) {
      nn_sequential(
        nn_linear(n_in, h1), nn_relu(), nn_dropout(p_drop),
        nn_linear(h1,   h2), nn_relu(), nn_dropout(p_drop),
        nn_linear(h2,   2)   # log_shape, log_rate
      )
    }

    train_mcd_gamma <- function(X_sc, y_obs, p_drop, lr,
                                h1 = 128, h2 = 64, epochs = 60, batch_size = 512) {
      net   <- make_mcd_gamma_net(n_input, h1, h2, p_drop)$to(device = device)
      optim <- optim_adam(net$parameters, lr = lr)
      dl    <- make_dl(X_sc, y_obs, batch_size)
      net$train()
      for (ep in seq_len(epochs)) {
        coro::loop(for (b in dl) {
          x_batch <- b[[1]]
          y_batch <- b[[2]]
          optim$zero_grad()
          params <- net(x_batch)
          loss   <- gamma_nll_loss(params, y_batch)
          loss$backward()
          optim$step()
        })
      }
      net  # keep in train mode for MC dropout inference
    }

    # MC dropout inference: average Gamma parameters over M passes, then quantiles
    mcd_gamma_quantiles <- function(net, X_sc, M_pass = 50) {
      X_sc <- as.matrix(X_sc)
      net$train()
      X_t   <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      acc_p <- matrix(0, nrow(X_sc), 2)
      with_no_grad({
        for (p in seq_len(M_pass)) {
          acc_p <- acc_p + as.matrix(net(X_t)$cpu())
        }
      })
      avg_params <- acc_p / M_pass
      shape <- exp(avg_params[, 1])
      rate  <- exp(avg_params[, 2])
      gamma_qgrid_fn(shape, rate, q_grid)
    }

    # ---- HP TUNING: MCD Gamma — CRPS-based (req. 2) ----
    cat("MCD: Tuning HPs under Gamma NLL using CRPS (6 configs)...\n")
    mcd_hp_grid <- expand.grid(
      p_drop = c(0.10, 0.20),
      lr     = c(1e-4, 5e-4, 1e-3),
      stringsAsFactors = FALSE
    )

    mcd_cv_crps <- sapply(seq_len(nrow(mcd_hp_grid)), function(i) {
      hp <- mcd_hp_grid[i, ]
      cat(sprintf("  [%d/6] dropout=%.2f lr=%.0e ... ", i, hp$p_drop, hp$lr))
      tryCatch({
        net <- train_mcd_gamma(X_tr_cv, y_tr_cv, hp$p_drop, hp$lr, epochs = 60)
        qm  <- mcd_gamma_quantiles(net, X_va_cv, M_pass = 20)
        cr  <- crps_from_qmat(y_va_cv, qm)
        cat(sprintf("CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat("ERROR\n"); Inf })
    })

    best_mcd_i  <- which.min(mcd_cv_crps)
    best_mcd_hp <- as.list(mcd_hp_grid[best_mcd_i, ])
    cat(sprintf("✓ Best MCD HP: dropout=%.2f lr=%.0e (CRPS=%.4f)\n",
                best_mcd_hp$p_drop, best_mcd_hp$lr, mcd_cv_crps[best_mcd_i]))

    set.seed(2026)
    cat("Training final MCD model (150 epochs)...\n")
    mc_dropout_net <- train_mcd_gamma(Xtr_ens_sc, y_tr,
                                      best_mcd_hp$p_drop, best_mcd_hp$lr,
                                      epochs = 150)

    M_dropout     <- 100
    cat(sprintf("Generating %d MC passes on test set...\n", M_dropout))
    mcd_qgrid_mat <- mcd_gamma_quantiles(mc_dropout_net, Xte_ens_sc, M_pass = M_dropout)
    mcd_qgrid_mat <- pmax(mcd_qgrid_mat, 1)

    mcd_q <- data.frame(q05 = mcd_qgrid_mat[, q05_idx],
                        q50 = mcd_qgrid_mat[, q50_idx],
                        q95 = mcd_qgrid_mat[, q95_idx])
    mcd_qgrid <- as.data.frame(mcd_qgrid_mat)
    colnames(mcd_qgrid) <- paste0("q", sprintf("%02d", round(100 * q_grid)))
    cat("MCD (Gamma) complete.\n\n")
  }, envir = env)
}

# ---- MC Dropout: 100 stochastic passes on the 39-level grid.
# Source: 39quantiles.R:660-667
stage_predict39_mcd <- function(env = .GlobalEnv) {
  local({
    cat("  MCD Gamma (100 passes)...")
    mc_dropout_net$train()
    acc_mcd <- matrix(0, nrow(Xte_ens_sc), 2)
    for (i in seq_len(100))
      with_no_grad({ acc_mcd <- acc_mcd + as.matrix(mc_dropout_net(X_te_t)$cpu()) })
    avg_mcd  <- acc_mcd / 100
    mcd_qmat <- pmax(gamma_qgrid_fn(exp(avg_mcd[, 1]), exp(avg_mcd[, 2]), q_grid), 1)
    cat(" ✓\n")
  }, envir = env)
}
