# ============================================================
# Deep Ensemble (Gamma).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Deep Ensemble (Gamma): tuning, 5 members, parameter-space averaging.
# Source: scripts.R:1084-1173
stage_train_ens <- function(env = .GlobalEnv) {
  local({
    cat("\n========================================\n")
    cat("MODEL 6: Deep Ensemble (Gamma parametric)\n")
    cat("========================================\n")

    # Network: outputs 2 Gamma params (log_shape, log_rate)
    make_ens_gamma_net <- function(n_in, h1, h2) {
      nn_sequential(
        nn_linear(n_in, h1), nn_relu(),
        nn_linear(h1,   h2), nn_relu(),
        nn_linear(h2,   2)   # log_shape, log_rate
      )
    }

    train_ens_gamma <- function(X_sc, y_obs, h1, h2, lr,
                                epochs = 60, batch_size = 512) {
      net   <- make_ens_gamma_net(n_input, h1, h2)$to(device = device)
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
      net$eval(); net
    }

    # ---- HP TUNING: ENS Gamma — CRPS-based (req. 2) ----
    cat("ENS: Tuning HPs under Gamma NLL using CRPS (12 configs, 1 CV fold)...\n")
    ens_hp_grid <- expand.grid(
      h1 = c(128, 256),
      h2 = c(32, 64),
      lr = c(1e-4, 5e-4, 1e-3),
      stringsAsFactors = FALSE
    )

    ens_cv_crps <- sapply(seq_len(nrow(ens_hp_grid)), function(i) {
      hp <- ens_hp_grid[i, ]
      cat(sprintf("  [%d/12] h1=%d h2=%d lr=%.0e ... ", i, hp$h1, hp$h2, hp$lr))
      tryCatch({
        net <- train_ens_gamma(X_tr_cv, y_tr_cv, hp$h1, hp$h2, hp$lr, epochs = 60)
        p   <- get_gamma_params(net, X_va_cv)
        qm  <- gamma_qgrid_fn(p$shape, p$rate, q_grid)
        cr  <- crps_from_qmat(y_va_cv, qm)
        cat(sprintf("CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat("ERROR\n"); Inf })
    })

    best_ens_i  <- which.min(ens_cv_crps)
    best_ens_hp <- as.list(ens_hp_grid[best_ens_i, ])
    cat(sprintf("✓ Best ENS HP: h1=%d h2=%d lr=%.0e (CRPS=%.4f)\n",
                best_ens_hp$h1, best_ens_hp$h2, best_ens_hp$lr, ens_cv_crps[best_ens_i]))

    # Train M members on full training set
    M_ens       <- 5
    ens_members <- vector("list", M_ens)
    cat(sprintf("Training %d ENS members (150 epochs each)...\n", M_ens))
    for (m in seq_len(M_ens)) {
      set.seed(2026 + m)
      cat(sprintf("  Member %d/%d\n", m, M_ens))
      ens_members[[m]] <- train_ens_gamma(Xtr_ens_sc, y_tr,
                                          best_ens_hp$h1, best_ens_hp$h2, best_ens_hp$lr,
                                          epochs = 150)
    }

    # Predict: average Gamma parameters across M members (parameter-space ensemble),
    # then convert to quantile grid via qgamma()
    cat("Computing ENS predictions (Gamma parameter averaging across members)...\n")
    ens_params_runs <- lapply(ens_members, function(net) get_gamma_params(net, Xte_ens_sc))
    ens_shape_avg   <- Reduce(`+`, lapply(ens_params_runs, `[[`, "shape")) / M_ens
    ens_rate_avg    <- Reduce(`+`, lapply(ens_params_runs, `[[`, "rate"))  / M_ens
    ens_qgrid_mat   <- gamma_qgrid_fn(ens_shape_avg, ens_rate_avg, q_grid)
    ens_qgrid_mat   <- pmax(ens_qgrid_mat, 1)

    ens_q <- data.frame(q05 = ens_qgrid_mat[, q05_idx],
                        q50 = ens_qgrid_mat[, q50_idx],
                        q95 = ens_qgrid_mat[, q95_idx])
    ens_qgrid <- as.data.frame(ens_qgrid_mat)
    colnames(ens_qgrid) <- paste0("q", sprintf("%02d", round(100 * q_grid)))

    # NLL: Gamma parametric (averaged parameters)
    nll_ens <- -mean(dgamma(y_te, shape = ens_shape_avg, rate = ens_rate_avg, log = TRUE))

    cat("ENS (Gamma) complete.\n\n")
  }, envir = env)
}

# ---- Deep Ensemble: 39-level prediction from loaded members.
# Source: 39quantiles.R:653-658
stage_predict39_ens <- function(env = .GlobalEnv) {
  local({
    cat("  ENS Gamma...")
    ens_params_runs <- lapply(ens_members, function(net) get_gamma_params(net, Xte_ens_sc))
    ens_shape_avg   <- Reduce(`+`, lapply(ens_params_runs, `[[`, "shape")) / M_ens
    ens_rate_avg    <- Reduce(`+`, lapply(ens_params_runs, `[[`, "rate"))  / M_ens
    ens_qmat        <- pmax(gamma_qgrid_fn(ens_shape_avg, ens_rate_avg, q_grid), 1)
    cat(" ✓\n")
  }, envir = env)
}
