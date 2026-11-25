# ============================================================
# Distributional Deep Neural Network (Gamma).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- DDNN (Gamma, log_mu/log_var): per-run tuning, 4 runs, qEns averaging.
# Source: scripts.R:1286-1478
stage_train_ddnn <- function(env = .GlobalEnv) {
  local({
    cat("========================================\n")
    cat("MODEL 8: DDNN (Gamma, Marcjasz et al. 2023)\n")
    cat("========================================\n")

    # ---- DDNN-specific Gamma helpers (log_mu / log_var parametrisation) ----
    gamma_nll_loss_ddnn <- function(params, y) {
      log_mu  <- params[, 1, drop = FALSE]
      log_var <- torch_clamp(params[, 2, drop = FALSE], min = -3, max = 4)
      mu    <- torch_exp(log_mu)
      v     <- torch_exp(log_var)
      shape <- mu^2 / v
      rate  <- mu   / v
      nll <- -(shape - 1) * torch_log(y) +
        y * rate -
        shape * torch_log(rate) +
        torch_lgamma(shape)
      torch_mean(nll)
    }

    get_gamma_params_ddnn <- function(net, X_sc) {
      X_sc <- as.matrix(X_sc)
      net$eval()
      Xv <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      with_no_grad({ p <- as.matrix(net(Xv)$cpu()) })
      net$train()
      log_mu  <- p[, 1]
      log_var <- pmin(pmax(p[, 2], -3), 4)
      mu  <- exp(log_mu)
      v   <- exp(log_var)
      list(shape = mu^2 / v, rate = mu / v)
    }

    # Network: outputs 2 params (log_mu, log_var)
    make_ddnn_net <- function(n_in, h1, h2, dropout_rate = 0.0) {
      nn_sequential(
        nn_linear(n_in, h1), nn_softplus(),
        if (dropout_rate > 0) nn_dropout(dropout_rate) else nn_identity(),
        nn_linear(h1,   h2), nn_softplus(),
        nn_linear(h2,    2)
      )
    }

    # Training with early stopping on val NLL
    train_ddnn <- function(X_tr, y_tr_in, X_va, y_va, h1, h2, lr,
                           dropout_rate = 0.0, l1_hidden = 0.0,
                           max_epochs = 200, patience = 20,
                           batch_size = 512, seed = 42) {
      set.seed(seed)
      torch_manual_seed(seed)
      net   <- make_ddnn_net(n_input, h1, h2, dropout_rate)$to(device = device)
      optim <- optim_adam(net$parameters, lr = lr)
      dl    <- make_dl(X_tr, y_tr_in, batch_size)

      Xva_t <- torch_tensor(as.matrix(X_va), dtype = torch_float(), device = device)
      yva_t <- torch_tensor(matrix(as.numeric(y_va), ncol = 1), dtype = torch_float(), device = device)

      best_val_nll <- Inf; best_state <- NULL; patience_cnt <- 0

      net$train()
      for (ep in seq_len(max_epochs)) {
        coro::loop(for (b in dl) {
          x_batch <- b[[1]]
          y_batch <- b[[2]]
          optim$zero_grad()
          params <- net(x_batch)
          nll    <- gamma_nll_loss_ddnn(params, y_batch)   # 使用新loss

          l1_pen <- torch_tensor(0.0, dtype = torch_float(), device = device)
          if (l1_hidden > 0) {
            for (nm in names(net$parameters)) {
              if (grepl("^[0-9]+\\.weight$", nm) && !grepl("^[56789]", nm)) {
                l1_pen <- l1_pen + l1_hidden * torch_sum(torch_abs(net$parameters[[nm]]))
              }
            }
          }
          loss <- nll + l1_pen
          loss$backward()
          optim$step()
        })

        net$eval()
        with_no_grad({
          va_params <- net(Xva_t)
          va_nll    <- as.numeric(gamma_nll_loss_ddnn(va_params, yva_t))
        })
        net$train()

        if (ep %% 50 == 0 || ep == 1) {
          cat(sprintf("    [ep %d] val NLL=%.4f (best=%.4f, patience=%d/%d)\n",
                      ep, va_nll, best_val_nll, patience_cnt, patience))
        }
        if (va_nll < best_val_nll - 1e-6) {
          best_val_nll <- va_nll
          best_state   <- lapply(net$state_dict(), function(p) p$clone())
          patience_cnt <- 0
        } else {
          patience_cnt <- patience_cnt + 1
          if (patience_cnt >= patience) {
            cat(sprintf("    Early stop at epoch %d (best val NLL=%.4f)\n", ep, best_val_nll))
            break
          }
        }
      }
      if (!is.null(best_state)) net$load_state_dict(best_state)
      net$eval()
      net
    }

    # ---- HP TUNING: DDNN — CRPS-based (req. 2) ----
    # Grid: h1 x h2 (two options each, h2 = h1/2) x lr x l1 = 8 configs per run
    cat("DDNN: Tuning HPs per run under Gamma NLL (log_mu/log_var), evaluated by CRPS...\n")
    ddnn_hp_search <- expand.grid(
      h1        = c(128, 256),
      h2        = c(64,  128),
      lr        = c(5e-4, 1e-3),
      l1_hidden = c(0.0, 1e-4),
      stringsAsFactors = FALSE
    )
    ddnn_hp_search <- ddnn_hp_search[ddnn_hp_search$h2 == ddnn_hp_search$h1 / 2, ]
    cat(sprintf("  DDNN HP grid: %d configs per run\n", nrow(ddnn_hp_search)))

    tune_ddnn_run <- function(run_id) {
      crps_run <- sapply(seq_len(nrow(ddnn_hp_search)), function(i) {
        hp <- ddnn_hp_search[i, ]
        tryCatch({
          net <- train_ddnn(X_tr_cv, y_tr_cv, X_va_cv, y_va_cv,
                            h1 = hp$h1, h2 = hp$h2, lr = hp$lr,
                            l1_hidden = hp$l1_hidden,
                            max_epochs = 150, patience = 20,
                            batch_size = 512, seed = 2026 + run_id * 31 + i)
          p  <- get_gamma_params_ddnn(net, X_va_cv)
          qm <- gamma_qgrid_fn(p$shape, p$rate, q_grid)
          cr <- crps_from_qmat(y_va_cv, qm)
          cr
        }, error = function(e) Inf)
      })
      best_i <- which.min(crps_run)
      cat(sprintf("  Run %d best: h1=%d h2=%d lr=%.0e l1=%.0e (CRPS=%.4f)\n",
                  run_id,
                  ddnn_hp_search$h1[best_i], ddnn_hp_search$h2[best_i],
                  ddnn_hp_search$lr[best_i], ddnn_hp_search$l1_hidden[best_i],
                  crps_run[best_i]))
      as.list(ddnn_hp_search[best_i, ])
    }

    N_runs_ddnn   <- 4
    ddnn_best_hps <- lapply(seq_len(N_runs_ddnn), tune_ddnn_run)

    # Final training: full training set with each run's best HP
    ddnn_runs <- vector("list", N_runs_ddnn)
    cat(sprintf("\nTraining %d DDNN runs on full train set with tuned HPs...\n", N_runs_ddnn))

    for (run in seq_len(N_runs_ddnn)) {
      hp <- ddnn_best_hps[[run]]
      cat(sprintf("\n--- DDNN Run %d/%d: h1=%d h2=%d lr=%.0e l1=%.0e ---\n",
                  run, N_runs_ddnn, hp$h1, hp$h2, hp$lr, hp$l1_hidden))
      set.seed(2026 + run * 31)
      ddnn_runs[[run]] <- list(
        net = train_ddnn(Xtr_ens_sc, y_tr,
                         X_va_cv, y_va_cv,
                         h1 = hp$h1, h2 = hp$h2, lr = hp$lr,
                         dropout_rate = 0.0, l1_hidden = hp$l1_hidden,
                         max_epochs = 300, patience = 30,
                         batch_size = 512, seed = 2026 + run * 31),
        best_hp = hp
      )
    }

    # qEns: horizontal quantile averaging across 4 runs (Gamma quantiles)
    cat("\nComputing DDNN qEns (horizontal quantile averaging across 4 runs)...\n")
    ddnn_qgrid_runs <- lapply(ddnn_runs, function(r) {
      p <- get_gamma_params_ddnn(r$net, Xte_ens_sc)
      gamma_qgrid_fn(p$shape, p$rate, q_grid)
    })
    ddnn_qgrid_mat <- Reduce(`+`, ddnn_qgrid_runs) / N_runs_ddnn
    ddnn_qgrid_mat <- pmax(ddnn_qgrid_mat, 1)

    ddnn_q <- data.frame(
      q05 = ddnn_qgrid_mat[, q05_idx],
      q50 = ddnn_qgrid_mat[, q50_idx],
      q95 = ddnn_qgrid_mat[, q95_idx]
    )
    ddnn_qgrid <- as.data.frame(ddnn_qgrid_mat)
    colnames(ddnn_qgrid) <- paste0("q", sprintf("%02d", round(100 * q_grid)))

    # NLL (Gamma parametric, averaged parameters)
    ddnn_params_all <- lapply(ddnn_runs, function(r) get_gamma_params_ddnn(r$net, Xte_ens_sc))
    ddnn_shape_avg  <- Reduce(`+`, lapply(ddnn_params_all, `[[`, "shape")) / N_runs_ddnn
    ddnn_rate_avg   <- Reduce(`+`, lapply(ddnn_params_all, `[[`, "rate"))  / N_runs_ddnn
    nll_ddnn        <- -mean(dgamma(y_te, shape = ddnn_shape_avg, rate = ddnn_rate_avg, log = TRUE))

    best_ddnn_hp <- lapply(ddnn_runs, `[[`, "best_hp")
    cat(sprintf("DDNN (Gamma, qEns, %d runs) complete.\n\n", N_runs_ddnn))
  }, envir = env)
}

# ---- DDNN: 39-level qEns prediction from loaded runs.
# Source: 39quantiles.R:669-675
stage_predict39_ddnn <- function(env = .GlobalEnv) {
  local({
    cat("  DDNN Gamma (4 runs)...")
    ddnn_qmat_runs <- lapply(ddnn_nets, function(net) {
      p <- get_gamma_params_ddnn(net, Xte_ens_sc)
      gamma_qgrid_fn(p$shape, p$rate, q_grid)
    })
    ddnn_qmat <- pmax(Reduce(`+`, ddnn_qmat_runs) / N_runs_ddnn, 1)
    cat(" ✓\n")
  }, envir = env)
}
