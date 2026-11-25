# ============================================================
# Recovery of torch models (retrain with saved hyper-parameters) and torch_save().
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Recovery: seed, 19-level q_grid, crps_from_qmat, quantile indices.
# Source: script_recovery.R:195-196; script_recovery.R:215-215; script_recovery.R:239-243
stage_recovery_helpers <- function(env = .GlobalEnv) {
  local({
    set.seed(2026)
    q_grid <- seq(0.05, 0.95, by = 0.05)

    crps_from_qmat <- function(y, qmat) crps_from_quantiles(y, as.matrix(qmat), q_grid)

    q05_idx <- which.min(abs(q_grid - 0.05))
    q50_idx <- which.min(abs(q_grid - 0.50))
    q95_idx <- which.min(abs(q_grid - 0.95))

    cat("✓ Shared helpers defined\n")
  }, envir = env)
}

# ---- Recovery: restore non-torch models and best hyper-parameters from the rds.
# Source: script_recovery.R:246-269
stage_recovery_load_rds <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # STEP 4: Load .rds — restore non-torch models + hyperparameters
    # ============================================================
    cat(sprintf("\nLoading RDS from: %s\n", rds_path))
    results <- readRDS(rds_path)

    # Non-torch models: restore directly (these work fine with saveRDS)
    drf_fit             <- results$drf_fit
    qrf_fit             <- results$qrf_fit
    rf_q_fit            <- results$rf_q_fit
    xgb_quantile_models <- results$xgb_quantile_models
    drf_base_full       <- results$edrf_members   # Ensemble DRF base learners

    # Hyperparameters: used to retrain torch models without tuning
    best_hyperparameters <- results$best_hyperparams

    # Summary / CRPS from original run (for sanity check later)
    summary_all <- results$summary
    crps_vec    <- results$crps

    cat("✓ Non-torch models and hyperparameters restored from RDS\n")
    cat("  Saved hyperparameters:\n")
    print(best_hyperparameters)
  }, envir = env)
}

# ---- Recovery: neural features, network definitions and train/predict functions.
# Source: script_recovery.R:272-541
stage_recovery_torch_infra <- function(env = .GlobalEnv) {
  local({
    # STEP 5: Rebuild shared torch infrastructure
    # ============================================================
    device <- torch_device(select_torch_device())
    cat("\nNeural models using device:", device$type, "\n")

    # One-hot encode + standardise (train-led, identical to original)
    Xtr_ens_mm <- model.matrix(~ . - 1, data = Xtr_all)
    Xte_ens_mm <- model.matrix(~ . - 1, data = Xte_all)
    common_nn  <- union(colnames(Xtr_ens_mm), colnames(Xte_ens_mm))
    Xtr_ens_mm <- Xtr_ens_mm[, common_nn, drop=FALSE]; Xtr_ens_mm[is.na(Xtr_ens_mm)] <- 0
    Xte_ens_mm <- Xte_ens_mm[, common_nn, drop=FALSE]; Xte_ens_mm[is.na(Xte_ens_mm)] <- 0

    x_mean     <- colMeans(Xtr_ens_mm)
    x_sd       <- apply(Xtr_ens_mm, 2, sd); x_sd[x_sd == 0] <- 1
    Xtr_ens_sc <- scale(Xtr_ens_mm, center = x_mean, scale = x_sd)
    Xte_ens_sc <- scale(Xte_ens_mm, center = x_mean, scale = x_sd)
    Xtr_ens_sc[Xtr_ens_sc >  10] <-  10; Xtr_ens_sc[Xtr_ens_sc < -10] <- -10
    Xte_ens_sc[Xte_ens_sc >  10] <-  10; Xte_ens_sc[Xte_ens_sc < -10] <- -10

    n_input     <- ncol(Xtr_ens_sc)
    n_quantiles <- length(q_grid)
    tau_t       <- torch_tensor(matrix(q_grid, nrow=1), dtype=torch_float(), device=device)

    # CV fold (needed by DDNN train_ddnn for early stopping val set)
    set.seed(2026)
    cv_folds_hp  <- caret::createFolds(y_tr, k=3, list=TRUE, returnTrain=FALSE)
    val_idx_nn   <- cv_folds_hp[[1]]
    train_idx_nn <- setdiff(seq_len(nrow(Xtr_ens_sc)), val_idx_nn)
    X_va_cv      <- Xtr_ens_sc[val_idx_nn,   , drop=FALSE]
    y_va_cv      <- y_tr[val_idx_nn]

    # ---- Shared torch helpers ----
    make_dl <- function(X_sc, y_obs, batch_size = 512) {
      X_t <- torch_tensor(as.matrix(X_sc), dtype=torch_float(), device=device)
      y_t <- torch_tensor(matrix(as.numeric(y_obs), ncol=1), dtype=torch_float(), device=device)
      ds  <- tensor_dataset(X_t, y_t)
      dataloader(ds, batch_size=batch_size, shuffle=TRUE)
    }

    gamma_nll_loss <- function(params, y_true) {
      log_shape <- params[, 1, drop=FALSE]; log_rate <- params[, 2, drop=FALSE]
      shape <- torch_exp(log_shape);        rate  <- torch_exp(log_rate)
      nll <- -(shape - 1) * torch_log(y_true) + y_true * rate -
        shape * torch_log(rate) + torch_lgamma(shape)
      torch_mean(nll)
    }

    get_gamma_params <- function(net, X_sc) {
      X_sc <- as.matrix(X_sc); net$eval()
      Xv <- torch_tensor(X_sc, dtype=torch_float(), device=device)
      with_no_grad({ p <- as.matrix(net(Xv)$cpu()) })
      net$train()
      list(shape = exp(p[,1]), rate = exp(p[,2]))
    }

    gamma_qgrid_fn <- function(shape, rate, qs)
      t(mapply(function(sh, ra) qgamma(qs, shape=sh, rate=ra), shape, rate))

    # ---- ENS architecture ----
    make_ens_gamma_net <- function(n_in, h1, h2) {
      nn_sequential(nn_linear(n_in, h1), nn_relu(),
                    nn_linear(h1,   h2), nn_relu(),
                    nn_linear(h2,   2))
    }

    train_ens_gamma <- function(X_sc, y_obs, h1, h2, lr, epochs=150, batch_size=512) {
      net   <- make_ens_gamma_net(n_input, h1, h2)$to(device=device)
      optim <- optim_adam(net$parameters, lr=lr)
      dl    <- make_dl(X_sc, y_obs, batch_size)
      net$train()
      for (ep in seq_len(epochs)) {
        coro::loop(for (b in dl) {
          optim$zero_grad()
          loss <- gamma_nll_loss(net(b[[1]]), b[[2]])
          loss$backward(); optim$step()
        })
      }
      net$eval(); net
    }

    # ---- MCD architecture ----
    make_mcd_gamma_net <- function(n_in, h1=128, h2=64, p_drop=0.1) {
      nn_sequential(nn_linear(n_in, h1), nn_relu(), nn_dropout(p_drop),
                    nn_linear(h1,   h2), nn_relu(), nn_dropout(p_drop),
                    nn_linear(h2,   2))
    }

    train_mcd_gamma <- function(X_sc, y_obs, p_drop, lr,
                                h1=128, h2=64, epochs=150, batch_size=512) {
      net   <- make_mcd_gamma_net(n_input, h1, h2, p_drop)$to(device=device)
      optim <- optim_adam(net$parameters, lr=lr)
      dl    <- make_dl(X_sc, y_obs, batch_size)
      net$train()
      for (ep in seq_len(epochs)) {
        coro::loop(for (b in dl) {
          optim$zero_grad()
          loss <- gamma_nll_loss(net(b[[1]]), b[[2]])
          loss$backward(); optim$step()
        })
      }
      net  # stays in train() mode for MC dropout at inference
    }

    mcd_gamma_quantiles <- function(net, X_sc, M_pass=100) {
      X_sc <- as.matrix(X_sc); net$train()
      X_t  <- torch_tensor(X_sc, dtype=torch_float(), device=device)
      acc  <- matrix(0, nrow(X_sc), 2)
      with_no_grad({ for (i in seq_len(M_pass)) acc <- acc + as.matrix(net(X_t)$cpu()) })
      avg  <- acc / M_pass
      gamma_qgrid_fn(exp(avg[,1]), exp(avg[,2]), q_grid)
    }

    # ---- DDNN architecture ----
    gamma_nll_loss_ddnn <- function(params, y) {
      log_mu  <- params[, 1, drop=FALSE]
      log_var <- torch_clamp(params[, 2, drop=FALSE], min=-3, max=4)
      mu <- torch_exp(log_mu); v <- torch_exp(log_var)
      shape <- mu^2 / v; rate <- mu / v
      nll <- -(shape - 1) * torch_log(y) + y * rate -
        shape * torch_log(rate) + torch_lgamma(shape)
      torch_mean(nll)
    }

    get_gamma_params_ddnn <- function(net, X_sc) {
      X_sc <- as.matrix(X_sc); net$eval()
      Xv <- torch_tensor(X_sc, dtype=torch_float(), device=device)
      with_no_grad({ p <- as.matrix(net(Xv)$cpu()) })
      net$train()
      log_mu <- p[,1]; log_var <- pmin(pmax(p[,2], -3), 4)
      mu <- exp(log_mu); v <- exp(log_var)
      list(shape = mu^2/v, rate = mu/v)
    }

    make_ddnn_net <- function(n_in, h1, h2, dropout_rate=0.0) {
      nn_sequential(nn_linear(n_in, h1), nn_softplus(),
                    if (dropout_rate > 0) nn_dropout(dropout_rate) else nn_identity(),
                    nn_linear(h1, h2), nn_softplus(),
                    nn_linear(h2, 2))
    }

    train_ddnn <- function(X_tr, y_tr_in, X_va, y_va, h1, h2, lr,
                           dropout_rate=0.0, l1_hidden=0.0,
                           max_epochs=300, patience=30, batch_size=512, seed=42) {
      set.seed(seed); torch_manual_seed(seed)
      net   <- make_ddnn_net(n_input, h1, h2, dropout_rate)$to(device=device)
      optim <- optim_adam(net$parameters, lr=lr)
      dl    <- make_dl(X_tr, y_tr_in, batch_size)
      Xva_t <- torch_tensor(as.matrix(X_va), dtype=torch_float(), device=device)
      yva_t <- torch_tensor(matrix(as.numeric(y_va), ncol=1), dtype=torch_float(), device=device)
      best_val_nll <- Inf; best_state <- NULL; patience_cnt <- 0
      net$train()
      for (ep in seq_len(max_epochs)) {
        coro::loop(for (b in dl) {
          optim$zero_grad()
          params  <- net(b[[1]])
          l1_pen  <- torch_tensor(0.0, dtype=torch_float(), device=device)
          if (l1_hidden > 0) {
            for (nm in names(net$parameters))
              if (grepl("^[0-9]+\\.weight$", nm) && !grepl("^[56789]", nm))
                l1_pen <- l1_pen + l1_hidden * torch_sum(torch_abs(net$parameters[[nm]]))
          }
          loss <- gamma_nll_loss_ddnn(params, b[[2]]) + l1_pen
          loss$backward(); optim$step()
        })
        net$eval()
        with_no_grad({ va_nll <- as.numeric(gamma_nll_loss_ddnn(net(Xva_t), yva_t)) })
        net$train()
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
      net$eval(); net
    }

    # ---- BNN architecture ----
    torch_log_scale_mixture_prior <- function(w, pi_mix, log_sigma1, log_sigma2) {
      half_log2pi <- 0.5 * log(2 * pi)
      log_p1 <- log(pi_mix)     - half_log2pi - log_sigma1 - 0.5*(w/exp(log_sigma1))^2
      log_p2 <- log(1-pi_mix)   - half_log2pi - log_sigma2 - 0.5*(w/exp(log_sigma2))^2
      m <- torch_maximum(log_p1, log_p2)
      torch_sum(m + torch_log(torch_exp(log_p1-m) + torch_exp(log_p2-m)))
    }

    bnn_linear <- nn_module(
      classname  = "bnn_linear_bbb",
      initialize = function(in_f, out_f, pi_mix, log_sigma1, log_sigma2) {
        self$in_f <- in_f; self$out_f <- out_f
        self$pi_mix <- pi_mix; self$log_sig1 <- log_sigma1; self$log_sig2 <- log_sigma2
        self$w_mu  <- nn_parameter(torch_randn(out_f, in_f) * 0.01)
        self$w_rho <- nn_parameter(torch_full(c(out_f, in_f), -3.0))
        self$b_mu  <- nn_parameter(torch_zeros(out_f))
        self$b_rho <- nn_parameter(torch_full(c(out_f), -3.0))
      },
      forward = function(x) {
        w_sigma <- torch_log1p(torch_exp(self$w_rho))
        b_sigma <- torch_log1p(torch_exp(self$b_rho))
        self$w_sample <- self$w_mu + w_sigma * torch_randn_like(w_sigma)
        self$b_sample <- self$b_mu + b_sigma * torch_randn_like(b_sigma)
        nnf_linear(x, self$w_sample, self$b_sample)
      },
      kl_term = function() {
        w_sigma <- torch_log1p(torch_exp(self$w_rho))
        b_sigma <- torch_log1p(torch_exp(self$b_rho))
        log_qw <- torch_sum(-0.5*log(2*pi) - torch_log(w_sigma) -
                              0.5*((self$w_sample - self$w_mu)/w_sigma)^2)
        log_qb <- torch_sum(-0.5*log(2*pi) - torch_log(b_sigma) -
                              0.5*((self$b_sample - self$b_mu)/b_sigma)^2)
        log_q  <- log_qw + log_qb
        log_pw <- torch_log_scale_mixture_prior(self$w_sample$detach(), self$pi_mix, self$log_sig1, self$log_sig2)
        log_pb <- torch_log_scale_mixture_prior(self$b_sample$detach(), self$pi_mix, self$log_sig1, self$log_sig2)
        log_q - (log_pw + log_pb)
      }
    )

    bnn_gamma_net <- nn_module(
      classname  = "bnn_gamma_bbb",
      initialize = function(n_in, h1, h2, pi_mix, log_sigma1, log_sigma2) {
        self$l1 <- bnn_linear(n_in, h1, pi_mix, log_sigma1, log_sigma2)
        self$l2 <- bnn_linear(h1,   h2, pi_mix, log_sigma1, log_sigma2)
        self$l3 <- bnn_linear(h2,    2, pi_mix, log_sigma1, log_sigma2)
      },
      forward = function(x) { x <- nnf_relu(self$l1(x)); x <- nnf_relu(self$l2(x)); self$l3(x) },
      kl      = function() { self$l1$kl_term() + self$l2$kl_term() + self$l3$kl_term() }
    )

    train_bnn_gamma <- function(X_sc, y_obs, h1, h2, lr,
                                pi_mix, log_sigma1, log_sigma2,
                                epochs=150, batch_size=128, seed=42) {
      set.seed(seed); torch_manual_seed(seed)
      N   <- nrow(X_sc)
      net <- bnn_gamma_net(n_input, h1, h2, pi_mix, log_sigma1, log_sigma2)$to(device=device)
      opt <- optim_adam(net$parameters, lr=lr)
      dl  <- make_dl(X_sc, y_obs, batch_size)
      n_batches  <- ceiling(N / batch_size)
      kl_weights <- sapply(seq_len(n_batches), function(i) (2^(n_batches-i))/(2^n_batches-1))
      net$train()
      for (ep in seq_len(epochs)) {
        batch_cnt <- 0L
        coro::loop(for (b in dl) {
          batch_cnt <- batch_cnt + 1L
          pi_i <- kl_weights[min(batch_cnt, n_batches)]
          opt$zero_grad()
          params <- net(b[[1]])
          loss   <- pi_i * net$kl() + gamma_nll_loss(params, b[[2]])
          loss$backward(); opt$step()
        })
        if (ep %% 50 == 0) cat(sprintf("    BNN Epoch %d/%d\n", ep, epochs))
      }
      net$eval(); net
    }

    bnn_gamma_quantiles <- function(net, X_sc, M_pass=100) {
      X_sc <- as.matrix(X_sc); net$train()
      X_t  <- torch_tensor(X_sc, dtype=torch_float(), device=device)
      acc  <- matrix(0, nrow(X_sc), 2)
      for (i in seq_len(M_pass)) with_no_grad({ acc <- acc + as.matrix(net(X_t)$cpu()) })
      avg  <- acc / M_pass
      gamma_qgrid_fn(exp(avg[,1]), exp(avg[,2]), q_grid)
    }

    cat("✓ Torch infrastructure and all model architectures defined\n")
  }, envir = env)
}

# ---- Recovery: re-predict DRF / QRF / RF / XGBoost.
# Source: script_recovery.R:547-590
stage_recovery_trees <- function(env = .GlobalEnv) {
  local({
    cat("\nRebuilding all_qgrid...\n")

    # ---- DRF ----
    cat("  DRF...")
    drf_pred      <- predict(drf_fit, Xte_tree)
    drf_qgrid_mat <- wquant_matrix(as.numeric(drf_pred$y), drf_pred$weights, q_grid)
    drf_qgrid     <- as.data.frame(drf_qgrid_mat)
    colnames(drf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    drf_q <- data.frame(q05=drf_qgrid_mat[,q05_idx], q50=drf_qgrid_mat[,q50_idx], q95=drf_qgrid_mat[,q95_idx])
    cat(" ✓\n")

    # ---- QRF ----
    cat("  QRF...")
    qrf_qgrid_mat <- as.matrix(predict(qrf_fit, Xte_all, what=q_grid))
    qrf_qgrid     <- as.data.frame(qrf_qgrid_mat)
    colnames(qrf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    qrf_q <- data.frame(q05=qrf_qgrid_mat[,q05_idx], q50=qrf_qgrid_mat[,q50_idx], q95=qrf_qgrid_mat[,q95_idx])
    cat(" ✓\n")

    # ---- Ranger RF ----
    cat("  Ranger RF...")
    rf_qgrid_mat <- as.matrix(predict(rf_q_fit, Xte_tree, type="quantiles", quantiles=q_grid)$predictions)
    rf_qgrid     <- as.data.frame(rf_qgrid_mat)
    colnames(rf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    rf_q <- data.frame(q05=rf_qgrid_mat[,q05_idx], q50=rf_qgrid_mat[,q50_idx], q95=rf_qgrid_mat[,q95_idx])
    cat(" ✓\n")

    # ---- XGBoost ----
    cat("  XGBoost...")
    Xtr_all_mm  <- model.matrix(~ . - 1, data = Xtr_all)
    Xte_all_mm  <- model.matrix(~ . - 1, data = Xte_all)
    common_cols <- union(colnames(Xtr_all_mm), colnames(Xte_all_mm))
    Xtr_all_mm  <- Xtr_all_mm[, common_cols, drop=FALSE]; Xtr_all_mm[is.na(Xtr_all_mm)] <- 0
    Xte_all_mm  <- Xte_all_mm[, common_cols, drop=FALSE]; Xte_all_mm[is.na(Xte_all_mm)] <- 0
    dtest_full  <- xgb.DMatrix(Xte_all_mm)
    xgb_qgrid_mat <- matrix(NA_real_, nrow(Xte_all_mm), length(q_grid))
    for (ki in seq_along(q_grid))
      xgb_qgrid_mat[, ki] <- predict(xgb_quantile_models[[ki]], dtest_full)
    xgb_qgrid_mat <- t(apply(xgb_qgrid_mat, 1, sort))
    xgb_qgrid_mat <- pmax(xgb_qgrid_mat, min(y_tr, na.rm=TRUE))
    xgb_qgrid     <- as.data.frame(xgb_qgrid_mat)
    colnames(xgb_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    xgb_q <- data.frame(q05=xgb_qgrid_mat[,q05_idx], q50=xgb_qgrid_mat[,q50_idx], q95=xgb_qgrid_mat[,q95_idx])
    cat(" ✓\n")
  }, envir = env)
}

# ---- Recovery: retrain 5 ENS members with saved HPs.
# Source: script_recovery.R:592-611
stage_recovery_ens <- function(env = .GlobalEnv) {
  local({
    # ---- ENS Gamma — retrain with saved HPs (no tuning) ----
    cat("  ENS Gamma (retraining 5 members with saved HPs)...\n")
    hp_ens  <- best_hyperparameters$ENS_Gamma   # list(h1, h2, lr)
    M_ens   <- 5
    ens_members <- vector("list", M_ens)
    for (m in seq_len(M_ens)) {
      set.seed(2026 + m)
      cat(sprintf("    Member %d/%d\n", m, M_ens))
      ens_members[[m]] <- train_ens_gamma(Xtr_ens_sc, y_tr,
                                          hp_ens$h1, hp_ens$h2, hp_ens$lr,
                                          epochs=150)
    }
    ens_params_runs <- lapply(ens_members, function(net) get_gamma_params(net, Xte_ens_sc))
    ens_shape_avg   <- Reduce(`+`, lapply(ens_params_runs, `[[`, "shape")) / M_ens
    ens_rate_avg    <- Reduce(`+`, lapply(ens_params_runs, `[[`, "rate"))  / M_ens
    ens_qgrid_mat   <- pmax(gamma_qgrid_fn(ens_shape_avg, ens_rate_avg, q_grid), 1)
    ens_qgrid       <- as.data.frame(ens_qgrid_mat)
    colnames(ens_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    ens_q <- data.frame(q05=ens_qgrid_mat[,q05_idx], q50=ens_qgrid_mat[,q50_idx], q95=ens_qgrid_mat[,q95_idx])
    cat("  ENS Gamma ✓\n")
  }, envir = env)
}

# ---- Recovery: retrain MCD with saved HPs.
# Source: script_recovery.R:613-624
stage_recovery_mcd <- function(env = .GlobalEnv) {
  local({
    # ---- MCD Gamma — retrain with saved HPs ----
    cat("  MCD Gamma (retraining with saved HPs)...\n")
    hp_mcd <- best_hyperparameters$MCD_Gamma    # list(p_drop, lr)
    set.seed(2026)
    mc_dropout_net <- train_mcd_gamma(Xtr_ens_sc, y_tr,
                                      hp_mcd$p_drop, hp_mcd$lr,
                                      epochs=150)
    mcd_qgrid_mat  <- pmax(mcd_gamma_quantiles(mc_dropout_net, Xte_ens_sc, M_pass=100), 1)
    mcd_qgrid      <- as.data.frame(mcd_qgrid_mat)
    colnames(mcd_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    mcd_q <- data.frame(q05=mcd_qgrid_mat[,q05_idx], q50=mcd_qgrid_mat[,q50_idx], q95=mcd_qgrid_mat[,q95_idx])
    cat("  MCD Gamma ✓\n")
  }, envir = env)
}

# ---- Recovery: retrain 4 DDNN runs with saved HPs.
# Source: script_recovery.R:626-653
stage_recovery_ddnn <- function(env = .GlobalEnv) {
  local({
    # ---- DDNN — retrain 4 runs with saved HPs ----
    cat("  DDNN (retraining 4 runs with saved HPs)...\n")
    hp_ddnn     <- best_hyperparameters$DDNN_Gamma   # list of 4 HP sets
    N_runs_ddnn <- 4
    ddnn_runs   <- vector("list", N_runs_ddnn)
    for (run in seq_len(N_runs_ddnn)) {
      hp <- hp_ddnn[[run]]
      cat(sprintf("    DDNN Run %d/%d: h1=%d h2=%d lr=%.0e l1=%.0e\n",
                  run, N_runs_ddnn, hp$h1, hp$h2, hp$lr, hp$l1_hidden))
      set.seed(2026 + run * 31)
      ddnn_runs[[run]] <- list(
        net = train_ddnn(Xtr_ens_sc, y_tr, X_va_cv, y_va_cv,
                         h1=hp$h1, h2=hp$h2, lr=hp$lr,
                         l1_hidden=hp$l1_hidden,
                         max_epochs=300, patience=30,
                         batch_size=512, seed=2026+run*31),
        best_hp = hp
      )
    }
    ddnn_qgrid_runs <- lapply(ddnn_runs, function(r) {
      p <- get_gamma_params_ddnn(r$net, Xte_ens_sc)
      gamma_qgrid_fn(p$shape, p$rate, q_grid)
    })
    ddnn_qgrid_mat  <- pmax(Reduce(`+`, ddnn_qgrid_runs) / N_runs_ddnn, 1)
    ddnn_qgrid      <- as.data.frame(ddnn_qgrid_mat)
    colnames(ddnn_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    ddnn_q <- data.frame(q05=ddnn_qgrid_mat[,q05_idx], q50=ddnn_qgrid_mat[,q50_idx], q95=ddnn_qgrid_mat[,q95_idx])
    cat("  DDNN ✓\n")
  }, envir = env)
}

# ---- Recovery: Ensemble DRF re-prediction.
# Source: script_recovery.R:655-710
stage_recovery_edrf <- function(env = .GlobalEnv) {
  local({
    # ---- Ensemble DRF — re-predict from saved base learners (no torch, works fine) ----
    cat("  Ensemble DRF...")
    library(quantreg)
    K_drf       <- length(drf_base_full)
    # Rebuild stacking split (same seed as original)
    set.seed(2026 + 900)
    n_tr     <- nrow(Xtr_tree)
    idx_set1 <- sample(n_tr, floor(n_tr/2), replace=FALSE)
    idx_set2 <- setdiff(seq_len(n_tr), idx_set1)
    X_set2   <- Xtr_tree[idx_set2, , drop=FALSE]; y_set2 <- y_tr[idx_set2]
    # Re-predict Set2 from base learners (retrained on full train in original)
    te_qpreds <- lapply(drf_base_full, function(fit) {
      pr <- predict(fit, Xte_tree)
      wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
    })
    # Rebuild QR combiners from Set2 predictions of base learners trained on Set1
    # (we need drf_base_s1 predictions — use drf_base_full as approximation,
    #  or rerun Set1 base learners with same seeds)
    drf_base_s1 <- lapply(seq_len(K_drf), function(k) {
      mtry_k <- if (k==1) floor(sqrt(ncol(Xtr_tree))) else floor(ncol(Xtr_tree)/3)
      drf(X=Xtr_tree[idx_set1,,drop=FALSE], Y=y_tr[idx_set1],
          num.trees=best_hyperparameters$Ensemble_DRF$num.trees,
          min.node.size=best_hyperparameters$Ensemble_DRF$min.node.size,
          mtry=mtry_k, sample.fraction=0.5, seed=2026+k*7)
    })
    set2_qpreds <- lapply(drf_base_s1, function(fit) {
      pr <- predict(fit, X_set2)
      wquant_matrix(as.numeric(pr$y), pr$weights, q_grid)
    })
    qr_combiners <- vector("list", length(q_grid))
    for (qi in seq_along(q_grid)) {
      tau_i   <- q_grid[qi]
      X_comb  <- do.call(cbind, lapply(set2_qpreds, function(qm) qm[,qi]))
      df_comb <- as.data.frame(X_comb)
      colnames(df_comb) <- paste0("bl", seq_len(K_drf))
      df_comb$y <- y_set2
      fmla <- as.formula(paste("y ~", paste(paste0("bl", seq_len(K_drf)), collapse=" + ")))
      qr_combiners[[qi]] <- tryCatch(rq(fmla, tau=tau_i, data=df_comb, method="br"),
                                     error=function(e) NULL)
    }
    edrf_qgrid_mat <- matrix(NA_real_, nrow(Xte_tree), length(q_grid))
    for (qi in seq_along(q_grid)) {
      X_te_comb  <- do.call(cbind, lapply(te_qpreds, function(qm) qm[,qi]))
      df_te_comb <- as.data.frame(X_te_comb)
      colnames(df_te_comb) <- paste0("bl", seq_len(K_drf))
      if (!is.null(qr_combiners[[qi]])) {
        edrf_qgrid_mat[,qi] <- pmax(predict(qr_combiners[[qi]], newdata=df_te_comb), 1)
      } else {
        edrf_qgrid_mat[,qi] <- rowMeans(X_te_comb)
      }
    }
    for (i in seq_len(nrow(edrf_qgrid_mat))) edrf_qgrid_mat[i,] <- sort(edrf_qgrid_mat[i,])
    edrf_qgrid <- as.data.frame(edrf_qgrid_mat)
    colnames(edrf_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    edrf_q <- data.frame(q05=edrf_qgrid_mat[,q05_idx], q50=edrf_qgrid_mat[,q50_idx], q95=edrf_qgrid_mat[,q95_idx])
    cat(" ✓\n")
  }, envir = env)
}

# ---- Recovery: retrain BNN with saved HPs.
# Source: script_recovery.R:712-728
stage_recovery_bnn <- function(env = .GlobalEnv) {
  local({
    # ---- BNN — retrain with saved HPs ----
    cat("  BNN (retraining with saved HPs)...\n")
    hp_bnn <- best_hyperparameters$BNN_Gamma    # list(h1, h2, lr, pi_mix, neg_log_s1, neg_log_s2)
    set.seed(2026 + 77)
    bnn_net <- train_bnn_gamma(
      Xtr_ens_sc, y_tr,
      h1=hp_bnn$h1, h2=hp_bnn$h2, lr=hp_bnn$lr,
      pi_mix     = hp_bnn$pi_mix,
      log_sigma1 = -hp_bnn$neg_log_s1,
      log_sigma2 = -hp_bnn$neg_log_s2,
      epochs=150, batch_size=128, seed=2026+77
    )
    bnn_qgrid_mat <- pmax(bnn_gamma_quantiles(bnn_net, Xte_ens_sc, M_pass=100), 1)
    bnn_qgrid     <- as.data.frame(bnn_qgrid_mat)
    colnames(bnn_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    bnn_q <- data.frame(q05=bnn_qgrid_mat[,q05_idx], q50=bnn_qgrid_mat[,q50_idx], q95=bnn_qgrid_mat[,q95_idx])
    cat("  BNN ✓\n")
  }, envir = env)
}

# ---- Recovery: retrain hybrid chain.
# Source: script_recovery.R:730-809
stage_recovery_hybrid <- function(env = .GlobalEnv) {
  local({
    # ---- Hybrid NN+DRF — retrain full chain ----
    cat("  Hybrid NN+DRF (retraining full hybrid chain)...\n")
    set.seed(2026)
    y_mean_ens <- mean(y_tr); y_sd_ens <- sd(y_tr)
    y_tr_sc    <- (y_tr - y_mean_ens) / y_sd_ens

    pinball_loss_hyb <- function(pred, y_true, tau) {
      err <- y_true - pred
      torch_mean(torch_where(err >= 0, tau * err, (tau - 1) * err))
    }

    nn_feat_extractor <- nn_sequential(
      nn_linear(n_input, 128), nn_relu(),
      nn_linear(128, 32), nn_relu()
    )
    nn_out_head <- nn_linear(32, n_quantiles)

    nn_feat_extractor$to(device=device)
    nn_out_head$to(device=device)

    opt_hyb <- optim_adam(c(nn_feat_extractor$parameters, nn_out_head$parameters), lr=5e-5)
    dl_hyb  <- make_dl(Xtr_ens_sc, y_tr_sc, batch_size=512)

    nn_feat_extractor$train()
    nn_out_head$train()

    for (epoch in seq_len(150)) {
      coro::loop(for (b in dl_hyb) {
        x_batch <- b[[1]]
        y_batch <- b[[2]]
        opt_hyb$zero_grad()
        loss <- pinball_loss_hyb(nn_out_head(nn_feat_extractor(x_batch)), y_batch, tau_t)
        loss$backward()
        opt_hyb$step()
      })
      if (epoch %% 50 == 0) {
        cat(sprintf("    Hybrid epoch %d, Loss: %.4f\n", epoch, as.numeric(loss$item())))
      }
    }

    nn_feat_extractor$eval()

    X_tr_t <- torch_tensor(as.matrix(Xtr_ens_sc), dtype=torch_float(), device=device)
    X_te_t <- torch_tensor(as.matrix(Xte_ens_sc), dtype=torch_float(), device=device)

    with_no_grad({
      neural_tr <- as.matrix(nn_feat_extractor(X_tr_t)$cpu())
      neural_te <- as.matrix(nn_feat_extractor(X_te_t)$cpu())
    })

    colnames(neural_tr) <- colnames(neural_te) <- paste0("nn_feat_", seq_len(ncol(neural_tr)))
    df_neural_tr <- as.data.frame(neural_tr)
    df_neural_te <- as.data.frame(neural_te)

    hp_drf <- best_hyperparameters$DRF

    hybrid_drf <- drf(
      X               = df_neural_tr,
      Y               = y_tr,
      num.trees       = best_hyperparameters$DRF$num.trees,
      min.node.size   = best_hyperparameters$DRF$min.node.size,
      mtry            = floor(sqrt(ncol(df_neural_tr))),
      sample.fraction = 0.5,
      seed            = 2026
    )
    hybrid_pred <- predict(hybrid_drf, newdata = df_neural_te)
    hybrid_qgrid_mat <- pmax(
      wquant_matrix(as.numeric(hybrid_pred$y), hybrid_pred$weights, q_grid),
      min(y_tr, na.rm = TRUE)
    )

    hybrid_qgrid <- as.data.frame(hybrid_qgrid_mat)
    colnames(hybrid_qgrid) <- paste0("q", sprintf("%02d", round(100*q_grid)))
    hybrid_q <- data.frame(
      q05 = hybrid_qgrid_mat[, q05_idx],
      q50 = hybrid_qgrid_mat[, q50_idx],
      q95 = hybrid_qgrid_mat[, q95_idx]
    )

    cat("  Hybrid NN+DRF ✓\n")
  }, envir = env)
}

# ---- Recovery: assemble all_q / all_qgrid.
# Source: script_recovery.R:812-852
stage_recovery_assemble <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # STEP 7: Assemble all_q and all_qgrid
    # ============================================================
    all_q <- list(
      DRF           = drf_q,
      QRF           = qrf_q,
      RF            = rf_q,
      XGB           = xgb_q,
      ENS_Gamma     = ens_q,
      MCD_Gamma     = mcd_q,
      DDNN_Gamma    = ddnn_q,
      Ensemble_DRF  = edrf_q,
      BNN_Gamma     = bnn_q,
      Hybrid_NN_DRF = hybrid_q
    )

    all_qgrid <- list(
      DRF           = drf_qgrid,
      QRF           = qrf_qgrid,
      RF            = rf_qgrid,
      XGB           = xgb_qgrid,
      ENS_Gamma     = ens_qgrid,
      MCD_Gamma     = mcd_qgrid,
      DDNN_Gamma    = ddnn_qgrid,
      Ensemble_DRF  = edrf_qgrid,
      BNN_Gamma     = bnn_qgrid,
      Hybrid_NN_DRF = hybrid_qgrid
    )

    model_names <- names(all_q)

    MODEL_COLOR <- c(
      Actual        = "gray20",
      DRF           = "#D65C5C", QRF          = "#5A9BD4", RF            = "#66BB66",
      XGB           = "#E8A838", ENS_Gamma    = "#9B59B6", MCD_Gamma     = "#E67E22",
      DDNN_Gamma    = "#1ABC9C", Ensemble_DRF = "#C0392B", BNN_Gamma     = "#2471A3",
      Hybrid_NN_DRF = "#6E2F1A"
    )

    cat("\n✓ all_q and all_qgrid assembled for all 10 models\n")
  }, envir = env)
}

# ---- Recovery: eval_df and CRPS sanity check against the saved values.
# Source: script_recovery.R:854-877
stage_recovery_check <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # STEP 8: Rebuild eval_df (base for all subgroup analyses)
    # ============================================================
    eval_df <- do.call(data.frame, lapply(names(all_q), function(m) {
      mn <- tolower(gsub("[^a-zA-Z0-9]","_", m))
      q  <- all_q[[m]]
      setNames(q, paste0(mn, c("_q05","_q50","_q95")))
    }))
    eval_df$y <- y_te

    # Quick CRPS sanity check vs saved values
    crps_score_full <- function(qmat, y) crps_from_quantiles(y, as.matrix(qmat), q_grid)
    crps_check <- sapply(model_names, function(m) crps_score_full(all_qgrid[[m]], y_te))
    cat("\nCRPS sanity check (should match saved values):\n")
    check_df <- data.frame(
      Model   = model_names,
      Saved   = round(crps_vec[model_names], 4),
      Rebuilt = round(crps_check,            4),
      Match   = abs(crps_vec[model_names] - crps_check) < 1e-3
    )
    print(check_df)

    cat("\n✓ Recovery complete — environment ready for analysis.\n")
    cat("  You can now run the Heterogeneity, Fairness, and Sensitivity sections.\n\n")
  }, envir = env)
}

# ---- Recovery: torch_save() of all neural models into torch_models/.
# Source: script_recovery.R:879-907
stage_recovery_save_torch <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # 保存 5 个 torch 模型（必须用 torch_save，不能用 saveRDS）
    # ============================================================

    torch_dir <- file.path(out_dir, "torch_models")
    if (!dir.exists(torch_dir)) dir.create(torch_dir, recursive = TRUE)

    # ENS Gamma — 5 个 member 分别保存
    for (m in seq_along(ens_members)) {
      torch_save(ens_members[[m]],
                 file.path(torch_dir, sprintf("ens_member_%d.pt", m)))
    }

    # MCD
    torch_save(mc_dropout_net,    file.path(torch_dir, "mcd_net.pt"))

    # DDNN — 4 个 run 分别保存（只保存 net，best_hp 已在 rds 里）
    for (run in seq_along(ddnn_runs)) {
      torch_save(ddnn_runs[[run]]$net,
                 file.path(torch_dir, sprintf("ddnn_run_%d.pt", run)))
    }

    # BNN
    torch_save(bnn_net,           file.path(torch_dir, "bnn_net.pt"))

    # Hybrid — feature extractor
    torch_save(nn_feat_extractor, file.path(torch_dir, "hybrid_feat_extractor.pt"))

    cat("✓ 所有 torch 模型已保存到:", torch_dir, "\n")
  }, envir = env)
}
