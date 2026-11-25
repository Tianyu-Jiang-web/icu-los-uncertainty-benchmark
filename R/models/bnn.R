# ============================================================
# Bayesian neural network (Bayes by Backprop, Gamma).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- BNN (Bayes-by-Backprop, Gamma): tuning, final net, Thompson sampling.
# Source: scripts.R:1602-1789
stage_train_bnn <- function(env = .GlobalEnv) {
  local({
    cat("========================================\n")
    cat("MODEL 10: BNN (Bayes by Backprop, Gamma parametric, req. 6)\n")
    cat("========================================\n")

    # ---- Scale mixture prior ----
    torch_log_scale_mixture_prior <- function(w, pi_mix, log_sigma1, log_sigma2) {
      half_log2pi <- 0.5 * log(2 * pi)
      log_p1 <- log(pi_mix)       - half_log2pi - log_sigma1 -
        0.5 * (w / exp(log_sigma1))^2
      log_p2 <- log(1 - pi_mix)   - half_log2pi - log_sigma2 -
        0.5 * (w / exp(log_sigma2))^2
      m      <- torch_maximum(log_p1, log_p2)
      torch_sum(m + torch_log(torch_exp(log_p1 - m) + torch_exp(log_p2 - m)))
    }

    # ---- Bayesian linear layer ----
    bnn_linear <- nn_module(
      classname  = "bnn_linear_bbb",
      initialize = function(in_f, out_f, pi_mix, log_sigma1, log_sigma2) {
        self$in_f      <- in_f
        self$out_f     <- out_f
        self$pi_mix    <- pi_mix
        self$log_sig1  <- log_sigma1
        self$log_sig2  <- log_sigma2
        self$w_mu   <- nn_parameter(torch_randn(out_f, in_f)  * 0.01)
        self$w_rho  <- nn_parameter(torch_full(c(out_f, in_f), -3.0))
        self$b_mu   <- nn_parameter(torch_zeros(out_f))
        self$b_rho  <- nn_parameter(torch_full(c(out_f), -3.0))
      },
      forward = function(x) {
        w_sigma        <- torch_log1p(torch_exp(self$w_rho))
        b_sigma        <- torch_log1p(torch_exp(self$b_rho))
        self$w_sample  <- self$w_mu + w_sigma * torch_randn_like(w_sigma)
        self$b_sample  <- self$b_mu + b_sigma * torch_randn_like(b_sigma)
        nnf_linear(x, self$w_sample, self$b_sample)
      },
      kl_term = function() {
        w_sigma <- torch_log1p(torch_exp(self$w_rho))
        b_sigma <- torch_log1p(torch_exp(self$b_rho))
        log_qw <- torch_sum(
          -0.5 * log(2 * pi) - torch_log(w_sigma) -
            0.5 * ((self$w_sample - self$w_mu) / w_sigma)^2
        )
        log_qb <- torch_sum(
          -0.5 * log(2 * pi) - torch_log(b_sigma) -
            0.5 * ((self$b_sample - self$b_mu) / b_sigma)^2
        )
        log_q  <- log_qw + log_qb
        log_pw <- torch_log_scale_mixture_prior(
          self$w_sample$detach(), self$pi_mix, self$log_sig1, self$log_sig2)
        log_pb <- torch_log_scale_mixture_prior(
          self$b_sample$detach(), self$pi_mix, self$log_sig1, self$log_sig2)
        log_q - (log_pw + log_pb)
      }
    )

    # ---- BNN net: outputs 2 Gamma params (log_shape, log_rate) ----
    bnn_gamma_net <- nn_module(
      classname  = "bnn_gamma_bbb",
      initialize = function(n_in, h1, h2, pi_mix, log_sigma1, log_sigma2) {
        self$l1 <- bnn_linear(n_in, h1, pi_mix, log_sigma1, log_sigma2)
        self$l2 <- bnn_linear(h1,   h2, pi_mix, log_sigma1, log_sigma2)
        self$l3 <- bnn_linear(h2,   2,  pi_mix, log_sigma1, log_sigma2)
      },
      forward = function(x) {
        x <- nnf_relu(self$l1(x))
        x <- nnf_relu(self$l2(x))
        self$l3(x)
      },
      kl = function() { self$l1$kl_term() + self$l2$kl_term() + self$l3$kl_term() }
    )

    # ---- Training with KL re-weighting (Eq. 9) + Gamma NLL ----
    train_bnn_gamma <- function(X_sc, y_obs, h1, h2, lr,
                                pi_mix, log_sigma1, log_sigma2,
                                epochs = 100, batch_size = 128, seed = 42) {
      set.seed(seed); torch_manual_seed(seed)
      N   <- nrow(X_sc)
      net <- bnn_gamma_net(n_input, h1, h2, pi_mix, log_sigma1, log_sigma2)$to(device = device)
      opt <- optim_adam(net$parameters, lr = lr)
      dl  <- make_dl(X_sc, y_obs, batch_size)

      n_batches  <- ceiling(N / batch_size)
      kl_weights <- sapply(seq_len(n_batches), function(i) {
        (2^(n_batches - i)) / (2^n_batches - 1)
      })

      net$train()
      for (ep in seq_len(epochs)) {
        batch_cnt <- 0L
        coro::loop(for (b in dl) {
          batch_cnt <- batch_cnt + 1L
          pi_i      <- kl_weights[min(batch_cnt, n_batches)]
          x_batch   <- b[[1]]
          y_batch   <- b[[2]]
          opt$zero_grad()
          params     <- net(x_batch)
          nll_loss   <- gamma_nll_loss(params, y_batch)
          kl_contrib <- net$kl()
          loss       <- pi_i * kl_contrib + nll_loss
          loss$backward()
          opt$step()
        })
        if (ep %% 25 == 0)
          cat(sprintf("    Epoch %d/%d\n", ep, epochs))
      }
      net$eval()
      net
    }

    # Thompson sampling: average Gamma params over M stochastic forward passes
    bnn_gamma_quantiles <- function(net, X_sc, M_pass = 100) {
      X_sc <- as.matrix(X_sc)
      net$train()  # keep stochasticity ON
      X_t   <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      acc_p <- matrix(0, nrow(X_sc), 2)
      for (p_i in seq_len(M_pass)) {
        with_no_grad({ pm <- as.matrix(net(X_t)$cpu()) })
        acc_p <- acc_p + pm
      }
      avg_p <- acc_p / M_pass
      shape <- exp(avg_p[, 1])
      rate  <- exp(avg_p[, 2])
      gamma_qgrid_fn(shape, rate, q_grid)
    }

    # ---- HP TUNING: BNN — CRPS-based (req. 2) ----
    cat("Tuning BNN hyperparameters (12 configs; CRPS-based, req. 2)...\n")
    bnn_hp_grid <- expand.grid(
      h1          = c(64, 128),
      h2          = c(32, 64),
      lr          = c(5e-4),
      pi_mix      = c(0.75),
      neg_log_s1  = c(0),
      neg_log_s2  = c(4, 6, 7),
      stringsAsFactors = FALSE
    )

    bnn_cv_crps <- sapply(seq_len(nrow(bnn_hp_grid)), function(i) {
      hp <- bnn_hp_grid[i, ]
      cat(sprintf("  [%d/%d] h1=%d h2=%d lr=%.4g pi=%.2f -log(s2)=%g ...",
                  i, nrow(bnn_hp_grid),
                  hp$h1, hp$h2, hp$lr, hp$pi_mix, hp$neg_log_s2))
      tryCatch({
        net <- train_bnn_gamma(
          X_tr_cv, y_tr_cv,
          h1 = hp$h1, h2 = hp$h2, lr = hp$lr,
          pi_mix      = hp$pi_mix,
          log_sigma1  = -hp$neg_log_s1,
          log_sigma2  = -hp$neg_log_s2,
          epochs = 50, batch_size = 128, seed = 2026 + i
        )
        qm  <- bnn_gamma_quantiles(net, X_va_cv, M_pass = 20)
        cr  <- crps_from_qmat(y_va_cv, qm)
        cat(sprintf(" CRPS=%.4f\n", cr)); cr
      }, error = function(e) { cat(sprintf(" ERROR: %s\n", conditionMessage(e))); Inf })
    })

    best_bnn_i  <- which.min(bnn_cv_crps)
    best_bnn_hp <- as.list(bnn_hp_grid[best_bnn_i, ])
    cat(sprintf("\n✓ Best BNN: h1=%d h2=%d lr=%.4g pi=%.2f -log(s2)=%g (CRPS=%.4f)\n",
                best_bnn_hp$h1, best_bnn_hp$h2, best_bnn_hp$lr,
                best_bnn_hp$pi_mix, best_bnn_hp$neg_log_s2,
                bnn_cv_crps[best_bnn_i]))

    # ---- Train final BNN ----
    set.seed(2026 + 77)
    cat("Training final BNN (150 epochs, full train set)...\n")
    bnn_net <- train_bnn_gamma(
      Xtr_ens_sc, y_tr,
      h1          = best_bnn_hp$h1,
      h2          = best_bnn_hp$h2,
      lr          = best_bnn_hp$lr,
      pi_mix      = best_bnn_hp$pi_mix,
      log_sigma1  = -best_bnn_hp$neg_log_s1,
      log_sigma2  = -best_bnn_hp$neg_log_s2,
      epochs = 150, batch_size = 128, seed = 2026 + 77
    )

    cat("BNN inference (M=100 stochastic forward passes)...\n")
    bnn_qgrid_mat  <- bnn_gamma_quantiles(bnn_net, Xte_ens_sc, M_pass = 100)
    bnn_qgrid_mat  <- pmax(bnn_qgrid_mat, 1)
    bnn_q          <- data.frame(q05 = bnn_qgrid_mat[, q05_idx],
                                 q50 = bnn_qgrid_mat[, q50_idx],
                                 q95 = bnn_qgrid_mat[, q95_idx])
    bnn_qgrid      <- as.data.frame(bnn_qgrid_mat)
    colnames(bnn_qgrid) <- paste0("q", sprintf("%02d", round(100 * q_grid)))
    cat("BNN (Gamma) complete.\n\n")
  }, envir = env)
}

# ---- BNN: 100 stochastic passes on the 39-level grid.
# Source: 39quantiles.R:724-731
stage_predict39_bnn <- function(env = .GlobalEnv) {
  local({
    cat("  BNN Gamma (100 passes)...")
    bnn_net$train()
    acc_bnn <- matrix(0, nrow(Xte_ens_sc), 2)
    for (i in seq_len(100))
      with_no_grad({ acc_bnn <- acc_bnn + as.matrix(bnn_net(X_te_t)$cpu()) })
    avg_bnn  <- acc_bnn / 100
    bnn_qmat <- pmax(gamma_qgrid_fn(exp(avg_bnn[, 1]), exp(avg_bnn[, 2]), q_grid), 1)
    cat(" ✓\n")
  }, envir = env)
}
