# ============================================================
# Shared torch infrastructure for the neural models.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Shared neural set-up: one-hot + standardised features, device, CV fold, pinball loss, dataloader builder.
# Source: scripts.R:967-1047
stage_train_torch_common <- function(env = .GlobalEnv) {
  local({
    library(torch)
    library(scoringutils)

    set.seed(2026)
    qs          <- c(0.05, 0.5, 0.95)
    q05_idx     <- which.min(abs(q_grid - 0.05))
    q50_idx     <- which.min(abs(q_grid - 0.50))
    q95_idx     <- which.min(abs(q_grid - 0.95))

    # One-hot encode + standardise X  (train-led, applied to test)
    Xtr_ens_mm  <- model.matrix(~ . - 1, data = Xtr_all)
    Xte_ens_mm  <- model.matrix(~ . - 1, data = Xte_all)
    # align columns
    common_nn   <- union(colnames(Xtr_ens_mm), colnames(Xte_ens_mm))
    Xtr_ens_mm  <- Xtr_ens_mm[, common_nn, drop = FALSE]; Xtr_ens_mm[is.na(Xtr_ens_mm)] <- 0
    Xte_ens_mm  <- Xte_ens_mm[, common_nn, drop = FALSE]; Xte_ens_mm[is.na(Xte_ens_mm)] <- 0

    x_mean      <- colMeans(Xtr_ens_mm)
    x_sd        <- apply(Xtr_ens_mm, 2, sd); x_sd[x_sd == 0] <- 1
    Xtr_ens_sc  <- scale(Xtr_ens_mm, center = x_mean, scale = x_sd)
    Xte_ens_sc  <- scale(Xte_ens_mm, center = x_mean, scale = x_sd)

    Xtr_ens_sc[Xtr_ens_sc > 10] <- 10
    Xtr_ens_sc[Xtr_ens_sc < -10] <- -10
    Xte_ens_sc[Xte_ens_sc > 10] <- 10
    Xte_ens_sc[Xte_ens_sc < -10] <- -10

    n_input     <- ncol(Xtr_ens_sc)
    n_quantiles <- length(q_grid)

    # Device
    device <- torch_device(select_torch_device())
    cat("Neural models using device:", device$type, "\n")

    # Shared CV fold
    val_idx_nn   <- cv_folds_hp[[1]]
    train_idx_nn <- setdiff(seq_len(nrow(Xtr_ens_sc)), val_idx_nn)
    X_tr_cv      <- Xtr_ens_sc[train_idx_nn, , drop = FALSE]
    y_tr_cv      <- y_tr[train_idx_nn]
    X_va_cv      <- Xtr_ens_sc[val_idx_nn,   , drop = FALSE]
    y_va_cv      <- y_tr[val_idx_nn]

    # ============================================================
    # SHARED: Pinball (quantile) loss — NON-PARAMETRIC
    # ============================================================
    pinball_loss_grid <- function(pred, y_true, tau_t) {
      err  <- y_true - pred
      loss <- torch_where(err >= 0, tau_t * err, (tau_t - 1) * err)
      torch_mean(loss)
    }

    # ============================================================
    # SHARED: Dataset / DataLoader builder
    # ============================================================
    make_dl <- function(X_sc, y_obs, batch_size = 512) {
      X_t <- torch_tensor(as.matrix(X_sc), dtype = torch_float(), device = device)
      y_t <- torch_tensor(matrix(as.numeric(y_obs), ncol = 1), dtype = torch_float(), device = device)
      ds  <- tensor_dataset(X_t, y_t)
      dataloader(ds, batch_size = batch_size, shuffle = TRUE)
    }

    # ============================================================
    # SHARED: Quantile head — output n_quantiles directly
    # ============================================================
    tau_t <- torch_tensor(matrix(q_grid, nrow = 1), dtype = torch_float(), device = device)

    # CRPS helper for non-parametric net (primary tuning metric, req. 2)
    qnet_crps <- function(net, X_sc, y_obs) {
      pm <- get_qgrid(net, X_sc)
      crps_from_qmat(y_obs, pm)
    }

    # Extract sorted quantile grid from non-parametric net
    get_qgrid <- function(net, X_sc) {
      net$eval()
      Xv <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      with_no_grad({ pm <- as.matrix(net(Xv)$cpu()) })
      net$train()
      t(apply(pm, 1, sort))
    }
  }, envir = env)
}

# ---- Gamma NLL loss, parameter extraction and Gamma quantile conversion shared by ENS/MCD/DDNN/BNN.
# Source: scripts.R:1059-1083
stage_train_gamma_helpers <- function(env = .GlobalEnv) {
  local({
    # ---- Shared Gamma helpers ----
    gamma_nll_loss <- function(params, y_true) {
      log_shape <- params[, 1, drop = FALSE]
      log_rate  <- params[, 2, drop = FALSE]
      shape <- torch_exp(log_shape)
      rate  <- torch_exp(log_rate)
      nll <- -(shape - 1) * torch_log(y_true) +
        y_true * rate -
        shape * torch_log(rate) +
        torch_lgamma(shape)
      torch_mean(nll)
    }

    get_gamma_params <- function(net, X_sc) {
      X_sc <- as.matrix(X_sc)
      net$eval()
      Xv <- torch_tensor(X_sc, dtype = torch_float(), device = device)
      with_no_grad({ p <- as.matrix(net(Xv)$cpu()) })
      net$train()
      list(shape = exp(p[, 1]), rate = exp(p[, 2]))
    }

    gamma_qgrid_fn <- function(shape, rate, qs) {
      t(mapply(function(sh, ra) qgamma(qs, shape = sh, rate = ra), shape, rate))
    }
  }, envir = env)
}
