# ============================================================
# Pure metric / quantile utilities (no dependence on global state).
# Verbatim from scripts.R (identical arithmetic to the copies in
# 39quantiles.R and script_recovery.R).
# ============================================================

pinball_fn_vec <- function(y, q, tau) {
  u <- y - q
  ifelse(u >= 0, tau * u, (tau - 1) * u)
}

# Compute per-observation CRPS then average (req. 7)
# qmat: [n x K] matrix of predicted quantiles
# taus: vector of K probability levels
# Returns: scalar mean CRPS
crps_from_quantiles <- function(y, qmat, taus) {
  ord   <- order(taus)
  taus  <- taus[ord]
  qmat  <- qmat[, ord, drop = FALSE]
  K     <- length(taus)
  
  # Trapezoidal integration weights
  w <- numeric(K)
  w[1]  <- (taus[2] - taus[1]) / 2
  w[K]  <- (taus[K] - taus[K-1]) / 2
  if (K > 2)
    for (j in 2:(K-1)) w[j] <- (taus[j+1] - taus[j-1]) / 2
  
  # Per-observation CRPS: 2 * sum_k w[k] * L_tau[k](y_i, q_ik)
  # [n x K] pinball matrix, then weighted row sums, then times 2
  n <- length(y)
  crps_per_obs <- numeric(n)
  for (j in seq_len(K)) {
    crps_per_obs <- crps_per_obs + w[j] * pinball_fn_vec(y, qmat[, j], taus[j])
  }
  crps_per_obs <- 2 * crps_per_obs  # per-observation CRPS
  
  mean(crps_per_obs, na.rm = TRUE)  # average CRPS (req. 7)
}

wquant_matrix <- function(y_train, W_sparse, probs) {
  ord <- order(y_train)
  y_sorted <- y_train[ord]
  W_sorted <- W_sparse[, ord, drop = FALSE]
  n_test <- nrow(W_sorted)
  n_probs <- length(probs)
  res <- matrix(NA, nrow = n_test, ncol = n_probs)
  W_r <- as(W_sorted, "RsparseMatrix")
  for (i in 1:n_test) {
    row_start <- W_r@p[i] + 1L
    row_end   <- W_r@p[i + 1L]
    if (row_start > row_end) next
    nz_idx <- W_r@j[row_start:row_end] + 1L
    nz_w   <- W_r@x[row_start:row_end]
    nz_y   <- y_sorted[nz_idx]
    cw <- cumsum(nz_w) / sum(nz_w)
    for (p_idx in 1:n_probs) {
      idx <- which(cw >= probs[p_idx])[1]
      res[i, p_idx] <- nz_y[idx]
    }
  }
  return(res)
}
