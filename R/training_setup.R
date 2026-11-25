# ============================================================
# Tuning set-up shared by all models (CV folds, 19-level grid).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- 3-fold CV folds for tuning and the 19-level training q_grid.
# Source: scripts.R:533-540
stage_train_setup_cv <- function(env = .GlobalEnv) {
  local({
    library(caret)

    # Create 3-fold CV splits (balance speed vs robustness)
    set.seed(2026)
    cv_folds_hp <- createFolds(y_tr, k = 3, list = TRUE, returnTrain = FALSE)

    # q_grid defined early — needed by CRPS tuning helpers
    q_grid <- seq(0.05, 0.95, by = 0.05)
  }, envir = env)
}

# ---- crps_from_qmat() bound to the current q_grid (pure CRPS/pinball/wquant live in R/utils_metrics.R).
# Source: scripts.R:585-587
stage_train_define_crps_helpers <- function(env = .GlobalEnv) {
  local({
    crps_from_qmat <- function(y, qmat) {
      crps_from_quantiles(y, as.matrix(qmat), q_grid)
    }
  }, envir = env)
}
