# ============================================================
# Package loading, one stage per pipeline (kept in the original order to preserve masking).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Libraries for the training pipeline (scripts.R header).
# Source: scripts.R:1-4
stage_packages_train <- function(env = .GlobalEnv) {
  local({
    library(dplyr)
    library(readr)
    library(lubridate)
    library(tidyr)
  }, envir = env)
}

# ---- Libraries for the torch-recovery pipeline (script_recovery.R STEP 1).
# Source: script_recovery.R:28-47
stage_packages_recovery <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    library(dplyr)
    library(readr)
    library(lubridate)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(scales)
    library(ggridges)
    library(caret)
    library(drf)
    library(quantregForest)
    library(ranger)
    library(xgboost)
    library(torch)
    library(scoringutils)
    library(coro)

    cat("✓ Libraries loaded\n")
  }, envir = env)
}

# ---- Libraries for the 39-quantile predict/visualise pipeline (39quantiles.R STEP 1).
# Source: 39quantiles.R:32-51
stage_packages_predict <- function(env = .GlobalEnv) {
  local({
    library(dplyr)
    library(readr)
    library(lubridate)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    library(scales)
    library(ggridges)
    library(caret)
    library(drf)
    library(quantregForest)
    library(ranger)
    library(xgboost)
    library(torch)
    library(scoringutils)
    library(coro)
    library(quantreg)
    library(ggforce)

    cat("✓ Libraries loaded\n")
  }, envir = env)
}
