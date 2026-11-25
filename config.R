# ============================================================
# Central configuration (paths, grids, switches)
# All settings can be overridden with environment variables so the
# project can be relocated / tested without editing code:
#   DRF_DATA_PATH     master CSV
#   DRF_OUTPUT_DIR    where this project writes results   (default: <root>/outputs)
#   DRF_MODELS_DIR    folder holding all_model_results.rds (+ torch_models/)
#   DRF_LEGACY_DIR    original (non-modular) project, used as fallback source
#   DRF_FORCE_CPU=1   never use the Apple "mps" device
# ============================================================

# ---- project root: directory of this file ----
.find_root <- function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f  <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f) > 0) return(dirname(normalizePath(f[1])))
  if (file.exists("config.R")) return(normalizePath("."))
  stop("Run from the project root (the folder containing config.R) or via Rscript.")
}
.root <- Sys.getenv("DRF_MODULAR_ROOT", unset = .find_root())

.legacy <- Sys.getenv("DRF_LEGACY_DIR",
                      unset = normalizePath(file.path(.root, "..", "drf-replication-ehr-v3"), mustWork = FALSE))
.output <- Sys.getenv("DRF_OUTPUT_DIR", unset = file.path(.root, "outputs"))

# models_dir: where all_model_results.rds lives.  Prefer this project's own
# training output, otherwise fall back to the original project's saved models.
.default_models <- if (file.exists(file.path(.output, "final", "all_model_results.rds"))) {
  file.path(.output, "final")
} else {
  file.path(.legacy, "figures", "final")
}
.models <- Sys.getenv("DRF_MODELS_DIR", unset = .default_models)

cfg <- list(
  root        = .root,
  data_path   = Sys.getenv("DRF_DATA_PATH",
                           unset = file.path(.legacy, "data", "FINAL_master_dataset_9Feb2026(in).csv")),
  output_dir  = .output,
  final_dir   = file.path(.output, "final"),                 # training / recovery outputs
  pred_dir    = file.path(.output, "final_39quantiles"),     # predict + visualise outputs
  models_dir  = .models,                                     # rds source
  seed        = 2026,
  q_grid_train   = seq(0.05, 0.95, by = 0.05),               # 19 levels (training / recovery)
  q_grid_predict = seq(0.025, 0.975, by = 0.025),            # 39 levels (final evaluation)
  run_crps_guide = TRUE                                      # draw the CRPS illustration at the end
)
# torch models: this project's own (after run_recovery.R) or next to the rds
cfg$torch_dir <- if (dir.exists(file.path(cfg$final_dir, "torch_models"))) {
  file.path(cfg$final_dir, "torch_models")
} else {
  file.path(cfg$models_dir, "torch_models")
}
rm(.root, .legacy, .output, .default_models, .models)
