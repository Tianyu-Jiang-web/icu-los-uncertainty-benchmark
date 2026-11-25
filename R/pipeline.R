# ============================================================
# Pipeline helpers: module loading, device selection, stage runner
# ============================================================

# Source every module under R/ (stage functions are only *defined* here;
# nothing is executed until run_stages()).
# plot_lancet_style.R is deliberately NOT sourced here: the original scripts
# source it at a specific point (inside a stage), which is preserved.
load_modules <- function(root = cfg$root) {
  files <- c("R/utils_metrics.R", "R/pipeline.R", "R/packages.R", "R/data_prep.R",
             "R/eda_plots.R", "R/table1.R", "R/training_setup.R",
             list.files(file.path(root, "R/models"), full.names = FALSE, pattern = "\\.R$") |>
               (\(x) file.path("R/models", x))(),
             "R/predict_setup.R",
             "R/analysis/evaluation.R", "R/analysis/heterogeneity.R", "R/analysis/fairness.R",
             "R/analysis/calibration.R", "R/analysis/sensitivity.R",
             "R/plots/heterogeneity_plots.R", "R/plots/fairness_plots.R",
             "R/plots/sensitivity_plots.R", "R/plots/crps_guide.R",
             "R/outputs.R", "R/recovery.R")
  for (f in unique(files)) source(file.path(root, f), local = FALSE)
  invisible(files)
}

# Same rule as the original scripts: MPS if available, else CPU.
# DRF_FORCE_CPU=1 forces CPU (useful for bit-comparable tests).
select_torch_device <- function() {
  if (identical(Sys.getenv("DRF_FORCE_CPU"), "1")) return("cpu")
  if (torch::backends_mps_is_available()) "mps" else "cpu"
}

# Run stages in order.  Optional env vars:
#   DRF_ONLY="stage_a,stage_b"  run only these
#   DRF_SKIP="stage_a,stage_b"  skip these
run_stages <- function(stages, env = .GlobalEnv) {
  only <- strsplit(Sys.getenv("DRF_ONLY", ""), ",")[[1]]
  skip <- strsplit(Sys.getenv("DRF_SKIP", ""), ",")[[1]]
  for (s in stages) {
    if (length(only) && !(s %in% only)) next
    if (s %in% skip) next
    fn <- get(s, mode = "function")
    cat(sprintf("\n>>>>>> [%s] %s\n", format(Sys.time(), "%H:%M:%S"), s))
    t0 <- Sys.time()
    fn(env)
    cat(sprintf("<<<<<< %s done in %.1f min\n", s, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  invisible(NULL)
}
