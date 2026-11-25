# ============================================================
# CRPS illustration figure.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- CRPS pedagogical illustration (needs q_grid, drf_qmat, y_te from the predict pipeline).
# Source: crps_guidance_lancet.R:6-133
stage_crps_guide <- function(env = .GlobalEnv) {
  local({
    library(ggplot2)
    library(dplyr)
    library(scales)

    # Load shared Lancet helpers (defines lancet_num, theme_lancet, etc.)
    source(file.path(cfg$root, "R", "plot_lancet_style.R"))

    stopifnot(exists("q_grid"), exists("drf_qmat"), exists("y_te"))
    if (!exists("out_dir")) out_dir <- getwd()

    # ------------------------------------------------------------
    # 1. Per-observation CRPS — same trapezoidal discretisation as
    #    crps_from_quantiles() in the main script
    # ------------------------------------------------------------
    pinball_vec <- function(y, q, tau) {
      u <- y - q; ifelse(u >= 0, tau * u, (tau - 1) * u)
    }
    crps_per_obs_fn <- function(qmat, y, taus) {
      ord  <- order(taus); taus <- taus[ord]; qmat <- qmat[, ord, drop = FALSE]
      K <- length(taus); w <- numeric(K)
      w[1] <- (taus[2] - taus[1]) / 2
      w[K] <- (taus[K] - taus[K - 1]) / 2
      if (K > 2) for (j in 2:(K - 1)) w[j] <- (taus[j + 1] - taus[j - 1]) / 2
      out <- numeric(length(y))
      for (j in seq_len(K))
        out <- out + w[j] * pinball_vec(y, qmat[, j], taus[j])
      2 * out
    }
    crps_i <- crps_per_obs_fn(drf_qmat, y_te, q_grid)

    # ------------------------------------------------------------
    # 2. Pick a typical test observation
    # ------------------------------------------------------------
    cand_idx <- which(y_te >= 2 & y_te <= 10)
    med_crps <- median(crps_i[cand_idx], na.rm = TRUE)
    pick_idx <- cand_idx[which.min(abs(crps_i[cand_idx] - med_crps))]

    cat(sprintf("Selected test obs #%d: y = %.2f days, CRPS = %.3f (median CRPS = %.3f)\n",
                pick_idx, y_te[pick_idx], crps_i[pick_idx], med_crps))

    y_obs     <- as.numeric(y_te[pick_idx])
    quantiles <- as.numeric(drf_qmat[pick_idx, ])

    # ------------------------------------------------------------
    # 3. Build predictive CDF by interpolating the 39 quantiles
    # ------------------------------------------------------------
    x_min <- max(min(quantiles) - 0.2 * diff(range(quantiles)), 0)
    x_max <- max(quantiles) + 0.2 * diff(range(quantiles))
    x_grid <- seq(x_min, x_max, length.out = 800)
    F_vals <- approx(
      x    = c(x_min, quantiles, x_max),
      y    = c(0,     q_grid,    1),
      xout = x_grid, rule = 2, ties = "ordered"
    )$y
    F_vals <- pmin(pmax(F_vals, 0), 1)
    H_vals <- as.numeric(x_grid >= y_obs)

    df_curves <- data.frame(x = x_grid, F = F_vals, H = H_vals)
    df_poly   <- data.frame(
      x = c(x_grid, rev(x_grid)),
      y = c(F_vals, rev(H_vals))
    )

    # ------------------------------------------------------------
    # 4. Figure — Lancet style, larger base type
    # ------------------------------------------------------------
    col_F <- "#00468B"   # Lancet deep blue — predictive CDF
    col_H <- "#AD002A"   # Lancet dark red  — Heaviside step

    p_crps_guide <- ggplot() +
      # shaded discrepancy region between F(x) and H(x)
      geom_polygon(
        data = df_poly,
        aes(x = x, y = y),
        fill = col_F, alpha = 0.12
      ) +
      # predictive CDF
      geom_line(
        data = df_curves,
        aes(x = x, y = F, colour = "Predictive CDF  F(x)"),
        linewidth = 1.0
      ) +
      # Heaviside step at y
      geom_line(
        data = df_curves,
        aes(x = x, y = H, colour = "Observed step  1{x \u2265 y}"),
        linewidth = 1.0
      ) +
      scale_colour_manual(
        values = c(
          "Predictive CDF  F(x)"        = col_F,
          "Observed step  1{x \u2265 y}" = col_H
        ),
        breaks = c("Predictive CDF  F(x)", "Observed step  1{x \u2265 y}")
      ) +
      scale_x_continuous(
        name   = "ICU length of stay (days)",
        labels = lancet_num(1),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      scale_y_continuous(
        name   = "Cumulative probability",
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0\u00b725", "0\u00b750", "0\u00b775", "1\u00b700"),
        limits = c(0, 1.02),
        expand = expansion(mult = c(0, 0.01))
      ) +
      theme_lancet(base_size = 16) +          # LARGER type for thesis
      theme(
        legend.key.width = unit(1.4, "cm")
      )

    print(p_crps_guide)

    # ------------------------------------------------------------
    # 5. Save — larger canvas so the larger text has room to breathe
    # ------------------------------------------------------------
    out_png <- file.path(out_dir, "crps_guide_illustration.png")
    out_pdf <- file.path(out_dir, "crps_guide_illustration.pdf")

    ggsave(out_png, plot = p_crps_guide,
           width = 7.5, height = 5.0, dpi = 600, bg = "white")
    ggsave(out_pdf, plot = p_crps_guide,
           width = 7.5, height = 5.0, bg = "white")

    cat(sprintf("\n\u2713 Saved: %s\n", out_png))
    cat(sprintf("\u2713 Saved: %s\n",   out_pdf))
  }, envir = env)
}
