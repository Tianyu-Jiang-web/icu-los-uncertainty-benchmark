# ============================================================
# Sensitivity figures (Lancet style).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Lancet-style sensitivity figures.
# Source: 39quantiles.R:2875-3026
stage_plot39_sensitivity <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # Lancet-style sensitivity plots
    # ============================================================
    # source("lancet_style.R")   # if not already loaded

    # Clean labels (used in both plots)
    model_label_map <- c(
      Hybrid_NN_DRF = "Hybrid",
      Ensemble_DRF  = "Ens DRF",
      ENS_Gamma     = "ENS",
      MCD_Gamma     = "MCD",
      DDNN_Gamma    = "DDNN",
      BNN_Gamma     = "BNN",
      DRF           = "DRF",
      RF            = "RF",
      QRF           = "QRF",
      XGB           = "XGBoost"
    )

    # Lancet 4-group palette (keys = clean labels)
    sens_colours <- c(
      "Hybrid"  = "#00468B",
      "Ens DRF" = "#5DA1DD",
      "ENS"     = "#AD002A",
      "MCD"     = "#ED0000",
      "DDNN"    = "#E18727",
      "BNN"     = "#FDAF91",
      "DRF"     = "#1B1919",
      "RF"      = "#5A5A5A",
      "QRF"     = "#8A8A8A",
      "XGBoost" = "#B8B8B8"
    )

    # ---- Relative degradation plot (Main text) -----------------
    # Order legend by 30% degradation (worst at top, most robust at bottom)
    legend_order <- sens_df %>%
      filter(miss_rate == 0.30) %>%
      arrange(desc(crps_pct_change)) %>%
      pull(model) %>%
      as.character()

    sens_delta_df <- sens_df %>%
      filter(miss_rate > 0) %>%
      mutate(
        model     = factor(model, levels = legend_order),
        model_lab = recode(as.character(model), !!!model_label_map),
        model_lab = factor(model_lab,
                           levels = unname(model_label_map[legend_order]))
      )

    p_sens_delta_pub <- ggplot(
      sens_delta_df,
      aes(x = miss_rate * 100, y = crps_pct_change,
          colour = model_lab, group = model_lab)
    ) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey35", linewidth = 0.5) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2.8) +
      scale_colour_manual(values = sens_colours, name = NULL) +
      scale_x_continuous(
        breaks = c(10, 20, 30),
        labels = c("10", "20", "30"),
        expand = expansion(mult = c(0.04, 0.06))
      ) +
      scale_y_continuous(
        labels = function(x) paste0(gsub("\\.", "\u00b7", sprintf("%g", x)), "%")
      ) +
      labs(
        x = "Additional missingness in test set (%)",
        y = "Change in CRPS from baseline"
      ) +
      theme_lancet(base_size = 14) +
      theme(
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),
        legend.position    = "right",
        legend.text        = element_text(size = 11, colour = "grey15"),
        legend.title       = element_blank(),
        legend.key.height  = unit(1.1, "lines"),
        legend.key.width   = unit(1.0, "cm"),
        legend.spacing.y   = unit(0.15, "cm"),
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text  = element_text(size = 11, colour = "grey20"),
        axis.title = element_text(size = 13, colour = "grey10"),
        axis.line  = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks = element_line(colour = "grey25", linewidth = 0.4),
        plot.margin = ggplot2::margin(t = 10, r = 14, b = 10, l = 12)
      )

    ggsave(file.path(out_dir, "sensitivity_missing_crps_pct_pub.pdf"),
           plot = p_sens_delta_pub, width = 9.0, height = 5.5, device = "pdf")
    ggsave(file.path(out_dir, "sensitivity_missing_crps_pct_pub.png"),
           plot = p_sens_delta_pub, width = 9.0, height = 5.5, dpi = 600, bg = "white")
    cat("\u2713 Saved: sensitivity_missing_crps_pct_pub.pdf / .png\n")

    # ---- Absolute CRPS plot (Supplementary) --------------------
    sens_plot_df <- sens_df %>%
      mutate(
        model     = factor(model, levels = legend_order),
        model_lab = recode(as.character(model), !!!model_label_map),
        model_lab = factor(model_lab,
                           levels = unname(model_label_map[legend_order]))
      )

    p_sens_abs <- ggplot(
      sens_plot_df,
      aes(x = miss_rate * 100, y = crps_mean,
          colour = model_lab, group = model_lab)
    ) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2.8) +
      scale_colour_manual(values = sens_colours, name = NULL) +
      scale_x_continuous(
        breaks = c(0, 10, 20, 30),
        expand = expansion(mult = c(0.02, 0.04))
      ) +
      scale_y_continuous(
        labels = function(x) gsub("\\.", "\u00b7", sprintf("%.2f", x))
      ) +
      labs(
        x = "Additional missingness in test set (%)",
        y = "CRPS"
      ) +
      theme_lancet(base_size = 14) +
      theme(
        plot.title         = element_blank(),
        legend.position    = "right",
        legend.text        = element_text(size = 11, colour = "grey15"),
        legend.title       = element_blank(),
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text  = element_text(size = 11, colour = "grey20"),
        axis.title = element_text(size = 13, colour = "grey10"),
        axis.line  = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks = element_line(colour = "grey25", linewidth = 0.4),
        plot.margin = ggplot2::margin(t = 10, r = 14, b = 10, l = 12)
      )

    ggsave(file.path(out_dir, "sensitivity_missing_crps_pub.pdf"),
           plot = p_sens_abs, width = 9.0, height = 5.5, device = "pdf")
    ggsave(file.path(out_dir, "sensitivity_missing_crps_pub.png"),
           plot = p_sens_abs, width = 9.0, height = 5.5, dpi = 600, bg = "white")
    cat("\u2713 Saved: sensitivity_missing_crps_pub.pdf / .png\n")

    write.csv(sens_df, file.path(out_dir, "sensitivity_missing_results.csv"), row.names = FALSE)
    cat("✓ Saved: sensitivity_missing_results.csv\n")
  }, envir = env)
}
