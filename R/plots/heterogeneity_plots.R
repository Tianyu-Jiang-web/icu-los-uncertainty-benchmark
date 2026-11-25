# ============================================================
# Heterogeneity figures (Lancet style).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Lancet-style heterogeneity heat-map (main figure).
# Source: 39quantiles.R:930-1229
stage_plot39_heterogeneity_heatmap <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # HETEROGENEITY HEATMAP — LANCET-STYLE MAIN FIGURE
    # Layout:
    #   LEFT  = brackets + category names (e.g. "Age group")
    #   RIGHT = subgroup level names (e.g. "<65", "CCU", "Sepsis")
    #   FILL  = within-subgroup CRPS z-score
    #   CELLS = raw CRPS (Lancet middle-dot decimals)
    # ============================================================

    # Recommended packages: ggplot2, dplyr, tidyr, scales, forcats
    # source("lancet_style.R")

    label_het <- function(df, category) {
      if (is.null(df)) return(NULL)
      df %>%
        rename(subgroup_level = Group) %>%
        mutate(category = category) %>%
        pivot_longer(
          cols      = -c(subgroup_level, category),
          names_to  = "model",
          values_to = "crps"
        ) %>%
        mutate(crps = as.numeric(crps))
    }

    het_long <- bind_rows(
      label_het(out_age_crps,         "Age group"),
      label_het(out_unit_crps,        "ICU care unit"),
      label_het(out_stroke_crps,      "Stroke"),
      label_het(out_sepsis_crps,      "Sepsis"),
      label_het(out_vent_crps,        "Ventilation"),
      label_het(out_vaso_crps,        "Vasopressors"),
      label_het(out_vasopressin_crps, "Vasopressin")
    )

    if (nrow(het_long) == 0) {
      cat("\u26a0 No heterogeneity data found \u2014 skipping heatmap.\n")
    } else {

      # ----------------------------------------------------------
      # 1. Clean labels
      # ----------------------------------------------------------
      het_long <- het_long %>%
        mutate(
          subgroup_level = dplyr::recode(
            subgroup_level,
            "Medical Intensive Care Unit (MICU)"               = "MICU",
            "Surgical Intensive Care Unit (SICU)"              = "SICU",
            "Cardiac Vascular Intensive Care Unit (CVICU)"     = "CVICU",
            "Coronary Care Unit (CCU)"                         = "CCU",
            "Trauma SICU (TSICU)"                              = "TSICU",
            "Neuro Surgical Intensive Care Unit (Neuro SICU)"  = "Neuro SICU",
            "Neuro Stepdown"                                   = "NSDU",
            "Medical/Surgical Intensive Care Unit (MICU/SICU)" = "MICU/SICU"
          )
        )

      # ----------------------------------------------------------
      # 2. Within-subgroup z-score
      # ----------------------------------------------------------
      het_long <- het_long %>%
        group_by(subgroup_level) %>%
        mutate(
          crps_z = (crps - mean(crps, na.rm = TRUE)) / sd(crps, na.rm = TRUE),
          crps_z = ifelse(is.nan(crps_z) | is.infinite(crps_z), 0, crps_z)
        ) %>%
        ungroup()

      # ----------------------------------------------------------
      # 3. Model order
      # ----------------------------------------------------------
      model_order <- c(
        "Hybrid_NN_DRF", "MCD_Gamma", "QRF", "XGB", "Ensemble_DRF",
        "ENS_Gamma", "DDNN_Gamma", "BNN_Gamma", "RF", "DRF"
      )
      model_order <- intersect(model_order, unique(het_long$model))

      model_labels <- c(
        "Hybrid_NN_DRF" = "Hybrid",
        "MCD_Gamma"     = "MCD",
        "QRF"           = "QRF",
        "XGB"           = "XGB",
        "Ensemble_DRF"  = "Ensemble DRF",
        "ENS_Gamma"     = "ENS",
        "DDNN_Gamma"    = "DDNN",
        "BNN_Gamma"     = "BNN",
        "RF"            = "RF",
        "DRF"           = "DRF"
      )

      # ----------------------------------------------------------
      # 4. Subgroup order
      # ----------------------------------------------------------
      unit_levels <- sort(unique(het_long$subgroup_level[het_long$category == "ICU care unit"]))

      subgroup_order <- c(
        "<65", "65-80", ">80",
        unit_levels,
        "No Stroke/Cerebro", "Stroke/Cerebro",
        "No Sepsis", "Sepsis",
        "No Ventilation", "Ventilated",
        "No Vasopressors", "Vasopressors",
        "No Vasopressin", "Vasopressin"
      )
      subgroup_order <- intersect(subgroup_order, unique(het_long$subgroup_level))

      category_order <- c(
        "Age group", "ICU care unit", "Stroke",
        "Sepsis", "Ventilation", "Vasopressors", "Vasopressin"
      )

      # ----------------------------------------------------------
      # 5. Final plotting data
      # ----------------------------------------------------------
      het_plot <- het_long %>%
        mutate(
          model          = factor(model, levels = model_order),
          model_label    = dplyr::recode(model, !!!model_labels),
          model_label    = factor(model_label, levels = unname(model_labels[model_order])),
          subgroup_level = factor(subgroup_level, levels = subgroup_order),
          category       = factor(category, levels = category_order),
          crps_label     = gsub("\\.", "\u00b7", sprintf("%.2f", crps))
        ) %>%
        filter(!is.na(model), !is.na(subgroup_level), !is.na(category))

      # ----------------------------------------------------------
      # 6. Lancet-style muted diverging palette
      # ----------------------------------------------------------
      heat_low  <- "#3B6EA8"
      heat_mid  <- "#F7F7F7"
      heat_high <- "#B85C4D"

      z_lim <- max(abs(het_plot$crps_z), na.rm = TRUE)
      z_lim <- ceiling(z_lim * 10) / 10

      het_plot <- het_plot %>%
        mutate(
          label_colour = ifelse(abs(crps_z) > 0.70 * z_lim, "white", "grey15")
        )

      # ----------------------------------------------------------
      # 7. Bracket geometry — drawn LEFT of the panel
      # ----------------------------------------------------------
      het_plot <- het_plot %>%
        mutate(y_num = as.integer(forcats::fct_rev(subgroup_level)))

      bracket_df <- het_plot %>%
        group_by(category) %>%
        summarise(
          y_lo  = min(y_num) - 0.5,
          y_hi  = max(y_num) + 0.5,
          y_mid = (y_lo + y_hi) / 2,
          .groups = "drop"
        )

      # X positions (data coords). All NEGATIVE so brackets sit LEFT of x=1.
      # coord_cartesian(clip = "off") lets us draw outside the panel.
      x_tip   <- 0.40    # bracket tips (closest to heatmap)
      x_stem  <- 0.05    # vertical stem of the bracket
      x_label <- -0.10   # category label (right-aligned to its x)

      bracket_lines <- bracket_df %>%
        rowwise() %>%
        do(data.frame(
          cat = .$category,
          x    = c(x_tip,   x_tip,   x_stem),
          xend = c(x_stem,  x_stem,  x_stem),
          y    = c(.$y_hi,  .$y_lo,  .$y_hi),
          yend = c(.$y_hi,  .$y_lo,  .$y_lo)
        )) %>%
        ungroup()

      bracket_stubs <- bracket_df %>%
        transmute(
          x    = x_stem,
          xend = x_label + 0.04,
          y    = y_mid,
          yend = y_mid
        )

      # ----------------------------------------------------------
      # 8. Figure — y-axis labels MOVED to RIGHT
      # ----------------------------------------------------------
      p_het_lancet <- ggplot(
        het_plot,
        aes(
          x    = model_label,
          y    = forcats::fct_rev(subgroup_level),
          fill = crps_z
        )
      ) +
        geom_tile(colour = "white", linewidth = 0.55) +
        geom_text(
          aes(label = crps_label, colour = label_colour),
          size = 5, show.legend = FALSE
        ) +
        scale_colour_identity() +

        # ---- Brackets + category names (LEFT side) ----
        geom_segment(
          data = bracket_lines,
          aes(x = x, xend = xend, y = y, yend = yend),
          inherit.aes = FALSE,
          colour = "grey25", linewidth = 0.45, lineend = "round"
        ) +
        geom_segment(
          data = bracket_stubs,
          aes(x = x, xend = xend, y = y, yend = yend),
          inherit.aes = FALSE,
          colour = "grey25", linewidth = 0.45, lineend = "round"
        ) +
        geom_text(
          data = bracket_df,
          aes(x = x_label, y = y_mid, label = category),
          inherit.aes = FALSE,
          hjust = 1, size = 5.5, colour = "grey10"
        ) +

        scale_fill_gradient2(
          low = heat_low, mid = heat_mid, high = heat_high,
          midpoint = 0,
          limits = c(-z_lim, z_lim), oob = scales::squish,
          name  = "CRPS z-score",
          guide = guide_colourbar(
            title.position = "top", title.hjust = 0.5,
            barwidth = unit(0.55, "cm"), barheight = unit(4.2, "cm"),
            ticks.colour = "grey35", frame.colour = "grey35"
          )
        ) +
        scale_x_discrete(position = "top", expand = expansion(add = 0)) +
        # KEY CHANGE: y-axis labels on the RIGHT
        scale_y_discrete(position = "right", expand = expansion(add = 0.25)) +

        coord_cartesian(clip = "off") +

        labs(x = NULL, y = NULL) +

        theme_lancet(base_size = 20) +
        theme(
          plot.title    = element_blank(),
          plot.subtitle = element_blank(),

          axis.text.x.top = element_text(
            size = 16, colour = "grey10",
            angle = 35, hjust = 0, vjust = 0, lineheight = 1.0
          ),
          # Subgroup labels now on RIGHT, left-aligned reading outward
          axis.text.y.right = element_text(
            size = 16, colour = "grey10",
            hjust = 0, margin = ggplot2::margin(l = 6)
          ),
          axis.line  = element_blank(),
          axis.ticks = element_blank(),

          strip.background = element_blank(),
          strip.text       = element_blank(),

          panel.border = element_blank(),
          panel.grid   = element_blank(),

          legend.position = "right",
          legend.title    = element_text(size = 12, colour = "grey10", hjust = 0.5),
          legend.text     = element_text(size = 11, colour = "grey20"),
          legend.margin   = ggplot2::margin(0, 0, 0, 14),

          # Wide LEFT margin for brackets + category names
          plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 150)
        )

      print(p_het_lancet)

      # ----------------------------------------------------------
      # 9. Save figure
      # ----------------------------------------------------------
      if (exists("out_dir") && exists("fig_w") && exists("fig_h")) {
        ggsave(file.path(out_dir, "heterogeneity_heatmap.pdf"),
               plot = p_het_lancet,
               width = fig_w + 1.0, height = fig_h, device = "pdf")
        ggsave(file.path(out_dir, "heterogeneity_heatmap.png"),
               plot = p_het_lancet,
               width = fig_w + 1.0, height = fig_h, dpi = 300)
        cat(sprintf("\u2713 Saved heatmap: heterogeneity_heatmap.pdf / .png\n"))
      }

      # ----------------------------------------------------------
      # 10. Save the merged long table for the supplement
      # ----------------------------------------------------------
      if (exists("out_dir")) {
        write.csv(
          het_plot %>%
            select(category, subgroup_level, model, crps, crps_z) %>%
            arrange(category, subgroup_level, model),
          file.path(out_dir, "heterogeneity_all_subgroups_combined.csv"),
          row.names = FALSE
        )
        cat("\u2713 Saved: heterogeneity_all_subgroups_combined.csv\n")
      }

    }
  }, envir = env)
}

# ---- Subgroup stability (SD of subgroup CRPS) bar chart.
# Source: 39quantiles.R:1230-1367
stage_plot39_heterogeneity_stability <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # SUBGROUP STABILITY BAR CHART — LANCET-STYLE
    # Metric = SD of subgroup-specific CRPS
    # Lower = more stable across subgroups
    # ============================================================

    source(file.path(cfg$root, "R", "plot_lancet_style.R"))  # if not already loaded

    stability_df <- het_long %>%
      group_by(model) %>%
      summarise(
        n_subgroups = n(),
        mean_crps   = mean(crps, na.rm = TRUE),
        sd_crps     = sd(crps, na.rm = TRUE),
        iqr_crps    = IQR(crps, na.rm = TRUE),
        min_crps    = min(crps, na.rm = TRUE),
        max_crps    = max(crps, na.rm = TRUE),
        range_crps  = max_crps - min_crps,
        .groups     = "drop"
      ) %>%
      mutate(
        model_label = dplyr::recode(
          model,
          DRF           = "DRF",
          Ensemble_DRF  = "Ensemble DRF",
          Hybrid_NN_DRF = "Hybrid",
          ENS_Gamma     = "ENS",
          MCD_Gamma     = "MCD",
          DDNN_Gamma    = "DDNN",
          BNN_Gamma     = "BNN",
          QRF           = "QRF",
          RF            = "RF",
          XGB           = "XGB"
        )
      ) %>%
      arrange(sd_crps) %>%
    mutate(
      model_label = factor(model_label, levels = rev(model_label)),
        is_best = model_label == "Hybrid",
        sd_label = gsub("\\.", "\u00b7", sprintf("%.3f", sd_crps))
      )

    # Colours: only highlight the best / most stable model
    highlight_col <- "#3B6EA8"
    other_col     <- "#D9D9D9"
    edge_col      <- "#8A8A8A"

    p_stability_sd_lancet <- ggplot(
      stability_df,
      aes(x = model_label, y = sd_crps)
    ) +
      geom_col(
        aes(fill = is_best),
        width = 0.68,
        colour = edge_col,
        linewidth = 0.25
      ) +
      geom_text(
        aes(label = sd_label),
        hjust = -0.18,
        size = 4.0,
        colour = "grey15"
      ) +
      scale_fill_manual(
        values = c(
          "TRUE" = highlight_col,
          "FALSE" = other_col
        ),
        guide = "none"
      ) +
      scale_y_continuous(
        limits = c(0, 0.44),
        breaks = seq(0, 0.4, by = 0.1),
        labels = function(x) gsub("\\.", "\u00b7", sprintf("%.1f", x)),
        expand = expansion(mult = c(0, 0))
      ) +
      labs(
        x = NULL,
        y = "SD of CRPS across subgroups"
      ) +
      coord_flip(clip = "off") +
      theme_lancet(base_size = 13) +
      theme(
        plot.title = element_blank(),
        plot.subtitle = element_blank(),

        axis.text.y = element_text(
          size = 11.5,
          colour = "grey10"
        ),
        axis.text.x = element_text(
          size = 10.5,
          colour = "grey20"
        ),
        axis.title.x = element_text(
          size = 12,
          colour = "grey10",
          margin = ggplot2::margin(t = 8)
        ),
        axis.title.y = element_blank(),

        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),

        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        plot.margin = ggplot2::margin(t = 8, r = 36, b = 8, l = 8)
      )

    print(p_stability_sd_lancet)

    ggsave(
      file.path(out_dir, "heterogeneity_stability_sd_lancet.pdf"),
      plot = p_stability_sd_lancet,
      width = 7.2,
      height = 4.5,
      device = "pdf"
    )

    ggsave(
      file.path(out_dir, "heterogeneity_stability_sd_lancet.png"),
      plot = p_stability_sd_lancet,
      width = 7.2,
      height = 4.5,
      dpi = 600
    )

    write.csv(
      stability_df,
      file.path(out_dir, "heterogeneity_stability_summary_lancet.csv"),
      row.names = FALSE
    )

    cat("✓ Saved: heterogeneity_stability_sd_lancet.pdf\n")
    cat("✓ Saved: heterogeneity_stability_sd_lancet.png\n")
    cat("✓ Saved: heterogeneity_stability_summary_lancet.csv\n")
  }, envir = env)
}
