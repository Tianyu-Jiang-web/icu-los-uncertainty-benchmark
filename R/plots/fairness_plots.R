# ============================================================
# Fairness figures (Lancet style).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Fairness disparity main figure.
# Source: 39quantiles.R:1418-1833
stage_plot39_fairness_disparity <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # FAIRNESS DISPARITY PLOT — MAIN FIGURE VERSION
    # Uses existing aggregated objects:
    #   out_gender_crps / out_race_crps / out_insurance_crps
    #
    # Main idea:
    #   delta_crps = subgroup CRPS - mean CRPS within the same
    #                attribute and model
    # So:
    #   delta > 0  : subgroup worse than that model's average for this attribute
    #   delta < 0  : subgroup better than that model's average for this attribute
    # ============================================================

    # ------------------------------------------------------------
    # 1. Helper: convert each fairness table to long format
    # ------------------------------------------------------------
    label_fair <- function(df, attribute) {
      if (is.null(df)) return(NULL)

      df %>%
        rename(subgroup = Group) %>%
        mutate(attribute = attribute) %>%
        pivot_longer(
          cols      = -c(subgroup, attribute),
          names_to  = "model",
          values_to = "crps"
        ) %>%
        mutate(
          subgroup = as.character(subgroup),
          crps     = as.numeric(crps)
        )
    }

    fair_long <- bind_rows(
      label_fair(out_gender_crps,    "Gender"),
      label_fair(out_race_crps,      "Ethnicity")
    )

    if (nrow(fair_long) == 0) {
      cat("⚠ No fairness subgroup tables found — skipping fairness plots.\n")
    } else {

      # ----------------------------------------------------------
      # 2. Model order and readable labels
      # ----------------------------------------------------------
      model_order <- c(
        "DRF", "Ensemble_DRF", "Hybrid_NN_DRF",
        "ENS_Gamma", "MCD_Gamma", "DDNN_Gamma", "BNN_Gamma",
        "QRF", "RF", "XGB"
      )
      model_order <- intersect(model_order, unique(fair_long$model))

      model_label_map <- c(
        DRF           = "DRF",
        Ensemble_DRF  = "Ensemble DRF",
        Hybrid_NN_DRF = "Hybrid",
        ENS_Gamma     = "ENS",
        MCD_Gamma     = "MCD",
        DDNN_Gamma    = "DDNN",
        BNN_Gamma     = "BNN",
        QRF           = "QRF",
        RF            = "RF",
        XGB           = "XGBoost"
      )

      # ----------------------------------------------------------
      # 3. Build subgroup counts from eval_df
      #    (natural counts, not forced CI-style approximation)
      # ----------------------------------------------------------
      count_gender <- NULL
      count_race <- NULL
      count_insurance <- NULL

      if ("gender_group" %in% names(eval_df)) {
        count_gender <- eval_df %>%
          filter(!is.na(gender_group)) %>%
          count(subgroup = as.character(gender_group), name = "n") %>%
          mutate(attribute = "Gender")
      }

      if ("race_group" %in% names(eval_df)) {
        count_race <- eval_df %>%
          filter(!is.na(race_group)) %>%
          count(subgroup = as.character(race_group), name = "n") %>%
          mutate(attribute = "Ethnicity")
      }

      if ("insurance_group" %in% names(eval_df)) {
        count_insurance <- eval_df %>%
          filter(!is.na(insurance_group)) %>%
          count(subgroup = as.character(insurance_group), name = "n") %>%
          mutate(attribute = "Insurance")
      }

      fair_counts <- bind_rows(count_gender, count_race, count_insurance)

      # 如果某些 count 没拿到，就至少保底给 NA
      # ----------------------------------------------------------
      # A. subgroup 标签不要旋转，改成更易读的短格式
      # ----------------------------------------------------------
      fair_long <- fair_long %>%
        left_join(fair_counts, by = c("attribute", "subgroup")) %>%
        mutate(
          subgroup_nlab = ifelse(
            is.na(n),
            subgroup,
            paste0(subgroup, " (n=", n, ")")
          )
        )
      # ----------------------------------------------------------
      # B. 顶刊风格、更稳妥的点图配色
      #    思路：色盲友好 + 同类模型同色系
      # ----------------------------------------------------------
      # Lancet-style 4-group palette:
      #   Hybrid (deep blue, best)  - Ensemble-DRF (mid blue)
      #   NN family (orange/red shades, 4 colours)
      #   Tree family incl. plain DRF (grey shades, 4 colours)
      model_colour_map <- c(
        "Hybrid" = "#00468B",   # deep blue  - best
        "Ensemble DRF"  = "#5DA1DD",   # mid blue

        "ENS"            = "#AD002A",   # darkest red
        "MCD"            = "#ED0000",   # red
        "DDNN"           = "#E18727",   # orange
        "BNN"            = "#FDAF91",   # salmon

        "DRF"            = "#1B1919",   # near-black (worst tree)
        "RF"             = "#5A5A5A",   # dark grey
        "QRF"            = "#8A8A8A",   # mid grey
        "XGBoost"        = "#B8B8B8"    # light grey
      )
      # ----------------------------------------------------------
      # 4. Order subgroup labels within each attribute
      #    Important: do NOT force bad factor conversion too early
      # ----------------------------------------------------------
      gender_order <- c("Female", "Male")
      race_order <- c("White", "Black", "Asian", "Hispanic", "Other")
      insurance_order <- c("Medicare", "Medicaid", "Private", "Government", "Self Pay", "Other")

      # helper: create ordered labels using subgroup values, then map to n-labels
      make_ordered_labels <- function(df_sub, preferred_order = NULL) {
        sg <- unique(df_sub$subgroup)
        if (!is.null(preferred_order)) {
          ordered_sg <- c(intersect(preferred_order, sg), setdiff(sg, preferred_order))
        } else {
          ordered_sg <- sg
        }

        df_sub %>%
          distinct(subgroup, subgroup_nlab) %>%
          filter(subgroup %in% ordered_sg) %>%
          arrange(match(subgroup, ordered_sg)) %>%
          pull(subgroup_nlab)
      }

      gender_levels <- fair_long %>%
        filter(attribute == "Gender", !is.na(subgroup)) %>%
        make_ordered_labels(preferred_order = gender_order)

      race_levels <- fair_long %>%
        filter(attribute == "Ethnicity", !is.na(subgroup)) %>%
        make_ordered_labels(preferred_order = race_order)

      insurance_levels <- fair_long %>%
        filter(attribute == "Insurance", !is.na(subgroup)) %>%
        make_ordered_labels(preferred_order = insurance_order)

      # ----------------------------------------------------------
      # 5. Compute disparity metric
      # ----------------------------------------------------------
      fair_long <- fair_long %>%
        group_by(attribute, model) %>%
        mutate(
          mean_crps_attr_model = mean(crps, na.rm = TRUE),
          delta_crps           = crps - mean_crps_attr_model,
          abs_delta_crps       = abs(delta_crps)
        ) %>%
        ungroup() %>%
        mutate(
          model       = factor(model, levels = model_order),
          model_label = recode(as.character(model), !!!model_label_map),
          model_label = factor(model_label, levels = unname(model_label_map[model_order])),
          attribute   = factor(attribute, levels = c("Gender", "Ethnicity"))
        )

      # ----------------------------------------------------------
      # 6. Create facet-specific y labels safely
      # ----------------------------------------------------------
      fair_plot_df <- bind_rows(
        fair_long %>%
          filter(attribute == "Gender") %>%
          mutate(subgroup_plot = factor(subgroup_nlab, levels = rev(gender_levels))),

        fair_long %>%
          filter(attribute == "Ethnicity") %>%
          mutate(subgroup_plot = factor(subgroup_nlab, levels = rev(race_levels))),

        fair_long %>%
          filter(attribute == "Insurance") %>%
          mutate(subgroup_plot = factor(subgroup_nlab, levels = rev(insurance_levels)))
      )

      # ----------------------------------------------------------
      # 7. Main plot
      # ----------------------------------------------------------
      # ----------------------------------------------------------
      # 7. Bracket geometry (replaces facet strip)
      # ----------------------------------------------------------
      fair_plot_df <- fair_plot_df %>%
        group_by(attribute) %>%
        mutate(y_num = as.integer(subgroup_plot)) %>%
        ungroup() %>%
        mutate(
          y_num_global = as.integer(factor(
            paste(attribute, subgroup_plot),
            levels = unique(paste(attribute, subgroup_plot))
          ))
        )

      # Re-anchor: a single global y order from top to bottom
      global_levels <- fair_plot_df %>%
        arrange(attribute, desc(as.integer(subgroup_plot))) %>%
        pull(subgroup_nlab) %>%
        unique()

      fair_plot_df <- fair_plot_df %>%
        mutate(
          subgroup_global = factor(subgroup_nlab, levels = rev(global_levels)),
          y_num = as.integer(subgroup_global)
        )

      bracket_df <- fair_plot_df %>%
        group_by(attribute) %>%
        summarise(
          y_lo  = min(y_num) - 0.5,
          y_hi  = max(y_num) + 0.5,
          y_mid = (y_lo + y_hi) / 2,
          .groups = "drop"
        )

      # X positions for brackets (data-axis units; x is delta_crps)
      x_range <- range(fair_plot_df$delta_crps, na.rm = TRUE)
      x_span  <- diff(x_range)
      x_tip   <- -0.104   # 紧贴 panel 左边缘外
      x_stem  <- -0.108
      x_label <- -0.115

      bracket_lines <- bracket_df %>%
        rowwise() %>%
        do(data.frame(
          cat = .$attribute,
          x    = c(x_tip,  x_tip,  x_stem),
          xend = c(x_stem, x_stem, x_stem),
          y    = c(.$y_hi, .$y_lo, .$y_hi),
          yend = c(.$y_hi, .$y_lo, .$y_lo)
        )) %>%
        ungroup()

      bracket_stubs <- bracket_df %>%
        transmute(
          x    = x_stem,
          xend = x_label + 0.005 * x_span,
          y    = y_mid,
          yend = y_mid
        )

      # ----------------------------------------------------------
      # 8. Main plot
      # ----------------------------------------------------------
      p_fair <- ggplot(
        fair_plot_df,
        aes(
          x = delta_crps,
          y = subgroup_global,
          colour = model_label
        )
      ) +
        # Zero reference line
        geom_vline(
          xintercept = 0,
          colour = "grey55",
          linewidth = 0.4,
          linetype = "dashed"
        ) +
        # Points (jittered vertically to avoid overlap)
        geom_point(
          size = 2.4,
          alpha = 0.95,
          position = position_jitter(width = 0, height = 0.18, seed = 42)
        ) +
        # Brackets (left)

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
          aes(x = x_label, y = y_mid, label = attribute),
          inherit.aes = FALSE,
          hjust = 1, size = 5.0, colour = "grey10"
        ) +

        scale_colour_manual(values = model_colour_map) +

        scale_x_continuous(
          labels = function(x) gsub("\\.", "\u00b7", sprintf("%.2f", x)),
          expand = expansion(mult = c(0.02, 0.04))
        ) +
        # Subgroup labels on the RIGHT
        scale_y_discrete(position = "right") +

        coord_cartesian(
          xlim = c(-0.10, 0.10),    # 锁住可视范围 = 数据实际范围
          clip = "off"
        ) +

        labs(
          x = expression(Delta * "CRPS (subgroup \u2212 attribute/model mean)"),
          y = NULL,
          colour = NULL
        ) +

        theme_lancet(base_size = 16) +
        theme(
          plot.title    = element_blank(),
          plot.subtitle = element_blank(),

          axis.text.x = element_text(size = 14, colour = "grey15"),
          axis.text.y.right = element_text(
            size = 14, colour = "grey15",
            hjust = 0, margin = ggplot2::margin(l = 6)
          ),
          axis.title.x = element_text(size = 15, colour = "grey10",
                                      margin = ggplot2::margin(t = 8)),
          axis.line.y  = element_blank(),
          axis.ticks.y = element_blank(),

          # Plan B: very faint horizontal grid for row anchoring
          panel.grid.major.y = element_line(colour = "grey95", linewidth = 0.3),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.border       = element_blank(),

          # No facet strips
          strip.background = element_blank(),
          strip.text       = element_blank(),

          legend.position    = "bottom",
          legend.text        = element_text(size = 13, colour = "grey15"),
          legend.title       = element_blank(),
          legend.key.width   = unit(0.7, "cm"),
          legend.spacing.x   = unit(0.25, "cm"),
          legend.box.margin  = ggplot2::margin(t = 4),

          # Wide LEFT margin for brackets + attribute labels
          plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 150)
        )

      print(p_fair)

      # ----------------------------------------------------------
      # 8. Save plot
      # ----------------------------------------------------------
      ggsave(
        file.path(out_dir, "fairness_disparity_plot.pdf"),
        plot   = p_fair,
        width  = 9.4,
        height = 4.5,
        device = "pdf"
      )

      ggsave(
        file.path(out_dir, "fairness_disparity_plot.png"),
        plot   = p_fair,
        width  = 9.4,
        height = 4.5,
        dpi    = 300
      )

      cat("✓ Saved: fairness_disparity_plot.pdf\n")
      cat("✓ Saved: fairness_disparity_plot.png\n")

      # ----------------------------------------------------------
      # 9. Save plotting data
      # ----------------------------------------------------------
      write.csv(
        fair_plot_df %>%
          select(
            attribute,
            subgroup,
            subgroup_nlab,
            n,
            model,
            crps,
            mean_crps_attr_model,
            delta_crps,
            abs_delta_crps
          ) %>%
          arrange(attribute, subgroup, model),
        file.path(out_dir, "fairness_disparity_long.csv"),
        row.names = FALSE
      )

      cat("✓ Saved: fairness_disparity_long.csv\n")
    }
  }, envir = env)
}
