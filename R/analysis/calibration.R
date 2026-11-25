# ============================================================
# Calibration diagnostics (PIT, coverage, density, conditional).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- PIT, transformed PIT, coverage calibration, conditional calibration - 19 levels.
# Source: scripts.R:2225-2575
stage_eval19_calibration <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # CALIBRATION PLOTS (req. 4: removed density distribution, WIS, KS, KL)
    # ============================================================
    library(ggplot2); library(scales)


    # PIT plot
    # ---- 对每个模型计算 PIT 值 ----
    pit_df <- bind_rows(lapply(names(all_qgrid), function(m) {
      qmat <- as.matrix(all_qgrid[[m]])

      pit_vals <- sapply(seq_len(nrow(qmat)), function(i) {
        q_vals <- qmat[i, ]
        y_i    <- y_te[i]

        if (any(is.na(q_vals)) || any(is.infinite(q_vals))) return(NA_real_)

        # 保留你原来的尾部截断逻辑
        if (y_i <= q_vals[1]) return(q_grid[1])
        if (y_i >= q_vals[length(q_vals)]) return(q_grid[length(q_vals)])

        idx <- findInterval(y_i, q_vals)

        if (idx < 1) idx <- 1
        if (idx >= length(q_vals)) idx <- length(q_vals) - 1

        p_lo <- q_grid[idx]
        p_hi <- q_grid[idx + 1]
        q_lo <- q_vals[idx]
        q_hi <- q_vals[idx + 1]

        if (q_hi == q_lo) return((p_lo + p_hi) / 2)

        p_lo + (y_i - q_lo) / (q_hi - q_lo) * (p_hi - p_lo)
      })

      data.frame(
        model = m,
        pit   = pit_vals
      )
    }))

    # ---- transformed PIT ----
    eps <- 1e-6

    pit_df <- pit_df %>%
      filter(!is.na(pit)) %>%
      mutate(
        pit_clipped = pmin(pmax(pit, eps), 1 - eps),
        x_trans     = qnorm(pit_clipped)
      )

    # ---- 统一模型顺序（可选）----
    pit_df$model <- factor(pit_df$model, levels = names(all_qgrid))

    # ============================================================
    # 图1：PIT histogram
    # ============================================================
    p_pit <- ggplot(pit_df, aes(x = pit)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = 20,
        fill = "grey85",
        color = "grey20",
        linewidth = 0.3
      ) +
      geom_hline(
        yintercept = 1,
        linetype = "dashed",
        color = "firebrick",
        linewidth = 0.5
      ) +
      facet_wrap(~model, ncol = 2) +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2),
        expand = c(0.01, 0.01)
      ) +
      labs(
        title = "PIT Histograms",
        subtitle = "A well-calibrated model should produce an approximately uniform PIT distribution",
        x = "PIT",
        y = "Density"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0.9, "lines"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")
      )

    p_pit

    ggsave(
      file.path(out_dir, "pit_histogram.png"),
      plot   = p_pit,
      width  = 14,
      height = ceiling(length(model_names) / 2) * 3.2,
      dpi    = 300
    )
    cat("✓ Saved: pit_histogram.png\n")

    # ============================================================
    # 图2：Transformed PIT histogram
    # ============================================================
    p_trans <- ggplot(pit_df, aes(x = x_trans)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins = 20,
        fill = "grey85",
        color = "grey20",
        linewidth = 0.3
      ) +
      stat_function(
        fun = dnorm,
        color = "black",
        linewidth = 0.6
      ) +
      facet_wrap(~model, ncol = 2) +
      scale_x_continuous(
        limits = c(-3.5, 3.5),
        breaks = seq(-3, 3, by = 1),
        expand = c(0.01, 0.01)
      ) +
      labs(
        title = "Transformed PIT Histograms",
        subtitle = "If PIT is uniform, the transformed values should approximately follow a standard normal distribution",
        x = "Transformed PIT  (qnorm(PIT))",
        y = "Density"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey30"),
        strip.background = element_rect(fill = "grey95", color = "grey70"),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.spacing = unit(0.9, "lines"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black")
      )

    p_trans

    ggsave(
      file.path(out_dir, "transformed_pit_histogram.png"),
      plot   = p_trans,
      width  = 14,
      height = ceiling(length(model_names) / 2) * 3.2,
      dpi    = 300
    )
    cat("✓ Saved: transformed_pit_histogram.png\n")




    # ============================================================
    # 3. Coverage Calibration Plot（单图，全局对比）
    #    原名 "Reliability Diagram" → 改名 "Coverage Calibration Plot"
    #
    #    MCE 嵌入图例标签（"DRF  [MCE=0.021]"），
    #    图例按 MCE 从小到大排序，兼顾全局对比与个体校准误差读取。
    # ============================================================
    rel_df <- bind_rows(lapply(names(all_qgrid), function(m) {
      qmat <- as.matrix(all_qgrid[[m]])
      emp_cov <- sapply(seq_along(alpha_levels), function(j) {
        mean(y_te <= qmat[, j], na.rm = TRUE)
      })
      data.frame(model = m, nominal = alpha_levels, empirical = emp_cov)
    }))

    # --- 计算 MCE ---
    rel_df <- rel_df %>%
      mutate(model = factor(model, levels = model_names))

    mce_df <- rel_df %>%
      group_by(model) %>%
      summarise(
        MCE = mean(abs(empirical - nominal), na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(MCE) %>%                                        # 按 MCE 升序
      mutate(
        legend_label = sprintf("%-14s [MCE=%.3f]", as.character(model), MCE)
      )

    # 把 legend_label 合并回 rel_df，并以 MCE 升序作为图例顺序
    rel_df <- rel_df %>%
      left_join(mce_df, by = "model") %>%
      mutate(
        model_leg = factor(legend_label, levels = mce_df$legend_label)
      )

    # MODEL_COLOR 需要映射到新的 legend_label key
    leg_colors <- setNames(MODEL_COLOR[as.character(mce_df$model)], mce_df$legend_label)

    p_cov_cal <- ggplot(rel_df, aes(x = nominal, y = empirical)) +

      # (i) Perfect calibration diagonal
      geom_abline(
        slope = 1, intercept = 0,
        linetype = "dashed", color = "grey40", linewidth = 0.6
      ) +

      # (iii) Model lines — 图例标签含 MCE
      geom_line(aes(color = model_leg), linewidth = 0.85) +
      geom_point(aes(color = model_leg), size = 1.8) +

      scale_color_manual(
        values = leg_colors,
        name   = "Model  [MCE ↑ = worse]"         # 图例标题提示排序含义
      ) +

      scale_x_continuous(
        breaks = seq(0.1, 0.9, by = 0.1),
        labels = scales::percent_format(accuracy = 1),
        expand = c(0.01, 0.01)
      ) +
      scale_y_continuous(
        breaks = seq(0.1, 0.9, by = 0.1),
        labels = scales::percent_format(accuracy = 1),
        expand = c(0.01, 0.01),
        limits = c(0, 1)
      ) +

      labs(
        title    = "Coverage Calibration Plot",
        subtitle = paste0(
          "Dashed diagonal = perfect calibration  |  ",
          "MCE = Mean Calibration Error (averaged over all quantile levels)  |  ",
          "Legend ordered best \u2192 worst"
        ),
        x = "Nominal quantile level",
        y = "Empirical coverage"
      ) +

      theme_bw(base_size = 12) +
      theme(
        plot.title      = element_text(face = "bold", size = 13),
        plot.subtitle   = element_text(size = 8.5, color = "grey35"),
        legend.position = "right",
        legend.title    = element_text(size = 9, face = "bold"),
        legend.text     = element_text(size = 8.5, family = "mono"),  # 等宽字体对齐
        legend.key.height = unit(1.1, "lines"),
        panel.grid.minor  = element_blank(),
        panel.grid.major  = element_line(color = "#EEEEEE"),
        axis.title        = element_text(size = 10)
      )

    print(p_cov_cal)

    ggsave(
      file.path(out_dir, "coverage_calibration_plot.png"),
      plot   = p_cov_cal,
      width  = 11,
      height = 7,
      dpi    = 300
    )
    cat("✓ Saved: coverage_calibration_plot.png\n")

    # ---- Conditional calibration by PI width quintile ----
    K <- 10  # 改成5组，对应 HRF 论文的 Q1-Q5

    cal_df <- bind_rows(lapply(names(all_q), function(m) {
      eval_df %>%
        transmute(model = m, y,
                  q05 = all_q[[m]]$q05,
                  q95 = all_q[[m]]$q95)
    })) %>%
      mutate(
        width   = q95 - q05,                          # ← 改：用PI宽度替代预测中位数
        covered = (y >= q05) & (y <= q95)
      ) %>%
      group_by(model) %>%
      mutate(bin = ntile(width, K)) %>%               # ← 改：按宽度分五分位
      group_by(model, bin) %>%
      summarise(
        x        = mean(width, na.rm = TRUE),          # ← 改：x轴是该bin的平均宽度
        coverage = mean(covered, na.rm = TRUE),
        n_bin    = n(),
        .groups  = "drop"
      ) %>%
      mutate(
        se    = sqrt(coverage * (1 - coverage) / n_bin),
        lower = pmax(0, coverage - 1.96 * se),
        upper = pmin(1, coverage + 1.96 * se),
        model = factor(model, levels = model_names),
        bin_label = factor(paste0("Q", bin),           # ← 新增：给x轴打Q1-Q5标签
                           levels = paste0("Q", 1:K))
      )

    p_cal <- ggplot(cal_df, aes(x = bin_label, y = coverage, group = model)) +
      geom_hline(yintercept = 0.90, linetype = "dashed", color = "red", alpha = 0.6) +
      geom_ribbon(aes(ymin = lower, ymax = upper, group = model),
                  fill = "steelblue", alpha = 0.15) +
      geom_line(color = "steelblue", linewidth = 0.8) +
      geom_point(color = "steelblue", size = 1.5, alpha = 0.8) +
      facet_wrap(~model, ncol = 2) +
      scale_y_continuous(labels = percent, limits = c(0.55, 1.0)) +
      labs(
        title    = "Conditional Calibration of 90% Prediction Intervals by Uncertainty Quintile",
        subtitle = "Q1 = narrowest intervals (low uncertainty)  |  Q5 = widest intervals (high uncertainty)",
        x        = "Uncertainty Quintile (PI Width)",   # ← 改
        y        = "Empirical Coverage (%)"
      ) +
      theme_bw(base_size = 11) +
      theme(
        strip.background = element_rect(fill = "#f0f0f0"),
        panel.grid.minor = element_blank(),
        plot.title       = element_text(face = "bold"),
        plot.subtitle    = element_text(size = 9, color = "grey40")
      )

    print(p_cal)
    ggsave(file.path(out_dir, "conditional_calibration_uncertainty_quintile.png"),
           plot = p_cal, width = 12, height = 18, dpi = 300)
    cat("✓ Saved: conditional_calibration_uncertainty_quintile.png\n")
    # ---- 配套表格：各模型在每个Decile的覆盖率 ----
    cal_table <- cal_df %>%
      mutate(
        ci_str    = sprintf("%.1f%% [%.1f%%, %.1f%%]",
                            coverage * 100,
                            lower    * 100,
                            upper    * 100)
      ) %>%
      select(model, bin_label, ci_str) %>%
      pivot_wider(names_from = bin_label, values_from = ci_str)

    # 打印到console
    cat("\n=== Conditional Calibration by Uncertainty Decile ===\n")
    cat("Format: Coverage% [95% CI lower%, upper%]\n\n")
    print(cal_table, n = Inf, width = Inf)

    # 保存为CSV
    write.csv(cal_table,
              file.path(out_dir, "conditional_calibration_decile_table.csv"),
              row.names = FALSE)
    cat("✓ Saved: conditional_calibration_decile_table.csv\n")
  }, envir = env)
}

# ---- PIT histograms and PIT statistics.
# Source: 39quantiles.R:1834-2037
stage_eval39_calibration_pit <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # STEP 13: CALIBRATION PLOTS — LANCET-STYLE PIT HISTOGRAMS
    # Layout: 2 rows x 5 cols, models ordered by overall CRPS,
    # 4-group palette (Hybrid / Ens-DRF / NN / Tree)
    # ============================================================
    cat("\n--- Calibration plots ---\n")

    # source("lancet_style.R")   # if not already loaded

    # ── 1. Compute PIT values ────────────────────────────────────
    pit_df <- bind_rows(lapply(names(all_qgrid), function(m) {
      qmat <- as.matrix(all_qgrid[[m]])
      pit_vals <- sapply(seq_len(nrow(qmat)), function(i) {
        q_vals <- qmat[i, ]; y_i <- y_te[i]
        if (any(is.na(q_vals)) || any(is.infinite(q_vals))) return(NA_real_)
        if (y_i <= q_vals[1]) return(q_grid[1])
        if (y_i >= q_vals[length(q_vals)]) return(q_grid[length(q_vals)])
        idx <- findInterval(y_i, q_vals)
        idx <- pmax(1, pmin(idx, length(q_vals) - 1))
        p_lo <- q_grid[idx]; p_hi <- q_grid[idx + 1]
        q_lo <- q_vals[idx]; q_hi <- q_vals[idx + 1]
        if (q_hi == q_lo) return((p_lo + p_hi) / 2)
        p_lo + (y_i - q_lo) / (q_hi - q_lo) * (p_hi - p_lo)
      })
      data.frame(model = m, pit = pit_vals)
    }))

    eps    <- 1e-6
    pit_df <- pit_df %>%
      filter(!is.na(pit)) %>%
      mutate(pit_clipped = pmin(pmax(pit, eps), 1 - eps),
             x_trans     = qnorm(pit_clipped))

    # Order models by PIT calibration quality:
    # KS statistic against Uniform(0,1) — smaller = better calibrated
    pit_model_order <- pit_df %>%
      group_by(model) %>%
      summarise(ks_stat = ks.test(pit, "punif")$statistic, .groups = "drop") %>%
      arrange(ks_stat) %>%
      pull(model) %>%
      as.character()

    cat("PIT calibration ranking (best -> worst):\n")
    print(pit_model_order)

    pit_model_labels <- c(
      "Hybrid_NN_DRF" = "Hybrid",
      "MCD_Gamma"     = "MCD",
      "QRF"           = "QRF",
      "XGB"           = "XGBoost",
      "Ensemble_DRF"  = "Ensemble DRF",
      "ENS_Gamma"     = "ENS",
      "DDNN_Gamma"    = "DDNN",
      "BNN_Gamma"     = "BNN",
      "RF"            = "RF",
      "DRF"           = "DRF"
    )

    # ── 3. Lancet 4-group palette (matches fairness/heatmap figures) ─
    pit_colours <- c(
      "Hybrid"       = "#00468B",   # deep blue   - best
      "Ensemble DRF" = "#5DA1DD",   # mid blue
      "ENS"          = "#AD002A",   # darkest red
      "MCD"          = "#ED0000",   # red
      "DDNN"         = "#E18727",   # orange
      "BNN"          = "#FDAF91",   # salmon
      "DRF"          = "#1B1919",   # near-black  - worst tree
      "RF"           = "#5A5A5A",   # dark grey
      "QRF"          = "#8A8A8A",   # mid grey
      "XGBoost"      = "#B8B8B8"    # light grey
    )

    # Apply ordering + readable labels to pit_df
    pit_df <- pit_df %>%
      filter(model %in% pit_model_order) %>%
      mutate(
        model       = factor(model, levels = pit_model_order),
        model_label = recode(as.character(model), !!!pit_model_labels),
        model_label = factor(model_label,
                             levels = unname(pit_model_labels[pit_model_order]))
      )

    # ── 4. PIT histogram figure ──────────────────────────────────
    p_pit <- ggplot(pit_df, aes(x = pit, fill = model_label)) +
      geom_histogram(
        aes(y = after_stat(density)),
        bins      = 20,
        colour    = "grey25",
        linewidth = 0.25,
        alpha     = 0.85
      ) +
      geom_hline(
        yintercept = 1,
        linetype   = "dashed",
        colour     = "grey35",
        linewidth  = 0.45
      ) +
      facet_wrap(~ model_label, nrow = 2, ncol = 5) +
      scale_fill_manual(values = pit_colours, guide = "none") +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0", "0\u00b725", "0\u00b750", "0\u00b775", "1\u00b700"),
        expand = c(0.01, 0.01)
      ) +
      scale_y_continuous(
        labels = function(x) gsub("\\.", "\u00b7", sprintf("%.1f", x)),
        expand = expansion(mult = c(0, 0.06))
      ) +
      labs(x = "PIT", y = "Density") +
      theme_lancet(base_size = 21) +
      theme(
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),

        # Subplot titles (model names)
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text       = element_text(
          size = 23, face = "plain", colour = "grey10",
          margin = ggplot2::margin(b = 4, t = 2)
        ),

        # Axes
        axis.text   = element_text(size = 15, colour = "grey20"),
        axis.title  = element_text(size = 23, colour = "grey10"),
        axis.line   = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks  = element_line(colour = "grey25", linewidth = 0.4),

        panel.spacing.x = unit(0.7, "lines"),
        panel.spacing.y = unit(1.0, "lines"),
        panel.grid      = element_blank(),

        plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 10)
      )

    print(p_pit)

    # ── 5. Save figure ───────────────────────────────────────────
    ggsave(
      file.path(out_dir, "pit_histogram.pdf"),
      plot = p_pit, width = 15, height = 6.5, device = "pdf"
    )
    ggsave(
      file.path(out_dir, "pit_histogram.png"),
      plot = p_pit, width = 15, height = 6.5, dpi = 600, bg = "white"
    )
    cat("\u2713 Saved: pit_histogram.pdf / .png\n")

    # ── 6. PIT CSV exports (unchanged) ────────────────────────────
    write.csv(
      pit_df %>% select(model, pit, x_trans) %>% arrange(model),
      file.path(out_dir, "pit_raw.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: pit_raw.csv\n")

    pit_hist_summary <- pit_df %>%
      group_by(model) %>%
      group_modify(~ {
        h <- hist(.x$pit, breaks = seq(0, 1, length.out = 21), plot = FALSE)
        data.frame(
          bin_mid   = h$mids,
          bin_lower = h$breaks[-length(h$breaks)],
          bin_upper = h$breaks[-1],
          count     = h$counts,
          density   = h$density
        )
      }) %>%
      ungroup()

    write.csv(
      pit_hist_summary,
      file.path(out_dir, "pit_histogram_summary.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: pit_histogram_summary.csv\n")

    pit_stats <- pit_df %>%
      group_by(model) %>%
      summarise(
        n          = n(),
        mean_pit   = mean(pit, na.rm = TRUE),
        sd_pit     = sd(pit, na.rm = TRUE),
        ks_stat    = ks.test(pit, "punif")$statistic,
        ks_pvalue  = ks.test(pit, "punif")$p.value,
        tail_ratio = mean(pit < 0.1 | pit > 0.9) / mean(pit >= 0.1 & pit <= 0.9) * (0.8 / 0.2),
        .groups    = "drop"
      ) %>%
      mutate(
        calibration_diagnosis = case_when(
          ks_pvalue > 0.05                        ~ "Well-calibrated",
          tail_ratio > 1.2                        ~ "Overconfident (U-shaped)",
          tail_ratio < 0.8                        ~ "Underconfident (hump-shaped)",
          TRUE                                    ~ "Miscalibrated (other)"
        )
      ) %>%
      arrange(ks_stat)

    write.csv(
      pit_stats,
      file.path(out_dir, "pit_stats_per_model.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: pit_stats_per_model.csv\n")
  }, envir = env)
}

# ---- Empirical-vs-nominal coverage calibration plot and MCE.
# Source: 39quantiles.R:2039-2210
stage_eval39_calibration_coverage <- function(env = .GlobalEnv) {
  local({
    # ── Coverage Calibration Plot — LANCET STYLE (4 representatives) ──
    # Story: two combination models (Hybrid, Ensemble-DRF) vs one tree
    #        baseline (QRF) vs one NN baseline (BNN)
    # source("lancet_style.R")   # if not already loaded

    # ── 1. Compute empirical coverage at each nominal level ──────
    alpha_levels <- q_grid
    rel_df <- bind_rows(lapply(names(all_qgrid), function(m) {
      qmat    <- as.matrix(all_qgrid[[m]])
      emp_cov <- sapply(seq_along(alpha_levels), function(j)
        mean(y_te <= qmat[, j], na.rm = TRUE))
      data.frame(model = m, nominal = alpha_levels, empirical = emp_cov)
    }))

    # ── 2. Per-model MCE for ALL models (CSV exports use this) ───
    mce_df_all <- rel_df %>%
      group_by(model) %>%
      summarise(MCE = mean(abs(empirical - nominal), na.rm = TRUE),
                .groups = "drop") %>%
      arrange(MCE)

    # ── 3. All 10 models — unified Lancet palette (consistent with all other figures)
    cov_label_map <- c(
      "Hybrid_NN_DRF" = "Hybrid",
      "Ensemble_DRF"  = "Ensemble DRF",
      "ENS_Gamma"     = "ENS",
      "MCD_Gamma"     = "MCD",
      "DDNN_Gamma"    = "DDNN",
      "BNN_Gamma"     = "BNN",
      "DRF"           = "DRF",
      "RF"            = "RF",
      "QRF"           = "QRF",
      "XGB"           = "XGBoost"
    )

    # Global Lancet palette — identical to CRPS / sensitivity figures
    cov_colours_10 <- c(
      "Hybrid"       = "#00468B",   # deep blue
      "Ensemble DRF" = "#5DA1DD",   # mid blue
      "ENS"          = "#AD002A",   # darkest red
      "MCD"          = "#ED0000",   # red
      "DDNN"         = "#E18727",   # orange
      "BNN"          = "#FDAF91",   # salmon
      "DRF"          = "#1B1919",   # near-black
      "RF"           = "#5A5A5A",   # dark grey
      "QRF"          = "#8A8A8A",   # mid grey
      "XGBoost"      = "#B8B8B8"    # light grey
    )

    # ── 4. Use all models + build legend labels with aligned MCE ──
    rel_sub <- rel_df %>%
      left_join(mce_df_all, by = "model") %>%
      mutate(model_clean = recode(as.character(model), !!!cov_label_map))

    # Sort by MCE ascending (best calibration on top of legend)
    mce_sub <- rel_sub %>%
      distinct(model, model_clean, MCE) %>%
      arrange(MCE)

    # Pad model name to fixed width so '[' aligns under monospace font
    name_width <- max(nchar(mce_sub$model_clean)) + 2
    mce_sub <- mce_sub %>%
      mutate(
        mce_str      = gsub("\\.", "\u00b7", sprintf("%.3f", MCE)),
        legend_label = sprintf(paste0("%-", name_width, "s [%s]"),
                               model_clean, mce_str)
      )

    rel_sub <- rel_sub %>%
      left_join(mce_sub %>% select(model_clean, legend_label), by = "model_clean") %>%
      mutate(model_leg = factor(legend_label, levels = mce_sub$legend_label))

    # Map colours by ordered legend label
    leg_colours <- setNames(
      cov_colours_10[mce_sub$model_clean],
      mce_sub$legend_label
    )

    # ── 5. Plot ──────────────────────────────────────────────────
    p_cov_cal <- ggplot(rel_sub,
                        aes(x = nominal, y = empirical,
                            group = model_leg, colour = model_leg)) +
      geom_abline(
        slope = 1, intercept = 0, linetype = "dashed",
        colour = "grey40", linewidth = 0.5
      ) +
      geom_line(linewidth = 1.0, alpha = 0.95, lineend = "round") +
      scale_colour_manual(values = leg_colours, name = "Model [MCE]") +
      scale_x_continuous(
        breaks = seq(0.1, 0.9, 0.1),
        labels = scales::percent_format(accuracy = 1),
        expand = c(0.01, 0.01),
        limits = c(0, 1)
      ) +
      scale_y_continuous(
        breaks = seq(0.1, 0.9, 0.1),
        labels = scales::percent_format(accuracy = 1),
        expand = c(0.01, 0.01),
        limits = c(0, 1)
      ) +
      labs(
        x = "Nominal quantile level",
        y = "Empirical coverage"
      ) +
      theme_lancet(base_size = 20) +
      theme(
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),

        legend.position   = "right",
        legend.text       = element_text(size = 14,
                                         family = "mono",
                                         colour = "grey15"),
        legend.title      = element_text(size = 17, colour = "grey15",
                                         face = "plain", family = "",
                                         margin = ggplot2::margin(b = 6)),
        legend.key.height = unit(1.2, "lines"),
        legend.key.width  = unit(1.2, "cm"),
        legend.spacing.y  = unit(0.2, "cm"),

        axis.text  = element_text(size = 14, colour = "grey20"),
        axis.title = element_text(size = 18, colour = "grey10"),
        axis.line  = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks = element_line(colour = "grey25", linewidth = 0.4),

        panel.grid = element_blank(),

        plot.margin = ggplot2::margin(t = 8, r = 12, b = 8, l = 10)
      )

    print(p_cov_cal)

    # ── 6. Save figure ───────────────────────────────────────────
    ggsave(
      file.path(out_dir, "coverage_calibration_plot.pdf"),
      plot = p_cov_cal, width = 12, height = 6.5, device = "pdf"
    )
    ggsave(
      file.path(out_dir, "coverage_calibration_plot.png"),
      plot = p_cov_cal, width = 12, height = 6.5, dpi = 600, bg = "white"
    )
    cat("\u2713 Saved: coverage_calibration_plot.pdf / .png\n")

    # ── 7. CSV exports (FULL 10 models, not just plotted 4) ──────
    write.csv(
      rel_df %>%
        left_join(mce_df_all, by = "model") %>%
        select(model, nominal, empirical, MCE) %>%
        arrange(model, nominal),
      file.path(out_dir, "coverage_calibration_long.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: coverage_calibration_long.csv\n")

    coverage_wide <- rel_df %>%
      select(model, nominal, empirical) %>%
      mutate(model = as.character(model)) %>%
      tidyr::pivot_wider(names_from = model, values_from = empirical) %>%
      arrange(nominal)
    write.csv(
      coverage_wide,
      file.path(out_dir, "coverage_calibration_wide.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: coverage_calibration_wide.csv\n")

    write.csv(
      mce_df_all %>% select(model, MCE) %>% arrange(MCE),
      file.path(out_dir, "coverage_mce_per_model.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: coverage_mce_per_model.csv\n")
  }, envir = env)
}

# ---- Average predictive densities vs actual LOS.
# Source: 39quantiles.R:2213-2538
stage_eval39_predictive_density <- function(env = .GlobalEnv) {
  local({
    # ── Predictive Density vs Actual LOS (robust smoothed per-patient) ──────


    x_grid <- seq(0, 21, length.out = 200)

    # ---- 小工具：梯形积分 ----
    trapz_area <- function(x, y) {
      idx <- is.finite(x) & is.finite(y)
      x <- x[idx]
      y <- y[idx]
      if (length(x) < 2) return(NA_real_)
      sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
    }

    # ---- 小工具：从单个病人的 quantiles 生成平滑密度 ----
    make_density_from_quantiles <- function(qs, q_grid, x_grid,
                                            spar_q = 0.55,
                                            spar_d = 0.40,
                                            min_slope = 1e-3) {
      # 1) 基本检查
      if (length(qs) != length(q_grid)) {
        stop("length(qs) != length(q_grid). 请检查 all_qgrid 的列数和 q_grid 是否一致。")
      }

      if (any(!is.finite(qs))) {
        return(rep(NA_real_, length(x_grid)))
      }

      # 2) 强制单调，避免 crossing
      qs_mon <- cummax(qs)

      # 如果完全常数，直接给一个很窄但有限的近似密度
      if (diff(range(qs_mon)) < 1e-8) {
        dens <- dnorm(x_grid, mean = qs_mon[1], sd = 0.15)
        area <- trapz_area(x_grid, dens)
        if (is.finite(area) && area > 0) dens <- dens / area
        return(dens)
      }

      tau <- q_grid

      # 3) 平滑 quantile curve Q(tau)
      fit_q <- tryCatch(
        smooth.spline(x = tau, y = qs_mon, spar = spar_q),
        error = function(e) NULL
      )

      if (is.null(fit_q)) {
        return(rep(NA_real_, length(x_grid)))
      }

      tau_fine <- seq(min(tau), max(tau), length.out = 400)
      q_smooth <- tryCatch(
        predict(fit_q, x = tau_fine)$y,
        error = function(e) rep(NA_real_, length(tau_fine))
      )

      if (all(!is.finite(q_smooth))) {
        return(rep(NA_real_, length(x_grid)))
      }

      q_smooth <- cummax(q_smooth)

      # 4) 导数 Q'(tau)
      dq_dtau <- diff(q_smooth) / diff(tau_fine)
      dq_dtau <- pmax(dq_dtau, min_slope)

      # 5) 密度 f(x)=1/Q'(tau)
      dens_seg <- 1 / dq_dtau
      x_mid <- (q_smooth[-1] + q_smooth[-length(q_smooth)]) / 2

      # 6) 去掉非法值
      keep <- is.finite(x_mid) & is.finite(dens_seg)
      x_mid <- x_mid[keep]
      dens_seg <- dens_seg[keep]

      if (length(x_mid) < 4) {
        return(rep(NA_real_, length(x_grid)))
      }

      # 7) x_mid 可能重复，先按相同 x 聚合
      df_tmp <- data.frame(x = x_mid, d = dens_seg) %>%
        group_by(x) %>%
        summarise(d = mean(d), .groups = "drop") %>%
        arrange(x)

      x_mid_u <- df_tmp$x
      dens_u  <- df_tmp$d

      if (length(x_mid_u) < 4) {
        # 不足以再 spline，就直接线性插值
        dens_interp <- approx(
          x = x_mid_u,
          y = dens_u,
          xout = x_grid,
          rule = 1,
          ties = mean
        )$y
        dens_interp[is.na(dens_interp)] <- 0
        dens_interp <- pmax(dens_interp, 0)

        area <- trapz_area(x_grid, dens_interp)
        if (is.finite(area) && area > 0) dens_interp <- dens_interp / area
        return(dens_interp)
      }

      # 8) 对 density 再轻微平滑（可能失败，所以 tryCatch）
      dens_smooth <- tryCatch({
        fit_d <- smooth.spline(x = x_mid_u, y = dens_u, spar = spar_d)
        predict(fit_d, x = x_mid_u)$y
      }, error = function(e) dens_u)

      dens_smooth <- pmax(dens_smooth, 0)

      # 9) 插值到统一 x_grid
      dens_interp <- approx(
        x = x_mid_u,
        y = dens_smooth,
        xout = x_grid,
        rule = 1,
        ties = mean
      )$y

      dens_interp[is.na(dens_interp)] <- 0
      dens_interp <- pmax(dens_interp, 0)

      # 10) 归一化
      area <- trapz_area(x_grid, dens_interp)
      if (is.finite(area) && area > 0) {
        dens_interp <- dens_interp / area
      }

      dens_interp
    }

    # ---- 主体：对每个模型求平均预测密度 ----
    dens_list <- lapply(names(all_qgrid), function(m) {

      qmat <- as.matrix(all_qgrid[[m]])

      # 再检查一次 q_grid 与列数
      if (ncol(qmat) != length(q_grid)) {
        stop(sprintf(
          "Model %s: ncol(qmat) = %d, but length(q_grid) = %d。两者必须一致。",
          m, ncol(qmat), length(q_grid)
        ))
      }

      dens_per_obs <- apply(qmat, 1, function(qs) {
        make_density_from_quantiles(
          qs     = qs,
          q_grid = q_grid,
          x_grid = x_grid,
          spar_q = 0.55,
          spar_d = 0.40,
          min_slope = 1e-3
        )
      })

      # apply 在某些极端情况下可能返回向量而不是矩阵，这里强制处理一下
      if (is.null(dim(dens_per_obs))) {
        dens_per_obs <- matrix(dens_per_obs, ncol = 1)
      }

      avg_dens <- rowMeans(dens_per_obs, na.rm = TRUE)

      # 如果某些 x 上全是 NA，rowMeans(na.rm=TRUE) 可能给 NaN
      avg_dens[!is.finite(avg_dens)] <- 0

      # 再归一化一次，确保群体平均密度也积分为1
      area_avg <- trapz_area(x_grid, avg_dens)
      if (is.finite(area_avg) && area_avg > 0) {
        avg_dens <- avg_dens / area_avg
      }

      data.frame(
        x = x_grid,
        density = avg_dens,
        model = m
      )
    })

    dens_df <- bind_rows(dens_list)

    # ---- 真实 LOS 密度 ----
    actual_dens <- density(y_te, from = 0, to = 21, n = 200, bw = "SJ")
    actual_df <- data.frame(
      x = actual_dens$x,
      density = actual_dens$y
    )

    # ---- Lancet style: clean labels + 4-group palette ----
    dens_label_map <- c(
      "Hybrid_NN_DRF" = "Hybrid",
      "MCD_Gamma"     = "MCD",
      "QRF"           = "QRF",
      "XGB"           = "XGBoost",
      "Ensemble_DRF"  = "Ensemble DRF",
      "ENS_Gamma"     = "ENS",
      "DDNN_Gamma"    = "DDNN",
      "BNN_Gamma"     = "BNN",
      "RF"            = "RF",
      "DRF"           = "DRF"
    )

    dens_colours <- c(
      "Hybrid"       = "#00468B",
      "Ensemble DRF" = "#5DA1DD",
      "ENS"          = "#AD002A",
      "MCD"          = "#ED0000",
      "DDNN"         = "#E18727",
      "BNN"          = "#FDAF91",
      "DRF"          = "#1B1919",
      "RF"           = "#5A5A5A",
      "QRF"          = "#8A8A8A",
      "XGBoost"      = "#B8B8B8"
    )

    # Model order = same CRPS ranking used in heatmap / PIT figures
    dens_panel_order <- c(
      "Hybrid_NN_DRF", "MCD_Gamma", "QRF", "XGB", "Ensemble_DRF",
      "ENS_Gamma", "DDNN_Gamma", "BNN_Gamma", "RF", "DRF"
    )
    dens_panel_order <- intersect(dens_panel_order, unique(dens_df$model))

    dens_df <- dens_df %>%
      mutate(
        model       = factor(model, levels = dens_panel_order),
        model_label = recode(as.character(model), !!!dens_label_map),
        model_label = factor(model_label,
                             levels = unname(dens_label_map[dens_panel_order]))
      )

    # ---- Plot ----
    p_dens <- ggplot() +
      # Empirical LOS (grey baseline, same in every panel)
      geom_area(
        data = actual_df,
        aes(x = x, y = density),
        fill = "grey75", alpha = 0.55
      ) +
      geom_line(
        data = actual_df,
        aes(x = x, y = density),
        colour = "grey35", linewidth = 0.55
      ) +
      # Model-specific predicted density
      geom_area(
        data = dens_df,
        aes(x = x, y = density, fill = model_label),
        alpha = 0.35
      ) +
      geom_line(
        data = dens_df,
        aes(x = x, y = density, colour = model_label),
        linewidth = 1.0
      ) +
      facet_wrap(~ model_label, nrow = 2, ncol = 5) +
      scale_fill_manual(values = dens_colours, guide = "none") +
      scale_colour_manual(values = dens_colours, guide = "none") +
      scale_x_continuous(
        limits = c(0, 21),
        breaks = c(0, 5, 10, 15, 20),
        expand = c(0.01, 0.01)
      ) +
      scale_y_continuous(
        labels = function(x) gsub("\\.", "\u00b7", sprintf("%.2f", x)),
        expand = expansion(mult = c(0, 0.05))
      ) +
      coord_cartesian(xlim = c(0, 21)) +
      labs(
        x = "Length of stay (days)",
        y = "Density"
      ) +
      theme_lancet(base_size = 16) +
      theme(
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),

        strip.background = element_rect(fill = "white", colour = NA),
        strip.text       = element_text(
          size = 15, face = "plain", colour = "grey10",
          margin = ggplot2::margin(b = 4, t = 2)
        ),

        axis.text  = element_text(size = 13, colour = "grey20"),
        axis.title = element_text(size = 15, colour = "grey10"),
        axis.line  = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks = element_line(colour = "grey25", linewidth = 0.4),

        panel.spacing.x = unit(0.7, "lines"),
        panel.spacing.y = unit(1.0, "lines"),
        panel.grid      = element_blank(),

        plot.margin = ggplot2::margin(t = 10, r = 14, b = 10, l = 12)
      )

    print(p_dens)

    ggsave(
      file.path(out_dir, "predictive_density_vs_actual.pdf"),
      plot = p_dens, width = 16, height = 7.5, device = "pdf"
    )
    ggsave(
      file.path(out_dir, "predictive_density_vs_actual.png"),
      plot = p_dens, width = 16, height = 7.5, dpi = 600, bg = "white"
    )
    cat("\u2713 Saved: predictive_density_vs_actual.pdf / .png\n")

    # ── Export density data for interactive HTML ────────────────
    # 1. Model average predicted densities (200 x-grid points × 10 models)
    write.csv(
      dens_df,   # already contains: x, density, model
      file.path(out_dir, "density_model_avg.csv"),
      row.names = FALSE
    )
    cat("✓ Saved: density_model_avg.csv\n")

    # 2. Actual LOS empirical density (200 points)
    write.csv(
      actual_df,  # already contains: x, density
      file.path(out_dir, "density_actual_los.csv"),
      row.names = FALSE
    )
    cat("✓ Saved: density_actual_los.csv\n")
  }, envir = env)
}

# ---- Conditional calibration by prediction-interval-width decile.
# Source: 39quantiles.R:2540-2694
stage_eval39_calibration_conditional <- function(env = .GlobalEnv) {
  local({
    # ── Conditional Calibration by PI-width decile — LANCET STYLE ──
    # source("lancet_style.R")   # if not already loaded

    K_bins <- 10
    cal_df <- bind_rows(lapply(names(all_q), function(m) {
      eval_df %>%
        transmute(model = m, y, q05 = all_q[[m]]$q05, q95 = all_q[[m]]$q95)
    })) %>%
      mutate(width = q95 - q05, covered = (y >= q05) & (y <= q95)) %>%
      group_by(model) %>%
      mutate(bin = ntile(width, K_bins)) %>%
      group_by(model, bin) %>%
      summarise(
        x         = mean(width, na.rm = TRUE),
        coverage  = mean(covered, na.rm = TRUE),
        n_bin     = n(),
        .groups   = "drop"
      ) %>%
      mutate(
        se        = sqrt(coverage * (1 - coverage) / n_bin),
        lower     = pmax(0, coverage - 1.96 * se),
        upper     = pmin(1, coverage + 1.96 * se),
        bin_label = factor(paste0("Q", bin), levels = paste0("Q", 1:K_bins))
      )

    # ── Clean model names + Lancet 4-group palette ──────────────
    cal_label_map <- c(
      "Hybrid_NN_DRF" = "Hybrid",
      "MCD_Gamma"     = "MCD",
      "QRF"           = "QRF",
      "XGB"           = "XGBoost",
      "Ensemble_DRF"  = "Ensemble DRF",
      "ENS_Gamma"     = "ENS",
      "DDNN_Gamma"    = "DDNN",
      "BNN_Gamma"     = "BNN",
      "RF"            = "RF",
      "DRF"           = "DRF"
    )

    cal_colours <- c(
      "Hybrid"       = "#00468B",
      "Ensemble DRF" = "#5DA1DD",
      "ENS"          = "#AD002A",
      "MCD"          = "#ED0000",
      "DDNN"         = "#E18727",
      "BNN"          = "#FDAF91",
      "DRF"          = "#1B1919",
      "RF"           = "#5A5A5A",
      "QRF"          = "#8A8A8A",
      "XGBoost"      = "#B8B8B8"
    )

    # Order models by overall CRPS ranking (matches other figures)
    cal_panel_order <- c(
      "Hybrid_NN_DRF", "MCD_Gamma", "QRF", "XGB", "Ensemble_DRF",
      "ENS_Gamma", "DDNN_Gamma", "BNN_Gamma", "RF", "DRF"
    )
    cal_panel_order <- intersect(cal_panel_order, unique(cal_df$model))

    cal_df <- cal_df %>%
      mutate(
        model       = factor(model, levels = cal_panel_order),
        model_label = recode(as.character(model), !!!cal_label_map),
        model_label = factor(model_label,
                             levels = unname(cal_label_map[cal_panel_order]))
      )

    # ── Figure ──────────────────────────────────────────────────
    p_cal <- ggplot(cal_df,
                    aes(x = bin_label, y = coverage,
                        group = model_label,
                        colour = model_label, fill = model_label)) +
      # 90% reference (perfect calibration target)
      geom_hline(yintercept = 0.90, linetype = "dashed",
                 colour = "grey35", linewidth = 0.6) +
      # Wilson 95% CI ribbon
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.25, colour = NA) +
      geom_line(linewidth = 1.0) +
      geom_point(size = 2.5, shape = 19) +

      scale_colour_manual(values = cal_colours, guide = "none") +
      scale_fill_manual(values = cal_colours, guide = "none") +

      facet_wrap(~ model_label, ncol = 2) +

      scale_y_continuous(
        labels = scales::percent_format(accuracy = 1),
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.1),
        expand = expansion(mult = c(0.01, 0.02))
      ) +
      labs(
        x = "Uncertainty decile (prediction interval width)",
        y = "Empirical coverage"
      ) +
      theme_lancet(base_size = 20) +
      theme(
        plot.title    = element_blank(),
        plot.subtitle = element_blank(),

        # Subplot titles (model names)
        strip.background = element_rect(fill = "white", colour = NA),
        strip.text = element_text(
          size = 20, face = "plain", colour = "grey10",
          margin = ggplot2::margin(b = 4, t = 2)
        ),

        # Faint anchoring grid lines around the 90% target
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),

        # Axes
        axis.text   = element_text(size = 17, colour = "grey20"),
        axis.title  = element_text(size = 21, colour = "grey10"),
        axis.line   = element_line(colour = "grey25", linewidth = 0.4),
        axis.ticks  = element_line(colour = "grey25", linewidth = 0.4),

        panel.spacing.x = unit(0.7, "lines"),
        panel.spacing.y = unit(1.0, "lines"),

        plot.margin = ggplot2::margin(t = 10, r = 14, b = 10, l = 12)
      )

    print(p_cal)

    # ── Save figure (large canvas: 5 rows x 2 cols, y from 0 needs space) ──
    ggsave(
      file.path(out_dir, "conditional_calibration_uncertainty_quintile.pdf"),
      plot = p_cal, width = 14.5, height = 21, device = "pdf"
    )
    ggsave(
      file.path(out_dir, "conditional_calibration_uncertainty_quintile.png"),
      plot = p_cal, width = 14.5, height = 21, dpi = 300, bg = "white"
    )
    cat("\u2713 Saved: conditional_calibration_uncertainty_quintile.pdf / .png\n")

    # ── Wide CSV with confidence intervals (unchanged) ──────────
    cal_table <- cal_df %>%
      mutate(ci_str = sprintf("%.1f%% [%.1f%%, %.1f%%]",
                              coverage * 100, lower * 100, upper * 100)) %>%
      select(model, bin_label, ci_str) %>%
      pivot_wider(names_from = bin_label, values_from = ci_str)

    cat("\n=== Conditional Calibration by Uncertainty Decile ===\n")
    cat("Format: Coverage% [95% CI lower%, upper%]\n\n")
    print(cal_table, n = Inf, width = Inf)

    write.csv(
      cal_table,
      file.path(out_dir, "conditional_calibration_decile_table.csv"),
      row.names = FALSE
    )
    cat("\u2713 Saved: conditional_calibration_decile_table.csv\n")
  }, envir = env)
}
