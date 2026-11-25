# ============================================================
# Exploratory plots of ICU LOS.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Exploratory ICU-LOS figures (los_fig1..6) saved to out_dir.
# Source: scripts.R:153-452
stage_eda_plots <- function(env = .GlobalEnv) {
  local({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(patchwork)
    library(scales)
    library(ggridges)

    # ── theme settings ────────────────────────────────────────────────────
    theme_los <- function() {
      theme_minimal(base_family = "sans") +
        theme(
          plot.title         = element_text(size = 13, face = "bold", color = "#1F3864", margin = margin(b = 6)),
          plot.subtitle      = element_text(size = 9,  color = "#595959", margin = margin(b = 8)),
          axis.title         = element_text(size = 9,  color = "#404040"),
          axis.text          = element_text(size = 8,  color = "#404040"),
          axis.text.x        = element_text(angle = 25, hjust = 1),
          strip.text         = element_text(size = 8,  face = "bold", color = "#1F3864"),
          panel.grid.major.y = element_line(color = "#E8E8E8", linewidth = 0.4),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position    = "none",
          plot.background    = element_rect(fill = "white", color = NA),
          panel.background   = element_rect(fill = "#FAFAFA", color = NA),
          plot.margin        = margin(12, 12, 8, 12)
        )
    }

    # color palette
    pal_cat <- c("#2E75B6","#ED7D31","#70AD47","#FFC000","#9B59B6",
                 "#E74C3C","#1ABC9C","#95A5A6","#F39C12","#2980B9",
                 "#8E44AD","#16A085","#D35400","#2C3E50","#27AE60","#C0392B")
    COL_1 <- "#2E75B6"
    COL_0 <- "#ED7D31"

    # full dataset for visualization
    df_vis <- X_tree %>% mutate(y = y)
    y_cap  <- quantile(df_vis$y, 0.99, na.rm = TRUE)
    df_vis <- df_vis %>% mutate(y_plot = pmin(y, y_cap))

    N_lab    <- scales::comma(nrow(df_vis))
    cap_note <- paste0("n = ", N_lab, "  |  LOS capped at 99th pct (",
                       round(y_cap, 1), " days)")

    # ════════════════════════════════════════════════════════════════════════════
    # HELPER 1: violin + boxplot
    # ════════════════════════════════════════════════════════════════════════════
    make_violin <- function(data, var, title) {
      d <- data %>%
        filter(!is.na(.data[[var]])) %>%
        mutate(grp = as.character(.data[[var]])) %>%
        group_by(grp) %>%
        mutate(n     = n(),
               med   = median(y_plot),
               label = paste0(grp, "\n(n=", scales::comma(n), ")")) %>%
        ungroup() %>%
        mutate(label = reorder(label, med))

      cols <- pal_cat[seq_len(n_distinct(d$label))]

      ggplot(d, aes(x = label, y = y_plot, fill = label, color = label)) +
        geom_violin(alpha = 0.30, linewidth = 0.3, trim = TRUE) +
        geom_boxplot(width = 0.18, alpha = 0.85, linewidth = 0.45,
                     outlier.size = 0.5, outlier.alpha = 0.25, color = "grey35") +
        stat_summary(fun = median, geom = "point", shape = 18, size = 2.8, color = "#1F3864") +
        scale_fill_manual(values = cols) +
        scale_color_manual(values = cols) +
        scale_y_continuous(labels = label_number(accuracy = 0.1),
                           expand = expansion(mult = c(0.02, 0.05))) +
        labs(title = title, x = NULL, y = "ICU LOS (days)") +
        theme_los()
    }

    # ════════════════════════════════════════════════════════════════════════════
    # HELPER 2: ridgeline for BINARY variables
    # ════════════════════════════════════════════════════════════════════════════
    make_ridge_binary <- function(data, vars, title, label_map = NULL, min_n = 50) {
      rows <- lapply(vars, function(v) {
        if (!v %in% names(data)) return(NULL)
        d  <- data %>% filter(!is.na(.data[[v]]))
        n1 <- sum(d[[v]] == 1, na.rm = TRUE)
        n0 <- sum(d[[v]] == 0, na.rm = TRUE)
        if (n1 < min_n || n0 < min_n) return(NULL)
        nm <- if (!is.null(label_map) && v %in% names(label_map)) label_map[[v]] else v
        bind_rows(
          d %>% filter(.data[[v]] == 1) %>%
            transmute(y_plot, variable = nm, gtype = "1",
                      n1 = n1, n0 = n0),
          d %>% filter(.data[[v]] == 0) %>%
            transmute(y_plot, variable = nm, gtype = "0",
                      n1 = n1, n0 = n0)
        )
      }) %>% bind_rows()

      if (is.null(rows) || nrow(rows) == 0) return(NULL)

      # Sort by gap (largest at top)
      var_order <- rows %>%
        group_by(variable, gtype) %>%
        summarise(med = median(y_plot), .groups = "drop") %>%
        pivot_wider(names_from = gtype, values_from = med) %>%
        mutate(gap = abs(`1` - `0`)) %>%
        arrange(gap) %>%
        pull(variable)

      rows <- rows %>% mutate(variable = factor(variable, levels = var_order))

      # y-axis labels with n
      n_labels <- rows %>%
        group_by(variable, gtype) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(variable) %>%
        summarise(
          lab = paste0(first(variable),
                       "  Yes=", scales::comma(n[gtype=="1"]),
                       " / No=", scales::comma(n[gtype=="0"])),
          .groups = "drop"
        )
      lev_map <- setNames(n_labels$lab, n_labels$variable)
      rows <- rows %>% mutate(var_lab = lev_map[as.character(variable)],
                              var_lab = factor(var_lab, levels = lev_map[levels(variable)]))

      ggplot(rows, aes(x = y_plot, y = var_lab, fill = gtype, color = gtype)) +
        geom_density_ridges(alpha = 0.45, linewidth = 0.35,
                            quantile_lines = TRUE, quantiles = 0.5,
                            scale = 0.88, rel_min_height = 0.01) +
        scale_fill_manual(values  = c("1" = COL_1, "0" = COL_0)) +
        scale_color_manual(values = c("1" = COL_1, "0" = COL_0)) +
        scale_x_continuous(labels = label_number(accuracy = 0.1),
                           expand = expansion(mult = c(0.01, 0.05))) +
        labs(title    = title,
             subtitle = paste0(cap_note,
                               "  |  Blue = Yes (=1)   Orange = No (=0)   ",
                               "vertical line = median  |  sorted by LOS gap"),
             x = "ICU LOS (days)", y = NULL) +
        theme_los() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.text.y = element_text(size = 7.5))
    }

    # ════════════════════════════════════════════════════════════════════════════
    # HELPER 3: ridgeline for MISSINGNESS
    # ════════════════════════════════════════════════════════════════════════════
    make_ridge_miss <- function(data, top_n = 20) {
      miss_vars <- grep("_miss$", names(data), value = TRUE)

      gaps <- lapply(miss_vars, function(v) {
        d      <- data %>% filter(!is.na(.data[[v]]))
        n_miss <- sum(d[[v]] == 1, na.rm = TRUE)
        n_obs  <- sum(d[[v]] == 0, na.rm = TRUE)
        if (n_miss < 30 || n_obs < 30) return(NULL)
        tibble(v = v, nm = sub("_miss$", "", v),
               gap  = abs(median(d$y_plot[d[[v]] == 1], na.rm = TRUE) -
                            median(d$y_plot[d[[v]] == 0], na.rm = TRUE)),
               n_miss = n_miss, n_obs = n_obs)
      }) %>% bind_rows() %>% arrange(desc(gap)) %>% slice_head(n = top_n)

      rows <- lapply(seq_len(nrow(gaps)), function(i) {
        v  <- gaps$v[i]; nm <- gaps$nm[i]
        n0 <- gaps$n_obs[i]; n1 <- gaps$n_miss[i]
        d  <- data %>% filter(!is.na(.data[[v]]))
        bind_rows(
          d %>% filter(.data[[v]] == 0) %>%
            transmute(y_plot, variable = nm, gtype = "obs",
                      lab = paste0(nm, "  obs=", scales::comma(n0),
                                   " / miss=", scales::comma(n1))),
          d %>% filter(.data[[v]] == 1) %>%
            transmute(y_plot, variable = nm, gtype = "mis",
                      lab = paste0(nm, "  obs=", scales::comma(n0),
                                   " / miss=", scales::comma(n1)))
        )
      }) %>% bind_rows()

      var_order <- gaps %>% arrange(gap) %>% pull(nm)
      lab_order <- lapply(var_order, function(nm) {
        rows %>% filter(variable == nm) %>% pull(lab) %>% first()
      }) %>% unlist()

      rows <- rows %>% mutate(lab = factor(lab, levels = lab_order))

      ggplot(rows, aes(x = y_plot, y = lab, fill = gtype, color = gtype)) +
        geom_density_ridges(alpha = 0.45, linewidth = 0.35,
                            quantile_lines = TRUE, quantiles = 0.5,
                            scale = 0.88, rel_min_height = 0.01) +
        scale_fill_manual(values  = c("obs" = COL_1, "mis" = COL_0)) +
        scale_color_manual(values = c("obs" = COL_1, "mis" = COL_0)) +
        scale_x_continuous(labels = label_number(accuracy = 0.1),
                           expand = expansion(mult = c(0.01, 0.05))) +
        labs(title    = paste0("ICU LOS Distribution by Missingness  (Top ", top_n, " by gap)"),
             subtitle = paste0(cap_note,
                               "  |  Blue = observed (=0)   Orange = missing (=1)   ",
                               "vertical line = median  |  sorted by LOS gap"),
             x = "ICU LOS (days)", y = NULL) +
        theme_los() +
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
              axis.text.y = element_text(size = 7.5))
    }

    # ════════════════════════════════════════════════════════════════════════════
    # BUILD FIGURES
    # ════════════════════════════════════════════════════════════════════════════

    # Fig1/2: violin
    p_gender    <- make_violin(df_vis, "gender",        "Gender")
    p_race      <- make_violin(df_vis, "race_clean",    "Race / Ethnicity")
    p_insurance <- make_violin(df_vis, "insurance",     "Insurance Type")
    p_marital   <- make_violin(df_vis, "marital_status","Marital Status")
    p_careunit  <- make_violin(df_vis, "first_careunit","First Care Unit")

    fig1 <- (p_gender | p_insurance) / (p_marital | p_race) +
      plot_annotation(
        title    = "ICU LOS by Categorical Variables",
        subtitle = cap_note,
        theme = theme(
          plot.title      = element_text(size=15, face="bold", color="#1F3864"),
          plot.subtitle   = element_text(size=9,  color="#595959"),
          plot.background = element_rect(fill="white", color=NA))
      )

    fig2 <- p_careunit +
      labs(title = "ICU LOS by First Care Unit", subtitle = cap_note)

    # Fig3: Comorbidities
    fig3 <- make_ridge_binary(df_vis,
                              vars = c("hypertension","copd_comorb","sepsis","diabetes","ckd","pad",
                                       "obesity","smoking","alcohol_use",
                                       "primary_af","primary_ami","primary_hf_cm",
                                       "primary_renal_failure","primary_copd_bronchitis","primary_stroke_cerebro"),
                              title = "ICU LOS Distribution by Comorbidity / Primary Diagnosis",
                              label_map = c(
                                hypertension             = "Hypertension",
                                copd_comorb              = "COPD",
                                sepsis                   = "Sepsis",
                                diabetes                 = "Diabetes",
                                ckd                      = "CKD",
                                pad                      = "PAD",
                                obesity                  = "Obesity",
                                smoking                  = "Smoking",
                                alcohol_use              = "Alcohol Use",
                                primary_af               = "Primary Dx: AF",
                                primary_ami              = "Primary Dx: AMI",
                                primary_hf_cm            = "Primary Dx: HF/CM",
                                primary_renal_failure    = "Primary Dx: Renal Failure",
                                primary_copd_bronchitis  = "Primary Dx: COPD/Bronchitis",
                                primary_stroke_cerebro   = "Primary Dx: Stroke/Cerebro"
                              )
    )

    # Fig4: Medications
    fig4 <- make_ridge_binary(df_vis,
                              vars = c("acei_arb","anticoagulants","antiplatelets",
                                       "beta_blockers","long_acting_bronchodilator","nephrotoxic","statins"),
                              title = "ICU LOS Distribution by Medication Use",
                              label_map = c(
                                acei_arb                   = "ACEi / ARB",
                                anticoagulants             = "Anticoagulants",
                                antiplatelets              = "Antiplatelets",
                                beta_blockers              = "Beta-Blockers",
                                long_acting_bronchodilator = "Long-Acting Bronchodilator",
                                nephrotoxic                = "Nephrotoxic Agents",
                                statins                    = "Statins"
                              )
    )

    # Fig5: Interventions
    fig5 <- make_ridge_binary(df_vis,
                              vars = c("vent_any","vasopressors","norepinephrine",
                                       "epinephrine","dopamine","vasopressin","CCU_flag","CVICU_flag"),
                              title = "ICU LOS Distribution by Interventions & Unit Type",
                              label_map = c(
                                vent_any       = "Mechanical Ventilation",
                                vasopressors   = "Any Vasopressor",
                                norepinephrine = "Norepinephrine",
                                epinephrine    = "Epinephrine",
                                dopamine       = "Dopamine",
                                vasopressin    = "Vasopressin",
                                CCU_flag       = "CCU",
                                CVICU_flag     = "CVICU"
                              )
    )

    # Fig6: Missingness
    fig6 <- make_ridge_miss(df_vis, top_n = 20)


    # ════════════════════════════════════════════════════════════════════════════
    # SAVE
    # ════════════════════════════════════════════════════════════════════════════
    ggsave(file.path(out_dir, "los_fig1_categorical.png"),   fig1, width=14, height=11, dpi=180, bg="white")
    ggsave(file.path(out_dir, "los_fig2_careunit.png"),      fig2, width=12, height=7,  dpi=180, bg="white")
    ggsave(file.path(out_dir, "los_fig3_comorbidities.png"), fig3, width=12, height=10, dpi=180, bg="white")
    ggsave(file.path(out_dir, "los_fig4_medications.png"),   fig4, width=12, height=7,  dpi=180, bg="white")
    ggsave(file.path(out_dir, "los_fig5_interventions.png"), fig5, width=12, height=7,  dpi=180, bg="white")
    ggsave(file.path(out_dir, "los_fig6_missingness.png"),   fig6, width=12, height=10, dpi=180, bg="white")

    cat("✓ los_fig1_categorical.png\n")
    cat("✓ los_fig2_careunit.png\n")
    cat("✓ los_fig3_comorbidities.png\n")
    cat("✓ los_fig4_medications.png\n")
    cat("✓ los_fig5_interventions.png\n")
    cat("✓ los_fig6_missingness.png\n")
  }, envir = env)
}
