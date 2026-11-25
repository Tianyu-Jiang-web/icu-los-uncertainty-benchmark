# ============================================================
# Baseline characteristics table (Table 1).
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Table 1 (gtsummary/gt) -> HTML + CSV in base_fig_dir.
# Source: 39quantiles.R:198-450
stage_table1 <- function(env = .GlobalEnv) {
  local({
    library(gt)
    library(gtsummary)

    # ---- ICU LOS strata (保留用于内部分析，不进主表) ----
    df_baseline <- df2 %>%
      mutate(
        los_group = case_when(
          y <  2           ~ "Short (<2 days)",
          y >= 2 & y < 7   ~ "Medium (2-6 days)",
          y >= 7           ~ "Long (>=7 days)",
          TRUE             ~ NA_character_
        ),
        los_group = factor(los_group,
                           levels = c("Short (<2 days)",
                                      "Medium (2-6 days)",
                                      "Long (>=7 days)"))
      )

    # ---- 变量列表（去掉 icu_sbp_last24h）----
    baseline_vars <- c(
      "y", "anchor_age", "gender", "race_clean", "insurance",
      "first_careunit", "sepsis", "primary_stroke_cerebro",
      "gcs_total_last24h", "icu_hr_last24h",
      "Creatinine", "Albumin", "vent_any", "vasopressors"
    )

    var_labels <- list(
      y                      ~ "ICU LOS, days",
      anchor_age             ~ "Age, years",
      gender                 ~ "Sex",
      race_clean             ~ "Race/ethnicity",
      insurance              ~ "Insurance",
      first_careunit         ~ "First care unit",
      sepsis                 ~ "Sepsis",
      primary_stroke_cerebro ~ "Stroke/cerebrovascular disease",
      gcs_total_last24h      ~ "GCS total",
      icu_hr_last24h         ~ "Heart rate",
      Creatinine             ~ "Creatinine",
      Albumin                ~ "Albumin",
      vent_any               ~ "Mechanical ventilation",
      vasopressors           ~ "Any vasopressor use"
    )

    # ---- gtsummary 主表 ----
    tbl_main <- df_baseline %>%
      select(all_of(baseline_vars)) %>%
      tbl_summary(
        label     = var_labels,
        statistic = list(
          all_continuous()  ~ "{median} ({p25}, {p75})",
          all_categorical() ~ "{n} ({p}%)"
        ),
        digits  = list(all_continuous() ~ 1, all_categorical() ~ c(0, 1)),
        missing = "no"
      ) %>%
      modify_header(
        label  ~ "**Variable**",
        stat_0 ~ "**Overall Cohort**\nN = {N}"
      ) %>%
      modify_footnote(
        stat_0 ~ "Continuous variables: median (IQR); Categorical variables: n (%)"
      )

    # ---- 计算 Missing (%) —— 强制去除names残留 ----
    miss_raw <- df_baseline %>%
      select(all_of(baseline_vars)) %>%
      summarise(across(everything(), ~ round(100 * mean(is.na(.)), 1))) %>%
      pivot_longer(everything(), names_to = "var", values_to = "pct") %>%
      mutate(miss_str = ifelse(pct == 0, "", sprintf("%.1f", pct)))

    # 用明确的 match() 取值，完全不依赖 named vector 的 [ ] 操作
    get_miss <- function(varname) {
      idx <- match(varname, miss_raw$var)
      if (is.na(idx)) return("—")
      miss_raw$miss_str[idx]
    }

    # ---- 共用：父变量缺失率映射 ----
    parent_missing <- setNames(
      sapply(c("y","anchor_age","gender","race_clean","insurance",
               "first_careunit","sepsis","primary_stroke_cerebro",
               "gcs_total_last24h","icu_hr_last24h",
               "Creatinine","Albumin","vent_any","vasopressors"),
             get_miss),
      c("ICU LOS, days","Age, years","Sex","Race/ethnicity","Insurance",
        "First care unit","Sepsis","Stroke/cerebrovascular disease",
        "GCS total","Heart rate","Creatinine","Albumin",
        "Mechanical ventilation","Any vasopressor use")
    )

    # 确认一下值是否正确
    print(parent_missing)

    # ---- 共用：父变量 -> 分组映射 ----
    parent_group <- c(
      "ICU LOS, days"                  = "Outcome",
      "Age, years"                     = "Demographics",
      "Sex"                            = "Demographics",
      "Race/ethnicity"                 = "Demographics",
      "Insurance"                      = "Demographics",
      "First care unit"                = "ICU context",
      "Sepsis"                         = "Clinical conditions",
      "Stroke/cerebrovascular disease" = "Clinical conditions",
      "GCS total"                      = "Clinical parameters",
      "Heart rate"                     = "Clinical parameters",
      "Creatinine"                     = "Laboratory values",
      "Albumin"                        = "Laboratory values",
      "Mechanical ventilation"         = "Interventions",
      "Any vasopressor use"            = "Interventions"
    )
    # ---- 预先计算每行的 missing 值（与 gt 内部行顺序一致）----
    gt_labels <- tbl_main %>%
      as_tibble() %>%
      setNames(c("Variable", "Overall Cohort")) %>%
      pull(Variable)

    missing_vec <- {
      result      <- rep(NA_character_, length(gt_labels))
      cur_miss    <- NA_character_
      for (i in seq_along(gt_labels)) {
        v <- gt_labels[i]
        if (v %in% names(parent_missing)) {
          cur_miss  <- unname(parent_missing[v])
          result[i] <- cur_miss
        }
        # level 行保持 NA，gt 里会显示为空
      }
      result
    }
    # ---- 预先算好每个组覆盖的行号范围 ----
    gt_row_groups <- local({
      labs <- gt_labels  # 之前已经算好的43个label向量
      parent_vars <- names(parent_group)

      # 找每个父变量在 gt_labels 里的位置
      parent_idx <- which(labs %in% parent_vars)

      # 每个父变量的行范围：从自己到下一个父变量前一行
      ranges <- lapply(seq_along(parent_idx), function(i) {
        start <- parent_idx[i]
        end   <- if (i < length(parent_idx)) parent_idx[i + 1] - 1 else length(labs)
        list(
          label = parent_group[labs[start]],
          rows  = start:end
        )
      })
      ranges
    })
    # ---- 转为 gt ----
    gt_final <- tbl_main %>%
      as_gt() %>%
      tab_header(
        title    = md("**Table 1. Baseline Characteristics of the Study Cohort**"),
        subtitle = md(sprintf("N = %d ICU admissions", nrow(df_baseline)))
      ) %>%
      cols_add(missing_col = missing_vec) %>%
      cols_label(missing_col = md("**Missing (%)**"))

    # 动态添加 row groups（从后往前声明，gt 会正序显示）
    # 先把相同 label 的组合并（Demographics 里有多个父变量要合成一个组）
    unique_groups <- unique(sapply(gt_row_groups, `[[`, "label"))

    for (grp_label in rev(unique_groups)) {
      # 找所有属于这个 label 的行号，合并成一个向量
      all_rows <- unlist(lapply(gt_row_groups, function(rg) {
        if (rg$label == grp_label) rg$rows else NULL
      }))
      all_rows <- sort(unique(all_rows))

      gt_final <- gt_final %>%
        tab_row_group(
          id    = grp_label,                            # 唯一ID
          label = md(paste0("**", grp_label, "**")),
          rows  = all_rows
        )
    }

    gt_final <- gt_final %>%
      tab_footnote(
        footnote  = "Missing values are excluded from percentage calculations for categorical variables.",
        locations = cells_column_labels(columns = missing_col)
      ) %>%
      tab_options(
        table.border.top.color            = "black",
        table.border.top.width            = px(2),
        table.border.bottom.color         = "black",
        table.border.bottom.width         = px(2),
        column_labels.border.top.color    = "black",
        column_labels.border.top.width    = px(2),
        column_labels.border.bottom.color = "black",
        column_labels.border.bottom.width = px(1),
        table_body.hlines.color           = "transparent",
        row_group.border.top.color        = "transparent",
        row_group.border.bottom.color     = "transparent",
        table.font.size                   = px(13),
        heading.title.font.size           = px(14),
        heading.subtitle.font.size        = px(12)
      ) %>%
      tab_style(
        style     = list(cell_text(weight = "bold"),
                         cell_fill(color = "#f2f2f2")),
        locations = cells_row_groups()
      ) %>%
      tab_style(
        style     = cell_text(align = "center"),
        locations = cells_column_labels(columns = c(stat_0, missing_col))
      ) %>%
      tab_style(
        style     = cell_text(align = "center"),
        locations = cells_body(columns = c(stat_0, missing_col))
      )

    baseline_html <- file.path(base_fig_dir, "Table1_baseline.html")
    gtsave(gt_final, baseline_html)
    cat(sprintf("✓ HTML saved: %s\n", baseline_html))

    # ---- 导出 CSV（含分组、level归属）----
    raw_tbl <- tbl_main %>%
      as_tibble() %>%
      setNames(c("Variable", "Overall Cohort"))

    parent_vars <- names(parent_group)

    # 逐行向下传播：遇到父变量记录当前组，非父变量行继承上一个父变量的组
    csv_df <- raw_tbl %>%
      mutate(
        is_parent  = Variable %in% parent_vars,
        Group      = NA_character_,
        Missing    = NA_character_
      )

    current_group   <- NA_character_
    current_missing <- NA_character_

    for (i in seq_len(nrow(csv_df))) {
      v <- csv_df$Variable[i]
      if (v %in% parent_vars) {
        current_group   <- parent_group[v]
        current_missing <- unname(parent_missing[v])
        csv_df$Group[i]   <- current_group
        csv_df$Missing[i] <- current_missing   # 父行显示缺失率
      } else {
        csv_df$Group[i]   <- current_group
        csv_df$Missing[i] <- ""               # level 行不显示缺失率
      }
    }

    csv_df <- csv_df %>%
      select(Group, Variable, `Overall Cohort`, `Missing (%)` = Missing)

    baseline_csv <- file.path(base_fig_dir, "Table1_baseline.csv")
    write.csv(csv_df, baseline_csv, row.names = FALSE)
    cat(sprintf("✓ CSV saved: %s\n", baseline_csv))
  }, envir = env)
}
