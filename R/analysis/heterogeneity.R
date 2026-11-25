# ============================================================
# Subgroup (heterogeneity) CRPS analysis.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Subgroup CRPS (vasopressor, ventilation, sepsis, stroke, unit, age) - 19 levels.
# Source: scripts.R:2050-2152
stage_eval19_heterogeneity <- function(env = .GlobalEnv) {
  local({
    # HETEROGENEITY ANALYSIS (subgroup CRPS)
    # Subgroups selected based on clinical prior knowledge:
    # interventions (vasopressin, vasopressors, ventilation),
    # diagnoses (sepsis, stroke), care unit, and age group.
    # LOS severity removed (outcome-based subgrouping = circular).
    # Missingness subgroups removed (data quality, not clinical heterogeneity).
    # ============================================================
    cat("\n========================================\n")
    cat("HETEROGENEITY ANALYSIS (CRPS-based subgroups)\n")
    cat("========================================\n\n")

    # ---- Subgroup 1: Vasopressin use ----
    if ("vasopressin" %in% names(Xte_tree)) {
      eval_df$vasopressin_group <- factor(
        Xte_tree$vasopressin, 0:1, c("No Vasopressin", "Vasopressin")
      )
      out_vasopressin_crps <- subgroup_crps_table(eval_df, "vasopressin_group")
      cat("\n=== CRPS by Vasopressin Subgroup ===\n"); print(out_vasopressin_crps)
    } else {
      out_vasopressin_crps <- NULL
      cat("\n(vasopressin not found in test features)\n")
    }

    # ---- Subgroup 2: Any Vasopressor use ----
    if ("vasopressors" %in% names(Xte_tree)) {
      eval_df$vaso_group <- factor(Xte_tree$vasopressors, 0:1, c("No Vasopressors", "Vasopressors"))
      out_vaso_crps <- subgroup_crps_table(eval_df, "vaso_group")
      cat("\n=== CRPS by Any Vasopressor Subgroup ===\n"); print(out_vaso_crps)
    } else {
      out_vaso_crps <- NULL
      cat("\n(vasopressors not found in test features)\n")
    }

    # ---- Subgroup 3: Mechanical Ventilation ----
    if ("vent_any" %in% names(Xte_tree)) {
      eval_df$vent_group <- factor(Xte_tree$vent_any, 0:1, c("No Ventilation", "Ventilated"))
      out_vent_crps <- subgroup_crps_table(eval_df, "vent_group")
      cat("\n=== CRPS by Mechanical Ventilation Subgroup ===\n"); print(out_vent_crps)
    } else {
      out_vent_crps <- NULL
      cat("\n(vent_any not found in test features)\n")
    }

    # ---- Subgroup 4: Sepsis ----
    if ("sepsis" %in% names(Xte_tree)) {
      eval_df$sepsis_group <- factor(Xte_tree$sepsis, 0:1, c("No Sepsis", "Sepsis"))
      out_sepsis_crps <- subgroup_crps_table(eval_df, "sepsis_group")
      cat("\n=== CRPS by Sepsis Subgroup ===\n"); print(out_sepsis_crps)
    } else {
      out_sepsis_crps <- NULL
      cat("\n(sepsis not found in test features)\n")
    }

    # ---- Subgroup 5: Primary Dx Stroke/Cerebrovascular ----
    if ("primary_stroke_cerebro" %in% names(Xte_tree)) {
      eval_df$stroke_group <- factor(
        Xte_tree$primary_stroke_cerebro, 0:1, c("No Stroke/Cerebro", "Stroke/Cerebro")
      )
      out_stroke_crps <- subgroup_crps_table(eval_df, "stroke_group")
      cat("\n=== CRPS by Stroke/Cerebrovascular Subgroup ===\n"); print(out_stroke_crps)
    } else {
      out_stroke_crps <- NULL
      cat("\n(primary_stroke_cerebro not found in test features)\n")
    }

    # ---- Subgroup 6: ICU Care Unit ----
    if ("first_careunit" %in% names(Xte_tree)) {
      unit_vals_te <- trimws(as.character(Xte_tree$first_careunit))
      unit_counts  <- sort(table(unit_vals_te), decreasing = TRUE)
      units_keep   <- names(unit_counts[unit_counts >= 200])
      cat(sprintf("\nCare units with n >= 200 in test set (%d units):\n", length(units_keep)))
      print(unit_counts[units_keep])
      eval_df$unit_group_raw <- ifelse(
        unit_vals_te %in% units_keep, unit_vals_te, NA_character_
      )
      eval_df$unit_group_raw <- factor(eval_df$unit_group_raw, levels = units_keep)
      out_unit_crps <- subgroup_crps_table(
        eval_df %>% filter(!is.na(unit_group_raw)), "unit_group_raw"
      )
      cat("\n=== CRPS by ICU Care Unit (n>=200 only) ===\n"); print(out_unit_crps)
    } else {
      out_unit_crps <- NULL
      cat("\n(first_careunit not found in test features)\n")
    }

    # ---- Subgroup 7: Age Group ----
    # Clinically motivated cut-points: <65 (non-elderly), 65-80 (elderly), >80 (oldest-old)
    # Older patients have greater physiological complexity and more variable LOS distributions
    if ("anchor_age" %in% names(Xte_tree)) {
      eval_df$age_group <- cut(
        Xte_tree$anchor_age,
        breaks = c(-Inf, 65, 80, Inf),
        labels = c("<65", "65-80", ">80"),
        right  = TRUE
      )
      age_counts <- table(eval_df$age_group)
      cat(sprintf("\nAge group distribution in test set:\n")); print(age_counts)
      out_age_crps <- subgroup_crps_table(eval_df, "age_group")
      cat("\n=== CRPS by Age Group Subgroup ===\n"); print(out_age_crps)
    } else {
      out_age_crps <- NULL
      cat("\n(anchor_age not found in test features)\n")
    }
  }, envir = env)
}

# ---- Subgroup CRPS tables (39 levels) + CSV export.
# Source: 39quantiles.R:858-929
stage_eval39_heterogeneity <- function(env = .GlobalEnv) {
  local({
    cat("\n========================================\n")
    cat("HETEROGENEITY ANALYSIS (CRPS-based subgroups)\n")
    cat("========================================\n\n")

    out_vasopressin_crps <- out_vaso_crps <- out_vent_crps <- NULL
    out_sepsis_crps <- out_stroke_crps <- out_unit_crps <- out_age_crps <- NULL

    if ("vasopressin" %in% names(Xte_tree)) {
      eval_df$vasopressin_group <- factor(Xte_tree$vasopressin, 0:1,
                                          c("No Vasopressin", "Vasopressin"))
      out_vasopressin_crps <- subgroup_crps_table(eval_df, "vasopressin_group")
      cat("\n=== CRPS by Vasopressin Subgroup ===\n"); print(out_vasopressin_crps)
    } else { cat("\n(vasopressin not found in test features)\n") }

    if ("vasopressors" %in% names(Xte_tree)) {
      eval_df$vaso_group <- factor(Xte_tree$vasopressors, 0:1,
                                   c("No Vasopressors", "Vasopressors"))
      out_vaso_crps <- subgroup_crps_table(eval_df, "vaso_group")
      cat("\n=== CRPS by Any Vasopressor Subgroup ===\n"); print(out_vaso_crps)
    } else { cat("\n(vasopressors not found in test features)\n") }

    if ("vent_any" %in% names(Xte_tree)) {
      eval_df$vent_group <- factor(Xte_tree$vent_any, 0:1,
                                   c("No Ventilation", "Ventilated"))
      out_vent_crps <- subgroup_crps_table(eval_df, "vent_group")
      cat("\n=== CRPS by Mechanical Ventilation Subgroup ===\n"); print(out_vent_crps)
    } else { cat("\n(vent_any not found in test features)\n") }

    if ("sepsis" %in% names(Xte_tree)) {
      eval_df$sepsis_group <- factor(Xte_tree$sepsis, 0:1, c("No Sepsis", "Sepsis"))
      out_sepsis_crps <- subgroup_crps_table(eval_df, "sepsis_group")
      cat("\n=== CRPS by Sepsis Subgroup ===\n"); print(out_sepsis_crps)
    } else { cat("\n(sepsis not found in test features)\n") }

    if ("primary_stroke_cerebro" %in% names(Xte_tree)) {
      eval_df$stroke_group <- factor(Xte_tree$primary_stroke_cerebro, 0:1,
                                     c("No Stroke/Cerebro", "Stroke/Cerebro"))
      out_stroke_crps <- subgroup_crps_table(eval_df, "stroke_group")
      cat("\n=== CRPS by Stroke/Cerebrovascular Subgroup ===\n"); print(out_stroke_crps)
    } else { cat("\n(primary_stroke_cerebro not found in test features)\n") }

    if ("first_careunit" %in% names(Xte_tree)) {
      unit_vals_te <- trimws(as.character(Xte_tree$first_careunit))
      unit_counts  <- sort(table(unit_vals_te), decreasing = TRUE)
      units_keep   <- names(unit_counts[unit_counts >= 200])
      cat(sprintf("\nCare units with n >= 200 in test set (%d units):\n", length(units_keep)))
      print(unit_counts[units_keep])
      eval_df$unit_group_raw <- factor(
        ifelse(unit_vals_te %in% units_keep, unit_vals_te, NA_character_),
        levels = units_keep)
      out_unit_crps <- subgroup_crps_table(
        eval_df %>% filter(!is.na(unit_group_raw)), "unit_group_raw")
      cat("\n=== CRPS by ICU Care Unit (n>=200) ===\n"); print(out_unit_crps)
    } else { cat("\n(first_careunit not found in test features)\n") }

    if ("anchor_age" %in% names(Xte_tree)) {
      eval_df$age_group <- cut(Xte_tree$anchor_age,
                               breaks = c(-Inf, 65, 80, Inf),
                               labels = c("<65", "65-80", ">80"), right = TRUE)
      cat("\nAge group distribution:\n"); print(table(eval_df$age_group))
      out_age_crps <- subgroup_crps_table(eval_df, "age_group")
      cat("\n=== CRPS by Age Group ===\n"); print(out_age_crps)
    } else { cat("\n(anchor_age not found in test features)\n") }

    if (!is.null(out_vasopressin_crps)) write.csv(out_vasopressin_crps, file.path(out_dir, "heterogeneity_vasopressin_crps.csv"),  row.names = FALSE)
    if (!is.null(out_vaso_crps))        write.csv(out_vaso_crps,        file.path(out_dir, "heterogeneity_vaso_crps.csv"),         row.names = FALSE)
    if (!is.null(out_vent_crps))        write.csv(out_vent_crps,        file.path(out_dir, "heterogeneity_vent_crps.csv"),          row.names = FALSE)
    if (!is.null(out_sepsis_crps))      write.csv(out_sepsis_crps,      file.path(out_dir, "heterogeneity_sepsis_crps.csv"),        row.names = FALSE)
    if (!is.null(out_stroke_crps))      write.csv(out_stroke_crps,      file.path(out_dir, "heterogeneity_stroke_crps.csv"),        row.names = FALSE)
    if (!is.null(out_unit_crps))        write.csv(out_unit_crps,        file.path(out_dir, "heterogeneity_unit_crps.csv"),          row.names = FALSE)
    if (!is.null(out_age_crps))         write.csv(out_age_crps,         file.path(out_dir, "heterogeneity_age_crps.csv"),           row.names = FALSE)
    cat("✓ Heterogeneity CSV tables saved\n")
  }, envir = env)
}
