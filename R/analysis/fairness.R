# ============================================================
# Fairness (protected attribute) CRPS analysis.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Subgroup CRPS by gender/race/insurance - 19 levels.
# Source: scripts.R:2154-2224
stage_eval19_fairness <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # FAIRNESS ANALYSIS (subgroup CRPS by protected attributes)
    # Examines whether UQ performance differs systematically
    # across sociodemographic groups (gender, race, insurance).
    # Uses the same CRPS subgroup framework as heterogeneity
    # analysis but interpreted through a health equity lens.
    # ============================================================
    cat("\n========================================\n")
    cat("FAIRNESS ANALYSIS (CRPS by protected attributes)\n")
    cat("========================================\n\n")

    # ---- Fairness 1: Gender ----
    if ("gender" %in% names(Xte_tree)) {
      gender_vals <- trimws(as.character(Xte_tree$gender))
      gender_vals[gender_vals == "" | gender_vals == "Unknown"] <- NA_character_
      valid_genders <- names(table(gender_vals)[table(gender_vals) >= 50])
      eval_df$gender_group <- factor(
        ifelse(gender_vals %in% valid_genders, gender_vals, NA_character_),
        levels = valid_genders
      )
      cat(sprintf("\nGender distribution in test set:\n"))
      print(table(eval_df$gender_group, useNA = "ifany"))
      out_gender_crps <- subgroup_crps_table(
        eval_df %>% filter(!is.na(gender_group)), "gender_group"
      )
      cat("\n=== CRPS by Gender ===\n"); print(out_gender_crps)
    } else {
      out_gender_crps <- NULL
      cat("\n(gender not found in test features)\n")
    }

    # ---- Fairness 2: Race / Ethnicity ----
    if ("race_clean" %in% names(Xte_tree)) {
      race_vals <- trimws(as.character(Xte_tree$race_clean))
      race_vals[race_vals == "" | race_vals == "Unknown"] <- NA_character_
      valid_races <- names(table(race_vals)[table(race_vals) >= 50])
      eval_df$race_group <- factor(
        ifelse(race_vals %in% valid_races, race_vals, NA_character_),
        levels = valid_races
      )
      cat(sprintf("\nRace/Ethnicity distribution in test set:\n"))
      print(table(eval_df$race_group, useNA = "ifany"))
      out_race_crps <- subgroup_crps_table(
        eval_df %>% filter(!is.na(race_group)), "race_group"
      )
      cat("\n=== CRPS by Race/Ethnicity ===\n"); print(out_race_crps)
    } else {
      out_race_crps <- NULL
      cat("\n(race_clean not found in test features)\n")
    }

    # ---- Fairness 3: Insurance Type ----
    # Insurance type is a proxy for socioeconomic status in the US healthcare context
    if ("insurance" %in% names(Xte_tree)) {
      ins_vals <- trimws(as.character(Xte_tree$insurance))
      ins_vals[ins_vals == "" | ins_vals == "Unknown"] <- NA_character_
      valid_ins <- names(table(ins_vals)[table(ins_vals) >= 50])
      eval_df$insurance_group <- factor(
        ifelse(ins_vals %in% valid_ins, ins_vals, NA_character_),
        levels = valid_ins
      )
      cat(sprintf("\nInsurance type distribution in test set:\n"))
      print(table(eval_df$insurance_group, useNA = "ifany"))
      out_insurance_crps <- subgroup_crps_table(
        eval_df %>% filter(!is.na(insurance_group)), "insurance_group"
      )
      cat("\n=== CRPS by Insurance Type ===\n"); print(out_insurance_crps)
    } else {
      out_insurance_crps <- NULL
      cat("\n(insurance not found in test features)\n")
    }
  }, envir = env)
}

# ---- Fairness subgroup CRPS tables (39 levels).
# Source: 39quantiles.R:1368-1417
stage_eval39_fairness <- function(env = .GlobalEnv) {
  local({
    # ============================================================
    # STEP 12: FAIRNESS ANALYSIS
    # ============================================================
    cat("\n========================================\n")
    cat("FAIRNESS ANALYSIS (CRPS by protected attributes)\n")
    cat("========================================\n\n")

    out_gender_crps <- out_race_crps <- out_insurance_crps <- NULL

    if ("gender" %in% names(Xte_tree)) {
      gender_vals   <- trimws(as.character(Xte_tree$gender))
      gender_vals[gender_vals == "" | gender_vals == "Unknown"] <- NA_character_
      valid_genders <- names(table(gender_vals)[table(gender_vals) >= 50])
      eval_df$gender_group <- factor(
        ifelse(gender_vals %in% valid_genders, gender_vals, NA_character_),
        levels = valid_genders)
      cat("\nGender distribution:\n"); print(table(eval_df$gender_group, useNA = "ifany"))
      out_gender_crps <- subgroup_crps_table(eval_df %>% filter(!is.na(gender_group)), "gender_group")
      cat("\n=== CRPS by Gender ===\n"); print(out_gender_crps)
    } else { cat("\n(gender not found in test features)\n") }

    if ("race_clean" %in% names(Xte_tree)) {
      race_vals   <- trimws(as.character(Xte_tree$race_clean))
      race_vals[race_vals == "" | race_vals == "Unknown"] <- NA_character_
      valid_races <- names(table(race_vals)[table(race_vals) >= 50])
      eval_df$race_group <- factor(
        ifelse(race_vals %in% valid_races, race_vals, NA_character_),
        levels = valid_races)
      cat("\nRace/Ethnicity distribution:\n"); print(table(eval_df$race_group, useNA = "ifany"))
      out_race_crps <- subgroup_crps_table(eval_df %>% filter(!is.na(race_group)), "race_group")
      cat("\n=== CRPS by Race/Ethnicity ===\n"); print(out_race_crps)
    } else { cat("\n(race_clean not found in test features)\n") }

    if ("insurance" %in% names(Xte_tree)) {
      ins_vals   <- trimws(as.character(Xte_tree$insurance))
      ins_vals[ins_vals == "" | ins_vals == "Unknown"] <- NA_character_
      valid_ins  <- names(table(ins_vals)[table(ins_vals) >= 50])
      eval_df$insurance_group <- factor(
        ifelse(ins_vals %in% valid_ins, ins_vals, NA_character_),
        levels = valid_ins)
      cat("\nInsurance type distribution:\n"); print(table(eval_df$insurance_group, useNA = "ifany"))
      out_insurance_crps <- subgroup_crps_table(eval_df %>% filter(!is.na(insurance_group)), "insurance_group")
      cat("\n=== CRPS by Insurance Type ===\n"); print(out_insurance_crps)
    } else { cat("\n(insurance not found in test features)\n") }

    if (!is.null(out_gender_crps))    write.csv(out_gender_crps,    file.path(out_dir, "fairness_gender_crps.csv"),    row.names = FALSE)
    if (!is.null(out_race_crps))      write.csv(out_race_crps,      file.path(out_dir, "fairness_race_crps.csv"),      row.names = FALSE)
    if (!is.null(out_insurance_crps)) write.csv(out_insurance_crps, file.path(out_dir, "fairness_insurance_crps.csv"), row.names = FALSE)
    cat("✓ Fairness CSV tables saved\n")
  }, envir = env)
}
