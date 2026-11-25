# ============================================================
# Data preparation shared by all pipelines.
# Stage functions run their (verbatim) original code inside `env`
# (default: global environment) so that data flows between stages
# exactly as in the original scripts.
# ============================================================

# ---- Read the master CSV, clean race/categoricals, add *_miss indicators, build df2 / df_imp / X_tree / X_all (full cohort, pre-split).
# Source: scripts.R:11-150
stage_prepare_features <- function(env = .GlobalEnv) {
  local({
    # read the dataset
    df <- read.csv(cfg$data_path) %>%
      filter(!is.na(icu_los_days))

    # ---- fix impossible RR=0 (treat as missing) ----
    if ("icu_rr_last24h_final" %in% names(df)) {
      df$icu_rr_last24h_final[df$icu_rr_last24h_final == 0] <- NA
    } else if ("icu_rr_last24h" %in% names(df)) {
      df$icu_rr_last24h[df$icu_rr_last24h == 0] <- NA
    }

    # variable mapping
    race_mapping <- c(
      "WHITE" = "White",
      "WHITE - OTHER EUROPEAN" = "White",
      "WHITE - EASTERN EUROPEAN" = "White",
      "WHITE - BRAZILIAN" = "White",
      "WHITE - RUSSIAN" = "White",
      "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER" = "White",
      "PORTUGUESE" = "White",

      "BLACK/AFRICAN AMERICAN" = "Black",
      "BLACK/CARIBBEAN ISLAND" = "Black",
      "BLACK/AFRICAN" = "Black",
      "BLACK/CAPE VERDEAN" = "Black",

      "UNKNOWN" = "Unknown",
      "PATIENT DECLINED TO ANSWER" = "Unknown",
      "UNABLE TO OBTAIN" = "Unknown",

      "OTHER" = "Other",
      "MULTIPLE RACE/ETHNICITY" = "Other",
      "AMERICAN INDIAN/ALASKA NATIVE" = "Other",

      "ASIAN" = "Asian",
      "ASIAN - SOUTH EAST ASIAN" = "Asian",
      "ASIAN - CHINESE" = "Asian",
      "ASIAN - ASIAN INDIAN" = "Asian",
      "ASIAN - KOREAN" = "Asian",

      "HISPANIC OR LATINO" = "Hispanic/Latino",
      "HISPANIC/LATINO - CUBAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - DOMINICAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - PUERTO RICAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - SALVADORAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - COLUMBIAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - GUATEMALAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - MEXICAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - CENTRAL AMERICAN" = "Hispanic/Latino",
      "HISPANIC/LATINO - HONDURAN" = "Hispanic/Latino",
      "SOUTH AMERICAN" = "Hispanic/Latino"
    )

    # create race clean variable and outcome variable
    df2 <- df %>%
      mutate(
        race = trimws(as.character(race)),
        race_clean = dplyr::recode(race, !!!race_mapping, .default = "Other"),
        race_clean = factor(race_clean,
                            levels = c("White", "Black", "Asian",
                                       "Hispanic/Latino", "Other", "Unknown")),
        y = icu_los_days
      )

    # build a dropping vector
    drop_cols <- c(
      "stay_id",
      "race",
      "admittime","dischtime","icu_intime_first","icu_outtime_first","intime",
      "primary_icd_long_title","primary_icd_code",
      "readmit_30d","readmit_90d","readmit_180d",
      "icu_los_days"
    )

    # drop the column based on created vector
    df2 <- df2 %>% select(-any_of(drop_cols))
    # ---- drop no-information / unwanted cols ----
    extra_drop <- c("arni", "mask_vent_or_intub", "icd_flag")
    df2 <- df2 %>% select(-any_of(extra_drop))

    names(df2) <- sub("_final$", "", names(df2))


    # ---- fix pseudo-missing: "" -> NA for categorical columns ----
    cat_cols_raw <- intersect(c("race_clean","marital_status","insurance","gender","first_careunit"), names(df2))

    df2 <- df2 %>%
      mutate(across(all_of(cat_cols_raw), ~ {
        x <- trimws(as.character(.x))
        x[x == ""] <- NA_character_
        x
      }))

    # transfer categories variables as factor

    cat_cols <- intersect(c("race_clean","marital_status","insurance","gender","first_careunit"), names(df2))
    df2 <- df2 %>% mutate(across(all_of(cat_cols), ~as.factor(trimws(.x))))

    # select all numeric variables except for los
    num_cols <- names(df2)[sapply(df2, is.numeric)]
    num_cols <- setdiff(num_cols, "y")

    # Create a new column named original_variable + "_miss"
    # Assign 1 if missing, 0 otherwise
    for (v in num_cols) {
      df2[[paste0(v, "_miss")]] <- as.integer(is.na(df2[[v]]))
    }

    # prepare a imputed dataset for QRF
    # extract and parse admission time from the original dataframe
    df_time <- df %>%
      mutate(admit_dt = dmy_hm(admittime)) %>%
      select(subject_id, hadm_id, admit_dt)

    # Merge admission datetime into df2
    df2 <- df2 %>%
      left_join(df_time, by = c("subject_id","hadm_id"))

    # Remove ID variables and their missingness indicators
    df2 <- df2 %>% select(-subject_id, -hadm_id, -subject_id_miss, -hadm_id_first_miss, -hadm_id_miss, -hadm_id_first)
    # create a new dataframe with median imputation applied on its numeric variables
    df_imp <- df2
    for (v in num_cols) {
      med <- median(df_imp[[v]], na.rm = TRUE)
      df_imp[[v]][is.na(df_imp[[v]])] <- med
    }
    # for categorical variables with missing values, assign level "unknown" to them
    for (v in cat_cols) {
      df_imp[[v]] <- addNA(df_imp[[v]])
      levels(df_imp[[v]])[is.na(levels(df_imp[[v]]))] <- "Unknown"
    }


    # finalize the dataset
    # using los as outcome variable
    y <- df2$y
    # X_tree dataset for DRF/tree model
    X_tree <- df2 %>% select(-y)
    # X_all dataset for QRF model
    X_all  <- df_imp %>% select(-y)
  }, envir = env)
}

# ---- 80:20 temporal split on admit_dt; train-led removal of all-zero *_miss and constant columns -> Xtr/Xte_tree, Xtr/Xte_all, y_tr, y_te.
# Source: scripts.R:460-527
stage_split_data <- function(env = .GlobalEnv) {
  local({
    # 1) data split (by admission time)
    cutoff <- quantile(df2$admit_dt, 0.8, na.rm = TRUE)
    idx_train <- which(df2$admit_dt <= cutoff)
    idx_test  <- which(df2$admit_dt >  cutoff)

    # 2) remove admit_dt from modeling tables (but AFTER using it for split)
    df2_model <- df2 %>% select(-admit_dt)
    df_imp2   <- df_imp %>% select(-any_of("admit_dt"))

    # 3) define outcome + feature tables
    y      <- df2_model$y
    X_tree <- df2_model %>% select(-y)
    X_all  <- df_imp2 %>% select(-y)

    # 4) split train/test
    Xtr_tree <- X_tree[idx_train, , drop = FALSE]
    Xte_tree <- X_tree[idx_test,  , drop = FALSE]
    y_tr     <- y[idx_train]
    y_te     <- y[idx_test]

    Xtr_all  <- X_all[idx_train, , drop = FALSE]
    Xte_all  <- X_all[idx_test,  , drop = FALSE]

    # ---------------------------------------------------------
    # 5) IMPORTANT: clean tree features AFTER split (train-led)
    # ---------------------------------------------------------
    miss_cols <- grep("_miss$", names(Xtr_tree), value = TRUE)

    miss_all_zero <- miss_cols[
      sapply(miss_cols, function(v) sum(Xtr_tree[[v]], na.rm = TRUE) == 0)
    ]

    const_cols <- names(Xtr_tree)[
      sapply(Xtr_tree, function(x) dplyr::n_distinct(x, na.rm = TRUE) <= 1)
    ]

    drop_cols_tree <- union(miss_all_zero, const_cols)

    Xtr_tree <- Xtr_tree %>% select(-any_of(drop_cols_tree))
    Xte_tree <- Xte_tree %>% select(-any_of(drop_cols_tree))

    cat("Tree feature cleanup:\n")
    cat("  Removed useless *_miss (all-zero in train):", length(miss_all_zero), "\n")
    cat("  Removed constant columns (in train):", length(const_cols), "\n")
    cat("  Xtr_tree cols:", ncol(Xtr_tree), " | Xte_tree cols:", ncol(Xte_tree), "\n")

    # ---------------------------------------------------------
    # 5b) clean all (QRF) features — train-led
    # ---------------------------------------------------------
    miss_cols_all <- grep("_miss$", names(Xtr_all), value = TRUE)

    miss_all_zero_all <- miss_cols_all[
      sapply(miss_cols_all, function(v) sum(Xtr_all[[v]], na.rm = TRUE) == 0)
    ]

    const_cols_all <- names(Xtr_all)[
      sapply(Xtr_all, function(x) dplyr::n_distinct(x, na.rm = TRUE) <= 1)
    ]

    drop_cols_all <- union(miss_all_zero_all, const_cols_all)

    Xtr_all <- Xtr_all %>% select(-any_of(drop_cols_all))
    Xte_all <- Xte_all %>% select(-any_of(drop_cols_all))

    cat("All (QRF) feature cleanup:\n")
    cat("  Removed useless *_miss (all-zero in train):", length(miss_all_zero_all), "\n")
    cat("  Removed constant columns (in train):", length(const_cols_all), "\n")
    cat("  Xtr_all cols:", ncol(Xtr_all), " | Xte_all cols:", ncol(Xte_all), "\n")
  }, envir = env)
}
