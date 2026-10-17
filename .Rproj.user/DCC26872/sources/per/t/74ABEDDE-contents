############################################################
## main_call.R: data loading and analysis execution
## Run from the Revised_code directory (or adjust the data paths below).
############################################################
#' Run the Long COVID DML Analysis
#'
#' Loads the N3C-derived analysis data, applies the selected missing-data and
#' mortality scenario, constructs cohort-specific exposure variables, selects
#' adjustment covariates, and runs cross-fitted DML for each comorbidity.
#'
#' @return No return value. Results and diagnostics are written below `results/`.
#' @keywords internal
NULL
## Resolve paths from this file so the script works when launched either from
## Revised_code or directly as R/main_call.R.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else if (!is.null(source_file) && nzchar(source_file)) {
  dirname(normalizePath(source_file))
} else {
  wd <- normalizePath(getwd())
  if (file.exists(file.path(wd, "dml.R"))) {
    wd
  } else if (file.exists(file.path(wd, "R", "dml.R"))) {
    file.path(wd, "R")
  } else {
    stop("Could not locate dml.R. Run from Revised_code or Revised_code/R.")
  }
}
project_root <- if (basename(script_dir) == "R") {
  normalizePath(file.path(script_dir, ".."))
} else {
  script_dir
}
setwd(project_root)
source(file.path(script_dir, "dml.R"))
#' Load the Analysis Data
#'
#' Reads the primary parquet data set and creates the `Dead` indicator used by
#' the sensitivity analyses.
############################################################
# Select the input with DML_DATA_PATH; the small simulated data set is the
# default for testing. Set DML_DATA_PATH explicitly for production analyses.
data_path <- Sys.getenv(
  "DML_DATA_PATH",
  unset = "./data/Simulated_data_dml.parquet"
)
# Data availability: All data used in this study is available through the N3C Enclave to approved users.
# See https://ncats.nih.gov/research/research-activities/n3c/overview for instructions on how to access the data.
#------ data loading ----------------------
############################################################
# source("./Main Code.R")
# write_parquet(Long_COVID_df[1:5E2,], "Simulated_data_dml.parquet")


#' Configure Sensitivity Scenarios
#'
#' Controls complete-case filtering and imputation/decedent inclusion.
#  ------ Initialize ----
Sensitivity.non.imputed.exc.deceased <- identical(
  Sys.getenv("DML_NON_IMPUTED_EXC_DECEASED", "TRUE"), "TRUE"
)
Sensitivity.imputed.inc.deceased <- identical(
  Sys.getenv("DML_IMPUTED_INC_DECEASED", "FALSE"), "TRUE"
)



# Long_COVID_df <- read_parquet("./data/Long_COVID_Secondary_Data_all_2026_05.parquet")
# Long_COVID_df <- read_parquet("./data/Long_COVID.parquet")
Long_COVID_df <- read_parquet(data_path)
Long_COVID_df[, Dead := 0L]
# Compute condition (vectorized, no copies)
Long_COVID_df[
  !is.na(death_date) &
    is.na(Long_COVID) &
    !is.na(COVID_19) &
    (death_date - COVID_19 <= 180),
  Dead := 1L
]

#' Define Candidate Adjustment Variables
#'
#' Collects demographic, vaccination, socioeconomic, insurance, and county-level
#' candidate covariates before correlation-based reduction.
Selected.variable <- c(
  "age",
  "sex",
  "race",
  "BMI",
  "total_number_of_COVID_vaccine_doses",
  "poverty_status",
  "health_insurance",
  colnames(Long_COVID_df)[c(79:102, 105:113)]
)
selected.ind <- which(colnames(Long_COVID_df) %in% Selected.variable)

#' Optional Multiple-Imputation Workflow
#'
#' Documents the `mice` workflow used to create the averaged imputed data object.
##Missing data
##--------------------------------------------------
# setDF(Long_COVID_df)
# temp <- Long_COVID_df[, selected.ind]
#
# # Perform multiple imputation by chained equations
# imp <- mice(temp, m = 5, method = "pmm", seed = 0)
#
# long_data <- complete(imp, action = "long", include = TRUE)
#
#
# avg_data <- long_data %>%
#   filter(.imp != 0) %>%
#   group_by(.id) %>%
#   dplyr::summarise(across(where(is.numeric), mean), .groups = "drop")




#' Apply the Selected Missing-Data Scenario
#'
#' Applies complete-case filtering or replaces covariates from the imputed data.
# ---- Load datasets ----
if(Sensitivity.non.imputed.exc.deceased){
  Long_COVID_df <- Long_COVID_df[
    complete.cases(Long_COVID_df[, ..selected.ind])
  ]
  Long_COVID_df <- Long_COVID_df %>% filter(Dead == 0)
}
if(Sensitivity.imputed.inc.deceased){
  # load("./data/avg_data_imputed_all.RData")
  imputed.ind <- which(colnames(Long_COVID_df) %in% colnames(avg_data))

  setDF(Long_COVID_df)
  Long_COVID_df[, imputed.ind] <- avg_data[, -c(1, ncol(avg_data))]

  rm(avg_data)
}

if((!Sensitivity.non.imputed.exc.deceased) &(!Sensitivity.imputed.inc.deceased)){
  # load("./data/avg_data_imputed_all.RData")
  imputed.ind <- which(colnames(Long_COVID_df) %in% colnames(avg_data))

  setDF(Long_COVID_df)
  Long_COVID_df[, imputed.ind] <- avg_data[, -c(1, ncol(avg_data))]
  Long_COVID_df <- Long_COVID_df %>% filter(Dead == 0)

  rm(avg_data)
}


#' Harmonize the Grouped Cancer Exposure
#'
#' Uses the earliest date among malignant and benign tumor components.
if (any(grouped_comorbidity_components("Cancer") %in% colnames(Long_COVID_df))) {
  Long_COVID_df$Cancer <- earliest_recorded_date(
    Long_COVID_df,
    grouped_comorbidity_components("Cancer")
  )
}


#' Reduce Correlated Covariates
#'
#' Retains a low-correlation adjustment set using `select_uncorrelated_vars()`.
ind.x <- which(colnames(Long_COVID_df) %in% c("age",
                                              "BMI",
                                              "total_number_of_COVID_vaccine_doses",
                                              "poverty_status",
                                              "health_insurance",
                                              colnames(Long_COVID_df)[c(79:102, 105:113)]
))
setDF(Long_COVID_df)
cor_mat <- cor(Long_COVID_df[, c(ind.x)], use = "complete.obs")

Selected.variable <- c("sex", "race", select_uncorrelated_vars(cor_mat, threshold = 0.7))




#' Define Individual and Grouped Comorbidities
#'
#' Lists the target conditions and grouped exposures used in the analyses.
indicator.Comorbidity <- c(
  "Hypertension",
  "Diabetes_without_complications",
  "Depression",
  "Chronic_lung_disease",
  "Diabetes_with_complications",
  "Kidney_disease",
  "Coronary_artery_disease",
  "Heart_failure",
  "Liver_disease",
  "Cerebrovascular_disease",
  "Peripheral_vascular_disease",
  "Myocardial_infarction",
  "Cardiomyopathy",
  "Dementia",
  "Malignant_tumor",
  "Benign_tumor"
)

selected.Comorbidity <- c(
  indicator.Comorbidity,
  "Diabetes",
  "Cardiovascular",
  "Neuro_psych_disorders",
  "Cancer",
  "Any_comorbidities"
)



############################################################
## Run by comorbidity
##########################################################
#' Configure the Comorbidity-Level DML Runs
#'
#' Sets readiness thresholds, cross-fitting folds, cohort selection, output
#' paths, and sensitivity-specific output directories before the execution loop.
min_treated_n_for_dml <- 10000
min_control_n_for_dml <- 10000
min_treated_events_for_dml <- 100
min_control_events_for_dml <- 100
min_exposure_prevalence_for_dml <- 1e-10
min_abs_ar_flag_for_dml <- 1e-5
setDT(Long_COVID_df)



# if saving detailed CV results is required
save.detailed.results <- TRUE

K <- 4  # Corresponding 5-folds cross validation
#------ 3.2 before ----
before <- 1
# Main procedures
if(before == 1){
  root.path <- "./results/Cohort_IV" # Cohort IV
}else if(before == 2){
  root.path <- "./results/Cohort_II" # Cohort II
}else if(before == 3){
  root.path <- "./results/Cohort_III" # Cohort III
}else if(before == 4){
  root.path <- "./results/Cohort_I" # Cohort I
}

# Sensitivity analysis with non-imputed data excluding decedents
if(Sensitivity.non.imputed.exc.deceased){
  save.detailed.results <- FALSE
  if(before == 1){
    root.path <- "./results/Cohort_IV_non_imputed_exc_deceased" # Cohort IV
  }else if(before == 2){
    root.path <- "./results/Cohort_II_non_imputed_exc_deceased" # Cohort II
  }else if(before == 3){
    root.path <- "./results/Cohort_III_non_imputed_exc_deceased" # Cohort III
  }else if(before == 4){
    root.path <- "./results/Cohort_I_non_imputed_exc_deceased" # Cohort I
  }
}

# Sensitivity analysis with imputed data including deceased individuals
if(Sensitivity.imputed.inc.deceased){
  save.detailed.results <- FALSE
  if(before == 1){
    root.path <- "./results/Cohort_IV_imputed_inc_deceased" # Cohort IV
  }else if(before == 2){
    root.path <- "./results/Cohort_II_imputed_inc_deceased" # Cohort II
  }else if(before == 3){
    root.path <- "./results/Cohort_III_imputed_inc_deceased" # Cohort III
  }else if(before == 4){
    root.path <- "./results/Cohort_I_imputed_inc_deceased" # Cohort I
  }
}





dir.create(root.path, recursive = TRUE, showWarnings = FALSE)
empirical.esti <- NULL
Long_COVID_df <- Long_COVID_df[!is.na(Long_COVID_df$COVID_19),]

drug_cols <- c(
  "min_Cardiomyopathy_start",
  "min_Cerebrovascular_disease_start",
  "min_Chronic_lung_disease_start",
  "min_Coronary_artery_disease_start",
  "min_Dementia_start",
  "min_Depression_start",
  "min_Diabetes_with_complications_start",
  "min_Diabetes_without_complications_start",
  "min_Heart_failure_start",
  "min_Hypertension_start",
  "min_Kidney_disease_start",
  "min_Liver_disease_start",
  "min_Malignant_tumor_start",
  "min_Benign_tumor_start",
  "min_Myocardial_infarction_start",
  "min_Peripheral_vascular_disease_start"
)
drug_cols <- drug_cols[drug_cols %in% colnames(Long_COVID_df)]



for (ci in selected.Comorbidity) {
  cat("\n========================================\n")
  cat("Running:", ci, "\n")
  cat("========================================\n")

  dml.df <- data.table::copy(Long_COVID_df)
  ind.d     <- which(colnames(dml.df) %in% ci)
  colname_d <- colnames(dml.df)[ind.d]
  if (length(colname_d) != 1) {
    cat("Skipping", ci, ": exposure column not found or ambiguous.\n")
    next
  }

  if(before %in% c(1, 3, 4)){
    # Match final-2.R's pre-COVID confounder construction: records with a
    # documented target condition use its date, while records without one use
    # the COVID-19 diagnosis date. Apply the same rule to Cohorts III-IV.
    cutoff_date <- ifelse(!is.na(dml.df[[colname_d]]),
                          dml.df[[colname_d]],
                          dml.df$COVID_19)
  }

  if(before == 2){
    analysis_end <- ifelse(
      is.na(dml.df$Long_COVID),
      dml.df$COVID_19 + 180,
      pmin(dml.df$Long_COVID - 90, dml.df$COVID_19 + 180)
    )

    # Match final-2.R: treated records use the target-condition date, while
    # untreated records use the end of the analysis window.
    cutoff_date <- ifelse(!is.na(dml.df[[colname_d]]),
                          dml.df[[colname_d]],
                          analysis_end)
  }
  dml.df[, (drug_cols) := lapply(.SD, function(x) {
    as.integer(!is.na(x) & !is.na(cutoff_date) & x < cutoff_date)
  }), .SDcols = drug_cols]

  setDF(dml.df)
  ind.ds <- get_other_comorbidity_cols(colnames(dml.df), ci, indicator.Comorbidity)

  dml.df$COVID_19 <- sanitize_date(dml.df$COVID_19)
  dml.df$Long_COVID <- sanitize_date(dml.df$Long_COVID)
  dml.df[[colname_d]] <- sanitize_date(dml.df[[colname_d]])

  sel_rule <- select_records_from_full_scenario(
    before = before,
    covid_date = dml.df$COVID_19,
    lc_date = dml.df$Long_COVID,
    exposure_date = dml.df[[colname_d]],
    has_other_comorbidity = has_any_recorded_comorbidity(dml.df, ind.ds)
  )

  dml.df$target.condition <- dml.df[[colname_d]]
  dml.df$D <- sel_rule$D
  eligible_rows <- !is.na(sel_rule$eligible) & sel_rule$eligible == 1
  dml.df <- dml.df[eligible_rows, , drop = FALSE]
  # `cutoff_date` was computed before eligibility filtering, so keep it aligned
  # with the rows retained in `dml.df` for the adjustment-date comparisons.
  cutoff_date <- cutoff_date[eligible_rows]


  if (ci %nin% "Any_comorbidities") {
    ind.ds.names <- colnames(dml.df)[ind.ds]

    if (length(ind.ds.names) > 0) {
      for (nm in ind.ds.names) {
        dml.df[[nm]] <- sanitize_date(dml.df[[nm]])
      }
    }

    ## Match final-2.R's pre-COVID rule for Cohorts I, III, and IV: construct other-
    ## comorbidity indicators relative to the target date when documented and
    ## relative to the COVID-19 date when the target condition is absent.
    if(before %in% c(1, 3, 4)){
      dml.df <- dml.df %>%
        mutate(across(
          all_of(ind.ds.names),
          ~ case_when(
            !is.na(.) & !is.na(target.condition) & !is.infinite(.) &
              (as.Date(.) < as.Date(target.condition)) ~ 1,
            !is.na(.) & is.na(target.condition) &
              (as.Date(.) < as.Date(COVID_19)) ~ 1,
            TRUE ~ 0
          )
        ))
    } else {
      dml.df <- dml.df %>%
        mutate(across(
          all_of(ind.ds.names),
          ~ case_when(
            !is.na(.) & !is.na(target.condition) & !is.infinite(.) &
              (as.Date(.) < as.Date(target.condition)) ~ 1,
            !is.na(.) & is.na(target.condition) & !is.na(Long_COVID) &
              (as.Date(.) < as.Date(Long_COVID))  ~ 1,
            !is.na(.) & is.na(target.condition) & is.na(Long_COVID)
            ~ 1,
            TRUE ~ 0
          )
        ))
    }

    if (length(ind.ds.names) > 0) {
      dml.df$comorbidity.count <- rowSums(dml.df[, ind.ds.names, drop = FALSE], na.rm = TRUE)
    } else {
      dml.df$comorbidity.count <- 0
    }
    dml.df$com.count <- ifelse(
      dml.df$comorbidity.count >= 2, 2,
      ifelse(dml.df$comorbidity.count == 1, 1, 0)
    )

    Selected.variable.new <- c(
      ind.ds.names,
      Selected.variable,
      drug_cols
    )

    strata <- interaction(
      dml.df$region,
      dml.df$Covid_19_Year,
      dml.df$com.count,
      drop = TRUE
    )
  } else {
    dml.df$comorbidity.count <- 0

    Selected.variable.new <- c(Selected.variable)

    strata <- interaction(
      dml.df$region,
      dml.df$Covid_19_Year,
      drop = TRUE
    )
  }

  setDT(dml.df)

  #
  dml.df <- dml.df[
    # !is.na(COVID_19) &
    #   !is.na(Long_COVID_binary) &
    !is.na(D)
  ]

  dml.df$D <- as.numeric(dml.df$D)
  dml.df$Long_COVID_binary <- as.numeric(dml.df$Long_COVID_binary)

  # Sample sizes
  n_total <- nrow(dml.df)
  n1 <- sum(dml.df$D == 1)
  n0 <- sum(dml.df$D == 0)
  events_D1 <- sum(dml.df$Long_COVID_binary[dml.df$D == 1], na.rm = TRUE)
  events_D0 <- sum(dml.df$Long_COVID_binary[dml.df$D == 0], na.rm = TRUE)
  exposure_prevalence <- if (n_total > 0) n1 / n_total else NA_real_

  # Risks
  R1 <- mean(dml.df$Long_COVID_binary[dml.df$D == 1])
  R0 <- mean(dml.df$Long_COVID_binary[dml.df$D == 0])

  # Effect measures
  RR <- R1 / R0
  AR <- R1 - R0
  dml_check <- assess_dml_readiness(
    n_total = n_total,
    n1 = n1,
    n0 = n0,
    events_d1 = events_D1,
    events_d0 = events_D0,
    exposure_prevalence = exposure_prevalence,
    ar = AR,
    min_treated_n = min_treated_n_for_dml,
    min_control_n = min_control_n_for_dml,
    min_treated_events = min_treated_events_for_dml,
    min_control_events = min_control_events_for_dml,
    min_exposure_prevalence = min_exposure_prevalence_for_dml,
    min_abs_ar_flag = min_abs_ar_flag_for_dml
  )

  # Output
  empirical.esti <- rbind(empirical.esti, data.frame(
    comorbidity = ci,
    n_total = n_total,
    n_D1 = n1,
    n_D0 = n0,
    events_D1 = events_D1,
    events_D0 = events_D0,
    exposure_prevalence = exposure_prevalence,
    event_rate_D1 = R1,
    event_rate_D0 = R0,
    R1 = R1,
    R0 = R0,
    RR = RR,
    AR = AR,
    dml_status = dml_check$status,
    dml_skip_reason = dml_check$skip_reason,
    dml_flag_reason = dml_check$flag_reason
  ))

  if (dml_check$skip) {
    cat("Skipping DML for", ci, "because", dml_check$skip_reason, "\n")
    next
  }
  if (nzchar(dml_check$flag_reason)) {
    cat("Flagging DML for", ci, "because", dml_check$flag_reason, "\n")
  }

  #
  Selected.variable.new <- unique(Selected.variable.new)
  Selected.variable.new <- Selected.variable.new[Selected.variable.new %in% colnames(dml.df)]

  #
  subgroup_cols <- c("race", "sex", "age_group", "com.count",  "Covid_19_Year", "UrbanRural", "region") # "age_group
  subgroup_cols <- subgroup_cols[subgroup_cols %in% colnames(dml.df)]

  # folds
  fold_id <- createFolds(strata[seq_len(nrow(dml.df))], k = K, list = TRUE, returnTrain = FALSE)

  #
  fit <- estimate_effects_by_subgroup_dr2_h2o_ensemble(
    dat = dml.df,
    y_col = "Long_COVID_binary",
    d_col = "D",
    x_cols = Selected.variable.new,
    s_cols = subgroup_cols,
    id_col = if ("ID" %in% colnames(dml.df)) "ID" else NULL,
    n_folds = K,
    fold_id = fold_id,
    trim = 0.01,
    seed = 12345,
    ntrees = 50,
    max_depth = 6,
    nfolds_ensemble = 5,
    ensemble_models = c("gbm", "drf", "glm"),#c("gbm", "drf", "glm"), #c("gbm", "drf", "dl", "glm")
    metalearner_algorithm = "glm",
    metalearner_params = list(family = "binomial"),
    balance_classes = FALSE,
    do_targeting = TRUE,
    h2o_max_mem_size = "50G",
    return_fold_level = TRUE,
    cleanup_h2o_each_fold = TRUE
  )

  # individual / fold level
  result_dt <- copy(fit$fold_level)
  result_dt[, Comorbidity := ci]

  if ("comorbidity.count" %in% names(dml.df)) {
    tmp_map <- dml.df[, .(row_index = seq_len(.N), comorbidity.count)]
    result_dt <- merge(result_dt, tmp_map, by = "row_index", all.x = TRUE)
  }


  overall.summary <- summarize_overall_from_fold_level(result_dt, do_targeting = TRUE)
  write_parquet(
    overall.summary,
    paste0(
      root.path, "/overall.", ci,
      ".Imputed_non_died_select_covariables_no_BMI_", K, ".parquet"
    )
  )
  summary_dt <- copy(fit$subgroup_results)
  summary_dt[, Comorbidity := ci]

  write_parquet(
    summary_dt,
    paste0(
      root.path, "/SUMMARY.", ci,
      ".Imputed_non_died_select_covariables_no_BMI_", K, ".parquet"
    )
  )


  if(save.detailed.results){

  diag_dir <- file.path(root.path, "diagnostics", ci)
  dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
  file_stub <- paste0(ci, ".Imputed_non_died_select_covariables_no_BMI_", K)

  diag_scores <- copy(fit$individual_scores)
  diag_scores[, `:=`(Comorbidity = ci, Question = before)]
  diag_overview <- build_propensity_summary(
    diag_scores,
    d_col = "D",
    y_col = "Long_COVID_binary",
    ps_col = "m_hat"
  )
  diag_by_treat <- build_propensity_by_treatment(
    diag_scores,
    d_col = "D",
    ps_col = "m_hat"
  )
  diag_bins <- build_propensity_bin_table(
    diag_scores,
    d_col = "D",
    y_col = "Long_COVID_binary",
    ps_col = "m_hat",
    n_bins = 10
  )
  diag_balance <- build_balance_table(
    diag_scores,
    covariate_cols = Selected.variable.new,
    d_col = "D",
    ps_col = "m_hat"
  )
  diag_subgroup_ps <- build_subgroup_positivity_table(
    diag_scores,
    subgroup_cols = subgroup_cols,
    d_col = "D",
    ps_col = "m_hat"
  )
  diag_model_perf <- build_model_performance_table(
    diag_scores,
    y_col = "Long_COVID_binary",
    d_col = "D"
  )
  diag_covariates <- data.table(
    covariate = Selected.variable.new,
    is_subgroup = Selected.variable.new %in% subgroup_cols
  )

  write_parquet(
    diag_scores,
    file.path(diag_dir, paste0("DIAG.SCORES.", file_stub, ".parquet"))
  )

  write_parquet(
    diag_overview,
    file.path(diag_dir, paste0("DIAG.OVERVIEW.", file_stub, ".parquet"))
  )
  write_parquet(
    diag_by_treat,
    file.path(diag_dir, paste0("DIAG.PS_BY_TREATMENT.", file_stub, ".parquet"))
  )
  write_parquet(
    diag_bins,
    file.path(diag_dir, paste0("DIAG.PS_BINS.", file_stub, ".parquet"))
  )
  write_parquet(
    diag_balance,
    file.path(diag_dir, paste0("DIAG.BALANCE.", file_stub, ".parquet"))
  )
  write_parquet(
    diag_subgroup_ps,
    file.path(diag_dir, paste0("DIAG.SUBGROUP_POSITIVITY.", file_stub, ".parquet"))
  )
  write_parquet(
    diag_model_perf,
    file.path(diag_dir, paste0("DIAG.MODEL_PERFORMANCE.", file_stub, ".parquet"))
  )
  write_parquet(
    diag_covariates,
    file.path(diag_dir, paste0("DIAG.COVARIATES.", file_stub, ".parquet"))
  )

  save_diagnostic_plots(
    score_dt = diag_scores,
    bin_dt = diag_bins,
    balance_dt = diag_balance,
    out_dir = diag_dir,
    ci_label = ci
  )



  write_parquet(
    result_dt,
    paste0(
      root.path, "/CV.", ci,
      ".Imputed_non_died_select_covariables_no_BMI_", K, ".parquet"
    )
  )

  rm(diag_scores)
  rm(result_dt)
  }


}


