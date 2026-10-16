############################################################
## dml.R: reusable DML helpers and estimator
############################################################
#' Long COVID Ensemble Double/Debiased Machine Learning Utilities
#'
#' Loads the packages and defines the reusable date, cohort-construction,
#' diagnostics, H2O ensemble, and cross-fitted DML functions used by
#' `main_call.R`. This file intentionally performs no data loading or model
#' fitting when sourced.
#'
#' @keywords internal
library(dplyr)
library(data.table)
library(h2o)
library(glmnet)
library(arrow)
library(caret)
library(Hmisc)
library(ggplot2)
library(mice)

############################################################
## 1. Helper functions
############################################################
#' Compute the Logit Transform
#'
#' Converts probabilities on the unit interval to the log-odds scale.
#'
#' @param p A numeric vector of probabilities.
#'
#' @return A numeric vector on the log-odds scale.
logit <- function(p) qlogis(p)
#' Compute the Inverse Logit Transform
#'
#' Converts values on the log-odds scale to probabilities.
#'
#' @param x A numeric vector on the log-odds scale.
#'
#' @return A numeric vector of probabilities.
inv_logit <- function(x) plogis(x)

#' Coerce Values to a Factor Safely
#'
#' Converts an object to a factor only when it is not already a factor.
#'
#' @param x An input vector.
#'
#' @return A factor vector.
.safe_factor <- function(x) {
  if (!is.factor(x)) factor(x) else x
}

#' Coerce a Binary Variable to Numeric Safely
#'
#' Converts a binary vector to numeric, handling factors first by converting
#' them through character values.
#'
#' @param x An input vector.
#'
#' @return A numeric vector.
.safe_binary_numeric <- function(x) {
  if (is.factor(x)) return(as.numeric(as.character(x)))
  as.numeric(x)
}

#' Apply a Logistic Fluctuation Update
#'
#' Fits a one-dimensional fluctuation model for targeted updating. When the
#' standard logistic regression is unstable, a ridge-penalized fallback is used.
#'
#' @param y A binary outcome vector.
#' @param H A clever covariate vector.
#' @param m_hat Initial predicted probabilities.
#' @param lambda Ridge penalty used in the fallback `glmnet` fit.
#' @param eps Small bound used to truncate probabilities away from 0 and 1.
#'
#' @return A list with the fluctuation coefficient, updated predictions, and
#'   an instability flag.
fluctuate <- function(y, H, m_hat, lambda = 1e-4, eps = 1e-6) {
  m_hat <- pmin(pmax(m_hat, eps), 1 - eps)
  logit_m <- logit(m_hat)

  fit <- tryCatch(
    glm(y ~ -1 + offset(logit_m) + H, family = binomial),
    warning = function(w) w,
    error   = function(e) e
  )

  unstable <- !inherits(fit, "glm") || any(!is.finite(coef(fit)))

  if (unstable) {
    X <- cbind(H = H, dummy = 0)
    fit_ridge <- glmnet(
      x = X,
      y = y,
      family = "binomial",
      offset = logit_m,
      alpha = 0,
      lambda = lambda,
      standardize = FALSE
    )
    epsilon <- as.numeric(coef(fit_ridge)["H", ])
  } else {
    epsilon <- as.numeric(coef(fit))
  }

  m_tilde <- inv_logit(logit_m + H * epsilon)

  list(
    epsilon = epsilon,
    m_tilde = m_tilde,
    unstable = unstable
  )
}

#' Convert Values to `Date`
#'
#' Standardizes common date-like inputs, including factors, numerics, `Date`,
#' and `POSIXt` objects.
#'
#' @param x A vector containing date-like values.
#'
#' @return A `Date` vector.
sanitize_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.factor(x)) x <- as.character(x)

  if (is.numeric(x)) {
    x[!is.finite(x)] <- NA_real_
    return(as.Date(x, origin = "1970-01-01"))
  }

  suppressWarnings(as.Date(x))
}

#' Select Eligible Records for a Time-Window Scenario
#'
#' Creates eligibility and exposure indicators for four study cohorts defined
#' relative to COVID-19 and long COVID dates.
#'
#' @param before Integer cohort code: `1` for Cohort IV, `2` for Cohort II,
#'   `3` for Cohort III, and `4` for Cohort I. Cohort IV otherwise matches
#'   Cohort III but removes the long-COVID-minus-90-days cutoff.
#' @param covid_date A vector of COVID-19 diagnosis dates.
#' @param lc_date A vector of long COVID dates.
#' @param exposure_date A vector of exposure dates.
#' @param has_other_comorbidity Logical vector indicating the presence of other
#'   comorbidities.
#' @return A list with `eligible` and `D` vectors.
select_records_from_full_scenario <- function(before,
                                              covid_date,
                                              lc_date,
                                              exposure_date,
                                              has_other_comorbidity = FALSE) {
  has_other_comorbidity <- rep_len(as.logical(has_other_comorbidity), length(exposure_date))

  has_covid <- !is.na(covid_date)
  has_lc <- !is.na(lc_date)
  has_e <- !is.na(exposure_date)
  q1_no_lc_window_end <- covid_date + 180
  q1_treated_cutoff <- covid_date
  q1_treated_cutoff[has_lc] <- pmin(covid_date[has_lc], lc_date[has_lc] - 90)
  cohort_ii_window_end <- covid_date + 180
  cohort_ii_lc_cutoff <- pmin(lc_date - 90, cohort_ii_window_end)

  e_before_covid <- has_e & has_covid & (exposure_date < q1_treated_cutoff)
  # Cohort IV removes the LC - 90 restriction used by Cohort III. Its treated
  # exposure window ends at COVID-19 diagnosis for every record.
  e_before_covid_without_lc90 <- has_e & has_covid & (exposure_date < covid_date)
  # Cohort I excludes LC cases with a target comorbidity documented from the
  # person-specific preexisting-comorbidity cutoff through before LC diagnosis.
  e_in_q1_lc_exclusion_window <- has_e & has_covid & has_lc &
    (exposure_date >= q1_treated_cutoff) &
    (exposure_date < lc_date)
  e_between_covid_q1_no_lc <- has_e & has_covid & !has_lc &
    (exposure_date >= covid_date) &
    (exposure_date <= q1_no_lc_window_end)
  e_in_cohort_ii_window <- has_e & has_covid &
    (
      (!has_lc & (exposure_date < cohort_ii_window_end)) |
        (has_lc & (exposure_date < cohort_ii_lc_cutoff))
    )

  eligible_cohort_i <- has_covid &
    !e_in_q1_lc_exclusion_window &
    !e_between_covid_q1_no_lc
  eligible_cohort_ii <- rep(TRUE, length(exposure_date))
  eligible_cohort_iii <- has_covid
  eligible_cohort_iv <- eligible_cohort_iii

  d_cohort_i <- as.integer(e_before_covid)
  d_cohort_ii <- as.integer(e_in_cohort_ii_window)
  d_cohort_iii <- as.integer(e_before_covid)
  d_cohort_iv <- as.integer(e_before_covid_without_lc90)

  if (before == 1) {
    list(eligible = as.integer(eligible_cohort_iv), D = d_cohort_iv)
  } else if (before == 2) {
    list(eligible = as.integer(eligible_cohort_ii), D = d_cohort_ii)
  } else if (before == 3) {
    list(eligible = as.integer(eligible_cohort_iii), D = d_cohort_iii)
  } else if (before == 4) {
    list(eligible = as.integer(eligible_cohort_i), D = d_cohort_i)
  } else {
    stop("before must be 1 (Cohort IV), 2 (Cohort II), 3 (Cohort III), or 4 (Cohort I).")
  }
}

#' Check Whether a DML Analysis Should Run
#'
#' Evaluates sample size, event count, prevalence, and effect-size thresholds
#' to decide whether a DML run should proceed, be flagged, or be skipped.
#'
#' @param n_total Total number of eligible records.
#' @param n1 Number of treated records.
#' @param n0 Number of control records.
#' @param events_d1 Number of events among treated records.
#' @param events_d0 Number of events among control records.
#' @param exposure_prevalence Prevalence of the exposure.
#' @param ar Empirical attributable risk.
#' @param min_treated_n Minimum treated sample size.
#' @param min_control_n Minimum control sample size.
#' @param min_treated_events Minimum treated event count.
#' @param min_control_events Minimum control event count.
#' @param min_exposure_prevalence Minimum exposure prevalence.
#' @param min_abs_ar_flag Minimum absolute attributable risk before flagging.
#'
#' @return A list describing whether to run, flag, or skip the analysis.
assess_dml_readiness <- function(n_total,
                                 n1,
                                 n0,
                                 events_d1,
                                 events_d0,
                                 exposure_prevalence,
                                 ar,
                                 min_treated_n = 50000,
                                 min_control_n = 50000,
                                 min_treated_events = 1000,
                                 min_control_events = 1000,
                                 min_exposure_prevalence = 0.01,
                                 min_abs_ar_flag = 0.001) {
  skip_reasons <- character(0)
  flags <- character(0)

  if (!is.finite(n_total) || n_total <= 0) {
    skip_reasons <- c(skip_reasons, "no eligible records")
  }
  if (!is.finite(n1) || n1 < min_treated_n) {
    skip_reasons <- c(skip_reasons, paste0("n_D1<", min_treated_n))
  }
  if (!is.finite(n0) || n0 < min_control_n) {
    skip_reasons <- c(skip_reasons, paste0("n_D0<", min_control_n))
  }
  if (!is.finite(events_d1) || events_d1 < min_treated_events) {
    skip_reasons <- c(skip_reasons, paste0("events_D1<", min_treated_events))
  }
  if (!is.finite(events_d0) || events_d0 < min_control_events) {
    skip_reasons <- c(skip_reasons, paste0("events_D0<", min_control_events))
  }
  if (!is.finite(exposure_prevalence) || exposure_prevalence < min_exposure_prevalence) {
    skip_reasons <- c(skip_reasons, paste0("exposure_prevalence<", min_exposure_prevalence))
  }
  if (!is.finite(ar) || is.na(ar)) {
    flags <- c(flags, "AR_missing")
  } else if (abs(ar) < min_abs_ar_flag) {
    flags <- c(flags, paste0("abs_AR<", min_abs_ar_flag))
  }

  list(
    skip = length(skip_reasons) > 0,
    status = if (length(skip_reasons) > 0) "skip" else if (length(flags) > 0) "flag" else "run",
    skip_reason = if (length(skip_reasons) > 0) paste(skip_reasons, collapse = "; ") else "",
    flag_reason = if (length(flags) > 0) paste(flags, collapse = "; ") else ""
  )
}

#' Get Other Comorbidity Column Indices
#'
#' Identifies comorbidity indicator columns excluding the components that belong
#' to the currently analyzed grouped comorbidity.
#'
#' @param df_colnames Column names from the data set.
#' @param ci Name of the comorbidity under analysis.
#' @param indicator_comorbidity Character vector of comorbidity indicator names.
#'
#' @return An integer vector of column indices.
get_other_comorbidity_cols <- function(df_colnames, ci, indicator_comorbidity) {
  exclude_cols <- grouped_comorbidity_components(ci)
  which(df_colnames %in% indicator_comorbidity[indicator_comorbidity %nin% exclude_cols])
}

#' Detect Any Recorded Comorbidity
#'
#' Returns a logical flag for whether each row has any recorded value among the
#' specified comorbidity columns.
#'
#' @param df A data frame.
#' @param col_idx Integer indices of columns to inspect.
#'
#' @return A logical vector.
has_any_recorded_comorbidity <- function(df, col_idx) {
  if (length(col_idx) == 0) return(rep(FALSE, nrow(df)))
  flags <- lapply(col_idx, function(j) {
    x <- df[[j]]
    if (is.factor(x)) x <- as.character(x)
    if (inherits(x, "Date") || inherits(x, "POSIXt")) return(!is.na(x))
    if (is.numeric(x)) return(!is.na(x) & is.finite(x) & x != 0)
    !is.na(x) & trimws(as.character(x)) != ""
  })
  Reduce(`|`, flags, init = rep(FALSE, nrow(df)))
}

#' Expand a Grouped Comorbidity Definition
#'
#' Maps grouped comorbidity labels to their component condition names.
#'
#' @param ci A comorbidity name.
#'
#' @return A character vector of component names.
grouped_comorbidity_components <- function(ci) {
  if (ci %in% "Diabetes") {
    return(c("Diabetes_with_complications", "Diabetes_without_complications"))
  }
  if (ci %in% "Cardiovascular") {
    return(c(
      "Heart_failure",
      "Coronary_artery_disease",
      "Peripheral_vascular_disease",
      "Cardiomyopathy",
      "Myocardial_infarction"
    ))
  }
  if (ci %in% "Neuro_psych_disorders") {
    return(c("Cerebrovascular_disease", "Dementia", "Depression"))
  }
  if (ci %in% "Cancer") {
    return(c("Malignant_tumor", "Benign_tumor"))
  }
  ci
}

#' Find the Earliest Recorded Date Across Columns
#'
#' Computes the row-wise minimum date across a set of candidate columns.
#'
#' @param df A data frame.
#' @param cols Character vector of column names to search.
#'
#' @return A `Date` vector containing the earliest recorded date per row.
earliest_recorded_date <- function(df, cols) {
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(rep(as.Date(NA), nrow(df)))
  dates <- lapply(cols, function(nm) sanitize_date(df[[nm]]))
  out <- do.call(pmin, c(dates, list(na.rm = TRUE)))
  out[!is.finite(out)] <- as.Date(NA)
  out
}

#' Initialize an H2O Cluster Safely
#'
#' Starts an H2O cluster only when one is not already available, then disables
#' progress output.
#'
#' @param max_mem_size Maximum memory size passed to `h2o::h2o.init()`.
#' @param nthreads Number of threads passed to `h2o::h2o.init()`.
#'
#' @return No return value. Called for its side effects.
init_h2o_safe <- function(max_mem_size = "100G", nthreads = 3) {
  if (!requireNamespace("h2o", quietly = TRUE)) {
    stop("Package 'h2o' is required.")
  }

  # Check whether a cluster is already available using a safe probe.
  is_up <- FALSE
  try({
    is_up <- h2o::h2o.clusterIsUp()
  }, silent = TRUE)

  if (!is_up) {
    h2o::h2o.init(max_mem_size = max_mem_size, nthreads = nthreads)
  }

  h2o::h2o.no_progress()
}

#' Fit a Base H2O Binary Model
#'
#' Trains one H2O binary classifier used as a base learner in an ensemble.
#'
#' @param train_h2o An H2O frame containing the training data.
#' @param x_cols Character vector of predictor column names.
#' @param y_col Outcome column name.
#' @param model_type Type of H2O model to fit.
#' @param ntrees Number of trees for tree-based models.
#' @param max_depth Maximum tree depth for tree-based models.
#' @param nfolds Number of cross-validation folds.
#' @param seed Random seed.
#' @param keep_cv_models Logical flag indicating whether to keep CV models.
#' @param balance_classes Logical flag indicating whether to balance classes.
#'
#' @return A fitted H2O model object.
fit_h2o_base_model <- function(train_h2o,
                               x_cols,
                               y_col,
                               model_type = c("gbm", "drf", "dl", "glm"),
                               ntrees = 200,
                               max_depth = 6,
                               nfolds = 3,
                               seed = 123,
                               keep_cv_models = FALSE,
                               balance_classes = TRUE) {
  model_type <- match.arg(model_type)

  if (model_type == "gbm") {
    fit <- h2o.gbm(
      x = x_cols,
      y = y_col,
      training_frame = train_h2o,
      distribution = "bernoulli",
      ntrees = ntrees,
      max_depth = max_depth,
      learn_rate = 0.05,
      nfolds = nfolds,
      keep_cross_validation_models = keep_cv_models,
      keep_cross_validation_predictions = TRUE,
      balance_classes = balance_classes,
      seed = seed
    )
  } else if (model_type == "drf") {
    fit <- h2o.randomForest(
      x = x_cols,
      y = y_col,
      training_frame = train_h2o,
      ntrees = ntrees,
      max_depth = max_depth,
      nfolds = nfolds,
      keep_cross_validation_models = keep_cv_models,
      keep_cross_validation_predictions = TRUE,
      balance_classes = balance_classes,
      binomial_double_trees = TRUE,
      seed = seed
    )
  } else if (model_type == "dl") {
    fit <- h2o.deeplearning(
      x = x_cols,
      y = y_col,
      training_frame = train_h2o,
      distribution = "bernoulli",
      hidden = c(32, 16),
      activation = "Rectifier",
      epochs = 50,
      nfolds = nfolds,
      keep_cross_validation_models = keep_cv_models,
      keep_cross_validation_predictions = TRUE,
      balance_classes = balance_classes,
      seed = seed
    )
  } else if (model_type == "glm") {
    fit <- h2o.glm(
      x = x_cols,
      y = y_col,
      training_frame = train_h2o,
      family = "binomial",
      nfolds = nfolds,
      keep_cross_validation_models = keep_cv_models,
      keep_cross_validation_predictions = TRUE,
      seed = seed
    )
  }

  fit
}

#' Fit a Stacked Binary H2O Ensemble
#'
#' Fits several base binary learners and combines them through an H2O stacked
#' ensemble. Small or degenerate samples fall back to mean predictions.
#'
#' @param train_df Training data frame.
#' @param test_df Test data frame.
#' @param y_col Outcome column name.
#' @param x_cols Character vector of predictor column names.
#' @param ensemble_models Character vector of base learner types.
#' @param metalearner_algorithm Name of the H2O metalearner algorithm.
#' @param metalearner_params Named list of metalearner settings.
#' @param ntrees Number of trees for tree-based base learners.
#' @param max_depth Maximum tree depth for tree-based base learners.
#' @param nfolds Number of cross-validation folds for each base learner.
#' @param seed Random seed.
#' @param balance_classes Logical flag indicating whether to balance classes.
#' @param min_rows_required Minimum number of rows needed before model fitting.
#' @param positive_min_count Minimum class count required for both classes.
#' @param trim Lower and upper truncation bound for predictions.
#'
#' @return A list containing predictions, fit mode, fallback mean prediction,
#'   and model identifiers.
fit_h2o_stacked_binary <- function(train_df,
                                   test_df,
                                   y_col,
                                   x_cols,
                                   ensemble_models = c("gbm", "drf", "glm"),
                                   metalearner_algorithm = "glm",
                                   metalearner_params = list(family = "binomial"),
                                   ntrees = 200,
                                   max_depth = 6,
                                   nfolds = 3,
                                   seed = 123,
                                   balance_classes = TRUE,
                                   min_rows_required = 30,
                                   positive_min_count = 2,
                                   trim = 0.01) {
  if (nrow(test_df) == 0) {
    return(list(
      pred = numeric(0),
      mode = "empty_test",
      mean_pred = NA_real_,
      model_ids = character(0)
    ))
  }

  if (nrow(train_df) == 0) {
    return(list(
      pred = rep(0.5, nrow(test_df)),
      mode = "mean_empty_train",
      mean_pred = 0.5,
      model_ids = character(0)
    ))
  }

  if (nrow(train_df) < min_rows_required) {
    mean_pred <- mean(.safe_binary_numeric(train_df[[y_col]]), na.rm = TRUE)
    if (!is.finite(mean_pred)) mean_pred <- 0.5
    mean_pred <- pmin(pmax(mean_pred, trim), 1 - trim)
    return(list(
      pred = rep(mean_pred, nrow(test_df)),
      mode = "mean_small_sample",
      mean_pred = mean_pred,
      model_ids = character(0)
    ))
  }

  y_num <- .safe_binary_numeric(train_df[[y_col]])
  tb <- table(y_num)

  if (length(tb) < 2 || min(tb) <= positive_min_count) {
    mean_pred <- mean(y_num, na.rm = TRUE)
    if (!is.finite(mean_pred)) mean_pred <- 0.5
    mean_pred <- pmin(pmax(mean_pred, trim), 1 - trim)
    return(list(
      pred = rep(mean_pred, nrow(test_df)),
      mode = "mean_single_class_or_too_few",
      mean_pred = mean_pred,
      model_ids = character(0)
    ))
  }

  train_df <- setDF(train_df)
  test_df  <- setDF(test_df)
  train_df[[y_col]] <- .safe_factor(train_df[[y_col]])

  train_h2o <- as.h2o(train_df)
  test_h2o  <- as.h2o(test_df)

  for (nm in x_cols) {
    if (is.character(train_df[[nm]]) || is.factor(train_df[[nm]])) {
      train_h2o[[nm]] <- as.factor(train_h2o[[nm]])
      test_h2o[[nm]]  <- as.factor(test_h2o[[nm]])
    }
  }
  train_h2o[[y_col]] <- as.factor(train_h2o[[y_col]])

  base_models <- list()
  model_ids <- character(0)

  for (j in seq_along(ensemble_models)) {
    mdl <- fit_h2o_base_model(
      train_h2o = train_h2o,
      x_cols = x_cols,
      y_col = y_col,
      model_type = ensemble_models[j],
      ntrees = ntrees,
      max_depth = max_depth,
      nfolds = nfolds,
      seed = seed + j,
      keep_cv_models = FALSE,
      balance_classes = balance_classes
    )
    base_models[[j]] <- mdl
    model_ids <- c(model_ids, mdl@model_id)
  }

  ens <- h2o.stackedEnsemble(
    x = x_cols,
    y = y_col,
    training_frame = train_h2o,
    base_models = model_ids,
    metalearner_algorithm = metalearner_algorithm,
    metalearner_params = metalearner_params,
    seed = seed + 999
  )

  pred_h2o <- h2o.predict(ens, test_h2o)
  pred <- as.vector(pred_h2o[, "p1"])
  pred <- pmin(pmax(pred, trim), 1 - trim)

  list(
    pred = pred,
    mode = "ensemble",
    mean_pred = NA_real_,
    model_ids = c(model_ids, ens@model_id)
  )
}

library(data.table)

#' Summarize Fold-Level Scores Overall
#'
#' Aggregates influence-function-based fold-level estimates into overall risks,
#' attributable risk, and relative risk summaries.
#'
#' @param fold_level A data frame or data table of fold-level predictions and
#'   influence-function values.
#' @param do_targeting Logical flag indicating whether targeted estimates should
#'   also be summarized.
#'
#' @return A one-row `data.table` of overall summary statistics.
summarize_overall_from_fold_level <- function(fold_level, do_targeting = TRUE) {
  summarize_phi <- function(phi1, phi0, phi_ar, prefix = "") {
    ok <- is.finite(phi1) & is.finite(phi0) & is.finite(phi_ar)
    n_eff <- sum(ok)
    if (n_eff == 0) {
      out <- list(
        n = 0L,
        risk1 = NA_real_,
        SE_risk1 = NA_real_,
        risk1_CI_low = NA_real_,
        risk1_CI_high = NA_real_,
        risk0 = NA_real_,
        SE_risk0 = NA_real_,
        risk0_CI_low = NA_real_,
        risk0_CI_high = NA_real_,
        AR = NA_real_,
        SE_AR = NA_real_,
        AR_CI_low = NA_real_,
        AR_CI_high = NA_real_,
        RR = NA_real_,
        SE_log_RR = NA_real_,
        RR_CI_low = NA_real_,
        RR_CI_high = NA_real_
      )
    } else {
      p1 <- phi1[ok]
      p0 <- phi0[ok]
      pa <- phi_ar[ok]
      mu1 <- mean(p1)
      mu0 <- mean(p0)
      ar  <- mean(pa)
      se_risk1 <- if (n_eff > 1) sd(p1) / sqrt(n_eff) else NA_real_
      se_risk0 <- if (n_eff > 1) sd(p0) / sqrt(n_eff) else NA_real_
      se_ar <- if (n_eff > 1) sd(pa) / sqrt(n_eff) else NA_real_
      rr <- if (mu0 <= 0) NA_real_ else mu1 / mu0

      if (is.na(mu1) || is.na(mu0) || mu1 <= 0 || mu0 <= 0 || n_eff <= 1) {
        se_log_rr <- NA_real_
        rr_ci_low <- NA_real_
        rr_ci_high <- NA_real_
      } else {
        # EIF-based delta-method for log(RR)
        eif_log_rr <- (p1 - mu1) / mu1 - (p0 - mu0) / mu0
        se_log_rr <- sd(eif_log_rr) / sqrt(n_eff)
        rr_ci_low <- exp(log(rr) - 1.96 * se_log_rr)
        rr_ci_high <- exp(log(rr) + 1.96 * se_log_rr)
      }

      out <- list(
        n = n_eff,
        risk1 = mu1,
        SE_risk1 = se_risk1,
        risk1_CI_low = ifelse(is.na(se_risk1), NA_real_, mu1 - 1.96 * se_risk1),
        risk1_CI_high = ifelse(is.na(se_risk1), NA_real_, mu1 + 1.96 * se_risk1),
        risk0 = mu0,
        SE_risk0 = se_risk0,
        risk0_CI_low = ifelse(is.na(se_risk0), NA_real_, mu0 - 1.96 * se_risk0),
        risk0_CI_high = ifelse(is.na(se_risk0), NA_real_, mu0 + 1.96 * se_risk0),
        AR = ar,
        SE_AR = se_ar,
        AR_CI_low = ifelse(is.na(se_ar), NA_real_, ar - 1.96 * se_ar),
        AR_CI_high = ifelse(is.na(se_ar), NA_real_, ar + 1.96 * se_ar),
        RR = rr,
        SE_log_RR = se_log_rr,
        RR_CI_low = rr_ci_low,
        RR_CI_high = rr_ci_high
      )
    }

    if (prefix == "") return(out)
    names(out) <- ifelse(names(out) == "n", "n", paste0(prefix, names(out)))
    out
  }

  base_sum <- summarize_phi(fold_level$phi1, fold_level$phi0, fold_level$phi_ar, prefix = "")
  out <- data.table::as.data.table(base_sum)

  if (do_targeting &&
      all(c("c_phi1", "c_phi0", "c_phi_ar") %in% names(fold_level))) {
    c_sum <- summarize_phi(fold_level$c_phi1, fold_level$c_phi0, fold_level$c_phi_ar, prefix = "c_")
    c_sum$n <- NULL
    out[, (names(c_sum)) := c_sum]
  }

  return(out)
}

#' Compute AUC Safely
#'
#' Calculates the area under the ROC curve and returns `NA` when the input does
#' not support a valid calculation.
#'
#' @param y Observed binary outcomes.
#' @param score Predicted scores or probabilities.
#'
#' @return A numeric scalar.
safe_auc <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  y <- y[ok]
  score <- score[ok]
  if (length(y) <= 1 || length(unique(y)) < 2) return(NA_real_)
  tryCatch(
    as.numeric(Hmisc::somers2(score, y)["C"]),
    error = function(e) NA_real_
  )
}

#' Compute the Brier Score Safely
#'
#' Calculates the mean squared error between observed binary outcomes and
#' predicted probabilities after removing invalid values.
#'
#' @param y Observed binary outcomes.
#' @param score Predicted probabilities.
#'
#' @return A numeric scalar.
safe_brier <- function(y, score) {
  ok <- is.finite(y) & is.finite(score)
  if (!any(ok)) return(NA_real_)
  mean((y[ok] - score[ok])^2)
}

#' Compute a Weighted Mean Safely
#'
#' Computes a weighted mean after filtering out non-finite values.
#'
#' @param x A numeric vector.
#' @param w A numeric vector of weights.
#'
#' @return A numeric scalar.
weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

#' Compute a Weighted Variance Safely
#'
#' Computes a normalized weighted variance after removing non-finite values.
#'
#' @param x A numeric vector.
#' @param w A numeric vector of weights.
#'
#' @return A numeric scalar.
weighted_var_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  if (sum(ok) <= 1) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  w <- w / sum(w)
  mu <- sum(w * x)
  sum(w * (x - mu)^2)
}

#' Compute a Standardized Mean Difference
#'
#' Calculates the standardized mean difference for a numeric covariate between
#' treated and control groups, optionally using weights.
#'
#' @param x A numeric covariate vector.
#' @param d A binary treatment indicator.
#' @param w Optional weights.
#'
#' @return A numeric scalar.
compute_smd_numeric <- function(x, d, w = NULL) {
  ok <- is.finite(x) & is.finite(d)
  if (!is.null(w)) ok <- ok & is.finite(w)
  x <- x[ok]
  d <- d[ok]
  if (!is.null(w)) w <- w[ok]
  if (length(x) <= 1 || sum(d == 1) == 0 || sum(d == 0) == 0) return(NA_real_)

  if (is.null(w)) {
    mu1 <- mean(x[d == 1])
    mu0 <- mean(x[d == 0])
    v1 <- var(x[d == 1])
    v0 <- var(x[d == 0])
  } else {
    mu1 <- weighted_mean_safe(x[d == 1], w[d == 1])
    mu0 <- weighted_mean_safe(x[d == 0], w[d == 0])
    v1 <- weighted_var_safe(x[d == 1], w[d == 1])
    v0 <- weighted_var_safe(x[d == 0], w[d == 0])
  }

  denom <- sqrt((v1 + v0) / 2)
  if (!is.finite(denom) || denom == 0) {
    if (isTRUE(all.equal(mu1, mu0))) return(0)
    return(NA_real_)
  }
  (mu1 - mu0) / denom
}

#' Build a Covariate Balance Table
#'
#' Computes unweighted and IPW-weighted standardized mean differences for
#' numeric and categorical covariates.
#'
#' @param dt A data table containing treatment, propensity scores, and
#'   covariates.
#' @param covariate_cols Character vector of covariate names.
#' @param d_col Treatment column name.
#' @param ps_col Propensity score column name.
#'
#' @return A `data.table` summarizing balance by covariate.
build_balance_table <- function(dt, covariate_cols, d_col = "D", ps_col = "m_hat") {
  if (length(covariate_cols) == 0) return(data.table())

  d <- dt[[d_col]]
  ps <- dt[[ps_col]]
  w_ate <- ifelse(is.finite(d) & d == 1, 1 / ps,
                  ifelse(is.finite(d) & d == 0, 1 / (1 - ps), NA_real_))

  balance_rows <- lapply(covariate_cols, function(v) {
    if (!v %in% names(dt)) return(NULL)
    x <- dt[[v]]

    if (inherits(x, "Date") || inherits(x, "POSIXt")) {
      x <- as.numeric(x)
    }

    if (is.character(x) || is.factor(x)) {
      x_chr <- as.character(x)
      levs <- sort(unique(x_chr[!is.na(x_chr)]))
      if (length(levs) == 0) return(NULL)
      rbindlist(lapply(levs, function(lv) {
        z <- as.numeric(x_chr == lv)
        data.table(
          covariate = v,
          level = lv,
          type = "categorical",
          smd_unweighted = compute_smd_numeric(z, d),
          smd_weighted = compute_smd_numeric(z, d, w_ate)
        )
      }), fill = TRUE)
    } else {
      z <- suppressWarnings(as.numeric(x))
      data.table(
        covariate = v,
        level = NA_character_,
        type = "numeric",
        smd_unweighted = compute_smd_numeric(z, d),
        smd_weighted = compute_smd_numeric(z, d, w_ate)
      )
    }
  })

  out <- rbindlist(balance_rows, fill = TRUE)
  if (nrow(out) == 0) return(out)
  out[, `:=`(
    abs_smd_unweighted = abs(smd_unweighted),
    abs_smd_weighted = abs(smd_weighted)
  )]
  out[order(-abs_smd_unweighted)]
}

#' Build an Overall Propensity Summary
#'
#' Produces high-level diagnostics for the propensity score distribution and
#' overlap in the analysis sample.
#'
#' @param dt A data table containing treatment, outcome, and propensity score
#'   columns.
#' @param d_col Treatment column name.
#' @param y_col Outcome column name.
#' @param ps_col Propensity score column name.
#'
#' @return A one-row `data.table` of propensity summary metrics.
build_propensity_summary <- function(dt, d_col = "D", y_col = "Long_COVID_binary", ps_col = "m_hat") {
  d <- dt[[d_col]]
  y <- dt[[y_col]]
  ps <- dt[[ps_col]]
  ok_ps <- is.finite(d) & is.finite(ps)
  if (sum(ok_ps) == 0) {
    return(data.table(
      n_all = nrow(dt),
      n_nonmissing_D = sum(is.finite(d)),
      n_nonmissing_Y = sum(is.finite(y)),
      treated_rate = mean(d[is.finite(d)]),
      outcome_rate = mean(y[is.finite(y)]),
      ps_mean = NA_real_,
      ps_sd = NA_real_,
      ps_min = NA_real_,
      ps_p01 = NA_real_,
      ps_p05 = NA_real_,
      ps_p50 = NA_real_,
      ps_p95 = NA_real_,
      ps_p99 = NA_real_,
      ps_max = NA_real_,
      pct_ps_lt_0_01 = NA_real_,
      pct_ps_lt_0_05 = NA_real_,
      pct_ps_gt_0_95 = NA_real_,
      pct_ps_gt_0_99 = NA_real_
    ))
  }

  data.table(
    n_all = nrow(dt),
    n_nonmissing_D = sum(is.finite(d)),
    n_nonmissing_Y = sum(is.finite(y)),
    treated_rate = mean(d[is.finite(d)]),
    outcome_rate = mean(y[is.finite(y)]),
    ps_mean = mean(ps[ok_ps]),
    ps_sd = sd(ps[ok_ps]),
    ps_min = min(ps[ok_ps]),
    ps_p01 = as.numeric(quantile(ps[ok_ps], 0.01)),
    ps_p05 = as.numeric(quantile(ps[ok_ps], 0.05)),
    ps_p50 = as.numeric(quantile(ps[ok_ps], 0.50)),
    ps_p95 = as.numeric(quantile(ps[ok_ps], 0.95)),
    ps_p99 = as.numeric(quantile(ps[ok_ps], 0.99)),
    ps_max = max(ps[ok_ps]),
    pct_ps_lt_0_01 = mean(ps[ok_ps] < 0.01),
    pct_ps_lt_0_05 = mean(ps[ok_ps] < 0.05),
    pct_ps_gt_0_95 = mean(ps[ok_ps] > 0.95),
    pct_ps_gt_0_99 = mean(ps[ok_ps] > 0.99)
  )
}

#' Summarize Propensity Scores by Treatment Group
#'
#' Computes distributional summaries of the propensity score separately for the
#' treated and control groups.
#'
#' @param dt A data table.
#' @param d_col Treatment column name.
#' @param ps_col Propensity score column name.
#'
#' @return A `data.table` grouped by treatment status.
build_propensity_by_treatment <- function(dt, d_col = "D", ps_col = "m_hat") {
  dt[is.finite(get(d_col)) & is.finite(get(ps_col)),
     .(
       n = .N,
       mean_ps = mean(get(ps_col)),
       sd_ps = sd(get(ps_col)),
       min_ps = min(get(ps_col)),
       p05_ps = as.numeric(quantile(get(ps_col), 0.05)),
       p25_ps = as.numeric(quantile(get(ps_col), 0.25)),
       p50_ps = as.numeric(quantile(get(ps_col), 0.50)),
       p75_ps = as.numeric(quantile(get(ps_col), 0.75)),
       p95_ps = as.numeric(quantile(get(ps_col), 0.95)),
       max_ps = max(get(ps_col))
     ),
     by = .(D = get(d_col))]
}

#' Build a Binned Propensity Calibration Table
#'
#' Groups observations by propensity score bins and summarizes treatment and
#' outcome rates within each bin.
#'
#' @param dt A data table.
#' @param d_col Treatment column name.
#' @param y_col Outcome column name.
#' @param ps_col Propensity score column name.
#' @param n_bins Number of propensity score bins.
#'
#' @return A `data.table` with one row per propensity score bin.
build_propensity_bin_table <- function(dt, d_col = "D", y_col = "Long_COVID_binary", ps_col = "m_hat", n_bins = 10) {
  tmp <- copy(dt[is.finite(get(d_col)) & is.finite(get(ps_col))])
  if (nrow(tmp) == 0) return(data.table())

  probs <- seq(0, 1, length.out = n_bins + 1)
  breaks <- unique(as.numeric(quantile(tmp[[ps_col]], probs = probs, na.rm = TRUE)))
  if (length(breaks) < 2) {
    min_ps <- min(tmp[[ps_col]], na.rm = TRUE)
    max_ps <- max(tmp[[ps_col]], na.rm = TRUE)
    if (min_ps == max_ps) max_ps <- min_ps + 1e-8
    breaks <- seq(min_ps, max_ps, length.out = 3)
  }

  tmp[, ps_bin := cut(get(ps_col), breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)]
  tmp[, .(
    n = .N,
    mean_ps = mean(get(ps_col)),
    observed_treat_rate = mean(get(d_col)),
    outcome_rate = mean(get(y_col), na.rm = TRUE)
  ), by = ps_bin]
}

#' Build Subgroup Positivity Diagnostics
#'
#' Summarizes treatment rates and propensity score support within each subgroup.
#'
#' @param dt A data table.
#' @param subgroup_cols Character vector of subgroup variable names.
#' @param d_col Treatment column name.
#' @param ps_col Propensity score column name.
#'
#' @return A `data.table` of subgroup-level positivity diagnostics.
build_subgroup_positivity_table <- function(dt, subgroup_cols, d_col = "D", ps_col = "m_hat") {
  subgroup_rows <- lapply(subgroup_cols, function(svar) {
    if (!svar %in% names(dt)) return(NULL)
    tmp <- copy(dt[!is.na(get(svar)) & is.finite(get(ps_col)) & is.finite(get(d_col))])
    if (nrow(tmp) == 0) return(NULL)
    res <- tmp[, .(
      n = .N,
      treated_rate = mean(get(d_col)),
      mean_ps = mean(get(ps_col)),
      min_ps = min(get(ps_col)),
      p05_ps = as.numeric(quantile(get(ps_col), 0.05)),
      p50_ps = as.numeric(quantile(get(ps_col), 0.50)),
      p95_ps = as.numeric(quantile(get(ps_col), 0.95)),
      max_ps = max(get(ps_col)),
      pct_ps_lt_0_05 = mean(get(ps_col) < 0.05),
      pct_ps_gt_0_95 = mean(get(ps_col) > 0.95)
    ), by = svar]
    setnames(res, svar, "subgroup_level")
    res[, subgroup_var := svar]
    res
  })

  out <- rbindlist(subgroup_rows, fill = TRUE)
  if (nrow(out) == 0) return(out)
  out[]
}

#' Build Model Performance Diagnostics
#'
#' Computes Brier scores and AUC values for propensity and outcome models using
#' available observed outcomes.
#'
#' @param dt A data table containing predictions and observed values.
#' @param y_col Outcome column name.
#' @param d_col Treatment column name.
#'
#' @return A `data.table` of model performance summaries.
build_model_performance_table <- function(dt, y_col = "Long_COVID_binary", d_col = "D") {
  out <- list(
    data.table(
      model = "propensity",
      subset = "all_nonmissing_D",
      n = sum(is.finite(dt[[d_col]]) & is.finite(dt$m_hat)),
      brier = safe_brier(dt[[d_col]], dt$m_hat),
      auc = safe_auc(dt[[d_col]], dt$m_hat)
    ),
    data.table(
      model = "outcome_g1",
      subset = "treated_nonmissing_Y",
      n = sum(dt[[d_col]] == 1 & is.finite(dt[[y_col]]) & is.finite(dt$g1_hat), na.rm = TRUE),
      brier = safe_brier(dt[[y_col]][dt[[d_col]] == 1], dt$g1_hat[dt[[d_col]] == 1]),
      auc = safe_auc(dt[[y_col]][dt[[d_col]] == 1], dt$g1_hat[dt[[d_col]] == 1])
    ),
    data.table(
      model = "outcome_g0",
      subset = "control_nonmissing_Y",
      n = sum(dt[[d_col]] == 0 & is.finite(dt[[y_col]]) & is.finite(dt$g0_hat), na.rm = TRUE),
      brier = safe_brier(dt[[y_col]][dt[[d_col]] == 0], dt$g0_hat[dt[[d_col]] == 0]),
      auc = safe_auc(dt[[y_col]][dt[[d_col]] == 0], dt$g0_hat[dt[[d_col]] == 0])
    )
  )

  if (all(c("c_g1_hat", "c_g0_hat") %in% names(dt))) {
    out <- c(out, list(
      data.table(
        model = "targeted_outcome_g1",
        subset = "treated_nonmissing_Y",
        n = sum(dt[[d_col]] == 1 & is.finite(dt[[y_col]]) & is.finite(dt$c_g1_hat), na.rm = TRUE),
        brier = safe_brier(dt[[y_col]][dt[[d_col]] == 1], dt$c_g1_hat[dt[[d_col]] == 1]),
        auc = safe_auc(dt[[y_col]][dt[[d_col]] == 1], dt$c_g1_hat[dt[[d_col]] == 1])
      ),
      data.table(
        model = "targeted_outcome_g0",
        subset = "control_nonmissing_Y",
        n = sum(dt[[d_col]] == 0 & is.finite(dt[[y_col]]) & is.finite(dt$c_g0_hat), na.rm = TRUE),
        brier = safe_brier(dt[[y_col]][dt[[d_col]] == 0], dt$c_g0_hat[dt[[d_col]] == 0]),
        auc = safe_auc(dt[[y_col]][dt[[d_col]] == 0], dt$c_g0_hat[dt[[d_col]] == 0])
      )
    ))
  }

  rbindlist(out, fill = TRUE)
}

#' Save Diagnostic Plots
#'
#' Writes propensity overlap, propensity calibration, and covariate balance
#' plots to disk when the necessary inputs are available.
#'
#' @param score_dt A data table containing treatment and propensity scores.
#' @param bin_dt A binned propensity summary table.
#' @param balance_dt A covariate balance table.
#' @param out_dir Output directory for image files.
#' @param ci_label Label used in plot titles.
#'
#' @return No return value. Called for its side effects.
save_diagnostic_plots <- function(score_dt, bin_dt, balance_dt, out_dir, ci_label) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ps_plot_dt <- score_dt[is.finite(D) & is.finite(m_hat)]
  if (nrow(ps_plot_dt) > 0) {
    p1 <- ggplot(ps_plot_dt, aes(x = m_hat, fill = factor(D), color = factor(D))) +
      geom_density(alpha = 0.2) +
      labs(
        title = paste("Propensity overlap -", ci_label),
        x = "Estimated propensity score",
        y = "Density",
        fill = "D",
        color = "D"
      ) +
      theme_minimal()
    ggsave(filename = file.path(out_dir, "propensity_overlap.png"), plot = p1, width = 8, height = 5, dpi = 300)
  }

  if (nrow(bin_dt) > 0) {
    p2 <- ggplot(bin_dt, aes(x = mean_ps, y = observed_treat_rate)) +
      geom_point(size = 2) +
      geom_line(group = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      labs(
        title = paste("Propensity calibration -", ci_label),
        x = "Mean predicted propensity",
        y = "Observed treatment rate"
      ) +
      theme_minimal()
    ggsave(filename = file.path(out_dir, "propensity_calibration.png"), plot = p2, width = 7, height = 5, dpi = 300)
  }

  if (nrow(balance_dt) > 0) {
    bal_plot <- copy(balance_dt)
    bal_plot[, label := ifelse(is.na(level), covariate, paste0(covariate, " = ", level))]
    bal_plot[, worst := pmax(abs_smd_unweighted, abs_smd_weighted, na.rm = TRUE)]
    bal_plot <- bal_plot[order(-worst)]
    bal_plot <- bal_plot[seq_len(min(nrow(bal_plot), 25))]
    bal_long <- melt(
      bal_plot[, .(label, smd_unweighted, smd_weighted)],
      id.vars = "label",
      variable.name = "stage",
      value.name = "smd"
    )
    bal_long[, stage := fifelse(stage == "smd_unweighted", "Unweighted", "IPW weighted")]
    p3 <- ggplot(bal_long, aes(x = smd, y = reorder(label, abs(smd)), color = stage)) +
      geom_point(size = 2) +
      geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray60") +
      labs(
        title = paste("Covariate balance -", ci_label),
        x = "Standardized mean difference",
        y = NULL,
        color = NULL
      ) +
      theme_minimal()
    ggsave(filename = file.path(out_dir, "covariate_balance.png"), plot = p3, width = 9, height = 7, dpi = 300)
  }
}

############################################################
## 2. Main estimation function
############################################################
#' Estimate Subgroup Effects with Cross-Fitted DML and H2O Ensembles
#'
#' Performs cross-fitted doubly robust estimation using stacked H2O binary
#' models for the propensity score and outcome regressions, with optional
#' targeting and subgroup summaries.
#'
#' @param dat Input analysis data.
#' @param y_col Outcome column name.
#' @param d_col Treatment column name.
#' @param x_cols Character vector of adjustment covariates.
#' @param s_cols Character vector of subgroup variables.
#' @param id_col Optional identifier column name.
#' @param n_folds Number of cross-fitting folds.
#' @param fold_id Optional predefined fold indices.
#' @param trim Lower and upper truncation bound for probabilities.
#' @param seed Random seed.
#' @param ntrees Number of trees for tree-based learners.
#' @param max_depth Maximum tree depth for tree-based learners.
#' @param nfolds_ensemble Number of folds used inside the H2O ensemble fits.
#' @param ensemble_models Character vector of base learner types.
#' @param metalearner_algorithm H2O metalearner algorithm.
#' @param metalearner_params Named list of metalearner settings.
#' @param balance_classes Logical flag indicating whether to balance classes in
#'   H2O learners.
#' @param do_targeting Logical flag indicating whether to perform targeted
#'   updating.
#' @param lambda_fluctuate Ridge penalty used in the fluctuation fallback.
#' @param h2o_max_mem_size Maximum H2O memory allocation.
#' @param return_fold_level Logical flag indicating whether to return fold-level
#'   records.
#' @param cleanup_h2o_each_fold Logical flag indicating whether to clear the H2O
#'   cluster after each fold.
#'
#' @return A list containing subgroup summaries, individual scores, and
#'   optionally fold-level results.
estimate_effects_by_subgroup_dr2_h2o_ensemble <- function(dat,
                                                          y_col,
                                                          d_col,
                                                          x_cols,
                                                          s_cols,
                                                          id_col = NULL,
                                                          n_folds = 5,
                                                          fold_id = NULL,
                                                          trim = 0.01,
                                                          seed = 123,
                                                          ntrees = 200,
                                                          max_depth = 6,
                                                          nfolds_ensemble = 3,
                                                          ensemble_models = c("gbm", "drf", "glm"), #c("gbm", "drf", "dl", "glm")
                                                          metalearner_algorithm = "glm",
                                                          metalearner_params = list(family = "binomial"),
                                                          balance_classes = TRUE,
                                                          do_targeting = TRUE,
                                                          lambda_fluctuate = 1e-4,
                                                          h2o_max_mem_size = "100G",
                                                          return_fold_level = TRUE,
                                                          cleanup_h2o_each_fold = TRUE) {
  set.seed(seed)
  # init_h2o_safe(max_mem_size = h2o_max_mem_size, nthreads = -1)

  dat.model <- data.table::copy(setDT(dat))

  if (!all(c(y_col, d_col, x_cols, s_cols) %in% names(dat.model))) {
    stop("Some specified columns are not in dat.model.")
  }

  if (!is.null(id_col) && !(id_col %in% names(dat.model))) {
    stop("id_col is not in dat.model.")
  }

  dat.model[[y_col]] <- .safe_binary_numeric(dat.model[[y_col]])
  dat.model[[d_col]] <- .safe_binary_numeric(dat.model[[d_col]])

  for (s in s_cols) {
    if (is.character(dat.model[[s]])) dat.model[[s]] <- factor(dat.model[[s]])
  }
  for (x in x_cols) {
    if (is.character(dat.model[[x]])) dat.model[[x]] <- factor(dat.model[[x]])
  }

  n <- nrow(dat.model)
  if (n == 0) stop("No rows in input data.")

  if (is.null(fold_id)) {
    fold_vec <- sample(rep(seq_len(n_folds), length.out = n))
    fold_id <- split(seq_len(n), fold_vec)
  } else {
    n_folds <- length(fold_id)
  }

  g1_hat <- rep(NA_real_, n)
  g0_hat <- rep(NA_real_, n)
  m_hat  <- rep(NA_real_, n)

  phi1   <- rep(NA_real_, n)
  phi0   <- rep(NA_real_, n)
  phi_ar <- rep(NA_real_, n)

  c_g1_hat <- rep(NA_real_, n)
  c_g0_hat <- rep(NA_real_, n)
  c_phi1 <- rep(NA_real_, n)
  c_phi0 <- rep(NA_real_, n)
  c_phi_ar <- rep(NA_real_, n)

  yhat_obs <- rep(NA_real_, n)
  D_obs    <- rep(NA_real_, n)
  Y_obs    <- rep(NA_real_, n)
  fold_used <- rep(NA_integer_, n)

  fold_records <- vector("list", n_folds)

  for (fold in seq_len(n_folds)) {
    init_h2o_safe(max_mem_size = h2o_max_mem_size, nthreads = 38)
    cat("Cross-fitting fold =", fold, "/", n_folds, "\n")

    test_idx  <- fold_id[[fold]]
    train_idx <- setdiff(seq_len(n), test_idx)

    train_dat <- dat.model[train_idx]
    test_dat  <- dat.model[test_idx]
    train_ps  <- train_dat[!is.na(get(d_col))]

    cat("  Training propensity model...\n")
    ps_fit <- fit_h2o_stacked_binary(
      train_df = as.data.frame(train_ps[, c(d_col, x_cols), with = FALSE]),
      test_df = as.data.frame(test_dat[, x_cols, with = FALSE]),
      y_col = d_col,
      x_cols = x_cols,
      ensemble_models = ensemble_models,
      metalearner_algorithm = metalearner_algorithm,
      metalearner_params = metalearner_params,
      ntrees = ntrees,
      max_depth = max_depth,
      nfolds = nfolds_ensemble,
      seed = seed + 1000 * fold,
      balance_classes = balance_classes,
      trim = trim
    )
    e_hat <- pmin(pmax(ps_fit$pred, trim), 1 - trim)

    cat("  Training outcome model for treated...\n")
    train_treat <- train_dat[get(d_col) == 1 & !is.na(get(y_col))]
    g1_fit <- fit_h2o_stacked_binary(
      train_df = as.data.frame(train_treat[, c(y_col, x_cols), with = FALSE]),
      test_df = as.data.frame(test_dat[, x_cols, with = FALSE]),
      y_col = y_col,
      x_cols = x_cols,
      ensemble_models = ensemble_models,
      metalearner_algorithm = metalearner_algorithm,
      metalearner_params = metalearner_params,
      ntrees = ntrees,
      max_depth = max_depth,
      nfolds = nfolds_ensemble,
      seed = seed + 2000 * fold,
      balance_classes = balance_classes,
      trim = trim
    )
    m1_hat <- pmin(pmax(g1_fit$pred, trim), 1 - trim)

    cat("  Training outcome model for control...\n")
    train_ctrl <- train_dat[get(d_col) == 0 & !is.na(get(y_col))]
    g0_fit <- fit_h2o_stacked_binary(
      train_df = as.data.frame(train_ctrl[, c(y_col, x_cols), with = FALSE]),
      test_df = as.data.frame(test_dat[, x_cols, with = FALSE]),
      y_col = y_col,
      x_cols = x_cols,
      ensemble_models = ensemble_models,
      metalearner_algorithm = metalearner_algorithm,
      metalearner_params = metalearner_params,
      ntrees = ntrees,
      max_depth = max_depth,
      nfolds = nfolds_ensemble,
      seed = seed + 3000 * fold,
      balance_classes = balance_classes,
      trim = trim
    )
    m0_hat <- pmin(pmax(g0_fit$pred, trim), 1 - trim)

    Y_test <- dat.model[[y_col]][test_idx]
    D_test <- dat.model[[d_col]][test_idx]
    valid_test <- is.finite(Y_test) & is.finite(D_test)

    phi1_fold <- rep(NA_real_, length(test_idx))
    phi0_fold <- rep(NA_real_, length(test_idx))
    phi1_fold[valid_test] <- m1_hat[valid_test] +
      D_test[valid_test] * (Y_test[valid_test] - m1_hat[valid_test]) / e_hat[valid_test]
    phi0_fold[valid_test] <- m0_hat[valid_test] +
      (1 - D_test[valid_test]) * (Y_test[valid_test] - m0_hat[valid_test]) / (1 - e_hat[valid_test])
    phi_ar_fold <- phi1_fold - phi0_fold

    g1_hat[test_idx] <- m1_hat
    g0_hat[test_idx] <- m0_hat
    m_hat[test_idx]  <- e_hat
    phi1[test_idx]   <- phi1_fold
    phi0[test_idx]   <- phi0_fold
    phi_ar[test_idx] <- phi_ar_fold

    yhat_obs[test_idx] <- ifelse(D_test == 1, m1_hat, m0_hat)
    D_obs[test_idx] <- D_test
    Y_obs[test_idx] <- Y_test
    fold_used[test_idx] <- fold

    if (do_targeting) {
      c_phi1_fold <- rep(NA_real_, length(test_idx))
      c_phi0_fold <- rep(NA_real_, length(test_idx))
      m1_tilde_all <- rep(NA_real_, length(test_idx))
      m0_tilde_all <- rep(NA_real_, length(test_idx))

      if (sum(valid_test) > 1) {
        Y_valid <- Y_test[valid_test]
        D_valid <- D_test[valid_test]
        e_valid <- e_hat[valid_test]
        m1_valid <- m1_hat[valid_test]
        m0_valid <- m0_hat[valid_test]
        H1_valid <- D_valid / e_valid
        H0_valid <- (1 - D_valid) / (1 - e_valid)

        res1 <- fluctuate(
          y = Y_valid,
          H = H1_valid,
          m_hat = m1_valid,
          lambda = lambda_fluctuate
        )
        res0 <- fluctuate(
          y = Y_valid,
          H = H0_valid,
          m_hat = m0_valid,
          lambda = lambda_fluctuate
        )

        m1_tilde_valid <- pmin(pmax(res1$m_tilde, trim), 1 - trim)
        m0_tilde_valid <- pmin(pmax(res0$m_tilde, trim), 1 - trim)

        c_phi1_fold[valid_test] <- m1_tilde_valid +
          D_valid * (Y_valid - m1_tilde_valid) / e_valid
        c_phi0_fold[valid_test] <- m0_tilde_valid +
          (1 - D_valid) * (Y_valid - m0_tilde_valid) / (1 - e_valid)
        m1_tilde_all[valid_test] <- m1_tilde_valid
        m0_tilde_all[valid_test] <- m0_tilde_valid
      }

      c_g1_hat[test_idx] <- m1_tilde_all
      c_g0_hat[test_idx] <- m0_tilde_all
      c_phi1[test_idx] <- c_phi1_fold
      c_phi0[test_idx] <- c_phi0_fold
      c_phi_ar[test_idx] <- c_phi1_fold - c_phi0_fold
    }

    fold_dt <- data.table(
      fold = fold,
      row_index = test_idx,
      g1_hat = m1_hat,
      g0_hat = m0_hat,
      m_hat  = e_hat,
      phi1   = phi1_fold,
      phi0   = phi0_fold,
      phi_ar = phi_ar_fold,
      y_obs  = Y_test,
      d_obs  = D_test,
      y_hat_obs = ifelse(D_test == 1, m1_hat, m0_hat)
    )

    if (do_targeting) {
      fold_dt[, `:=`(
        c_g1_hat = c_g1_hat[test_idx],
        c_g0_hat = c_g0_hat[test_idx],
        c_phi1 = c_phi1[test_idx],
        c_phi0 = c_phi0[test_idx],
        c_phi_ar = c_phi_ar[test_idx]
      )]
    }

    if (!is.null(id_col)) {
      fold_dt[, (id_col) := dat.model[[id_col]][test_idx]]
    }

    for (ss in s_cols) {
      fold_dt[, (ss) := dat.model[[ss]][test_idx]]
    }

    fold_records[[fold]] <- fold_dt

    if (cleanup_h2o_each_fold) {
      h2o.removeAll()
      gc()
    }
  }

  dat_out <- data.table::copy(data.table::as.data.table(dat))
  dat_out[, `:=`(
    g1_hat = g1_hat,
    g0_hat = g0_hat,
    m_hat  = m_hat,
    phi1   = phi1,
    phi0   = phi0,
    phi_ar = phi_ar,
    fold   = fold_used
  )]

  if (do_targeting) {
    dat_out[, `:=`(
      c_g1_hat = c_g1_hat,
      c_g0_hat = c_g0_hat,
      c_phi1 = c_phi1,
      c_phi0 = c_phi0,
      c_phi_ar = c_phi_ar
    )]
  }


  summarize_phi_group <- function(phi1, phi0, phi_ar, prefix = "") {
    ok <- is.finite(phi1) & is.finite(phi0) & is.finite(phi_ar)
    n_eff <- sum(ok)

    if (n_eff == 0) {
      out <- list(
        n = 0L,
        risk1 = NA_real_,
        SE_risk1 = NA_real_,
        risk1_CI_low = NA_real_,
        risk1_CI_high = NA_real_,
        risk0 = NA_real_,
        SE_risk0 = NA_real_,
        risk0_CI_low = NA_real_,
        risk0_CI_high = NA_real_,
        AR = NA_real_,
        SE_AR = NA_real_,
        AR_CI_low = NA_real_,
        AR_CI_high = NA_real_,
        RR = NA_real_,
        SE_log_RR = NA_real_,
        RR_CI_low = NA_real_,
        RR_CI_high = NA_real_
      )
    } else {
      p1 <- phi1[ok]
      p0 <- phi0[ok]
      pa <- phi_ar[ok]
      mu1 <- mean(p1)
      mu0 <- mean(p0)
      ar <- mean(pa)
      se_risk1 <- if (n_eff > 1) sd(p1) / sqrt(n_eff) else NA_real_
      se_risk0 <- if (n_eff > 1) sd(p0) / sqrt(n_eff) else NA_real_
      se_ar <- if (n_eff > 1) sd(pa) / sqrt(n_eff) else NA_real_
      rr <- if (mu0 <= 0) NA_real_ else mu1 / mu0

      if (is.na(mu1) || is.na(mu0) || mu1 <= 0 || mu0 <= 0 || n_eff <= 1) {
        se_log_rr <- NA_real_
        rr_ci_low <- NA_real_
        rr_ci_high <- NA_real_
      } else {
        eif_log_rr <- (p1 - mu1) / mu1 - (p0 - mu0) / mu0
        se_log_rr <- sd(eif_log_rr) / sqrt(n_eff)
        rr_ci_low <- exp(log(rr) - 1.96 * se_log_rr)
        rr_ci_high <- exp(log(rr) + 1.96 * se_log_rr)
      }

      out <- list(
        n = n_eff,
        risk1 = mu1,
        SE_risk1 = se_risk1,
        risk1_CI_low = ifelse(is.na(se_risk1), NA_real_, mu1 - 1.96 * se_risk1),
        risk1_CI_high = ifelse(is.na(se_risk1), NA_real_, mu1 + 1.96 * se_risk1),
        risk0 = mu0,
        SE_risk0 = se_risk0,
        risk0_CI_low = ifelse(is.na(se_risk0), NA_real_, mu0 - 1.96 * se_risk0),
        risk0_CI_high = ifelse(is.na(se_risk0), NA_real_, mu0 + 1.96 * se_risk0),
        AR = ar,
        SE_AR = se_ar,
        AR_CI_low = ifelse(is.na(se_ar), NA_real_, ar - 1.96 * se_ar),
        AR_CI_high = ifelse(is.na(se_ar), NA_real_, ar + 1.96 * se_ar),
        RR = rr,
        SE_log_RR = se_log_rr,
        RR_CI_low = rr_ci_low,
        RR_CI_high = rr_ci_high
      )
    }

    if (prefix == "") return(out)
    names(out) <- ifelse(names(out) == "n", "n", paste0(prefix, names(out)))
    out
  }

  has_target_cols <- do_targeting && all(c("c_phi1", "c_phi0", "c_phi_ar") %in% colnames(dat_out))

  subgroup_list <- lapply(s_cols, function(svar) {
    res <- dat_out[, {
      out <- summarize_phi_group(phi1, phi0, phi_ar, prefix = "")
      if (has_target_cols) {
        c_out <- summarize_phi_group(c_phi1, c_phi0, c_phi_ar, prefix = "c_")
        c_out$n <- NULL
        out <- c(out, c_out)
      }
      out
    }, by = svar]

    res[, subgroup_var := svar]
    setnames(res, svar, "subgroup_level")

    return(res)
  })

  subgroup_results <- rbindlist(subgroup_list, fill = TRUE)
  fold_level <- if (return_fold_level) rbindlist(fold_records, fill = TRUE) else NULL




  list(
    subgroup_results = subgroup_results,
    individual_scores = dat_out,
    fold_level = fold_level
  )
}










#' Select a Low-Correlation Subset of Variables
#'
#' Greedily selects variables from a correlation matrix so that no retained pair
#' exceeds the chosen absolute-correlation threshold. When several variables are
#' highly correlated, the function keeps the one with the smallest average
#' absolute correlation to the remaining variables.
#'
#' @param cor_mat A square correlation matrix with named rows and columns.
#' @param threshold Absolute-correlation threshold used to define high
#'   correlation.
#'
#' @return A character vector of selected variable names.
select_uncorrelated_vars <- function(cor_mat, threshold = 0.5) {
  # Get all variable names.
  vars <- colnames(cor_mat)

  # Build an adjacency matrix for pairs above the threshold.
  adj_mat <- abs(cor_mat) >= threshold
  diag(adj_mat) <- FALSE  # Remove self-links on the diagonal.

  selected <- c()
  remaining <- vars

  while (length(remaining) > 0) {
    # Start from one of the remaining variables.
    v <- remaining[1]

    # Compute the average absolute correlation for the current variable.
    ind <- which(vars %in% v)
    cor.v <- mean(abs(cor_mat[ind, -ind]))  # Exclude the variable itself.

    # Find variables that are highly correlated with the current one.
    high_corr_vars <- vars[adj_mat[v, ]]

    if (length(high_corr_vars) > 0) {
      # Keep only high-correlation variables that are still remaining.
      high_corr_vars <- intersect(high_corr_vars, remaining)

      if (length(high_corr_vars) > 0) {
        cor.h <- sapply(high_corr_vars, function(hv) {
          ind.h <- which(vars %in% hv)
          mean(abs(cor_mat[ind.h, -ind.h]))
        })

        # Replace `v` when another candidate has lower average correlation.
        if (min(cor.h) < cor.v) {
          v <- high_corr_vars[which.min(cor.h)]
        }
      }
    }

    # Retain the chosen variable.
    selected <- c(selected, v)
    # Remove the chosen variable and its high-correlation neighbors.
    remaining <- setdiff(remaining, c(v, high_corr_vars))
  }

  selected
}
