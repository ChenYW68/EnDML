# rm(list = ls())
############################################################
## 0. Packages
############################################################
library(dplyr)
library(data.table)
library(h2o)
library(glmnet)
# library(arrow)
# library(caret)
library(dplyr)
library(data.table)
library(ROSE)
library(h2o)
library(Hmisc)
library(ggplot2)
library(arrow)
# Cohort definitions and date utilities synchronized with LC_risk_ensemble_DML.R
`%nin%` <- Negate(`%in%`)
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
    list(eligible = as.integer(eligible_cohort_i), D = d_cohort_i)
  } else if (before == 2) {
    list(eligible = as.integer(eligible_cohort_ii), D = d_cohort_ii)
  } else if (before == 3) {
    list(eligible = as.integer(eligible_cohort_iii), D = d_cohort_iii)
  } else if (before == 4) {
    list(eligible = as.integer(eligible_cohort_iv), D = d_cohort_iv)
  } else {
    stop("before must be 1 (Cohort I), 2 (Cohort II), 3 (Cohort III), or 4 (Cohort IV).")
  }
}

get_other_comorbidity_cols <- function(df_colnames, ci, indicator_comorbidity) {
  exclude_cols <- grouped_comorbidity_components(ci)
  which(df_colnames %in% indicator_comorbidity[indicator_comorbidity %nin% exclude_cols])
}

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

earliest_recorded_date <- function(df, cols) {
  cols <- cols[cols %in% colnames(df)]
  if (length(cols) == 0) return(rep(as.Date(NA), nrow(df)))
  dates <- lapply(cols, function(nm) sanitize_date(df[[nm]]))
  out <- do.call(pmin, c(dates, list(na.rm = TRUE)))
  out[!is.finite(out)] <- as.Date(NA)
  out
}

Long_COVID_df <- read_parquet("./data/Long_COVID_Secondary_Data_all_2026_05.parquet")


# Initialize column
Long_COVID_df[, Dead := 0L]

# Compute condition (vectorized, no copies)
Long_COVID_df[
  !is.na(death_date) &
    is.na(Long_COVID) &
    !is.na(COVID_19) &
    (death_date - COVID_19 <= 180),
  Dead := 1L
]
# Long_COVID_df <- data.table::copy(Long_COVID_Secondary_Data_all)

# idx <- !is.na(Long_COVID_Secondary_Data_all$death_date) &
#   is.na(Long_COVID_Secondary_Data_all$Long_COVID) &
#   !is.na(Long_COVID_Secondary_Data_all$COVID_19) &
#   (as.Date(Long_COVID_Secondary_Data_all$death_date) -
#      Long_COVID_Secondary_Data_all$COVID_19 <= 90)
#
# Long_COVID_Secondary_Data_all$Dead[idx] <- 1L
#
# Long_COVID_df <- Long_COVID_Secondary_Data_all %>%
#   filter(age >= 18 | is.na(age))

# Selected.variable <- c(
#   "age",
#   "sex",
#   "race",
#   "BMI",
#   "total_number_of_COVID_vaccine_doses",
#   "poverty_status",
#   "health_insurance",
#   colnames(Long_COVID_df)[c(79:102, 105:113)]
# )
#
# selected.ind <- which(colnames(Long_COVID_df) %in% Selected.variable)


# # #--------------------------------------------------
# setDF(Long_COVID_df)
# temp <- Long_COVID_df[, selected.ind]
#
#
# # library(DataExplorer)
# # DataExplorer::plot_missing(temp)
#
# #Missing data
# library(mice)
#
# # Suppose df is your dataframe with missing values
# # Perform multiple imputation by chained equations
# imp <- mice(temp, m = 5, method = "pmm", seed = 0)
#
# # # Extract the completed dataset
# # df_imputed <- complete(imp, 1)
#
# long_data <- complete(imp, action = "long", include = TRUE)
#
# # long 格式里 .imp 表示第几个插补
# library(dplyr)
# avg_data <- long_data %>%
#   filter(.imp != 0) %>%  # 去掉原始数据
#   group_by(.id) %>%
#   dplyr::summarise(across(where(is.numeric), mean), .groups = "drop")
#
# save(avg_data, file = "./data/avg_data_imputed_all.RData")


load("./data/avg_data_imputed_all.RData")

imputed.ind <- which(colnames(Long_COVID_df) %in% colnames(avg_data))

setDF(Long_COVID_df)
Long_COVID_df[, imputed.ind] <- avg_data[, -c(1, ncol(avg_data))]
Long_COVID_df <- Long_COVID_df %>% filter(Dead == 0)

if (any(grouped_comorbidity_components("Cancer") %in% colnames(Long_COVID_df))) {
  Long_COVID_df$Cancer <- earliest_recorded_date(
    Long_COVID_df,
    grouped_comorbidity_components("Cancer")
  )
}


rm(avg_data)

ind.x <- which(colnames(Long_COVID_df) %in% c("age",
                                              "BMI",
                                              "total_number_of_COVID_vaccine_doses",
                                              "poverty_status",
                                              "health_insurance",
                                              colnames(Long_COVID_df)[c(79:102, 105:113)]
))
setDF(Long_COVID_df)
cor_mat <- cor(Long_COVID_df[, c(ind.x)], use = "complete.obs")
select_uncorrelated_vars <- function(cor_mat, threshold = 0.5) {
  # 获取所有变量名
  vars <- colnames(cor_mat)

  # 创建邻接矩阵：相互相关性超过 threshold 的连边
  adj_mat <- abs(cor_mat) >= threshold
  diag(adj_mat) <- FALSE  # 去除对角线

  selected <- c()
  remaining <- vars

  while (length(remaining) > 0) {
    # 任选一个当前剩下的变量
    v <- remaining[1]

    # v 的平均相关度
    ind <- which(vars %in% v)
    cor.v <- mean(abs(cor_mat[ind, -ind]))  # 排除自身

    # 找到与 v 高相关的变量
    high_corr_vars <- vars[adj_mat[v, ]]

    if (length(high_corr_vars) > 0) {
      # 只考虑 remaining 中的相关变量
      high_corr_vars <- intersect(high_corr_vars, remaining)

      if (length(high_corr_vars) > 0) {
        cor.h <- sapply(high_corr_vars, function(hv) {
          ind.h <- which(vars %in% hv)
          mean(abs(cor_mat[ind.h, -ind.h]))
        })

        # 如果某个高相关变量比 v 更“重要”，替换
        if (min(cor.h) < cor.v) {
          v <- high_corr_vars[which.min(cor.h)]
        }
      }
    }

    # 保留 v
    selected <- c(selected, v)
    # 从剩余变量中去掉 v 及其高相关变量
    remaining <- setdiff(remaining, c(v, high_corr_vars))
  }

  selected
}

Selected.variable <- c("sex",
                       "race",
                       select_uncorrelated_vars(cor_mat, threshold = 0.7))




indicator.Comorbidity <- c(
  "Hypertension",
  "Diabetes_without_complications",
  "Depression",
  "Chronic_lung_disease",
  "Diabetes_with_complications",
  "Kidney_disease",  #6
  "Coronary_artery_disease",
  "Heart_failure",
  "Liver_disease",
  "Cerebrovascular_disease", #10
  "Peripheral_vascular_disease",
  "Myocardial_infarction",
  "Cardiomyopathy",
  "Dementia",
  "Malignant_tumor", #15
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
## 4. Run by comorbidity
##    Keep your data-preparation logic, remove old model code
############################################################
K <- 4
# init_h2o_safe(max_mem_size = "20G")
setDT(Long_COVID_df)
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



comorbidity_order <- c(
  "Any_comorbidities",
  "Hypertension",
  "Chronic_lung_disease",
  "Kidney_disease",
  "Liver_disease",
  "Neuro_psych_disorders",
  "Depression",
  "Cerebrovascular_disease",
  "Dementia",
  "Diabetes",
  "Diabetes_without_complications",
  "Diabetes_with_complications",
  "Cardiovascular",
  "Coronary_artery_disease",
  "Heart_failure",
  "Peripheral_vascular_disease",
  "Myocardial_infarction",
  "Cardiomyopathy",
  "Cancer",
  "Malignant_tumor",
  "Benign_tumor"
)

comorbidity_labels <- c(
  Any_comorbidities = "Any comorbidities",
  Hypertension = "Hypertension",
  Chronic_lung_disease = "Chronic lung disease",
  Kidney_disease = "Kidney disease",
  Liver_disease = "Liver disease",
  Neuro_psych_disorders = "Neuro psych disorders",
  Depression = "Depression",
  Cerebrovascular_disease = "Cerebrovascular disease",
  Dementia = "Dementia",
  Diabetes = "Diabetes",
  Diabetes_without_complications = "Diabetes without complications",
  Diabetes_with_complications = "Diabetes with complications",
  Cardiovascular = "Cardiovascular",
  Coronary_artery_disease = "Coronary artery disease",
  Heart_failure = "Heart failure",
  Peripheral_vascular_disease = "Peripheral vascular disease",
  Myocardial_infarction = "Myocardial infarction",
  Cardiomyopathy = "Cardiomyopathy",
  Cancer = "Cancer",
  Malignant_tumor = "Malignant tumor",
  Benign_tumor = "Benign tumor"
)

`%nin%` <- Negate(`%in%`)

build_empirical_esti <- function(before, data = Long_COVID_df, verbose = TRUE) {
  if (!before %in% 1:4 || length(before) != 1L || is.na(before)) {
    stop("before must be 1, 2, 3, or 4.", call. = FALSE)
  }
  needed <- c("COVID_19", "Long_COVID", "Long_COVID_binary")
  missing <- setdiff(needed, names(data))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

  results <- vector("list", length(comorbidity_order))
  for (j in seq_along(comorbidity_order)) {
    ci <- comorbidity_order[j]
    if (verbose) message("Cohort ", before, ": ", ci)
    if (!ci %in% names(data)) {
      warning("Missing exposure column: ", ci)
      next
    }
    df <- data.table::copy(data.table::as.data.table(data))
    df[, COVID_19 := sanitize_date(COVID_19)]
    df[, Long_COVID := sanitize_date(Long_COVID)]
    df[, (ci) := sanitize_date(get(ci))]
    df <- df[!is.na(COVID_19)]

    other_idx <- get_other_comorbidity_cols(names(df), ci, indicator.Comorbidity)
    sel <- select_records_from_full_scenario(
      before = before,
      covid_date = df$COVID_19,
      lc_date = df$Long_COVID,
      exposure_date = df[[ci]],
      has_other_comorbidity = has_any_recorded_comorbidity(df, other_idx)
    )
    eligible <- !is.na(sel$eligible) & sel$eligible == 1L
    df <- df[eligible]
    df[, D := sel$D[eligible]]
    df <- df[!is.na(D)]
    y <- if (is.factor(df$Long_COVID_binary)) {
      suppressWarnings(as.numeric(as.character(df$Long_COVID_binary)))
    } else {
      suppressWarnings(as.numeric(df$Long_COVID_binary))
    }
    valid_y <- !is.na(y) & is.finite(y) & y %in% c(0, 1)
    # Only observed binary outcomes contribute to empirical risks.
    df <- df[valid_y]
    y <- y[valid_y]
    n_total <- nrow(df)
    n_D1 <- sum(df$D == 1L)
    n_D0 <- sum(df$D == 0L)
    events_D1 <- sum(y[df$D == 1L])
    events_D0 <- sum(y[df$D == 0L])
    R1 <- if (n_D1 > 0) events_D1 / n_D1 else NA_real_
    R0 <- if (n_D0 > 0) events_D0 / n_D0 else NA_real_
    RR <- if (is.finite(R0) && R0 > 0) R1 / R0 else NA_real_
    AR <- R1 - R0
    results[[j]] <- data.frame(
      comorbidity = ci,
      comorbidity.label = unname(comorbidity_labels[ci]),
      cohort = paste0("Cohort_", c("I", "II", "III", "IV")[before]),
      n_total = n_total, n_D1 = n_D1, n_D0 = n_D0,
      events_D1 = events_D1, events_D0 = events_D0,
      exposure_prevalence = if (n_total > 0) n_D1 / n_total else NA_real_,
      event_rate_D1 = R1, event_rate_D0 = R0,
      R1 = R1, R0 = R0, RR = RR, AR = AR
    )
  }
  dplyr::bind_rows(results)
}

if (exists("Long_COVID_df") && exists("indicator.Comorbidity")) {
  empirical.by.cohort <- lapply(1:4, build_empirical_esti)
  names(empirical.by.cohort) <- paste0("Cohort_", c("I", "II", "III", "IV"))
  empirical.esti.before1 <- empirical.by.cohort[[1]]
  empirical.esti.before2 <- empirical.by.cohort[[2]]
  empirical.esti.before3 <- empirical.by.cohort[[3]]
  empirical.esti.before4 <- empirical.by.cohort[[4]]
  empirical.esti <- empirical.esti.before2
  for (nm in names(empirical.by.cohort)) {
    out_dir <- file.path("results", nm)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(empirical.by.cohort[[nm]],
                     file.path(out_dir, "empirical_estimates.csv"), row.names = FALSE)
  }
}

# Cohort-specific DML summaries and empirical-risk comparison (1--4).
read_overall_summary <- function(root.path, scenario) {
  rows <- lapply(comorbidity_order, function(ci) {
    path <- file.path(root.path, paste0(
      "overall.", ci, ".Imputed_non_died_select_covariables_no_BMI_", K, ".parquet"))
    if (!file.exists(path)) return(NULL)
    dat <- as.data.frame(arrow::read_parquet(path))
    dat$Comorbidity <- ci
    dat$Scenario <- scenario
    dat
  })
  out <- dplyr::bind_rows(rows)
  if (!nrow(out)) stop("No overall parquet files found in ", root.path, call. = FALSE)
  out
}

read_empirical_fallback <- function(cohort_number) {
  path <- file.path("results", paste0("Cohort_", c("I", "II", "III", "IV")[cohort_number]),
                    "empirical_estimates.csv")
  if (!file.exists(path)) return(NULL)
  utils::read.csv(path, check.names = FALSE)
}

require_columns <- function(dat, columns, object_name) {
  missing_columns <- setdiff(columns, names(dat))
  if (length(missing_columns)) {
    stop(object_name, " is missing required columns: ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
}

standardize_empirical_risks <- function(dat, suffix) {
  require_columns(dat, c("comorbidity", "R0", "R1"), paste0("empirical.esti.", suffix))
  out <- data.frame(comorbidity = dat$comorbidity,
                    R0 = sprintf("%.3f", dat$R0 * 100),
                    R1 = sprintf("%.3f", dat$R1 * 100), check.names = FALSE)
  names(out)[-1] <- paste0(c("R0", "R1"), " of empirical.esti using before = ", suffix)
  out
}

format_risk_ci <- function(est, low, high, digits = 3, scale = 100) {
  out <- rep(NA_character_, length(est))
  complete <- is.finite(est) & is.finite(low) & is.finite(high)
  out[complete] <- sprintf(paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
                           est[complete] * scale, low[complete] * scale, high[complete] * scale)
  out
}

standardize_summary_risks <- function(dat, suffix) {
  require_columns(dat, c("Comorbidity", "c_risk0", "c_risk0_CI_low", "c_risk0_CI_high",
                         "c_risk1", "c_risk1_CI_low", "c_risk1_CI_high"), paste0("Summary.", suffix))
  out <- data.frame(comorbidity = dat$Comorbidity,
                    R0 = format_risk_ci(dat$c_risk0, dat$c_risk0_CI_low, dat$c_risk0_CI_high),
                    R1 = format_risk_ci(dat$c_risk1, dat$c_risk1_CI_low, dat$c_risk1_CI_high),
                    check.names = FALSE)
  names(out)[-1] <- paste0(c("R0", "R1"), " of Summary.", suffix)
  out
}

# Accept named lists indexed by cohort number; preserve R0 columns before R1 columns.
build_r0_r1_table <- function(empirical_by_cohort, summary_by_cohort) {
  ids <- as.character(1:4)
  if (!all(ids %in% names(empirical_by_cohort)) ||
      !all(ids %in% names(summary_by_cohort))) {
    stop("Both input lists must be named '1', '2', '3', '4'.", call. = FALSE)
  }
  out <- data.frame(comorbidity = comorbidity_order,
                    comorbidity.label = unname(comorbidity_labels[comorbidity_order]),
                    check.names = FALSE)
  for (i in ids) {
    a <- standardize_empirical_risks(empirical_by_cohort[[i]], i)
    b <- standardize_summary_risks(summary_by_cohort[[i]], i)
    if (anyDuplicated(a$comorbidity) || anyDuplicated(b$comorbidity)) {
      stop("Duplicate comorbidity rows in Cohort ", i, call. = FALSE)
    }
    out <- merge(out, a, by = "comorbidity", all.x = TRUE, sort = FALSE)
    out <- merge(out, b, by = "comorbidity", all.x = TRUE, sort = FALSE)
  }
  out <- out[match(comorbidity_order, out$comorbidity), , drop = FALSE]
  rownames(out) <- NULL
  risk_columns <- unlist(lapply(c("R0", "R1"), function(risk) {
    unlist(lapply(ids, function(i) c(
      paste0(risk, " of empirical.esti using before = ", i),
      paste0(risk, " of Summary.", i))))
  }))
  out[, c("comorbidity.label", risk_columns), drop = FALSE]
}

# Do not use legacy 'Before'/'Combined' objects or CSVs: they can refer to
# different cohorts and silently contaminate a four-cohort comparison.
empirical_for_table <- setNames(lapply(1:4, function(i) {
  object_name <- paste0("empirical.esti.before", i)
  if (exists(object_name, inherits = FALSE)) get(object_name, inherits = FALSE)
  else read_empirical_fallback(i)
}), as.character(1:4))

summary_for_table <- setNames(lapply(1:4, function(i) {
  path <- file.path("results", paste0("Cohort_", c("I", "II", "III", "IV")[i]))
  if (!dir.exists(path)) return(NULL)
  tryCatch(read_overall_summary(path, as.character(i)), error = function(e) {
    message("Cohort ", i, ": ", conditionMessage(e))
    NULL
  })
}), as.character(1:4))

# Expose explicitly numbered summary objects for interactive use.
for (i in 1:4) {
  if (!is.null(summary_for_table[[i]])) {
    assign(paste0("Summary.", i), summary_for_table[[i]])
  }
}

if (all(vapply(empirical_for_table, function(x) !is.null(x), logical(1))) &&
    all(vapply(summary_for_table, function(x) !is.null(x), logical(1)))) {
  R0_R1_table <- build_r0_r1_table(empirical_for_table, summary_for_table)
  dir.create(file.path("Results", "table"), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(R0_R1_table,
                   file.path("Results", "table", "R0_R1_empirical_summary_comparison_1_4.csv"),
                   row.names = FALSE)
  print(R0_R1_table)
} else {
  message("Four-cohort R0/R1 table not written: missing empirical estimates or DML summaries for: ",
          paste(names(empirical_for_table)[vapply(empirical_for_table, is.null, logical(1))], collapse = ", "),
          "; missing summaries: ",
          paste(names(summary_for_table)[vapply(summary_for_table, is.null, logical(1))], collapse = ", "))
}
