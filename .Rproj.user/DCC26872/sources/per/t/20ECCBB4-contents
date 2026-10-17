suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(maps)
  library(RColorBrewer)
})

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

selected.Comorbidity <- comorbidity_order

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

K <- 4
`%nin%` <- Negate(`%in%`)
table_dir <- file.path("Results", "table")
before_root_miss <- file.path("results", "Cohort_IV_non_imputed_exc_deceased")
combined_root_miss <- file.path("results", "Cohort_II_non_imputed_exc_deceased")
before_root_died <- file.path("results", "Cohort_IV_imputed_inc_deceased")
combined_root_died <- file.path("results", "Cohort_II_imputed_inc_deceased")

read_overall_summary <- function(root.path, scenario) {
  ## This reproduces the user's Summary object for one scenario:
  ## Summary <- rbind(Summary, overall.<Comorbidity> parquet)
  Summary <- lapply(selected.Comorbidity, function(ci) {
    path <- file.path(
      root.path,
      paste0("overall.", ci, ".Imputed_non_died_select_covariables_no_BMI_", K, ".parquet")
    )
    if (!file.exists(path)) return(NULL)
    dat <- as.data.frame(read_parquet(path))
    dat$Comorbidity <- ci
    dat$Scenario <- scenario
    dat
  })
  Summary <- bind_rows(Summary)
  if (nrow(Summary) == 0) {
    stop("No overall parquet files found under ", normalizePath(root.path, mustWork = FALSE), call. = FALSE)
  }
  Summary
}

require_targeted_columns <- function(dat) {
  required <- c(
    "c_AR", "c_AR_CI_low", "c_AR_CI_high",
    "c_RR", "c_RR_CI_low", "c_RR_CI_high",
    "Comorbidity"
  )
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0) {
    stop(
      "Missing required targeted DML columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
  dat
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_est_ci <- function(est, low, high, digits = 3, scale = 1, suffix = "") {
  est <- est * scale
  low <- low * scale
  high <- high * scale
  ifelse(
    is.na(est) | is.na(low) | is.na(high),
    "",
    paste0(
      fmt_num(est, digits), suffix,
      " (", fmt_num(low, digits), suffix, ", ", fmt_num(high, digits), suffix, ")"
    )
  )
}

format_sensitivity_summary <- function(dat, source_name) {
  dat %>%
    require_targeted_columns() %>%
    mutate(
      Comorbidity = factor(Comorbidity, levels = comorbidity_order),
      c_AR_fmt = fmt_est_ci(
        c_AR, c_AR_CI_low, c_AR_CI_high,
        digits = 3, scale = 100
      ),
      c_RR_fmt = fmt_est_ci(
        c_RR, c_RR_CI_low, c_RR_CI_high,
        digits = 3, scale = 1
      )
    ) %>%
    filter(!is.na(Comorbidity)) %>%
    transmute(
      Comorbidity,
      !!paste0("c_AR (%) of ", source_name, " (95% CI)") := c_AR_fmt,
      !!paste0("c_RR of ", source_name, " (95% CI)") := c_RR_fmt
    )
}

build_sensitivity_table <- function(
    Summary.Before_miss,
    Summary.Combined_miss,
    Summary.Before_died,
    Summary.Combined_died
) {
  table_order <- tibble(
    Comorbidity = factor(comorbidity_order, levels = comorbidity_order),
    order_id = seq_along(comorbidity_order),
    comorbidity = unname(comorbidity_labels[comorbidity_order])
  )

  table_order %>%
    left_join(format_sensitivity_summary(Summary.Before_miss, "Summary.Before_miss"), by = "Comorbidity") %>%
    left_join(format_sensitivity_summary(Summary.Combined_miss, "Summary.Combined_miss"), by = "Comorbidity") %>%
    left_join(format_sensitivity_summary(Summary.Before_died, "Summary.Before_died"), by = "Comorbidity") %>%
    left_join(format_sensitivity_summary(Summary.Combined_died, "Summary.Combined_died"), by = "Comorbidity") %>%
    arrange(order_id) %>%
    select(
      comorbidity,
      `c_AR (%) of Summary.Before_miss (95% CI)`,
      `c_AR (%) of Summary.Combined_miss (95% CI)`,
      `c_AR (%) of Summary.Before_died (95% CI)`,
      `c_AR (%) of Summary.Combined_died (95% CI)`,
      `c_RR of Summary.Before_miss (95% CI)`,
      `c_RR of Summary.Combined_miss (95% CI)`,
      `c_RR of Summary.Before_died (95% CI)`,
      `c_RR of Summary.Combined_died (95% CI)`
    )
}

main <- function() {
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

  Summary.Before_miss <- read_overall_summary(before_root_miss, "Before")
  # colnames(Summary.Before_miss)
  # [1] "n"               "risk1"           "SE_risk1"        "risk1_CI_low"
  # [5] "risk1_CI_high"   "risk0"           "SE_risk0"        "risk0_CI_low"
  # [9] "risk0_CI_high"   "AR"              "SE_AR"           "AR_CI_low"
  # [13] "AR_CI_high"      "RR"              "SE_log_RR"       "RR_CI_low"
  # [17] "RR_CI_high"      "c_risk1"         "c_SE_risk1"      "c_risk1_CI_low"
  # [21] "c_risk1_CI_high" "c_risk0"         "c_SE_risk0"      "c_risk0_CI_low"
  # [25] "c_risk0_CI_high" "c_AR"            "c_SE_AR"         "c_AR_CI_low"
  # [29] "c_AR_CI_high"    "c_RR"            "c_SE_log_RR"     "c_RR_CI_low"
  # [33] "c_RR_CI_high"    "Comorbidity"     "Scenario"

  Summary.Combined_miss <- read_overall_summary(combined_root_miss, "Combined")
  Summary.Before_died <- read_overall_summary(before_root_died, "Before")
  Summary.Combined_died <- read_overall_summary(combined_root_died, "Combined")

  sensitivity_table <- build_sensitivity_table(
    Summary.Before_miss,
    Summary.Combined_miss,
    Summary.Before_died,
    Summary.Combined_died
  )

  write.csv(
    sensitivity_table,
    file.path(table_dir, "Sensitivity_table.csv"),
    row.names = FALSE
  )
  # saveRDS(sensitivity_table, file.path(table_dir, "Sensitivity_table.rds"))

  sensitivity_table
}

if (!identical(Sys.getenv("SENSITIVITY_SKIP_MAIN"), "true")) {
  sensitivity_table <- main()
  print(sensitivity_table, n = Inf)
}
