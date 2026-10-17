############################################################
## Build Table 1 and Figures 1-4 from DML parquet results
## Table 1 and Figures 1-3 compare Cohorts I-IV.
## Figure 4 compares Cohorts I, III, and IV.
## Figures S4_1-S4_3 compare Cohorts I-IV.
##
## Inputs:
##   results/Cohort_I/overall.*.parquet and SUMMARY.*.parquet
##   results/Cohort_II/overall.*.parquet and SUMMARY.*.parquet
##   results/Cohort_III/overall.*.parquet and SUMMARY.*.parquet
##   results/Cohort_IV/overall.*.parquet and SUMMARY.*.parquet
##
## Outputs:
##   Results/table/Table1_DML_targeted.csv
##   Results/table/Table1_DML_targeted.rds
##   Results/figure/Fig1_multimorbidity_combined.png
##   Results/figure/FigS1.1_aRR_multimorbidity_cohorts.png
##   Results/figure/FigS1.2_aR0_multimorbidity_cohorts.png
##   Results/figure/FigS1.3_aR1_multimorbidity_cohorts.png
##   Results/figure/FigS1.1.0_aAR_multimorbidity_other_cohorts.png
##   Results/figure/FigS1.1.1_aRR_multimorbidity_other_cohorts.png
##   Results/figure/FigS1.2.1_aR0_multimorbidity_other_cohorts.png
##   Results/figure/FigS1.3.1_aR1_multimorbidity_other_cohorts.png
##   Results/figure/Fig2_subpopulations_combined.png
##   Results/figure/FigS2_subpopulations_other_combined.png
##   Results/figure/Fig3_ARs.png
##   Results/figure/FigS3_aRR.png
##   Results/figure/FigS3_aR0.png
##   Results/figure/FigS3_aR1.png
##   Results/figure/FigS3_other_aAR.png
##   Results/figure/FigS3_other_aRR.png
##   Results/figure/FigS3_other_aR0.png
##   Results/figure/FigS3_other_aR1.png
##   Results/figure/Fig4_covid_year_aAR_combined.png
##   Results/figure/FigS4_1_covid_year_aR0_combined.png
##   Results/figure/FigS4_2_covid_year_aR1_combined.png
##   Results/figure/FigS4_3_covid_year_aRR_combined.png
############################################################

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(maps)
  library(RColorBrewer)
})

K <- 4
`%nin%` <- Negate(`%in%`)
cohort_i_root <- file.path("results", "Cohort_IV")
cohort_ii_root <- file.path("results", "Cohort_II")
cohort_iii_root <- file.path("results", "Cohort_III")
cohort_iv_root <- file.path("results", "Cohort_I")

table_dir <- file.path("new_Results", "table")
figure_dir <- file.path("new_Results", "figure")

dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

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

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Missing required file: ", normalizePath(path, mustWork = FALSE), call. = FALSE)
  }
}

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

read_sub_summary <- function(root.path, scenario) {
  ## This reproduces the user's Sub.summary object for one scenario:
  ## Sub.summary <- rbind(Sub.summary, SUMMARY.<Comorbidity> parquet)
  Sub.summary <- lapply(selected.Comorbidity, function(ci) {
    path <- file.path(
      root.path,
      paste0("SUMMARY.", ci, ".Imputed_non_died_select_covariables_no_BMI_", K, ".parquet")
    )
    if (!file.exists(path)) return(NULL)
    dat <- as.data.frame(read_parquet(path))
    if (!"Comorbidity" %in% names(dat)) dat$Comorbidity <- ci
    dat$Scenario <- scenario
    dat
  })
  Sub.summary <- bind_rows(Sub.summary)
  if (nrow(Sub.summary) == 0) {
    stop("No subgroup parquet files found under ", normalizePath(root.path, mustWork = FALSE), call. = FALSE)
  }
  Sub.summary
}

require_targeted_columns <- function(dat) {
  required <- c(
    "c_risk0", "c_risk0_CI_low", "c_risk0_CI_high",
    "c_risk1", "c_risk1_CI_low", "c_risk1_CI_high",
    "c_AR", "c_AR_CI_low", "c_AR_CI_high",
    "c_RR", "c_RR_CI_low", "c_RR_CI_high"
  )
  missing_required <- setdiff(required, names(dat))
  if (length(missing_required) > 0) {
    stop(
      "Missing built-in targeted DML columns: ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
  dat
}

e_value_rr <- function(rr) {
  rr <- as.numeric(rr)
  out <- rep(NA_real_, length(rr))
  valid <- is.finite(rr) & rr > 0
  rr2 <- ifelse(valid & rr < 1, 1 / rr, rr)
  out[valid] <- rr2[valid] + sqrt(rr2[valid] * (rr2[valid] - 1))
  out
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(x, format = "f", digits = digits))
}

fmt_est_ci <- function(est, low, high, digits = 3, scale = 1) {
  est <- est * scale
  low <- low * scale
  high <- high * scale
  ifelse(
    is.na(est) | is.na(low) | is.na(high),
    "",
    paste0(fmt_num(est, digits), " (", fmt_num(low, digits), ", ", fmt_num(high, digits), ")")
  )
}

label_comorbidity <- function(x) {
  x <- as.character(x)
  out <- unname(comorbidity_labels[x])
  out[is.na(out)] <- str_to_sentence(gsub("_", " ", x[is.na(out)]))
  out
}

prepare_metric_long <- function(dat) {
  if (!"Scenario" %in% names(dat)) dat$Scenario <- NA_character_
  dat <- require_targeted_columns(dat)

  bind_rows(
    dat %>% transmute(Scenario, Comorbidity, subgroup_var, subgroup_level, Metric = "aR0 (%)",
                      estimate = c_risk0 * 100, low = c_risk0_CI_low * 100, high = c_risk0_CI_high * 100),
    dat %>% transmute(Scenario, Comorbidity, subgroup_var, subgroup_level, Metric = "aR1 (%)",
                      estimate = c_risk1 * 100, low = c_risk1_CI_low * 100, high = c_risk1_CI_high * 100),
    dat %>% transmute(Scenario, Comorbidity, subgroup_var, subgroup_level, Metric = "aAR (%)",
                      estimate = c_AR * 100, low = c_AR_CI_low * 100, high = c_AR_CI_high * 100),
    dat %>% transmute(Scenario, Comorbidity, subgroup_var, subgroup_level, Metric = "aRR",
                      estimate = c_RR, low = c_RR_CI_low, high = c_RR_CI_high)
  ) %>%
    mutate(
      Comorbidity = factor(Comorbidity, levels = comorbidity_order),
      Comorbidity_label = factor(label_comorbidity(Comorbidity), levels = label_comorbidity(comorbidity_order)),
      Metric = factor(Metric, levels = c("aR0 (%)", "aR1 (%)", "aAR (%)", "aRR"))
    )
}

build_table1 <- function(Summary) {
  ## T1 uses only overall Summary columns:
  ## c_risk0, c_risk1, c_AR, c_RR and their built-in CI columns.
  ## It does not use Sub.summary.
  table_order <- tibble(
    Comorbidity = factor(comorbidity_order, levels = comorbidity_order),
    order_id = seq_along(comorbidity_order),
    Preexisting_comorbidities = unname(comorbidity_labels[comorbidity_order])
  )

  Summary %>%
    require_targeted_columns() %>%
    mutate(
      Comorbidity = factor(Comorbidity, levels = comorbidity_order),
      aAR = fmt_est_ci(c_AR, c_AR_CI_low, c_AR_CI_high, digits = 3, scale = 100),
      aRR = fmt_est_ci(c_RR, c_RR_CI_low, c_RR_CI_high, digits = 3, scale = 1),
      E_value = fmt_num(e_value_rr(c_RR), digits = 3)
    ) %>%
    filter(!is.na(Comorbidity)) %>%
    select(Comorbidity, Scenario, aAR, aRR, E_value) %>%
    pivot_wider(
      names_from = Scenario,
      values_from = c(aAR, aRR, E_value),
      names_glue = "{.value}_{Scenario}"
    ) %>%
    right_join(table_order, by = "Comorbidity") %>%
    arrange(order_id) %>%
    transmute(
      `Preexisting comorbidities` = Preexisting_comorbidities,
      `aARs (%;95% CI) Cohort I` = aAR_CohortI,
      `aARs (%;95% CI) Cohort II` = aAR_CohortII,
      `aARs (%;95% CI) Cohort III` = aAR_CohortIII,
      `aARs (%;95% CI) Cohort IV` = aAR_CohortIV,
      `aRR (95% CI) Cohort I` = aRR_CohortI,
      `aRR (95% CI) Cohort II` = aRR_CohortII,
      `aRR (95% CI) Cohort III` = aRR_CohortIII,
      `aRR (95% CI) Cohort IV` = aRR_CohortIV,
      `E value Cohort I` = E_value_CohortI,
      `E value Cohort II` = E_value_CohortII,
      `E value Cohort III` = E_value_CohortIII,
      `E value Cohort IV` = E_value_CohortIV
    )
}

base_plot_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey92", color = "grey70"),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

fig1_main_comorbidities <- c(
  "Hypertension",
  "Chronic_lung_disease",
  "Kidney_disease",
  "Liver_disease",
  "Neuro_psych_disorders",
  "Diabetes",
  "Cardiovascular",
  "Cancer"
)

fig3_main_comorbidities <- c(
  "Hypertension",
  "Chronic_lung_disease",
  "Kidney_disease",
  "Liver_disease",
  "Neuro_psych_disorders",
  "Diabetes",
  "Cardiomyopathy",
  "Cancer"
)

fig3_four_row_panel_levels <- function() {
  block1 <- c(
    "Hypertension",
    "Chronic_lung_disease",
    "Kidney_disease",
    "Liver_disease"
  )
  block2 <- c(
    "Neuro_psych_disorders",
    "Diabetes",
    "Cardiomyopathy",
    "Cancer"
  )

  unlist(lapply(list(block1, block2), function(block) {
    unlist(lapply(unname(cohort_labels), function(cohort) {
      paste(cohort, label_comorbidity(block), sep = "\n")
    }), use.names = FALSE)
  }), use.names = FALSE)
}

cohort_labels <- c(
  CohortI = "Cohort I",
  CohortII = "Cohort II",
  CohortIII = "Cohort III",
  CohortIV = "Cohort IV"
)
cohort_colors <- c(
  "Cohort I" = "black",
  "Cohort II" = "green",
  "Cohort III" = "red",
  "Cohort IV" = "blue"
)
cohort_shapes <- c(
  "Cohort I" = 16,
  "Cohort II" = 18,
  "Cohort III" = 17,
  "Cohort IV" = 15
)
cohort_linetypes <- c(
  "Cohort I" = "solid",
  "Cohort II" = "44",
  "Cohort III" = "22",
  "Cohort IV" = "42"
)
fig4_cohort_names <- c("Cohort I", "Cohort III", "Cohort IV")

label_cohort <- function(x) {
  x <- as.character(x)
  out <- unname(cohort_labels[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

plot_multimorbidity_metric_by_cohort <- function(sub_summary,
                                                 metric,
                                                 selected_comorbidities = fig1_main_comorbidities,
                                                 include_selected = TRUE) {
  comorbidity_levels <- if (include_selected) {
    selected_comorbidities
  } else {
    setdiff(comorbidity_order, c(selected_comorbidities, "Any_comorbidities"))
  }

  plot_df <- prepare_metric_long(sub_summary) %>%
    filter(subgroup_var == "com.count", Metric == metric) %>%
    filter(as.character(Comorbidity) %in% comorbidity_levels) %>%
    mutate(
      Cohort = factor(label_cohort(Scenario), levels = unname(cohort_labels)),
      subgroup_level = case_when(
        as.character(subgroup_level) %in% c("0", "0.0") ~ "0",
        as.character(subgroup_level) %in% c("1", "1.0") ~ "1",
        as.character(subgroup_level) %in% c("2", "2.0", ">=2", ">= 2") ~ ">=2",
        TRUE ~ as.character(subgroup_level)
      ),
      subgroup_level = factor(subgroup_level, levels = c("0", "1", ">=2")),
      Comorbidity_label = factor(
        str_wrap(label_comorbidity(as.character(Comorbidity)), width = 13),
        levels = str_wrap(label_comorbidity(comorbidity_levels), width = 13)
      )
    ) %>%
    filter(!is.na(subgroup_level), !is.na(Cohort), !is.na(Comorbidity_label))

  cohort_position <- position_dodge(width = 0.45)

  # p <- ggplot(plot_df, aes(x = subgroup_level, y = estimate, color = Cohort, group = Cohort)) +
  #   geom_errorbar(aes(ymin = low, ymax = high), width = 0.5, position = cohort_position, na.rm = TRUE) +
  #   geom_point(aes(shape = Cohort), size = 2, position = cohort_position, na.rm = TRUE) +
  #   facet_grid(. ~ Comorbidity_label, scales = "free_y") +
  #   labs(x = "Number of other medical conditions before the targeted comorbidity", y = metric) +
  #   base_plot_theme(base_size = 12) +
  #   theme(
  #     axis.text.x = element_text(size = 20),
  #     strip.text.x = element_text(size = 16),
  #     strip.text.y = element_blank(),
  #     legend.title = element_blank(),
  #     plot.margin = margin(
  #       t = 5,
  #       r = 5,
  #       b = 7,
  #       l = 5
  #     )
  #   )


  p <- ggplot(plot_df, aes(x = subgroup_level, y = estimate, color = Cohort, group = Cohort)) +
    geom_errorbar(
      aes(ymin = low, ymax = high),
      width = 0.5,
      position = cohort_position,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_point(
      aes(shape = Cohort, fill = Cohort),
      size = 2.8,
      stroke = 0.7,
      position = cohort_position,
      na.rm = TRUE
    ) +
    facet_grid(. ~ Comorbidity_label, scales = "free_y") +
    scale_color_manual(values = cohort_colors) +
    scale_fill_manual(values = cohort_colors) +
    scale_shape_manual(values = cohort_shapes) +
    labs(
      x = "Number of other medical conditions before the targeted comorbidity",
      y = metric
    ) +
    base_plot_theme(base_size = 20) +
    guides(
      fill = "none",
      shape = "none",
      color = guide_legend(
        override.aes = list(
          shape = unname(cohort_shapes),
          fill = unname(cohort_colors),
          size = 4,
          stroke = 0.7
        )
      )
    ) +
    theme(
      legend.position = "top",

      legend.margin = margin(
        b = -10
      ),
      # legend.justification = c(0.5, 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 16),
      strip.text.x = element_text(size = 16),
      # strip.text.y = element_blank(),
      legend.title = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 7, l = 5)
    )
  if (metric == "aRR") {
    p <- p + geom_hline(yintercept = 1, linetype = 2, color = "grey55")
  }

  p
}

plot_fig1_multimorbidity <- function(sub_summary,
                                     selected_comorbidities = fig1_main_comorbidities,
                                     include_selected = TRUE) {
  plot_df <- prepare_metric_long(sub_summary) %>%
    filter(subgroup_var == "com.count") %>%
    filter(if (include_selected) {
      as.character(Comorbidity) %in% selected_comorbidities
    } else {
      !(as.character(Comorbidity) %in% selected_comorbidities)
    }) %>%
    mutate(
      subgroup_level = case_when(
        as.character(subgroup_level) %in% c("0", "0.0") ~ "0",
        as.character(subgroup_level) %in% c("1", "1.0") ~ "1",
        as.character(subgroup_level) %in% c("2", "2.0", ">=2", ">= 2") ~ ">=2",
        TRUE ~ as.character(subgroup_level)
      ),
      subgroup_level = factor(subgroup_level, levels = c("0", "1", ">=2")),
      Comorbidity_label = factor(
        str_wrap(as.character(Comorbidity_label), width = 13),
        levels = str_wrap(label_comorbidity(comorbidity_order), width = 13)
      )
    ) %>%
    filter(!is.na(subgroup_level))

  ggplot(plot_df, aes(x = subgroup_level, y = estimate)) +
    geom_hline(data = data.frame(Metric = factor("aRR", levels = levels(plot_df$Metric))),
               aes(yintercept = 1), linetype = 2, color = "grey55", inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.18, color = "grey25", na.rm = TRUE) +
    geom_point(size = 1.1, shape = 21, fill = "grey20", color = "black", na.rm = TRUE) +
    facet_grid(Metric ~ Comorbidity_label, scales = "free_y") +
    labs(x = "Number of other medical conditions before the targeted comorbidity", y = NULL) +
    base_plot_theme(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12),
      strip.text.x = element_text(size = 11),
      strip.text.y = element_text(size = 12),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 7,
        l = 5
      )
    )
}

radar_ci_panel <- function(df,
                           group_col = "Subgroup_level",
                           color_col = NULL,
                           linetype_col = NULL,
                           value_col = "estimate",
                           low_col = "low",
                           high_col = "high",
                           title = "",
                           metric_label = "",
                           panel_label = NULL,
                           var_order = label_comorbidity(comorbidity_order),
                           palette = NULL,
                           seg = 4) {

  df <- df %>%
    mutate(
      Comorbidity_label = as.character(Comorbidity_label),
      Comorbidity_label = factor(Comorbidity_label, levels = var_order)
    ) %>%
    filter(!is.na(Comorbidity_label))

  groups <- unique(as.character(df[[group_col]]))
  groups <- groups[!is.na(groups)]
  legend_labels <- groups
  if (length(groups) == 0) {
    plot.new()
    title(main = title, cex.main = 0.9)
    text(0.5, 0.5, "No data")
    return(invisible(NULL))
  }

  K <- length(var_order)
  thetas <- seq(0, 2 * pi, length.out = K + 1)[-(K + 1)]
  thetas <- rev(thetas + pi / 2)
  value_range <- c(df[[value_col]], df[[low_col]], df[[high_col]])
  max_val <- max(value_range[is.finite(value_range)], na.rm = TRUE)
  if (!is.finite(max_val) || max_val <= 0) max_val <- 1
  max_val <- max_val * 1.08
  group_meta <- tibble(group = groups)
  if (!is.null(color_col)) {
    color_meta <- df %>%
      transmute(
        group = as.character(.data[[group_col]]),
        color_key = as.character(.data[[color_col]])
      ) %>%
      distinct(group, .keep_all = TRUE)
    group_meta <- left_join(group_meta, color_meta, by = "group")
  }
  if (!is.null(linetype_col)) {
    linetype_meta <- df %>%
      transmute(
        group = as.character(.data[[group_col]]),
        linetype_key = as.character(.data[[linetype_col]])
      ) %>%
      distinct(group, .keep_all = TRUE)
    group_meta <- left_join(group_meta, linetype_meta, by = "group")
  }

  if (is.null(palette)) {
    cols <- grDevices::hcl.colors(length(groups), "Dark 3")
  } else if (!is.null(color_col) && !is.null(names(palette))) {
    cols <- unname(palette[group_meta$color_key])
  } else {
    cols <- rep_len(palette, length(groups))
  }
  cols[is.na(cols)] <- "black"

  if (!is.null(linetype_col)) {
    linetype_levels <- unique(group_meta$linetype_key)
    linetypes <- seq_along(linetype_levels)
    names(linetypes) <- linetype_levels
    group_linetypes <- unname(linetypes[group_meta$linetype_key])
  } else {
    group_linetypes <- rep(1, length(groups))
  }

  polar_xy <- function(vals) {
    vals <- pmax(vals, 0)
    r <- vals / max_val
    data.frame(x = r * cos(thetas), y = r * sin(thetas))
  }

  plot.new()
  plot.window(xlim = c(-1.1, 1.1), ylim = c(-1.5, 1.1), asp = 1)
  title_text <- if (is.null(panel_label)) title else paste(panel_label, title)
  text(-2.2, 1.55, labels = title_text, cex = 1.5, font = 2, adj = c(0, 1))

  for (s in seq_len(seg)) {
    r <- s / seg
    polygon(r * cos(c(thetas, thetas[1])),
            r * sin(c(thetas, thetas[1])),
            border = "grey82", col = NA, lwd = 0.6)
  }
  segments(0, 0, cos(thetas), sin(thetas), col = "grey86", lwd = 0.5)

  for (i in seq_along(groups)) {
    d <- df %>%
      filter(.data[[group_col]] == groups[i]) %>%
      arrange(Comorbidity_label)
    d <- right_join(
      tibble(Comorbidity_label = factor(var_order, levels = var_order)),
      d,
      by = "Comorbidity_label"
    )
    lo_vals <- d[[low_col]]
    hi_vals <- d[[high_col]]
    mean_vals <- d[[value_col]]
    lo_vals[!is.finite(lo_vals)] <- NA_real_
    hi_vals[!is.finite(hi_vals)] <- NA_real_
    mean_vals[!is.finite(mean_vals)] <- NA_real_
    if (any(is.na(mean_vals))) next

    lo_vals[is.na(lo_vals)] <- mean_vals[is.na(lo_vals)]
    hi_vals[is.na(hi_vals)] <- mean_vals[is.na(hi_vals)]
    xy_hi <- polar_xy(hi_vals)
    xy_lo <- polar_xy(lo_vals)
    polygon(c(xy_hi$x, xy_hi$x[1], rev(c(xy_lo$x, xy_lo$x[1]))),
            c(xy_hi$y, xy_hi$y[1], rev(c(xy_lo$y, xy_lo$y[1]))),
            border = NA, col = adjustcolor(cols[i], alpha.f = 0.18))
    xy_mean <- polar_xy(mean_vals)
    lines(c(xy_mean$x, xy_mean$x[1]), c(xy_mean$y, xy_mean$y[1]),
          col = cols[i], lty = group_linetypes[i], lwd = 1.6)
  }

  axis_vals <- seq(0, max_val, length.out = seg + 1)
  label_digits <- ifelse(max_val >= 10, 1, 2)
  ax_x <- cos(thetas[1])
  ax_y <- sin(thetas[1])
  text(0, 0, labels = "0", cex = 0.42, col = "grey35", pos = 3)
  for (s in seq_len(seg)) {
    r <- axis_vals[s + 1] / max_val
    text(
      x = ax_x * r,
      y = ax_y * r,
      labels = formatC(axis_vals[s + 1], format = "f", digits = label_digits),
      cex = 0.42,
      col = "grey35",
      pos = 3
    )
  }
  label_text <- str_wrap(var_order, width = 11)
  text(1.17 * cos(thetas), 1.17 * sin(thetas), labels = label_text, cex = 1.4)
  legend(x = -2.4, y = 1.3, legend = legend_labels, col = cols,
         lty = group_linetypes, lwd = 2,
         bty = "n", horiz = FALSE, cex = 1.0)
  invisible(NULL)
}
plot_radar_subpopulations <- function(sub_summary,
                                      out_path,
                                      metrics,
                                      title,
                                      selected_comorbidities = fig1_main_comorbidities,
                                      include_selected = TRUE) {
  options(bitmapType = "cairo")
  race_levels <- c(
    "1" = "Black or African American",
    "2" = "White",
    "5" = "Asian"
  )
  sex_levels <- c("1" = "Male", "2" = "Female")

  plot_df <- prepare_metric_long(sub_summary) %>%
    filter(if (include_selected) {
      as.character(Comorbidity) %in% selected_comorbidities
    } else {
      !(as.character(Comorbidity) %in% selected_comorbidities)
    }) %>%
    filter(subgroup_var %in% c("race", "sex", "age_group", "UrbanRural")) %>%
    mutate(level_chr = as.character(subgroup_level)) %>%
    mutate(level_chr = sub("\\.0$", "", level_chr)) %>%
    filter(
      (subgroup_var == "race" & level_chr %in% names(race_levels)) |
        (subgroup_var == "sex" & level_chr %in% names(sex_levels)) |
        (subgroup_var == "age_group" & level_chr %in% c("< 65", ">= 65")) |
        (subgroup_var == "UrbanRural" & level_chr %in% c("Urban", "Rural"))
    ) %>%
    mutate(
      Cohort = factor(label_cohort(Scenario), levels = unname(cohort_labels)),
      Subgroup = case_when(
        subgroup_var == "race" ~ "Race",
        subgroup_var == "sex" ~ "Sex",
        subgroup_var == "age_group" ~ "Age group",
        subgroup_var == "UrbanRural" ~ "Urban/rural",
        TRUE ~ subgroup_var
      ),
      Level = case_when(
        subgroup_var == "race" ~ unname(race_levels[level_chr]),
        subgroup_var == "sex" ~ unname(sex_levels[level_chr]),
        TRUE ~ level_chr
      ),
      Curve = paste(Cohort, Level, sep = " - "),
      Subgroup_level = paste0(Subgroup, ": ", Level),
      Comorbidity_label = str_wrap(as.character(Comorbidity_label), width = 22)
    ) %>%
    mutate(
      Metric = factor(Metric, levels = c("aR0 (%)", "aR1 (%)", "aAR (%)", "aRR")),
      Subgroup = factor(Subgroup, levels = c("Race", "Sex", "Age group", "Urban/rural"))
    ) %>%
    filter(!is.na(Cohort)) %>%
    filter(Metric %in% metrics) %>%
    mutate(Metric = factor(as.character(Metric), levels = metrics))

  radar_comorbidity_order <- if (include_selected) {
    selected_comorbidities
  } else {
    setdiff(comorbidity_order, selected_comorbidities)
  }
  radar_var_order <- label_comorbidity(radar_comorbidity_order)

  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  n_rows <- length(levels(plot_df$Subgroup))
  n_cols <- length(levels(plot_df$Metric))
  png(out_path, width = 1120 * n_cols, height = 880 * n_rows, res = 300)
  op <- par(mfrow = c(n_rows, n_cols), mar = c(1.8, 5, 1.35, 2), oma = c(0.25, 0.2, 1.7, 0.2), xpd = NA)
  on.exit({
    par(op)
    dev.off()
  }, add = TRUE)

  panel_index <- 0L
  for (sg in levels(plot_df$Subgroup)) {
    for (mt in levels(plot_df$Metric)) {
      panel_index <- panel_index + 1
      panel_df <- plot_df %>% filter(Subgroup == sg, Metric == mt)
      radar_ci_panel(
        panel_df,
        group_col = "Curve",
        color_col = "Cohort",
        linetype_col = "Level",
        title = paste(sg, mt, sep = " - "),
        metric_label = mt,
        panel_label = paste0("(", letters[panel_index], ")"),
        var_order = radar_var_order,
        palette = cohort_colors
      )
    }
  }
  # mtext(title, outer = TRUE, cex = 1.1, font = 2)
  invisible(out_path)
}

plot_fig2_subpopulations <- function(sub_summary, out_path) {
  plot_radar_subpopulations(
    sub_summary = sub_summary,
    out_path = out_path,
    metrics = c("aAR (%)", "aRR"),
    title = "Figure 2. Shaded radar charts of targeted aAR and aRR across subpopulations",
    selected_comorbidities = fig1_main_comorbidities,
    include_selected = TRUE
  )
}

plot_figS2_subpopulations_other <- function(sub_summary, out_path) {
  plot_radar_subpopulations(
    sub_summary = sub_summary,
    out_path = out_path,
    metrics = c("aAR (%)", "aRR"),
    title = "Figure S2. Shaded radar charts of targeted aAR and aRR for other comorbidities",
    selected_comorbidities = fig1_main_comorbidities,
    include_selected = FALSE
  )
}

state_region_lookup <- function() {
  tibble(
    state_name = state.name,
    census_region = as.character(state.region)
  ) %>%
    mutate(
      region = case_when(
        census_region == "West" ~ "West",
        census_region == "North Central" ~ "Midwest",
        census_region == "South" ~ "South",
        census_region == "Northeast" ~ "Northeast",
        TRUE ~ NA_character_
      )
    )
}

plot_region_map_metric <- function(sub_summary,
                                   metric,
                                   fill_label,
                                   selected_comorbidities = fig1_main_comorbidities,
                                   include_selected = TRUE,
                                   map_layout = c("cohort_rows", "fig3_four_rows")) {
  map_layout <- match.arg(map_layout)
  comorbidity_levels <- if (include_selected) {
    selected_comorbidities
  } else {
    setdiff(comorbidity_order, c(selected_comorbidities, "Any_comorbidities"))
  }

  metric_df <- prepare_metric_long(sub_summary) %>%
    filter(subgroup_var == "region", Metric == metric) %>%
    filter(as.character(Comorbidity) %in% comorbidity_levels) %>%
    mutate(region = as.character(subgroup_level)) %>%
    filter(region %in% c("West", "Midwest", "South", "Northeast")) %>%
    mutate(
      Cohort = factor(label_cohort(Scenario), levels = unname(cohort_labels)),
      Comorbidity_label = factor(
        str_wrap(as.character(Comorbidity_label), width = 16),
        levels = str_wrap(label_comorbidity(comorbidity_levels), width = 16)
      ),
      Fig3_panel = factor(
        paste(label_cohort(Scenario), label_comorbidity(as.character(Comorbidity)), sep = "\n"),
        levels = fig3_four_row_panel_levels()
      ),
      label = paste0(
        sprintf(ifelse(metric == "aRR", "%.3f", "%.2f"), estimate),
        "\n[",
        sprintf(ifelse(metric == "aRR", "%.3f", "%.2f"), low),
        ", ",
        sprintf(ifelse(metric == "aRR", "%.3f", "%.2f"), high),
        "]"
      )
    ) %>%
    filter(!is.na(Cohort), !is.na(Comorbidity_label))

  if (map_layout == "fig3_four_rows") {
    metric_df <- metric_df %>% filter(!is.na(Fig3_panel))
  }

  states_map <- map_data("state") %>%
    as_tibble() %>%
    rename(state_map_name = region) %>%
    mutate(state_name = tools::toTitleCase(state_map_name)) %>%
    left_join(state_region_lookup(), by = "state_name") %>%
    left_join(metric_df, by = "region", relationship = "many-to-many") %>%
    filter(!is.na(Comorbidity_label))

  label_df <- states_map %>%
    group_by(region, Cohort, Comorbidity_label, Fig3_panel) %>%
    summarise(
      long = mean(range(long), na.rm = TRUE),
      lat = mean(range(lat), na.rm = TRUE),
      estimate = first(estimate),
      label = first(label),
      .groups = "drop"
    ) %>%
    filter(is.finite(estimate)) %>%
    mutate(
      lat = case_when(
        region == "Midwest" ~ lat + 4.3,
        region == "Northeast" ~ lat - 1.0,
        TRUE ~ lat
      ),
      long = ifelse(region == "Northeast", long - 3.0, long)
    )

  p <- ggplot(states_map, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = estimate, color = Cohort), linewidth = 0.18) +
    geom_text(
      data = label_df,
      inherit.aes = FALSE,
      aes(x = long, y = lat, label = label),
      size = 4,
      lineheight = 0.9,
      color = "black"
    ) +
    coord_fixed(1.3)

  if (map_layout == "fig3_four_rows") {
    p <- p + facet_wrap(~Fig3_panel, ncol = 4)
  } else {
    p <- p + facet_grid(Cohort ~ Comorbidity_label)
  }

  p +
    scale_color_manual(values = cohort_colors, guide = "none") +
    scale_fill_gradientn(
      colours = colorRampPalette(rev(brewer.pal(11, "Spectral")))(100),
      name = fill_label,
      na.value = "white",
      guide = guide_colorbar(
        direction = "horizontal",
        # title.position = "top",
        barwidth = unit(18, "cm"),
        barheight = unit(0.5, "cm")
      )
    ) +
    theme_void(base_size = 16) +
    theme(
      # legend.position = "inside",
      legend.position = "bottom",
      legend.background = element_rect(
        fill = alpha("white", 0.7),
        colour = NA
      ),
      legend.title = element_text(size = 20,
                                  margin = margin(b = 25)),
      legend.text = element_text(size = 18),
      strip.text = element_text(size = 14, face = "bold"),
      panel.spacing = unit(1, "lines"),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background  = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "white", colour = NA)
    )
}


plot_fig3_region <- function(sub_summary,
                             selected_comorbidities = fig3_main_comorbidities,
                             include_selected = TRUE,
                             map_layout = "fig3_four_rows") {
  plot_region_map_metric(
    sub_summary = sub_summary,
    metric = "aAR (%)",
    fill_label = "aAR (%)",
    selected_comorbidities = selected_comorbidities,
    include_selected = include_selected,
    map_layout = map_layout
  )
}

plot_fig4_covid_year <- function(sub_summary,
                                 metric = "aAR (%)",
                                 y_label = "Adjusted attributable risk (%)",
                                 cohort_names = unname(cohort_labels)) {
  op <- par(mar = c(5.8, 5, 1.35, 2), oma = c(0.25, 0.2, 1.7, 0.2), xpd = NA)
  plot_colors <- cohort_colors[cohort_names]
  plot_shapes <- cohort_shapes[cohort_names]
  plot_linetypes <- cohort_linetypes[cohort_names]
  sub_summary <- sub_summary[sub_summary$Comorbidity %nin% c("Any_comorbidities"),]
  plot_df <- prepare_metric_long(sub_summary) %>%
    filter(subgroup_var == "Covid_19_Year", Metric == metric) %>%
    mutate(Year = as.integer(as.character(subgroup_level))) %>%
    filter(Year %in% 2020:2024) %>%
    mutate(
      Cohort = factor(label_cohort(Scenario), levels = cohort_names),
      Comorbidity_label = factor(
        str_wrap(as.character(Comorbidity_label), width = 22),
        levels = str_wrap(label_comorbidity(comorbidity_order), width = 22)
      )
    ) %>%
    filter(!is.na(Cohort)) %>%
    arrange(Comorbidity_label, Cohort, Year)

  reference_line <- ifelse(metric == "aRR", 1, 0)

  ggplot(
    plot_df,
    aes(
      x = Year,
      y = estimate,
      color = Cohort,
      fill = Cohort,
      linetype = Cohort,
      group = Cohort
    )
  ) +
    geom_hline(yintercept = reference_line, linetype = 2, color = "grey60") +
    geom_ribbon(aes(ymin = low, ymax = high), alpha = 0.18, color = NA, na.rm = TRUE) +
    geom_line(linewidth = 1, na.rm = TRUE) +
    geom_point(
      aes(shape = Cohort),
      size = 3.5,
      stroke = 0.5,
      fill = "white",
      na.rm = TRUE
    ) +
    # geom_point(shape = 21, color = "white", size = 1.4, stroke = 0.25, na.rm = TRUE) +
    facet_wrap(~Comorbidity_label, ncol = 5) +
    scale_x_continuous(breaks = 2020:2024) +
    scale_color_manual(values = plot_colors) +
    scale_fill_manual(values = plot_colors) +
    scale_linetype_manual(values = plot_linetypes) +
    labs(x = "COVID-19 diagnosis year", y = y_label) +
    scale_shape_manual(
      values = plot_shapes
    ) +
    base_plot_theme(base_size = 12) +
    guides(
      fill = "none",
      shape = "none",
      linetype = "none",
      color = guide_legend(
        override.aes = list(
          shape = unname(plot_shapes),
          color = unname(plot_colors),
          fill = unname(plot_colors),
          linetype = unname(plot_linetypes),
          linewidth = 1.2,
          size = 3,
          stroke = 1
        )
      )
    ) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.text = element_text(size = 16),
      legend.title = element_blank(),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key.width = unit(5.2, "cm"),
      strip.text = element_text(size = 14, face = "bold"),)
}

stop_if_missing(cohort_i_root)
stop_if_missing(cohort_ii_root)
stop_if_missing(cohort_iii_root)
stop_if_missing(cohort_iv_root)

Summary.CohortI <- read_overall_summary(cohort_i_root, "CohortI")
Summary.CohortII <- read_overall_summary(cohort_ii_root, "CohortII")
Summary.CohortIII <- read_overall_summary(cohort_iii_root, "CohortIII")
Summary.CohortIV <- read_overall_summary(cohort_iv_root, "CohortIV")
Summary <- bind_rows(
  Summary.CohortIV,
  Summary.CohortII,
  Summary.CohortIII,
  Summary.CohortI
)

Sub.summary.CohortI <- read_sub_summary(cohort_i_root, "CohortI")
Sub.summary.CohortII <- read_sub_summary(cohort_ii_root, "CohortII")
Sub.summary.CohortIII <- read_sub_summary(cohort_iii_root, "CohortIII")
Sub.summary.CohortIV <- read_sub_summary(cohort_iv_root, "CohortIV")
Sub.summary.Cohorts <- bind_rows(
  Sub.summary.CohortI,
  Sub.summary.CohortII,
  Sub.summary.CohortIII,
  Sub.summary.CohortIV
)
Sub.summary.Fig4 <- bind_rows(
  Sub.summary.CohortI,
  Sub.summary.CohortIII,
  Sub.summary.CohortIV
)
Sub.summary <- Sub.summary.Cohorts

table1 <- build_table1(Summary)




write.csv(table1, file.path(table_dir, "Table1_DML_targeted_I_IV.csv"), row.names = FALSE)
# saveRDS(table1, file.path(table_dir, "Table1_DML_targeted.rds"))

fig1 <- plot_multimorbidity_metric_by_cohort(
  Sub.summary.Cohorts,
  metric = "aAR (%)",
  selected_comorbidities = fig1_main_comorbidities,
  include_selected = TRUE
)
# ggsave(file.path(figure_dir, "Fig1_multimorbidity_combined.png"),
#        fig1, width = 18, height = 4, dpi = 300)
#
# figS1_1 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aRR",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = TRUE
# )
# ggsave(file.path(figure_dir, "FigS1.1_aRR_multimorbidity_cohorts.png"),
#        figS1_1, width = 18, height = 4, dpi = 300)
#
# figS1_2 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aR0 (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = TRUE
# )
# ggsave(file.path(figure_dir, "FigS1.2_aR0_multimorbidity_cohorts.png"),
#        figS1_2, width = 18, height = 4, dpi = 300)
#
# figS1_3 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aR1 (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = TRUE
# )
# ggsave(file.path(figure_dir, "FigS1.3_aR1_multimorbidity_cohorts.png"),
#        figS1_3, width = 18, height = 4, dpi = 300)
#
# figS1_1_0 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aAR (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS1.1.0_aAR_multimorbidity_other_cohorts.png"),
#        figS1_1_0, width = 23, height = 4, dpi = 300)
#
# figS1_1_1 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aRR",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS1.1.1_aRR_multimorbidity_other_cohorts.png"),
#        figS1_1_1, width = 23, height = 4, dpi = 300)
#
# figS1_2_1 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aR0 (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS1.2.1_aR0_multimorbidity_other_cohorts.png"),
#        figS1_2_1, width = 23, height = 4, dpi = 300)
#
# figS1_3_1 <- plot_multimorbidity_metric_by_cohort(
#   Sub.summary.Cohorts,
#   metric = "aR1 (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS1.3.1_aR1_multimorbidity_other_cohorts.png"),
#        figS1_3_1, width = 23, height = 4, dpi = 300)
#
# plot_fig2_subpopulations(
#   Sub.summary,
#   file.path(figure_dir, "Fig2_subpopulations_combined.png")
# )
#
# plot_figS2_subpopulations_other(
#   Sub.summary,
#   file.path(figure_dir, "FigS2_subpopulations_other_combined.png")
# )
#
# fig3 <- plot_fig3_region(
#   Sub.summary.Cohorts,
#   selected_comorbidities = fig3_main_comorbidities,
#   include_selected = TRUE,
#   map_layout = "fig3_four_rows"
# )
# ggsave(file.path(figure_dir, "Fig3_ARs_I_II_III_IV.png"),
#        fig3, width = 14, height = 28, dpi = 300)
#
# figS3_aRR <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aRR",
#   fill_label = "aRR",
#   selected_comorbidities = fig3_main_comorbidities,
#   include_selected = TRUE,
#   map_layout = "fig3_four_rows"
# )
# ggsave(file.path(figure_dir, "FigS3_aRR_I_II_III_IV.png"),
#        figS3_aRR, width = 14, height = 28, dpi = 300)
#
# figS3_aR0 <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aR0 (%)",
#   fill_label = "aR0 (%)",
#   selected_comorbidities = fig3_main_comorbidities,
#   include_selected = TRUE,
#   map_layout = "fig3_four_rows"
# )
# ggsave(file.path(figure_dir, "FigS3_aR0_I_II_III_IV.png"),
#        figS3_aR0, width = 14, height = 28, dpi = 300)
#
# figS3_aR1 <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aR1 (%)",
#   fill_label = "aR1 (%)",
#   selected_comorbidities = fig3_main_comorbidities,
#   include_selected = TRUE,
#   map_layout = "fig3_four_rows"
# )
# ggsave(file.path(figure_dir, "FigS3_aR1_I_II_III_IV.png"),
#        figS3_aR1, width = 14, height = 28, dpi = 300)
#
# figS3_other_aAR <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aAR (%)",
#   fill_label = "aAR (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS3_other_aAR_I_II_III_IV.png"),
#        figS3_other_aAR, width = 22, height = 16, dpi = 300)
#
# figS3_other_aRR <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aRR",
#   fill_label = "aRR",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS3_other_aRR_I_II_III_IV.png"),
#        figS3_other_aRR, width = 22, height = 16, dpi = 300)
#
# figS3_other_aR0 <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aR0 (%)",
#   fill_label = "aR0 (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS3_other_aR0_I_II_III_IV.png"),
#        figS3_other_aR0, width = 22, height = 16, dpi = 300)
#
# figS3_other_aR1 <- plot_region_map_metric(
#   sub_summary = Sub.summary.Cohorts,
#   metric = "aR1 (%)",
#   fill_label = "aR1 (%)",
#   selected_comorbidities = fig1_main_comorbidities,
#   include_selected = FALSE
# )
# ggsave(file.path(figure_dir, "FigS3_other_aR1_I_II_III_IV.png"),
#        figS3_other_aR1, width = 22, height = 16, dpi = 300)

# fig4 <- plot_fig4_covid_year(
#   Sub.summary.Fig4,
#   cohort_names = fig4_cohort_names
# )
# ggsave(file.path(figure_dir, "Fig4_covid_year_aAR_Cohorts_I_III_IV.png"),
#        fig4, width = 14, height = 11, dpi = 300)

figS4_1 <- plot_fig4_covid_year(
  Sub.summary.Cohorts,
  metric = "aR0 (%)",
  y_label = "Adjusted R0 (%)"
)
ggsave(file.path(figure_dir, "FigS12_aR0_year_Cohorts_I_II_III_IV.png"),
       figS4_1, width = 14, height = 11, dpi = 300)

figS4_2 <- plot_fig4_covid_year(
  Sub.summary.Cohorts,
  metric = "aR1 (%)",
  y_label = "Adjusted R1 (%)"
)
ggsave(file.path(figure_dir, "FigS13_aR1_year_Cohorts_I_II_III_IV.png"),
       figS4_2, width = 14, height = 11, dpi = 300)

figS4_3 <- plot_fig4_covid_year(
  Sub.summary.Cohorts,
  metric = "aAR (%)",
  y_label = "Adjusted attributable risks (%)"
)
ggsave(file.path(figure_dir, "FigS14_aAR_year_Cohorts_I_II_III_IV.png"),
       figS4_3, width = 14, height = 11, dpi = 300)

figS4_3 <- plot_fig4_covid_year(
  Sub.summary.Cohorts,
  metric = "aRR",
  y_label = "Adjusted relative risks"
)
ggsave(file.path(figure_dir, "FigS15_aRR_year_Cohorts_I_II_III_IV.png"),
       figS4_3, width = 14, height = 11, dpi = 300)

message("Done. Outputs saved to:")
message("  ", normalizePath(table_dir, mustWork = FALSE))
message("  ", normalizePath(figure_dir, mustWork = FALSE))
