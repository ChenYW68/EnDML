############################################################
## Propensity-score overlap / positivity diagnostics
##
## Inputs produced by LC_risk_ensemble_DML.R:
##   results/<cohort>/diagnostics/<comorbidity>/
##     DIAG.SCORES.*.parquet
##     DIAG.SUBGROUP_POSITIVITY.*.parquet
##
## This script does not use CV.*.parquet and does not refit models.
############################################################

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

## User-adjustable settings ---------------------------------
output_dir <- file.path("new_Results", "propensity_overlap")
max_density_n_per_treatment <- 20000L
central_support_quantiles <- c(0.01, 0.99)
histogram_breaks <- seq(0, 1, by = 0.01)
minimum_subgroup_n <- 50L
set.seed(20260925)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cohort_levels <- c("Cohort I", "Cohort II", "Cohort III", "Cohort IV")
cohort_colors <- c(
  "Cohort I" = "black",
  "Cohort II" = "green",
  "Cohort III" = "red",
  "Cohort IV" = "blue"
)
status_levels <- c("No major concern", "Concern", "Serious concern", "Failure")
status_colors <- c(
  "No major concern" = "#2E8B57",
  "Concern" = "#F2B134",
  "Serious concern" = "#D95F02",
  "Failure" = "#B2182B"
)

candidate_results_roots <- unique(c(
  "results",
  file.path("Results", "results"),
  file.path("..", "results")
))
candidate_results_roots <- candidate_results_roots[dir.exists(candidate_results_roots)]

if (length(candidate_results_roots) == 0) {
  stop(
    "No results directory found. Run from the project root or the Results directory.",
    call. = FALSE
  )
}

score_files <- unique(unlist(lapply(candidate_results_roots, function(root) {
  list.files(
    root,
    pattern = "^DIAG[.]SCORES[.].*[.]parquet$",
    recursive = TRUE,
    full.names = TRUE
  )
})))

if (length(score_files) == 0) {
  stop(
    "No DIAG.SCORES.*.parquet files found under: ",
    paste(candidate_results_roots, collapse = ", "),
    call. = FALSE
  )
}

score_files <- score_files[!duplicated(normalizePath(score_files, mustWork = TRUE))]

cohort_from_path <- function(path, question = NA_integer_) {
  clean_path <- gsub("\\\\", "/", normalizePath(path, mustWork = FALSE))
  cohort_dir <- sub("/diagnostics/.*$", "", clean_path)
  cohort_dir <- basename(cohort_dir)

  known <- c(
    Cohort_I = "Cohort IV",
    Cohort_II = "Cohort II",
    Cohort_III = "Cohort III",
    Cohort_IV = "Cohort I"#,
    # new4_before = "Cohort I",
    # new4_combined = "Cohort II",
    # new4_before_retained_unexposed = "Cohort III"
  )
  label <- unname(known[cohort_dir])
  if (!is.na(label)) return(label)

  question_map <- c(
    `1` = "Cohort I",
    `2` = "Cohort II",
    `3` = "Cohort III",
    `4` = "Cohort IV"
  )
  question_label <- unname(question_map[as.character(question[1])])
  if (!is.na(question_label)) return(question_label)
  cohort_dir
}

pretty_comorbidity <- function(x) {
  str_to_sentence(gsub("_", " ", x))
}

safe_quantile <- function(x, probability) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probability, names = FALSE, na.rm = TRUE))
}

histogram_overlap <- function(ps0, ps1) {
  if (length(ps0) == 0 || length(ps1) == 0) return(NA_real_)
  h0 <- hist(ps0, breaks = histogram_breaks, plot = FALSE)$counts
  h1 <- hist(ps1, breaks = histogram_breaks, plot = FALSE)$counts
  sum(pmin(h0 / sum(h0), h1 / sum(h1)))
}

classify_overall <- function(n0,
                             n1,
                             common_width,
                             overlap_coefficient,
                             pct_outside_common,
                             pct_extreme,
                             ess_fraction) {
  if (n0 == 0 || n1 == 0 || !is.finite(common_width) || common_width <= 0) {
    return("Failure")
  }
  if (
    overlap_coefficient < 0.10 ||
      pct_outside_common > 0.20 ||
      pct_extreme > 0.10 ||
      ess_fraction < 0.10
  ) {
    return("Serious concern")
  }
  if (
    overlap_coefficient < 0.30 ||
      pct_outside_common > 0.10 ||
      pct_extreme > 0.05 ||
      ess_fraction < 0.25
  ) {
    return("Concern")
  }
  "No major concern"
}

summarize_score_file <- function(path) {
  dat <- as.data.frame(read_parquet(
    path,
    col_select = c("D", "m_hat", "Comorbidity", "Question")
  )) %>%
    filter(is.finite(D), D %in% c(0, 1), is.finite(m_hat)) %>%
    mutate(m_hat = pmin(pmax(as.numeric(m_hat), 0), 1))

  if (nrow(dat) == 0) return(NULL)

  comorbidity <- as.character(dat$Comorbidity[which(!is.na(dat$Comorbidity))[1]])
  if (is.na(comorbidity) || !nzchar(comorbidity)) {
    comorbidity <- basename(dirname(path))
  }
  cohort <- cohort_from_path(path, dat$Question)

  ps0 <- dat$m_hat[dat$D == 0]
  ps1 <- dat$m_hat[dat$D == 1]
  n0 <- length(ps0)
  n1 <- length(ps1)

  q0_low <- safe_quantile(ps0, central_support_quantiles[1])
  q0_high <- safe_quantile(ps0, central_support_quantiles[2])
  q1_low <- safe_quantile(ps1, central_support_quantiles[1])
  q1_high <- safe_quantile(ps1, central_support_quantiles[2])
  common_low <- max(q0_low, q1_low, na.rm = TRUE)
  common_high <- min(q0_high, q1_high, na.rm = TRUE)
  common_width <- max(0, common_high - common_low)
  union_width <- max(q0_high, q1_high, na.rm = TRUE) -
    min(q0_low, q1_low, na.rm = TRUE)
  central_support_ratio <- ifelse(union_width > 0, common_width / union_width, NA_real_)

  ps <- dat$m_hat
  pct_outside_common <- if (common_width > 0) {
    mean(ps < common_low | ps > common_high)
  } else {
    1
  }
  pct_extreme <- mean(ps < 0.01 | ps > 0.99)
  pct_near_extreme <- mean(ps < 0.05 | ps > 0.95)
  overlap_coefficient <- histogram_overlap(ps0, ps1)

  bounded_ps <- pmin(pmax(ps, 1e-6), 1 - 1e-6)
  weights <- dat$D / bounded_ps + (1 - dat$D) / (1 - bounded_ps)
  ess <- sum(weights)^2 / sum(weights^2)
  ess_fraction <- ess / length(weights)

  summary_row <- tibble(
    Cohort = cohort,
    Comorbidity = comorbidity,
    Comorbidity_label = pretty_comorbidity(comorbidity),
    n = nrow(dat),
    n_control = n0,
    n_treated = n1,
    treated_rate = mean(dat$D),
    control_p01 = q0_low,
    control_p05 = safe_quantile(ps0, 0.05),
    control_p50 = safe_quantile(ps0, 0.50),
    control_p95 = safe_quantile(ps0, 0.95),
    control_p99 = q0_high,
    treated_p01 = q1_low,
    treated_p05 = safe_quantile(ps1, 0.05),
    treated_p50 = safe_quantile(ps1, 0.50),
    treated_p95 = safe_quantile(ps1, 0.95),
    treated_p99 = q1_high,
    common_support_low = common_low,
    common_support_high = common_high,
    common_support_width = common_width,
    central_support_ratio = central_support_ratio,
    histogram_overlap_coefficient = overlap_coefficient,
    pct_outside_common_support = pct_outside_common,
    pct_ps_below_0_01_or_above_0_99 = pct_extreme,
    pct_ps_below_0_05_or_above_0_95 = pct_near_extreme,
    ipw_p99 = safe_quantile(weights, 0.99),
    ipw_max = max(weights),
    ipw_ess = ess,
    ipw_ess_fraction = ess_fraction
  )
  summary_row$positivity_status <- classify_overall(
    n0 = n0,
    n1 = n1,
    common_width = common_width,
    overlap_coefficient = overlap_coefficient,
    pct_outside_common = pct_outside_common,
    pct_extreme = pct_extreme,
    ess_fraction = ess_fraction
  )

  sampled <- dat %>%
    group_by(D) %>%
    group_modify(~ slice_sample(
      .x,
      n = min(nrow(.x), max_density_n_per_treatment)
    )) %>%
    ungroup() %>%
    transmute(
      Cohort = cohort,
      Comorbidity = comorbidity,
      Comorbidity_label = pretty_comorbidity(comorbidity),
      D,
      Treatment = factor(D, levels = c(0, 1), labels = c("Control (D=0)", "Treated (D=1)")),
      m_hat
    )

  list(summary = summary_row, sample = sampled)
}

message("Reading ", length(score_files), " propensity-score diagnostic files...")
results <- lapply(score_files, summarize_score_file)
results <- results[!vapply(results, is.null, logical(1))]

if (length(results) == 0) {
  stop("The diagnostic files contained no finite D/m_hat observations.", call. = FALSE)
}

overall_summary <- bind_rows(lapply(results, `[[`, "summary")) %>%
  mutate(
    Cohort = factor(Cohort, levels = unique(c(cohort_levels, as.character(Cohort)))),
    positivity_status = factor(positivity_status, levels = status_levels)
  ) %>%
  arrange(Cohort, Comorbidity_label)

density_sample <- bind_rows(lapply(results, `[[`, "sample")) %>%
  mutate(Cohort = factor(Cohort, levels = levels(overall_summary$Cohort)))

write.csv(
  overall_summary,
  file.path(output_dir, "positivity_overlap_diagnostics.csv"),
  row.names = FALSE
)
# saveRDS(overall_summary, file.path(output_dir, "positivity_overlap_diagnostics.rds"))

## Simple D=0 versus D=1 density plots ----------------------
## This reproduces the propensity-overlap plot written by
## LC_risk_ensemble_DML.R, using a bounded plotting sample when a diagnostic
## dataset is very large. Numerical diagnostics above always use all records.
individual_plot_dir <- file.path(output_dir, "individual_propensity_overlap")
dir.create(individual_plot_dir, recursive = TRUE, showWarnings = FALSE)

individual_groups <- density_sample %>%
  distinct(Cohort, Comorbidity, Comorbidity_label)

# for (i in seq_len(nrow(individual_groups))) {
#   cohort_i <- as.character(individual_groups$Cohort[i])
#   comorbidity_i <- individual_groups$Comorbidity[i]
#   comorbidity_label_i <- individual_groups$Comorbidity_label[i]
#   ps_plot_dt <- density_sample %>%
#     filter(
#       as.character(Cohort) == cohort_i,
#       Comorbidity == comorbidity_i,
#       is.finite(D),
#       is.finite(m_hat)
#     )
#
#   if (nrow(ps_plot_dt) > 0) {
#     p1 <- ggplot(
#       ps_plot_dt,
#       aes(
#         x = m_hat,
#         y = after_stat(scaled),
#         fill = factor(D),
#         color = factor(D)
#       )
#     ) +
#       geom_density(alpha = 0.2) +
#       labs(
#         title = paste("Propensity overlap -", cohort_i, "-", comorbidity_label_i),
#         x = "Estimated propensity score",
#         y = "Standardized density",
#         fill = "D",
#         color = "D"
#       ) +
#       theme_minimal()
#
#     cohort_dir <- gsub("[^A-Za-z0-9]+", "_", cohort_i)
#     comorbidity_dir <- gsub("[^A-Za-z0-9]+", "_", comorbidity_i)
#     out_dir <- file.path(individual_plot_dir, cohort_dir, comorbidity_dir)
#     dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
#     ggsave(
#       filename = file.path(out_dir, "propensity_overlap.png"),
#       plot = p1,
#       width = 8,
#       height = 5,
#       dpi = 300
#     )
#   }
# }

## Overall status heatmap ------------------------------------
heatmap_plot <- ggplot(
  overall_summary %>%
    mutate(Comorbidity_label = factor(
      Comorbidity_label,
      levels = sort(unique(Comorbidity_label))
    )),
  aes(x = Cohort, y = Comorbidity_label, fill = positivity_status)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = sprintf("%.2f", histogram_overlap_coefficient)),
    size = 3.2
  ) +
  scale_fill_manual(values = status_colors, drop = FALSE) +
  labs(
    title = "Propensity-score overlap diagnostic",
    subtitle = "Cell labels are histogram overlap coefficients (0=no overlap; 1=identical distributions)",
    x = NULL,
    y = NULL,
    fill = "Diagnostic status"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "bottom")

ggsave(
  file.path(output_dir, "positivity_status_heatmap.png"),
  heatmap_plot,
  width = 11,
  height = max(7, 0.35 * length(unique(overall_summary$Comorbidity))),
  dpi = 300
)

## Cohort-specific overlap figures --------------------------
for (cohort in unique(as.character(overall_summary$Cohort))) {
  cohort_density <- density_sample %>% filter(as.character(Cohort) == cohort)
  cohort_summary <- overall_summary %>% filter(as.character(Cohort) == cohort)
  if (nrow(cohort_density) == 0) next

  density_plot <- ggplot() +
    geom_rect(
      data = cohort_summary,
      aes(
        xmin = common_support_low,
        xmax = common_support_high,
        ymin = -Inf,
        ymax = Inf
      ),
      inherit.aes = FALSE,
      fill = "grey80",
      alpha = 0.35
    ) +
    geom_density(
      data = cohort_density,
      aes(
        x = m_hat,
        y = after_stat(scaled),
        color = Treatment,
        fill = Treatment
      ),
      alpha = 0.18,
      linewidth = 0.7,
      adjust = 1
    ) +
    geom_vline(xintercept = c(0.05, 0.95), linetype = "dashed", color = "grey45") +
    facet_wrap(~Comorbidity_label, scales = "free_y", ncol = 4) +
    scale_color_manual(values = c("Control (D=0)" = "#2166AC", "Treated (D=1)" = "#B2182B")) +
    scale_fill_manual(values = c("Control (D=0)" = "#2166AC", "Treated (D=1)" = "#B2182B")) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = paste("Propensity-score overlap -", cohort),
      subtitle = "Grey area: common central (1st-99th percentile) support; dashed lines: 0.05 and 0.95",
      x = "Cross-fitted estimated propensity score",
      y = "Standardized density",
      color = NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top", panel.grid.minor = element_blank())

  safe_cohort <- gsub("[^A-Za-z0-9]+", "_", cohort)
  ggsave(
    file.path(output_dir, paste0("propensity_overlap_density_", safe_cohort, ".png")),
    density_plot,
    width = 14,
    height = max(8, 2.4 * ceiling(n_distinct(cohort_density$Comorbidity) / 4)),
    dpi = 300
  )

  interval_data <- bind_rows(
    cohort_summary %>% transmute(
      Comorbidity_label,
      Treatment = "Control (D=0)",
      low = control_p01,
      median = control_p50,
      high = control_p99
    ),
    cohort_summary %>% transmute(
      Comorbidity_label,
      Treatment = "Treated (D=1)",
      low = treated_p01,
      median = treated_p50,
      high = treated_p99
    )
  )

  interval_plot <- ggplot(
    interval_data,
    aes(y = reorder(Comorbidity_label, median), color = Treatment)
  ) +
    geom_linerange(
      aes(xmin = low, xmax = high),
      position = position_dodge(width = 0.55),
      linewidth = 1
    ) +
    # geom_point(
    #   aes(x = median),
    #   position = position_dodge(width = 0.55),
    #   size = 2
    # ) +#2166AC/#B2182B
    scale_color_manual(values = c("Control (D=0)" = "gray80",
                                  "Treated (D=1)" = "gray30")) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = paste("Central propensity-score support -", cohort),
      subtitle = "Lines show 1st-99th percentiles",#; points show medians
      x = "Cross-fitted estimated propensity score",
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top",
          panel.grid.minor = element_blank(),
          legend.key.width = unit(3.8, "cm")
          )

  ggsave(
    file.path(output_dir, paste0("FigS18_propensity_common_support_", safe_cohort, ".png")),
    interval_plot,
    width = 10,
    height = max(7, 0.35 * n_distinct(interval_data$Comorbidity_label)),
    dpi = 300
  )
}

## Subgroup positivity summaries ----------------------------
subgroup_files <- sub("DIAG[.]SCORES[.]", "DIAG.SUBGROUP_POSITIVITY.", score_files)
subgroup_files <- subgroup_files[file.exists(subgroup_files)]

if (length(subgroup_files) > 0) {
  subgroup_rows <- lapply(subgroup_files, function(path) {
    dat <- as.data.frame(read_parquet(path))
    if (nrow(dat) == 0) return(NULL)
    cohort <- cohort_from_path(path)
    comorbidity <- basename(dirname(path))
    dat %>%
      mutate(
        Cohort = cohort,
        Comorbidity = comorbidity,
        Comorbidity_label = pretty_comorbidity(comorbidity),
        subgroup_level = as.character(subgroup_level),
        positivity_status = case_when(
          n < minimum_subgroup_n ~ "Too small to assess",
          treated_rate <= 0 | treated_rate >= 1 ~ "Failure",
          treated_rate < 0.01 | treated_rate > 0.99 |
            pct_ps_lt_0_05 > 0.50 | pct_ps_gt_0_95 > 0.50 ~ "Serious concern",
          treated_rate < 0.05 | treated_rate > 0.95 |
            pct_ps_lt_0_05 > 0.20 | pct_ps_gt_0_95 > 0.20 ~ "Concern",
          TRUE ~ "No major concern"
        )
      )
  })
  subgroup_summary <- bind_rows(subgroup_rows) %>%
    arrange(Cohort, Comorbidity_label, subgroup_var, subgroup_level)

  write.csv(
    subgroup_summary,
    file.path(output_dir, "subgroup_positivity_diagnostics.csv"),
    row.names = FALSE
  )
  saveRDS(subgroup_summary, file.path(output_dir, "subgroup_positivity_diagnostics.rds"))

  subgroup_status_counts <- subgroup_summary %>%
    count(Cohort, positivity_status, name = "n_subgroups")
  write.csv(
    subgroup_status_counts,
    file.path(output_dir, "subgroup_positivity_status_counts.csv"),
    row.names = FALSE
  )
}

overall_status_counts <- overall_summary %>%
  count(Cohort, positivity_status, name = "n_comorbidities")
report_lines <- c(
  "Propensity-score overlap / empirical positivity assessment",
  "=========================================================",
  "",
  "Important: estimated-score diagnostics cannot prove the theoretical positivity assumption.",
  "They identify empirical overlap problems in the observed analysis sample.",
  "",
  "Overall status counts:",
  capture.output(print(as.data.frame(overall_status_counts), row.names = FALSE)),
  "",
  "Overall classification thresholds:",
  "- Failure: one treatment group is absent or central 1st-99th percentile support does not overlap.",
  "- Serious concern: overlap coefficient <0.10, >20% outside common support, >10% extreme scores, or IPW ESS fraction <0.10.",
  "- Concern: overlap coefficient <0.30, >10% outside common support, >5% extreme scores, or IPW ESS fraction <0.25.",
  "- No major concern: none of the preceding flags.",
  ""
)
if (exists("subgroup_status_counts")) {
  report_lines <- c(
    report_lines,
    "Subgroup status counts:",
    capture.output(print(as.data.frame(subgroup_status_counts), row.names = FALSE)),
    "",
    paste0("Subgroups with n < ", minimum_subgroup_n, " are labeled too small to assess.")
  )
}
writeLines(report_lines, file.path(output_dir, "positivity_assessment.txt"))

cat("\nPositivity/overlap diagnostics completed.\n")
cat("Files written to:", normalizePath(output_dir, mustWork = FALSE), "\n")
cat("\nInterpretation note:\n")
cat("These diagnostics can identify empirical positivity concerns, but they cannot prove the positivity assumption.\n")
cat("Inspect density plots, common-support intervals, extreme-score proportions, IPW ESS, and subgroup failures together.\n")
