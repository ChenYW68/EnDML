
library(data.table)
library(arrow)
library(ggplot2)

# ============================================================
# 1. Paths and configuration
# ============================================================

# Use the same root.path and K as your modeling code.
# Example:
# root.path <- "D:/Your_Project/Results"
# K <- 5
cohort <- "Cohort_IV"
root.path <- paste0("./results/", cohort)
stopifnot(
  exists("root.path"),
  exists("K")
)

diag_root <- file.path(
  root.path,
  "diagnostics"
)

output_dir <- file.path(
  root.path,
  "Covariate_Balance_Figures"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# Display all covariate-level comparisons.
# Set to 25L to reproduce the original top-25 selection.
max_covariates <- Inf

# SMD reference threshold
smd_threshold <- 0.1

# Point styling
unweighted_color <- "grey50"
weighted_color   <- "black"

unweighted_shape <- 16  # Filled circle
weighted_shape   <- 17  # Filled triangle

# ============================================================
# 2. Original comorbidity identifiers
# ============================================================

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

# ============================================================
# 3. Required plotting order
# ============================================================

comorbidity_ids <- c(
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

stopifnot(
  length(comorbidity_ids) == 21L,
  setequal(
    comorbidity_ids,
    selected.Comorbidity
  )
)

comorbidity_labels <- gsub(
  "_",
  " ",
  comorbidity_ids
)

panel_labels <- LETTERS[
  seq_along(comorbidity_ids)
]

# ============================================================
# 4. Locate diagnostic balance file
# ============================================================

find_balance_file <- function(ci) {

  file_stub <- paste0(
    ci,
    ".Imputed_non_died_select_covariables_no_BMI_",
    K
  )

  file_path <- file.path(
    diag_root,
    ci,
    paste0(
      "DIAG.BALANCE.",
      file_stub,
      ".parquet"
    )
  )

  if (!file.exists(file_path)) {

    warning(
      "Missing balance file for ", ci,
      "\nExpected: ", file_path
    )

    return(NULL)
  }

  file_path
}

# ============================================================
# 5. Read balance data
# ============================================================

read_balance_data <- function(ci, label) {

  file_path <- find_balance_file(ci)

  if (is.null(file_path)) {
    return(NULL)
  }

  dt <- as.data.table(
    arrow::read_parquet(file_path)
  )

  # Exact column names from your build_balance_table()
  required_cols <- c(
    "covariate",
    "level",
    "smd_unweighted",
    "smd_weighted"
  )

  missing_cols <- setdiff(
    required_cols,
    names(dt)
  )

  if (length(missing_cols) > 0L) {

    stop(
      "Missing columns for ", ci, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  result <- dt[, .(
    covariate = as.character(covariate),
    level = as.character(level),
    Unweighted = as.numeric(smd_unweighted),
    Weighted = as.numeric(smd_weighted)
  )]

  # Create a unique label for each categorical level
  result[, covariate_label := fifelse(
    is.na(level),
    covariate,
    paste0(covariate, " = ", level)
  )]

  result[, Comorbidity := label]

  # Remove invalid rows
  result <- result[
    !is.na(covariate_label) &
      is.finite(Unweighted) &
      is.finite(Weighted)
  ]

  # Remove exact duplicate records only
  result <- unique(result)

  if (nrow(result) == 0L) {

    warning(
      "No valid balance observations for: ",
      ci
    )

    return(NULL)
  }

  # Check uniqueness of covariate-level labels
  if (anyDuplicated(result$covariate_label)) {

    stop(
      "Multiple distinct SMD estimates for the same ",
      "covariate-level combination: ", ci
    )
  }

  # Same ordering criterion as your original code
  result[, worst := pmax(
    abs(Unweighted),
    abs(Weighted)
  )]

  setorder(result, -worst)

  return(result)
}

# ============================================================
# 6. Covariate-balance plot function
# ============================================================

plot_covariate_balance <- function(
    balance_dt,
    panel,
    title,
    threshold = 0.1,
    max_covariates = Inf
) {

  dt <- copy(balance_dt)

  # Sort by maximum absolute SMD
  setorder(dt, -worst)

  # Optional top-N selection
  if (is.finite(max_covariates)) {
    dt <- head(
      dt,
      as.integer(max_covariates)
    )
  }

  # Highest imbalance appears at the top
  dt[, covariate_label := factor(
    covariate_label,
    levels = rev(
      as.character(covariate_label)
    )
  )]

  # Convert to long format
  plot_dt <- melt(
    dt,
    id.vars = c(
      "covariate_label",
      "Comorbidity"
    ),
    measure.vars = c(
      "Unweighted",
      "Weighted"
    ),
    variable.name = "Weighting",
    value.name = "SMD"
  )

  plot_dt[, Weighting := factor(
    Weighting,
    levels = c(
      "Unweighted",
      "Weighted"
    )
  )]

  # Independent x-axis range for each figure
  smd_range <- range(
    c(
      plot_dt$SMD,
      -threshold,
      threshold
    ),
    finite = TRUE
  )

  x_padding <- max(
    0.04,
    0.05 * diff(smd_range)
  )

  x_limits <- smd_range +
    c(-x_padding, x_padding)

  # ==========================================================
  # Plot
  # ==========================================================

  p <- ggplot(
    plot_dt,
    aes(
      x = SMD,
      y = covariate_label,
      shape = Weighting,
      color = Weighting
    )
  ) +

    # Reference lines
    geom_vline(
      xintercept = c(
        -threshold,
        threshold
      ),
      linetype = "dashed",
      color = "grey60",
      linewidth = 0.6
    ) +

    # Unweighted: gray circles
    # Weighted: black triangles
    geom_point(
      size = 2.8,
      stroke = 0
    ) +

    scale_shape_manual(
      values = c(
        "Unweighted" = unweighted_shape,
        "Weighted" = weighted_shape
      )
    ) +

    scale_color_manual(
      values = c(
        "Unweighted" = unweighted_color,
        "Weighted" = weighted_color
      )
    ) +

    scale_x_continuous(
      breaks = scales::breaks_pretty(n = 5)
    ) +

    coord_cartesian(
      xlim = x_limits
    ) +

    labs(
      title = paste0(
        "(", panel, ") ",
        title
      ),
      x = "Standardized mean difference",
      y = NULL,
      color = NULL,
      shape = NULL
    ) +

    theme_bw(base_size = 13) +

    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      panel.border = element_rect(
        color = "grey75",
        linewidth = 0.5
      ),

      # Y-axis tick labels: bold
      axis.text.y = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),

      # X-axis tick labels: bold
      axis.text.x = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),

      # X-axis title
      axis.title.x = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),

      # Figure title
      plot.title = element_text(
        size = 17,
        face = "bold",
        hjust = 0
      ),

      # Legend position
      legend.position = "inside",
      legend.position.inside = c(0.84, 0.18),
      legend.justification = c(0.5, 0.5),

      # Legend labels: larger and bold
      legend.text = element_text(
        size = 14,
        face = "bold",
        color = "black"
      ),

      # Legend symbols
      legend.key.size = grid::unit(0.7, "cm"),

      legend.background = element_blank(),
      legend.key = element_blank(),

      plot.margin = margin(8, 12, 8, 8)
    )

  return(p)
}

# ============================================================
# 7. Generate individual figures
# ============================================================

all_plots <- list()
all_balance <- list()

plot_status <- data.table(
  Panel = paste0("(", panel_labels, ")"),
  Comorbidity = comorbidity_labels,
  Status = "Not processed"
)

for (i in seq_along(comorbidity_ids)) {

  ci <- comorbidity_ids[i]

  label <- comorbidity_labels[i]

  panel <- panel_labels[i]

  cat(
    "\nProcessing (", panel, ") ",
    label, "\n",
    sep = ""
  )

  # Read diagnostic balance data
  dt <- read_balance_data(
    ci = ci,
    label = label
  )

  if (is.null(dt)) {

    plot_status[i, Status := "Skipped"]

    next
  }

  all_balance[[panel]] <- dt

  # Create plot
  p <- plot_covariate_balance(
    balance_dt = dt,
    panel = panel,
    title = label,
    threshold = smd_threshold,
    max_covariates = max_covariates
  )

  all_plots[[panel]] <- p

  # Output name, e.g. B_Hypertension
  filename <- paste0(
    cohort,
    "_",
    panel,
    "_",
    ci
  )

  # Adjust height to number of plotted entries
  n_cov <- if (is.finite(max_covariates)) {
    min(nrow(dt), max_covariates)
  } else {
    nrow(dt)
  }

  plot_height <- max(
    6,
    1.8 + 0.23 * n_cov
  )

  # Save individual PDF
  # pdf_file <- file.path(
  #   output_dir,
  #   paste0(filename, ".pdf")
  # )

  # ggsave(
  #   filename = pdf_file,
  #   plot = p,
  #   width = 10,
  #   height = plot_height,
  #   units = "in",
  #   device = grDevices::cairo_pdf,
  #   limitsize = FALSE
  # )

  # Save individual PNG
  png_file <- file.path(
    output_dir,
    paste0(filename, ".png")
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = 10,
    height = plot_height,
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )

  plot_status[i, Status := "Saved"]

  cat(
    "Saved: ",
    filename,
    ".pdf and .png\n",
    sep = ""
  )
}

# ============================================================
# 8. Save balance data and status
# ============================================================

if (length(all_balance) > 0L) {

  combined_balance <- rbindlist(
    all_balance,
    use.names = TRUE,
    fill = TRUE
  )

  fwrite(
    combined_balance,
    file.path(
      output_dir,
      "All_Comorbidities_Balance_Data.csv"
    )
  )
}

fwrite(
  plot_status,
  file.path(
    output_dir,
    "Figure_Generation_Status.csv"
  )
)

# ============================================================
# 9. Completion summary
# ============================================================

cat("\n========================================\n")
cat("Covariate-balance figure generation done\n")
cat("========================================\n")

print(plot_status)

cat(
  "\nNumber of saved figures: ",
  sum(plot_status$Status == "Saved"),
  " / ",
  length(comorbidity_ids),
  "\n",
  sep = ""
)

cat(
  "\nOutput directory:\n",
  output_dir,
  "\n"
)
