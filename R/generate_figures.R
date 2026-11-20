###############################################################################
# Reproduce manuscript and supplementary figures from Source Data.xlsx
#
# Input workbook:
#   Source Data.xlsx, located in the same directory as this script.
#
# Output:
#   generated_figures/ in the same directory. Each figure is written as a
#   vector PDF.
###############################################################################

required_packages <- c(
  "readxl", "dplyr", "tidyr", "ggplot2", "stringr", "maps",
  "patchwork", "scales"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(maps)
  library(patchwork)
  library(scales)
  library(ggtext)
})

# Locate files relative to the script so the code works from any working
# directory. When sourced interactively, fall back to the current directory.
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
f_index <- which(all_args %in% c("-f", "--file"))
script_path <- if (length(file_arg) == 1) {
  sub("^--file=", "", file_arg)
} else if (length(f_index) > 0 && f_index[1] < length(all_args)) {
  all_args[f_index[1] + 1]
} else {
  NA_character_
}
script_dir <- if (!is.na(script_path)) dirname(normalizePath(script_path)) else normalizePath(getwd())


candidate_paths <- c(
  file.path(script_dir, "Source Data.xlsx"),
  file.path(script_dir, "Figure and data", "Source Data.xlsx")
)
workbook_path <- candidate_paths[file.exists(candidate_paths)][1]


output_root <- if (basename(script_dir) == "Figure and data") dirname(script_dir) else script_dir
output_dir <- file.path(output_root, "generated_figures")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (is.na(workbook_path) || !file.exists(workbook_path)) {
  stop(
    "Source Data.xlsx was not found. Looked in:\n  ",
    paste(candidate_paths, collapse = "\n  "),
    call. = FALSE
  )
}


dry_run <- identical(Sys.getenv("GENERATE_FIGURES_DRY_RUN"), "1")


fig2_comorbidities <- c(
  "Hypertension", "Chronic_lung_disease",
  "Kidney_disease", "Liver_disease"
)

s2_comorbidities <- c(
  "Neuro_psych_disorders", "Depression", "Cerebrovascular_disease",
  "Dementia", "Diabetes",
  "Diabetes_without_complications", "Diabetes_with_complications"
)

s3_comorbidities <- c(
  "Cardiovascular", "Coronary_artery_disease", "Heart_failure",
  "Peripheral_vascular_disease", "Myocardial_infarction",
  "Cardiomyopathy", "Cancer", "Malignant_tumor", "Benign_tumor"
)

all_comorbidities <- c(
  "Hypertension", "Chronic_lung_disease", "Kidney_disease",
  "Liver_disease", "Neuro_psych_disorders", "Depression",
  "Cerebrovascular_disease", "Dementia", "Diabetes",
  "Diabetes_without_complications", "Diabetes_with_complications",
  "Cardiovascular", "Coronary_artery_disease", "Heart_failure",
  "Peripheral_vascular_disease", "Myocardial_infarction",
  "Cardiomyopathy", "Cancer", "Malignant_tumor", "Benign_tumor"
)


main_comorbidities <- c(
  fig2_comorbidities,
  "Neuro_psych_disorders", "Diabetes", "Cardiovascular", "Cancer"
)
other_comorbidities <- setdiff(all_comorbidities, main_comorbidities)

label_comorbidity <- function(x) {
  labels <- str_replace_all(x, "_", " ")
  str_to_sentence(labels)
}

cohort_colors <- c(
  "Cohort I" = "#222222",
  "Cohort II" = "#777777",
  "Cohort III" = "#2B6CB0",
  "Cohort IV" = "#B45309"
)
cohort_shapes <- c(
  "Cohort I" = 16,
  "Cohort II" = 1,
  "Cohort III" = 17,
  "Cohort IV" = 15
)
cohort_linetypes <- c(
  "Cohort I" = "solid",
  "Cohort II" = "dashed",
  "Cohort III" = "dotdash",
  "Cohort IV" = "longdash"
)

base_theme <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "grey94", color = "grey75"),
      strip.text = element_text(face = "bold"),
      legend.position = "top"
    )
}


read_figure_sheet <- function(sheet_name) {
  raw <- read_excel(workbook_path, sheet = sheet_name, col_names = FALSE)
  header_row <- which(as.character(raw[[1]]) == "Cohort")[1]
  if (is.na(header_row)) stop("Could not locate the data header in ", sheet_name, call. = FALSE)
  header <- as.character(unlist(raw[header_row, ], use.names = FALSE))
  dat <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  names(dat) <- header
  names(dat)[1:5] <- c(
    "cohort", "comorbidity", "subgroup_var", "subgroup_level", "n"
  )
  dat %>%
    mutate(
      cohort = factor(
        cohort,
        levels = c("Cohort I", "Cohort II", "Cohort III", "Cohort IV")
      ),
      comorbidity = as.character(comorbidity),
      subgroup_var = as.character(subgroup_var),
      subgroup_level = as.character(subgroup_level),
      n = as.numeric(n)
    )
}


metric_columns <- function(dat, metric) {
  escaped <- stringr::str_replace_all(metric, "([.])", "\\\\\\1")
  list(
    estimate = grep(paste0("^", escaped, ".* estimate"), names(dat), value = TRUE),
    low = grep(paste0("^", escaped, ".* 95% CI low"), names(dat), value = TRUE),
    high = grep(paste0("^", escaped, ".* 95% CI high"), names(dat), value = TRUE)
  )
}

metric_data <- function(dat, metric) {
  cols <- metric_columns(dat, metric)
  if (any(lengths(cols) != 1)) {
    stop("Could not uniquely identify columns for metric ", metric, call. = FALSE)
  }
  dat %>%
    transmute(
      cohort, comorbidity, subgroup_var, subgroup_level, n,
      metric = metric,
      estimate = as.numeric(.data[[cols$estimate]]),
      low = as.numeric(.data[[cols$low]]),
      high = as.numeric(.data[[cols$high]])
    )
}

save_figure <- function(plot, stem, width, height) {
  if (dry_run) return(invisible(NULL))
  # ggsave(
  #   file.path(output_dir, paste0(stem, ".png")), plot,
  #   width = width, height = height, units = "in", dpi = 300, bg = "white"
  # )
  ggsave(
    file.path(output_dir, paste0(stem, ".pdf")), plot,
    width = width, height = height, units = "in", device = cairo_pdf,
    bg = "white"
  )
}

###############################################################################
###############################################################################
radar_subgroups <- c("race", "UrbanRural", "age_group", "sex")
radar_legend_titles <- c(race = "Race", UrbanRural = "Area", age_group = "Age", sex = "Sex")
radar_comorbidity_order <- c(
  "Any_comorbidities", "Hypertension", "Chronic_lung_disease", "Kidney_disease",
  "Liver_disease", "Neuro_psych_disorders", "Depression", "Cerebrovascular_disease",
  "Dementia", "Diabetes", "Diabetes_without_complications", "Diabetes_with_complications",
  "Cardiovascular", "Coronary_artery_disease", "Heart_failure", "Peripheral_vascular_disease",
  "Myocardial_infarction", "Cardiomyopathy", "Cancer", "Malignant_tumor", "Benign_tumor"
)

radar_pretty_names <- function(x) {
  x <- gsub("_", "\n", x)
  x <- gsub("Cerebrovascular\ndisease", "Cerebro-\nvascular\ndisease", x)
  gsub("Cardiovascular", "Cardiov-\nascular", x)
}

radar_theta_offset <- c(
  "Liver\ndisease" = 0.10, "Neuro\npsych\ndisorders" = 0.12,
  "Depression" = 0.16, "Cerebro-\nvascular\ndisease" = 0.18,
  "Dementia" = 0.20, "Diabetes" = 0.28,
  "Diabetes\nwithout\ncomplications" = 0.22,
  "Diabetes\nwith\ncomplications" = 0.08
)
radar_group_headers <- c("Neuro\npsych\ndisorders", "Diabetes", "Cardiov-\nascular", "Cancer")

radar_palette <- function(n) {
  if (n == 1) return("#377EB8")
  if (n == 2) return(c("#377EB8", "#D7301F"))
  colorRampPalette(c("#377EB8", "#636363", "#D7301F"))(n)
}

radar_group_order <- function(dat, subgroup_name) {
  if (subgroup_name == "UrbanRural") return(c("Urban", "Rural"))
  dat %>%
    filter(cohort == "Cohort I", subgroup_var == subgroup_name) %>%
    group_by(subgroup_level) %>%
    summarise(value = mean(estimate, na.rm = TRUE), .groups = "drop") %>%
    arrange(value) %>%
    pull(subgroup_level)
}

radar_shared_scale <- function(dat, subgroup_name, group_order, metric_name) {
  selected <- dat %>%
    filter(subgroup_var == subgroup_name, subgroup_level %in% group_order)
  shift <- if (metric_name == "aAR") max(0, -min(selected$low, na.rm = TRUE)) else 0
  list(max = (max(selected$high, na.rm = TRUE) + shift) * 1.05, shift = shift)
}

draw_radar <- function(
    dat,
    subgroup_name,
    group_order,
    scale_info,
    legend_title,
    legend_cex = 1.6,
    legend_x_shift = 0.03,
    legend_y_spacing = 1.4,
    band_width_in = 0.20,
    band_height_in = 0.065
) {

  # ============================================================
  # 1. Prepare data
  # ============================================================

  d <- dat %>%
    filter(
      subgroup_var == subgroup_name,
      subgroup_level %in% group_order,
      !is.na(estimate),
      !is.na(low),
      !is.na(high)
    ) %>%
    mutate(
      variable = factor(
        radar_pretty_names(comorbidity),
        levels = radar_pretty_names(radar_comorbidity_order)
      ),
      group = factor(
        subgroup_level,
        levels = group_order
      )
    )

  wide <- function(column) {
    tapply(
      d[[column]],
      list(d$group, d$variable),
      mean
    ) + scale_info$shift
  }

  estimates <- wide("estimate")
  lower <- wide("low")
  upper <- wide("high")

  categories <- colnames(estimates)
  category_count <- length(categories)
  group_names <- rownames(estimates)
  group_count <- nrow(estimates)

  # ============================================================
  # 2. Legend labels using plotmath
  # ============================================================

  # Keep original group_names for data matching.
  # Only change the displayed legend labels.

  legend_labels <- as.expression(
    lapply(group_names, function(x) {

      x_clean <- gsub("\\s+", "", trimws(x))

      if (x_clean %in% c(">=65", "≥65")) {

        # Mathematical greater-than-or-equal symbol
        bquote("" >= 65)

      } else if (x_clean == "<65") {

        # Mathematical less-than symbol
        bquote("" < 65)

      } else {

        # Other groups remain unchanged
        x

      }
    })
  )

  # ============================================================
  # 3. Accessible colours
  # ============================================================

  radar_colors <- list(
    race = c(
      "asian" = "#0072B2",
      "white" = "black",
      "black" = "#D55E00"
    ),
    urbanicity = c(
      "urban" = "#0072B2",
      "rural" = "#D55E00"
    ),
    area = c(
      "urban" = "#0072B2",
      "rural" = "#D55E00"
    ),
    age = c(
      ">=65" = "#0072B2",
      "<65" = "#D55E00"
    ),
    sex = c(
      "male" = "#0072B2",
      "female" = "#D55E00"
    )
  )

  default_colors <- c(
    "#0072B2",
    "#D55E00",
    "black",
    "#CC79A7",
    "#009E73",
    "#E69F00"
  )

  subgroup_key <- tolower(trimws(subgroup_name))
  group_keys <- tolower(trimws(group_names))

  if (
    subgroup_key %in% names(radar_colors) &&
    all(group_keys %in% names(radar_colors[[subgroup_key]]))
  ) {

    palette <- unname(
      radar_colors[[subgroup_key]][group_keys]
    )

  } else {

    palette <- rep(
      default_colors,
      length.out = group_count
    )
  }

  # ============================================================
  # 4. Line types
  # ============================================================

  # Race:
  # Asian = solid
  # White = dotted
  # Black = dashed

  if (subgroup_key == "race") {

    race_linetypes <- c(
      "asian" = 1,
      "white" = 3,
      "black" = 2
    )

    line_types <- unname(
      race_linetypes[group_keys]
    )

    if (anyNA(line_types)) {
      warning(
        "Unmatched race labels; using default line types."
      )

      line_types <- rep(
        c(1, 2, 3, 4, 5, 6),
        length.out = group_count
      )
    }

  } else {

    line_types <- rep(
      c(1, 2, 3, 4, 5, 6),
      length.out = group_count
    )
  }

  # ============================================================
  # 5. Radar geometry
  # ============================================================

  angles <- rev(
    seq(
      0,
      2 * pi,
      length.out = category_count + 1
    )[-(category_count + 1)] + pi / 2
  )

  offsets <- radar_theta_offset[categories]

  angles <- angles +
    ifelse(is.na(offsets), 0, offsets)

  closed_angles <- c(angles, angles[1])

  to_xy <- function(values, theta) {
    cbind(
      values / scale_info$max * cos(theta),
      values / scale_info$max * sin(theta)
    )
  }

  # ============================================================
  # 6. Initialize plot
  # ============================================================

  par(
    mar = c(0, 0, 0, 0),
    xaxs = "i",
    yaxs = "i"
  )

  plot.new()

  plot.window(
    xlim = c(-1.45, 1.45),
    ylim = c(-1.45, 1.45),
    asp = 1
  )

  old_xpd <- par(xpd = NA)
  on.exit(par(old_xpd))

  # ============================================================
  # 7. Radar grid
  # ============================================================

  for (ring in 1:4) {

    grid_xy <- to_xy(
      rep(
        scale_info$max * ring / 4,
        category_count
      ),
      angles
    )

    polygon(
      grid_xy[, 1],
      grid_xy[, 2],
      border = "grey70",
      lwd = 1,
      col = NA
    )
  }

  segments(
    0,
    0,
    cos(angles),
    sin(angles),
    col = "grey70",
    lwd = 1
  )

  # ============================================================
  # 8. Confidence interval ribbons
  # ============================================================

  for (i in seq_len(group_count)) {

    high_xy <- to_xy(
      c(upper[i, ], upper[i, 1]),
      closed_angles
    )

    low_xy <- to_xy(
      c(lower[i, ], lower[i, 1]),
      closed_angles
    )

    polygon(
      c(
        high_xy[, 1],
        rev(low_xy[, 1])
      ),
      c(
        high_xy[, 2],
        rev(low_xy[, 2])
      ),
      border = NA,
      col = adjustcolor(
        palette[i],
        alpha.f = 0.20
      )
    )
  }

  # ============================================================
  # 9. Mean estimate lines
  # ============================================================

  for (i in seq_len(group_count)) {

    point_xy <- to_xy(
      c(estimates[i, ], estimates[i, 1]),
      closed_angles
    )

    lines(
      point_xy[, 1],
      point_xy[, 2],
      col = palette[i],
      lwd = 2,
      lty = line_types[i]
    )
  }

  # ============================================================
  # 10. Radial tick labels
  # ============================================================

  tick_values <- seq(
    0,
    scale_info$max,
    length.out = 5
  )

  shown_values <- tick_values - scale_info$shift

  digits <- if (
    max(abs(shown_values)) >= 10
  ) {
    0
  } else {
    3
  }

  tick_labels <- vapply(
    shown_values,
    function(value) {
      if (value == 0) {
        "0"
      } else {
        format(
          round(value, digits),
          nsmall = digits
        )
      }
    },
    character(1)
  )

  radius <- tick_values / scale_info$max - 0.10

  text(
    cos(angles[1]) * radius,
    sin(angles[1]) * radius,
    tick_labels,
    cex = 0.85,
    col = "black",
    pos = 3,
    font = 2
  )

  # ============================================================
  # 11. Comorbidity labels
  # ============================================================

  is_header <- categories %in% radar_group_headers

  text(
    1.25 * cos(angles),
    1.25 * sin(angles),
    categories,
    cex = ifelse(is_header, 1.3, 1.1),
    font = ifelse(is_header, 2, 1)
  )

  # ============================================================
  # 12. Custom legend: CI shading + mean line
  # ============================================================

  lg <- legend(
    "topright",
    inset = c(0.05, 0),
    legend = legend_labels,
    title = legend_title,
    cex = legend_cex,
    text.font = 2,
    y.intersp = legend_y_spacing,
    x.intersp = 0.6,
    seg.len = 1.2,
    pch = rep(NA_integer_, group_count),
    lty = rep(0, group_count),
    bty = "n",
    plot = FALSE
  )

  inch_to_x <- function(x) {
    grconvertX(x, from = "inches", to = "user") -
      grconvertX(0, from = "inches", to = "user")
  }

  inch_to_y <- function(y) {
    abs(
      grconvertY(y, from = "inches", to = "user") -
        grconvertY(0, from = "inches", to = "user")
    )
  }

  band_width <- inch_to_x(
    band_width_in * legend_cex
  )

  band_height <- inch_to_y(
    band_height_in * legend_cex
  )

  x_shift <- legend_x_shift * diff(par("usr")[1:2])

  text_x <- lg$text$x + x_shift
  text_y <- lg$text$y

  gap <- inch_to_x(0.04 * legend_cex)

  title_y <- max(text_y) +
    inch_to_y(0.17 * legend_cex)

  title_x <- text_x[1] - gap - band_width

  text(
    x = title_x,
    y = title_y,
    labels = legend_title,
    adj = c(0, 0.5),
    font = 2,
    cex = legend_cex,
    xpd = NA
  )

  # ============================================================
  # 13. Draw CI bands, lines, and labels
  # ============================================================

  for (i in seq_along(group_names)) {

    x_right <- text_x[i] - gap
    x_left <- x_right - band_width
    y_mid <- text_y[i]

    # Shaded 95% CI band
    rect(
      xleft = x_left,
      ybottom = y_mid - band_height / 2,
      xright = x_right,
      ytop = y_mid + band_height / 2,
      col = adjustcolor(
        palette[i],
        alpha.f = 0.20
      ),
      border = NA,
      xpd = NA
    )

    # Mean estimate line
    segments(
      x0 = x_left,
      y0 = y_mid,
      x1 = x_right,
      y1 = y_mid,
      col = palette[i],
      lwd = 2,
      lty = line_types[i],
      xpd = NA
    )

    # Legend label
    text(
      x = text_x[i],
      y = y_mid,
      labels = legend_labels[i],
      adj = c(0, 0.5),
      font = 2,
      cex = legend_cex,
      xpd = NA
    )
  }
}

generate_radar_figure <- function(sheet_name, metric_name, output_stem,
                                  cex = c(1.6, 1, 2.2),
                                  legend_cex = 1.2,
                                  width = 15, height = 24.48) {
  radar_data <- metric_data(read_figure_sheet(sheet_name), metric_name)
  settings <- lapply(radar_subgroups, function(subgroup_name) {
    order <- radar_group_order(radar_data, subgroup_name)
    list(
      order = order,
      scale = radar_shared_scale(radar_data, subgroup_name, order, metric_name)
    )
  })
  names(settings) <- radar_subgroups
  if (dry_run) return(invisible(NULL))

  draw_panel <- function() {
    layout(matrix(1:10, ncol = 2, byrow = TRUE), heights = c(0.48, rep(6, 4)))
    for (title in paste0(c("Cohort I (", "Cohort II ("), metric_name, ")")) {
      par(mar = c(0, 0, 0, 0))
      plot.new()
      text(0.5, 0.5, title, font = 2, cex = cex[1])
    }
    for (i in seq_along(radar_subgroups)) {
      for (j in 1:2) {
        subgroup_name <- radar_subgroups[i]
        cohort_name <- c("Cohort I", "Cohort II")[j]
        setting <- settings[[subgroup_name]]
        par(cex = cex[2])
        draw_radar(
          radar_data %>% filter(cohort == cohort_name), subgroup_name,
          setting$order, setting$scale, radar_legend_titles[[subgroup_name]],
          legend_cex = legend_cex
        )
        text(
          grconvertX(0.02, "npc", "user"), grconvertY(0.98, "npc", "user"),
          paste0("(", letters[2 * (i - 1) + j], ")"),
          adj = c(0, 1), font = 2, cex = cex[3], xpd = NA
        )
      }
    }
  }

  pdf(file.path(output_dir, paste0(output_stem, ".pdf")),
      width = width, height = height, bg = "white", onefile = TRUE) #210 x 297 mm = 8.27 × 11.69in
  draw_panel()
  dev.off()
  # png(file.path(output_dir, paste0(output_stem, ".png")),
  #     width = 3000, height = 4896, res = 200, bg = "white")
  # draw_panel()
  # dev.off()
}

# Figure 3: adjusted attributable risk (aAR).
generate_radar_figure("Figure 3", "aAR", "Figure_3_subpopulations",
                      cex = c(1.4, 0.49, 2.2),
                      legend_cex = 1.2,
                      width = 7.0, height = 11.0)
# Figure S4: adjusted baseline risk (aR0).
generate_radar_figure("Figure S4", "aR0", "Figure_S4_aR0_subpopulations",
                      cex = c(1.4, 0.49, 2.2),
                      legend_cex = 1.2,
                      width = 7.0, height = 11.0)
# Figure S5: adjusted exposed risk (aR1).
generate_radar_figure("Figure S5", "aR1", "Figure_S5_aR1_subpopulations",
                      cex = c(1.4, 0.49, 2.2),
                      legend_cex = 1.2,
                      width = 7.0, height = 11.0)
# Figure S6: adjusted relative risk (aRR).
generate_radar_figure("Figure S6", "aRR", "Figure_S6_aRR_subpopulations",
                      cex = c(1.4, 0.49, 2.2),
                      legend_cex = 1.2,
                      width = 7.0, height = 11.0)

###############################################################################
###############################################################################

make_multimorbidity_figure <- function(
    sheet_name,
    stem,
    comorbidity_order,
    tag_hjust = c(4.3, 6.2, 4.5, 6.2),
    tag_vjust = c(2.8, 0.8, 0.8, 0.8),
    width = NULL,
    height = 10,
    legend_point_size = 3.5,
    legend_ci_linewidth = 0.6,
    legend_ci_height = 0.8,
    point_sizes = c(1.0, 1.5, 2.1),
    ci_thresholds = c(0.05, 0.10),
    top_margin_a = 0,
    bottom_margin_d = 0
) {

  raw <- read_figure_sheet(sheet_name)


  point_shapes <- cohort_shapes
  point_shapes[point_shapes == 1] <- 21

  # ==========================================================
  # 1. Generate four metric panels
  # ==========================================================

  plots <- lapply(
    c("aR0", "aR1", "aAR", "aRR"),
    function(mt) {

      dat <- metric_data(raw, mt) %>%
        mutate(
          comorbidity = factor(
            comorbidity,
            levels = comorbidity_order
          ),
          subgroup_level = factor(
            subgroup_level,
            levels = c("0", "1", "2"),
            labels = c("0", "1", ">=2")
          )
        ) %>%
        group_by(comorbidity) %>%
        mutate(
          ci_width = pmax(high - low, 0),

          panel_range = max(high, na.rm = TRUE) -
            min(low, na.rm = TRUE),

          relative_ci_width = if_else(
            is.finite(panel_range) & panel_range > 0,
            ci_width / panel_range,
            0
          ),

          point_size = case_when(
            relative_ci_width < ci_thresholds[1] ~ point_sizes[1],
            relative_ci_width < ci_thresholds[2] ~ point_sizes[2],
            TRUE ~ point_sizes[3]
          )
        ) %>%
        ungroup()

      cohorts_in_plot <- intersect(
        names(point_shapes),
        unique(as.character(dat$cohort))
      )

      pos <- position_dodge(width = 0.45)

      # ======================================================
      # 2. Custom legend: vertical CI + point
      # ======================================================

      vertical_ci_key <- function(data, params, size) {

        key_colour <- as.character(data$colour[1])

        idx <- match(
          key_colour,
          unname(cohort_colors[cohorts_in_plot])
        )

        if (is.na(idx)) idx <- 1L

        cohort_name <- cohorts_in_plot[idx]
        key_shape <- unname(point_shapes[cohort_name])

        y_low <- (1 - legend_ci_height) / 2
        y_high <- (1 + legend_ci_height) / 2

        key_fill <- if (key_shape %in% 21:25) {
          "white"
        } else {
          key_colour
        }

        grid::grobTree(

          # Vertical CI line
          grid::segmentsGrob(
            x0 = grid::unit(0.5, "npc"),
            x1 = grid::unit(0.5, "npc"),
            y0 = grid::unit(y_low, "npc"),
            y1 = grid::unit(y_high, "npc"),
            gp = grid::gpar(
              col = key_colour,
              lwd = legend_ci_linewidth * ggplot2::.pt
            )
          ),

          # Upper and lower CI caps
          grid::segmentsGrob(
            x0 = grid::unit(c(0.35, 0.35), "npc"),
            x1 = grid::unit(c(0.65, 0.65), "npc"),
            y0 = grid::unit(c(y_low, y_high), "npc"),
            y1 = grid::unit(c(y_low, y_high), "npc"),
            gp = grid::gpar(
              col = key_colour,
              lwd = legend_ci_linewidth * ggplot2::.pt
            )
          ),

          # Point drawn last
          grid::pointsGrob(
            x = grid::unit(0.5, "npc"),
            y = grid::unit(0.5, "npc"),
            pch = key_shape,
            size = grid::unit(
              legend_point_size,
              "mm"
            ),
            gp = grid::gpar(
              col = key_colour,
              fill = key_fill,
              lwd = 0.8 * ggplot2::.pt
            )
          )
        )
      }

      # ======================================================
      # 3. Draw plot
      # ======================================================

      p <- ggplot(
        dat,
        aes(
          x = subgroup_level,
          y = estimate,
          color = cohort,
          group = cohort
        )
      )

      # Draw reference line first
      if (mt == "aRR") {
        p <- p +
          geom_hline(
            yintercept = 1,
            linetype = 2,
            color = "grey60",
            linewidth = 0.5
          )
      }

      p <- p +

        # Confidence intervals
        geom_errorbar(
          aes(ymin = low, ymax = high),
          width = 0.25,
          linewidth = 0.5,
          position = pos,
          show.legend = FALSE
        ) +

        # Automatically sized points
        geom_point(
          aes(
            shape = cohort,
            size = point_size
          ),
          stroke = 0.6,
          fill = "white",
          position = pos,
          key_glyph = vertical_ci_key,
          show.legend = TRUE
        ) +

        scale_size_identity(guide = "none") +

        facet_wrap(
          ~comorbidity,
          scales = "free_y",
          nrow = 1,
          labeller = as_labeller(
            function(x) {
              str_wrap(
                label_comorbidity(x),
                15
              )
            }
          )
        ) +

        scale_color_manual(
          values = cohort_colors,
          breaks = cohorts_in_plot,
          name = NULL
        ) +

        scale_shape_manual(
          values = point_shapes,
          breaks = cohorts_in_plot,
          guide = "none"
        ) +

        guides(
          color = guide_legend(
            override.aes = list(
              size = legend_point_size
            )
          )
        ) +

        labs(
          x = "Number of other medical conditions",
          y = switch(
            mt,
            aR0 = bquote(aR[0] ~ "(%)"),
            aR1 = bquote(aR[1] ~ "(%)"),
            aAR = "aAR (%)",
            aRR = "aRR"
          ),
          color = NULL,
          shape = NULL
        ) +

        base_theme(13) +

        theme(
          axis.text.x = element_text(size = 12),

          strip.text = element_text(
            size = 13,
            face = "bold",
            colour = "black"
          ),

          legend.position = "top",

          legend.key.height = grid::unit(
            1.1, "cm"
          ),

          legend.key.width = grid::unit(
            0.8, "cm"
          )
        )

      p
    }
  )

  # ==========================================================
  # 4. Legend only in panel (a)
  # ==========================================================

  plots[[1]] <- plots[[1]] +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      legend.key.height = grid::unit(1.1, "cm"),
      legend.key.width = grid::unit(0.8, "cm")
    )

  for (i in 2:4) {
    plots[[i]] <- plots[[i]] +
      theme(legend.position = "none")
  }

  # ==========================================================
  # 5. X-axis title only in panel (d)
  # ==========================================================

  for (i in 1:3) {
    plots[[i]] <- plots[[i]] +
      labs(x = NULL) +
      theme(axis.title.x = element_blank())
  }

  # ==========================================================
  # 6. Panel tags
  # ==========================================================

  tag_labels <- c("(a)", "(b)", "(c)", "(d)")
  tag_x <- c(0.08, 0.12, 0.08, 0.12)
  tag_y <- c(0.94, 0.98, 0.98, 0.98)

  plots <- lapply(
    seq_along(plots),
    function(i) {
      plots[[i]] +
        labs(tag = tag_labels[i]) +
        theme(
          plot.tag.position = c(
            tag_x[i],
            tag_y[i]
          ),
          plot.tag = element_text(
            size = 18,
            face = "bold",
            hjust = tag_hjust[i],
            vjust = tag_vjust[i]
          ),
          plot.margin = ggplot2::margin(
            t = 10,
            r = 5,
            b = 10,
            l = 5,
            unit = "pt"
          )
        )
    }
  )

  # ==========================================================
  # 7. Adjust ONLY the two outer margins
  # ==========================================================

  # Panel (a): reduce top margin only
  plots[[1]] <- plots[[1]] +
    theme(
      plot.margin = ggplot2::margin(
        t = top_margin_a,
        r = 5,
        b = 10,
        l = 5,
        unit = "pt"
      )
    )

  # Panel (d): reduce bottom margin only
  plots[[4]] <- plots[[4]] +
    theme(
      plot.margin = ggplot2::margin(
        t = 10,
        r = 5,
        b = bottom_margin_d,
        l = 5,
        unit = "pt"
      )
    )

  # ==========================================================
  # 8. Combine four panels
  # ==========================================================

  combined <- patchwork::wrap_plots(
    plots,
    ncol = 1
  )

  # ==========================================================
  # 9. Figure dimensions
  # ==========================================================

  if (is.null(width)) {
    width <- if (length(comorbidity_order) <= 4) {
      8.0
    } else if (length(comorbidity_order) <= 8) {
      16
    } else {
      22
    }
  }

  save_figure(
    combined,
    stem,
    width,
    height
  )

  invisible(combined)
}


# ============================================================
# Figure 2
# ============================================================

make_multimorbidity_figure(
  "Figure 2",
  "Figure_2_multimorbidity",
  fig2_comorbidities,
  tag_hjust = c(2.2, 3.1, 2.2, 3.0),
  tag_vjust = c(4.0, 0.8, 0.8, 0.6),
  width = 8.0,
  height = 10,
  legend_point_size = 3.5,
  legend_ci_linewidth = 0.6,
  legend_ci_height = 0.5,
  point_sizes = c(0.1, 0.5, 2.0),
  ci_thresholds = c(0.05, 0.10),
  top_margin_a = 0,
  bottom_margin_d = 0
)


## Figure S2 (seven comorbidities): neurological, psychiatric, diabetes.
make_multimorbidity_figure(
  "Figure S2", "Figure_S2_multimorbidity", s2_comorbidities,
  tag_hjust = c(5.0, 7.0, 5.0, 7.0),
  tag_vjust = c(4.0, 0.8, 0.8, 0.8),
  point_sizes = c(0.1, 0.5, 2.0),
  width = 18, height = 10
)

## Figure S3 (nine comorbidities): cardiovascular + cancer groups.
make_multimorbidity_figure(
  "Figure S3", "Figure_S3_multimorbidity_other", s3_comorbidities,
  tag_hjust = c(5.0, 7, 5.0, 7.0),
  tag_vjust = c(4, 0.8, 0.8, 0.8),
  point_sizes = c(0.1, 0.5, 2.0),
  width = 18, height = 10
)


###############################################################################
###############################################################################
state_region_lookup <- tibble(
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

state_polygons <- map_data("state") %>%
  as_tibble() %>%
  rename(state_map_name = region) %>%
  mutate(state_name = tools::toTitleCase(state_map_name)) %>%
  left_join(state_region_lookup, by = "state_name")

library("ggtext")
make_region_map <- function(sheet_name, metric_name, stem,
                            width = 16, height = 12,
                            size = c(5, 20, 20),
                            text = c("<span style='font-size:18pt'><b>",
                                     "<span style='font-size:18pt'>"
                            ),
                            vjust = 0.5,
                            barheight  = 1,
                            legend_breaks = ggplot2::waiver(),
                            legend_labels = ggplot2::waiver()) {
  number_format <- if (metric_name == "aRR") "%.3f" else "%.2f"

  dat <- metric_data(read_figure_sheet(sheet_name), metric_name) %>%
    mutate(
      region = subgroup_level,
      comorbidity_label = str_wrap(label_comorbidity(comorbidity), 40),
      panel = paste0(
        text[1],
        cohort, "</b></span><br>",
        text[2],
        gsub("\n", "<br>", comorbidity_label, fixed = TRUE),
        "</span>"
      ),
      label = paste0(
        sprintf(number_format, estimate),
        "\n[",
        sprintf(number_format, low),
        ", ",
        sprintf(number_format, high),
        "]"
      )
    )

  cohort_order <- intersect(
    c("Cohort I", "Cohort II", "Cohort III", "Cohort IV"),
    unique(as.character(dat$cohort))
  )
  comorbidity_order <- unique(as.character(dat$comorbidity))
  ncol_panels <- 4
  panel_levels <- unlist(
    lapply(
      split(
        comorbidity_order,
        (seq_along(comorbidity_order) - 1) %/% ncol_panels
      ),
      function(chunk_com) {
        unlist(lapply(cohort_order, function(coh) {
          dat %>%
            filter(cohort == coh, comorbidity %in% chunk_com) %>%
            distinct(comorbidity, panel) %>%
            arrange(match(comorbidity, chunk_com)) %>%
            pull(panel)
        }))
      }
    ),
    use.names = FALSE
  )
  dat <- dat %>% mutate(panel = factor(panel, levels = panel_levels))

  map_df <- state_polygons %>%
    left_join(dat, by = "region", relationship = "many-to-many") %>%
    filter(!is.na(panel))

  label_df <- map_df %>%
    group_by(region, panel) %>%
    summarise(
      long = mean(range(long, na.rm = TRUE)),
      lat = mean(range(lat, na.rm = TRUE)),
      estimate = first(estimate),
      label = first(label),
      .groups = "drop"
    ) %>%
    mutate(
      lat = case_when(
        region == "Midwest" ~ lat + 4.3,
        region == "Northeast" ~ lat - 1.0,
        TRUE ~ lat
      ),
      long = ifelse(region == "Northeast", long - 3.0, long)
    )

  p <- ggplot(map_df, aes(long, lat, group = group)) +
    geom_polygon(aes(fill = estimate), color = "white", linewidth = 0.5) +
    geom_text(
      data = label_df,
      inherit.aes = FALSE,
      aes(long, lat, label = label),
      size = size[1],
      lineheight = 1.2
    ) +
    # coord_fixed(1.5,
    #             expand = FALSE) +
    facet_wrap(~panel, ncol = 4) +
    scale_fill_gradientn(
      colors = colorRampPalette(
        c("#D7E98B", "#83C8C3", "#A6A4D4")
      )(100),
      name = switch(
        metric_name,
        aR0 = bquote(aR[0] ~ "(%) "),
        aR1 = bquote(aR[1] ~ "(%) "),
        aAR = "aAR (%) ",
        aRR = "aRR "
      ),
      breaks = legend_breaks,
      labels = legend_labels,
      guide = guide_colorbar(
        direction = "horizontal",
        title.hjust = 0.5,
        barwidth = unit(12, "cm"),
        barheight = grid::unit(barheight, "cm")  # Legend height
      )
    ) +
    theme_void(base_size = size[2]) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(
        size = size[3],
        vjust = vjust,
        margin = ggplot2::margin(t = -30, unit = "pt")
      ),
      strip.text = element_markdown(lineheight = 1.15, halign = 0.5),
      panel.spacing = unit(0.7, "lines"),
      # PDF page margins
      plot.margin = ggplot2::margin(
        t = 1,
        r = 1,
        b = 1,
        l = 1,
        unit = "mm"
      )
    )

  save_figure(p, stem, width = width, height = height)
}
###############################################################################
# Figure 4: aAR across census regions for the eight primary groups.
###############################################################################
make_region_map("Figure 4", "aAR", "Figure_4_aAR_regions",
                width = 8.1, height = 8.2,
                size = c(3.0, 14, 13),
                text = c("<span style='font-size:12pt'><b>",
                         "<span style='font-size:12pt'>"
                ),
                vjust = 0.4,
                barheight  = 0.3,
                legend_breaks = c(0.1, 0.5, 1.0, 1.6),
                legend_labels = c("0.1", "0.5", "1.0", "1.6"))

###############################################################################
# Figure 5: aRR across census regions for the eight primary groups.
###############################################################################
make_region_map("Figure 5", "aRR", "Figure_5_aRR_regions",
                width = 8.1, height = 8.2,
                size = c(3.0, 14, 13),
                text = c("<span style='font-size:12pt'><b>",
                         "<span style='font-size:12pt'>"
                ),
                vjust = 0.4,
                barheight  = 0.3,
                legend_breaks = c(1.05, 1.55, 2.05, 2.55),
                legend_labels = c("1.05", "1.55", "2.05", "2.55"))

###############################################################################
# Figure S7: two regional maps for aR0 and aR1. They are saved as two separate
###############################################################################
make_region_map("Figure S7", "aR0", "Figure_S7a_aR0_regions", 16, 12)
make_region_map("Figure S7", "aR1", "Figure_S7b_aR1_regions", 16, 12)

###############################################################################
# Figures S8-S11: regional aR0, aR1, aAR, and aRR for the remaining groups.
###############################################################################
make_region_map("Figure S8", "aR0", "Figure_S8_aR0_regions_other", 16, 17)
make_region_map("Figure S9", "aR1", "Figure_S9_aR1_regions_other", 16, 17)
make_region_map("Figure S10", "aAR", "Figure_S10_aAR_regions_other", 16, 17)
make_region_map("Figure S11", "aRR", "Figure_S11_aRR_regions_other", 16, 17)

###############################################################################
# Temporal-trend helper for Figures 6-7 and S12-S15
###############################################################################
make_temporal_figure <- function(
    sheet_name,
    metric_name,
    stem,
    cohorts,
    width = 16,
    height = 11,
    size = c(3, 16, 14, 14),
    legend_point_size = 4,
    legend_linewidth = 1.2,
    legend_ci_alpha = 0.15
) {

  dat <- metric_data(
    read_figure_sheet(sheet_name),
    metric_name
  ) %>%
    filter(as.character(cohort) %in% cohorts) %>%
    mutate(
      cohort = factor(
        as.character(cohort),
        levels = cohorts
      ),
      year = as.integer(subgroup_level),
      comorbidity = factor(
        comorbidity,
        levels = all_comorbidities
      )
    )

  reference_line <- ifelse(metric_name == "aRR", 1, 0)

  # ----------------------------------------------------------
  # Cohort styling
  # ----------------------------------------------------------

  if (length(cohorts) <= 2) {

    color_values <- c(
      "Cohort I" = "grey40",
      "Cohort II" = "grey70"
    )

    fill_values <- c(
      "Cohort I" = "grey10",
      "Cohort II" = "grey40"
    )

    # Shape 21 supports white fill
    shape_values <- c(
      "Cohort I" = 16,
      "Cohort II" = 21
    )

    linetype_values <- c(
      "Cohort I" = "solid",
      "Cohort II" = "solid"
    )

    legend_key_width <- 2.2

  } else {

    color_values <- c(
      "Cohort I" = "black",
      "Cohort II" = "#2CA02C",
      "Cohort III" = "#D55E00",
      "Cohort IV" = "#0072B2"
    )

    fill_values <- color_values

    shape_values <- c(
      "Cohort I" = 16,
      "Cohort II" = 16,
      "Cohort III" = 17,
      "Cohort IV" = 15
    )

    linetype_values <- c(
      "Cohort I" = "solid",
      "Cohort II" = "dashed",
      "Cohort III" = "dotted",
      "Cohort IV" = "longdash"
    )

    legend_key_width <- 2.8
  }

  # ----------------------------------------------------------
  # Axis label
  # ----------------------------------------------------------

  y_label <- switch(
    metric_name,
    aR0 = bquote("Adjusted" ~ R[0] ~ "(%)"),
    aR1 = bquote("Adjusted" ~ R[1] ~ "(%)"),
    aAR = "Adjusted attributable risks (%)",
    aRR = "Adjusted relative risks"
  )

  # ----------------------------------------------------------
  # Custom legend key
  # ----------------------------------------------------------

  # ggplot2 supplies one row per legend key.
  # Use the mapped colour to identify the cohort style.

  temporal_key <- function(data, params, size) {

    key_color <- as.character(data$colour[1])

    idx <- match(
      key_color,
      unname(color_values[cohorts])
    )

    if (is.na(idx)) idx <- 1L

    cohort_name <- cohorts[idx]

    key_shape <- unname(shape_values[cohort_name])
    key_linetype <- unname(linetype_values[cohort_name])
    key_fill <- unname(fill_values[cohort_name])

    # White fill for hollow circles
    point_fill <- if (key_shape %in% 21:25) {
      "white"
    } else {
      key_fill
    }

    grid::grobTree(

      # CI ribbon background
      grid::rectGrob(
        x = 0.5,
        y = 0.5,
        width = grid::unit(0.95, "npc"),
        height = grid::unit(0.65, "npc"),
        gp = grid::gpar(
          fill = scales::alpha(
            key_fill,
            legend_ci_alpha
          ),
          col = NA
        )
      ),

      # Connecting line
      grid::segmentsGrob(
        x0 = grid::unit(0.05, "npc"),
        x1 = grid::unit(0.95, "npc"),
        y0 = grid::unit(0.5, "npc"),
        y1 = grid::unit(0.5, "npc"),
        gp = grid::gpar(
          col = key_color,
          lwd = legend_linewidth * ggplot2::.pt,
          lty = key_linetype
        )
      ),

      # Point drawn AFTER line
      # White fill masks the line inside hollow circles
      grid::pointsGrob(
        x = grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        pch = key_shape,
        size = grid::unit(
          legend_point_size,
          "mm"
        ),
        gp = grid::gpar(
          col = key_color,
          fill = point_fill,
          lwd = 0.8
        )
      )
    )
  }

  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------

  p <- ggplot(
    dat,
    aes(
      x = year,
      y = estimate,
      color = cohort,
      fill = cohort,
      group = cohort,
      linetype = cohort
    )
  ) +

    # Reference line
    geom_hline(
      yintercept = reference_line,
      linetype = 2,
      color = "grey60"
    ) +

    # CI ribbon
    geom_ribbon(
      aes(
        ymin = low,
        ymax = high
      ),
      alpha = 0.15,
      color = NA,
      na.rm = TRUE,
      show.legend = FALSE
    ) +

    # Trend lines
    geom_line(
      linewidth = 0.9,
      na.rm = TRUE,
      key_glyph = temporal_key
    ) +

    # Points are drawn AFTER the lines
    # Shape 21 with white fill masks the line underneath
    geom_point(
      aes(shape = cohort),
      size = size[1],
      stroke = 0.6,
      fill = "white",
      na.rm = TRUE,
      show.legend = FALSE
    ) +

    facet_wrap(
      ~comorbidity,
      ncol = 5,
      scales = "fixed",
      labeller = as_labeller(
        function(x) {
          str_wrap(
            label_comorbidity(x),
            18
          )
        }
      )
    ) +

    scale_x_continuous(
      breaks = 2020:2024
    ) +

    scale_color_manual(
      values = color_values,
      breaks = cohorts,
      name = NULL
    ) +

    scale_fill_manual(
      values = fill_values,
      breaks = cohorts,
      guide = "none"
    ) +

    scale_shape_manual(
      values = shape_values,
      breaks = cohorts,
      guide = "none"
    ) +

    scale_linetype_manual(
      values = linetype_values,
      breaks = cohorts,
      guide = "none"
    ) +

    labs(
      x = "COVID-19 diagnosis year",
      y = y_label
    ) +

    base_theme(12) +

    theme(
      axis.title = element_text(
        size = size[2]
      ),

      axis.text = element_text(
        size = size[2]-3
      ),

      axis.text.x = element_text(
        angle = 35,
        hjust = 1
      ),

      legend.position = "top",

      legend.text = element_text(
        size = size[4]
      ),

      legend.title = element_blank(),

      legend.key = element_blank(),

      legend.background = element_blank(),

      legend.box.background = element_blank(),

      legend.key.width = grid::unit(
        legend_key_width,
        "cm"
      ),

      legend.key.height = grid::unit(
        max(0.5, legend_point_size / 10),
        "cm"
      ),

      strip.text = element_text(
        size = size[3],
        face = "bold"
      )
    )

  # ----------------------------------------------------------
  # Save
  # ----------------------------------------------------------

  save_figure(
    p,
    stem,
    width = width,
    height = height
  )

  invisible(p)
}

# ============================================================
# Figure 6: aAR temporal trends
# ============================================================

make_temporal_figure(
  "Figure 6",
  "aAR",
  "Figure_6_aAR_year",
  c("Cohort I", "Cohort II"),
  width = 8.1,
  height = 8.1,
  size = c(2.5, 12, 12, 12),
  legend_point_size = 4,
  legend_linewidth = 1.2,
  legend_ci_alpha = 0.15
)



###############################################################################
# Figure 7: temporal trends in aRR for the two primary cohorts.
###############################################################################
make_temporal_figure("Figure 7", "aRR", "Figure_7_aRR_year", c("Cohort I", "Cohort II"),
                     width = 8.1,
                     height = 8.1,
                     size = c(2.5, 12, 12, 12),
                     legend_point_size = 4,
                     legend_linewidth = 1.2,
                     legend_ci_alpha = 0.15)

###############################################################################
# Figures S12-S15: four-cohort temporal sensitivity analyses.
###############################################################################
four_cohorts <- c("Cohort I", "Cohort II", "Cohort III", "Cohort IV")
make_temporal_figure("Figure S12", "aR0", "Figure_S12_aR0_year", four_cohorts)
make_temporal_figure("Figure S13", "aR1", "Figure_S13_aR1_year", four_cohorts)
make_temporal_figure("Figure S14", "aAR", "Figure_S14_aAR_year", four_cohorts)
make_temporal_figure("Figure S15", "aRR", "Figure_S15_aRR_year", four_cohorts)

###############################################################################
# Figures S16-S19: Covariate balance for Cohorts I-IV
###############################################################################
balance_comorbidity_ids <- c(
  "Any_comorbidities", "Hypertension", "Chronic_lung_disease",
  "Kidney_disease", "Liver_disease", "Neuro_psych_disorders",
  "Depression", "Cerebrovascular_disease", "Dementia", "Diabetes",
  "Diabetes_without_complications", "Diabetes_with_complications",
  "Cardiovascular", "Coronary_artery_disease", "Heart_failure",
  "Peripheral_vascular_disease", "Myocardial_infarction",
  "Cardiomyopathy", "Cancer", "Malignant_tumor", "Benign_tumor"
)
balance_panel_letters <- LETTERS[seq_along(balance_comorbidity_ids)]
balance_figures <- c(
  "S16" = "Cohort I",
  "S17" = "Cohort II",
  "S18" = "Cohort III",
  "S19" = "Cohort IV"
)
smd_threshold <- 0.1


read_balance_sheet <- function(sheet_name) {
  raw <- suppressMessages(
    read_excel(workbook_path, sheet = sheet_name, col_names = FALSE)
  )
  header_row <- which(as.character(raw[[1]]) == "Cohort")[1]
  if (is.na(header_row)) {
    stop("Could not locate the data header in ", sheet_name, call. = FALSE)
  }
  header <- as.character(unlist(raw[header_row, ], use.names = FALSE))
  dat <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  names(dat) <- header
  dat %>%
    mutate(
      Cohort = as.character(Cohort),
      Comorbidity = as.character(Comorbidity),
      covariate = as.character(covariate),
      level = as.character(level),
      type = as.character(type),
      smd_unweighted = as.numeric(smd_unweighted),
      smd_weighted = as.numeric(smd_weighted)
    ) %>%
    mutate(
      covariate_label = ifelse(
        is.na(level) | level == "" | level == "NA",
        covariate,
        paste0(covariate, " = ", level)
      ),
      worst = pmax(abs(smd_unweighted), abs(smd_weighted), na.rm = TRUE)
    ) %>%
    filter(is.finite(smd_unweighted), is.finite(smd_weighted)) %>%
    arrange(desc(worst))
}

draw_balance_plot <- function(dat, panel_letter, title, threshold = smd_threshold) {
  long_dat <- dat %>%
    mutate(
      covariate_label = factor(
        covariate_label,
        levels = rev(unique(covariate_label))
      )
    ) %>%
    tidyr::pivot_longer(
      cols = c("smd_unweighted", "smd_weighted"),
      names_to = "Weighting",
      values_to = "SMD"
    ) %>%
    mutate(
      Weighting = factor(
        ifelse(Weighting == "smd_unweighted", "Unweighted", "Weighted"),
        levels = c("Unweighted", "Weighted")
      )
    )

  smd_range <- range(c(long_dat$SMD, -threshold, threshold), finite = TRUE)
  x_padding <- max(0.04, 0.05 * diff(smd_range))
  x_limits <- smd_range + c(-x_padding, x_padding)

  ggplot(
    long_dat,
    aes(SMD, covariate_label, shape = Weighting, color = Weighting)
  ) +
    geom_vline(
      xintercept = c(-threshold, threshold),
      linetype = "dashed", color = "grey60", linewidth = 0.6
    ) +
    geom_point(size = 2.8, stroke = 0) +
    scale_shape_manual(
      values = c("Unweighted" = 16, "Weighted" = 17)
    ) +
    scale_color_manual(
      values = c("Unweighted" = "grey50", "Weighted" = "black")
    ) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 5)) +
    coord_cartesian(xlim = x_limits) +
    labs(
      title = paste0("(", panel_letter, ") ", title),
      x = "Standardized mean difference",
      y = NULL, color = NULL, shape = NULL
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "grey75", linewidth = 0.5),
      axis.text.y = element_text(size = 10, face = "bold", color = "black"),
      axis.text.x = element_text(size = 12, face = "bold", color = "black"),
      axis.title.x = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(size = 17, face = "bold", hjust = 0),
      legend.position = "inside",
      legend.position.inside = c(0.84, 0.18),
      legend.justification = c(0.5, 0.5),
      legend.text = element_text(size = 14, face = "bold", color = "black"),
      legend.key.size = grid::unit(0.7, "cm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      plot.margin = margin(8, 12, 8, 8)
    )
}

make_balance_figure <- function(figure_code) {
  for (j in seq_along(balance_comorbidity_ids)) {
    ci <- balance_comorbidity_ids[j]
    letter <- balance_panel_letters[j]
    sheet_name <- paste0("Figure ", figure_code, letter)
    dat <- read_balance_sheet(sheet_name)
    if (nrow(dat) == 0L) next
    plot_height <- max(6, 1.8 + 0.23 * nrow(dat))
    p <- draw_balance_plot(dat, tolower(letter), label_comorbidity(ci))
    save_figure(
      p,
      paste0("Figure_", figure_code, letter, "_covariate_balance"),
      10, plot_height
    )
  }
}

for (figure_code in names(balance_figures)) make_balance_figure(figure_code)

###############################################################################
# Figures S20-S21: Propensity-score central-support intervals
###############################################################################
propensity_figures <- list(
  list(code = "S20A", sheet = "Figure S20A", cohort = "Cohort I"),
  list(code = "S20B", sheet = "Figure S20B", cohort = "Cohort II"),
  list(code = "S21A", sheet = "Figure S21A", cohort = "Cohort III"),
  list(code = "S21B", sheet = "Figure S21B", cohort = "Cohort IV")
)

read_propensity_sheet <- function(sheet_name) {
  raw <- suppressMessages(
    read_excel(workbook_path, sheet = sheet_name, col_names = FALSE)
  )
  header_row <- which(as.character(raw[[1]]) == "Cohort")[1]
  if (is.na(header_row)) {
    stop("Could not locate the data header in ", sheet_name, call. = FALSE)
  }
  header <- as.character(unlist(raw[header_row, ], use.names = FALSE))
  dat <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  names(dat) <- header
  dat %>%
    mutate(
      Cohort = as.character(Cohort),
      Comorbidity_label = as.character(Comorbidity_label),
      Treatment = factor(
        as.character(Treatment),
        levels = c("Control (D=0)", "Treated (D=1)"),
        labels = c("Unexposed (X=0)", "Exposed (X=1)")
      ),
      low = as.numeric(low),
      median = as.numeric(median),
      high = as.numeric(high)
    ) %>%
    filter(!is.na(Treatment))
}

draw_propensity_support <- function(dat, cohort_label) {
  p <- ggplot(
    dat,
    aes(
      y = reorder(Comorbidity_label, median),
      color = Treatment
    )
  ) +
    geom_linerange(
      aes(xmin = low, xmax = high),
      position = position_dodge(width = 0.55),
      linewidth = 1
    ) +
    scale_color_manual(
      values = c(
        "Unexposed (X=0)" = "gray80",
        "Exposed (X=1)" = "gray30"
      )
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = paste("Central propensity-score support -", cohort_label),
      subtitle = "Lines show 1st-99th percentiles",
      x = "Cross-fitted estimated propensity score",
      y = NULL,
      color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      legend.key.width = grid::unit(3.8, "cm")
    )
  p
}

make_propensity_figure <- function(spec) {
  dat <- read_propensity_sheet(spec$sheet)
  if (nrow(dat) == 0L) return(invisible(NULL))
  n_comorbidities <- length(unique(dat$Comorbidity_label))
  plot_height <- max(7, 0.35 * n_comorbidities)
  p <- draw_propensity_support(dat, spec$cohort)
  save_figure(
    p,
    paste0("Figure_", spec$code, "_propensity_common_support"),
    10, plot_height
  )
}

for (spec in propensity_figures) make_propensity_figure(spec)

message(
  if (dry_run) "Dry run completed: all worksheets and metric columns validated."
  else paste("Done. Figures were saved to", normalizePath(output_dir))
)
