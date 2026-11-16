# Comorbidity Exposure-Window Definitions and Multidimensional Disparities in Long COVID Risk
This Github page provides code and simulated data for the manuscript: Chen Y, Chen Z, Yang G, et al. [Comorbidity Exposure Windows and Multidimensional Disparities in Long COVID Risk in a US National Cohort](https://www.medrxiv.org/content/10.64898/2026.07.11.26357794v1). medRxiv, 2026: 2026.07. 11.26357794.

This paper addresses two key scientific questions: (1) whether pre-infection comorbidities alone are sufficient to assess the risk of developing long COVID (LC); and (2) which comorbidities have risk estimates that are most sensitive to the definitions of comorbidity exposure windows. To address these questions, we employ ensemble double/debiased machine learning (EnDML) for binary exposures and outcomes to estimate the heterogeneous risks of developing LC associated with preexisting comorbidities.


# Software and R-package dependencies

The code has been tested using R version 4.4.0, with the required R packages and Java/H2O dependencies installed. The tested software versions include arrow (21.0.0.1), data.table (1.16.0), h2o (3.46.0.7), glmnet (4.1-9), caret (6.0-94), Hmisc (5.1-3), ggplot2 (4.0.3), mice (3.18.0), and Java (23.0.1, 64-bit). A properly configured Java environment is required to support the H2O machine learning framework. Compatibility with other software versions has not been verified.

## 1. Demo

Load the provided simulated dataset using the arrow package in R:

 ```r
 library(arrow)
 dat <- read_parquet("data/Simulated_data_dml.parquet")
 dim(dat)
 ```
The default `main_call.R` configuration uses the simulated parquet file and the complete-case, non-imputed testing scenario. Run:

```text
Rscript R/main_call.R
```
The script runs the selected cohort and writes results below `results/`. 


## 2. Instructions for use

To analyze a data set with the same schema as the demonstration data:

1. Place the parquet file under `data/`, or provide its path through `DML_DATA_PATH`.
2. Confirm that the data contain the COVID-19 date, Long COVID date, target comorbidity dates, covariates, and subgroup variables expected by `main_call.R`.
3. Set the cohort code in `R/main_call.R`:
   ```r
   before <- 1  # Cohort IV; outputs use the Cohort_I_* naming convention
   ```
with 
   ```text
   1 = Cohort IV
   2 = Cohort II
   3 = Cohort III
   4 = Cohort I
   ```

## 3. Expected output

The output files include predicted outcomes for the test sets and estimated propensity scores.

# Data

This folder reproduces every main-text and supplementary figure of the Long COVID manuscript from a single curated workbook.

| File | Role |
| --- | --- |
| `Simulated_data_dml.parquet` | Simulated data. |
| `Source Data.xlsx` | Curated numeric source data, one worksheet per figure. The deposited Source Data file for the paper. |
| `generate_figures.R` | Reproducibility script. Reads `Source Data.xlsx` and emits every figure as a vector PDF and a high-resolution PNG. |
| `generated_figures/` | Output directory written by `generate_figures.R`. Not shipped. |

---

## `Source Data.xlsx`

### Structure

108 worksheets in ascending order:

- **Main text — 6 sheets.** `Figure 2`, `Figure 3`, `Figure 4`, `Figure 5`, `Figure 6`, `Figure 7`.
- **Supplementary — 102 sheets.** `Figure S2` … `Figure S15`, then the lettered propensity-diagnostic blocks `Figure S16A` … `Figure S21B`.

Every worksheet uses the same five-row preamble + data layout:

| Row | Content |
| --- | --- |
| 1 | Blank spacer. |
| 2 | Title cell `"<sheet name> source data"` (e.g. `Figure 2 source data`). |
| 3 | One-line description of what the figure shows. |
| 4 | Column header row, starting with the literal cell `Cohort`. |
| 5+ | Data rows. First five columns are always `Cohort`, `Comorbidity`, `Subgroup variable`, `Subgroup level`, `n`; remaining columns hold the metric estimates and 95 % CIs for that figure (`aR0`, `aR1`, `aAR`, `aRR`, plus `*_lo` / `*_hi` where applicable). |

### Visual style

Two palettes distinguish main-text from supplementary data:

| Band | Main text (`Figure 2–7`) | Supplementary (`Figure S…`) |
| --- | --- | --- |
| Title row | fill `#1F4E79`, white bold 14 pt | fill `#595959`, white bold 14 pt |
| Subtitle row | fill `#D9E1F2`, italic `#1F3864` 11 pt | fill `#F2F2F2`, italic `#404040` 11 pt |
| Header row | fill `#2E75B6`, white bold 11 pt, bottom border | fill `#808080`, white bold 11 pt, bottom border |
| Data rows | alternating `#DEEBF7` / white, 10 pt | alternating `#EDEDED` / white, 10 pt |
| Sheet tab | default | gray `#808080` |

Each sheet has `freeze_panes = A5` so the title, subtitle and header stay
pinned while scrolling.

### Main-text sheets

| Sheet | Figure | Description |
| --- | --- | --- |
| `Figure 2` | Fig. 2 | Targeted risks by number of additional comorbidities for the four primary comorbidity groups (Hypertension, Chronic lung disease, Kidney disease, Liver disease). |
| `Figure 3` | Fig. 3 | Adjusted attributable risks and adjusted relative risks across race, sex, age, and urban/rural subpopulations (radar panels). |
| `Figure 4` | Fig. 4 | Adjusted attributable risks across U.S. census regions for the eight primary comorbidity groups. |
| `Figure 5` | Fig. 5 | Adjusted relative risks across U.S. census regions. |
| `Figure 6` | Fig. 6 | Temporal trends in adjusted attributable risks from 2020 through 2024. |
| `Figure 7` | Fig. 7 | Temporal trends in adjusted relative risks from 2020 through 2024. |

### Supplementary sheets

| Sheet(s) | Figure | Description |
| --- | --- | --- |
| `Figure S2` | S2 | Targeted risks by additional-comorbidity count for seven neurological, psychiatric, and diabetes groups. |
| `Figure S3` | S3 | Targeted risks by additional-comorbidity count for nine cardiovascular and cancer groups. |
| `Figure S4`–`S6` | S4–S6 | Adjusted baseline risk (aR0), exposed risk (aR1), and relative risk (aRR) across all subpopulations. |
| `Figure S7` | S7 | Adjusted baseline and exposed risks across U.S. census regions for the eight primary groups. |
| `Figure S8`–`S11` | S8–S11 | Regional aR0 / aR1 / aAR / aRR for the remaining twelve comorbidities. |
| `Figure S12`–`S15` | S12–S15 | Temporal trends (2020–2024) for aR0 / aR1 / aAR / aRR. |
| `Figure S16A`–`S16U` | S16 | Covariate balance under Cohort I; one panel per comorbidity (unweighted vs IPW-weighted). |
| `Figure S17A`–`S17U` | S17 | Covariate balance under Cohort II. |
| `Figure S18A`–`S18U` | S18 | Covariate balance under Cohort III. |
| `Figure S19A`–`S19U` | S19 | Covariate balance under Cohort IV. |
| `Figure S20A`, `S20B` | S20 | Central propensity-score support (1st–99th percentiles) under Cohorts I and II. |
| `Figure S21A`, `S21B` | S21 | Central propensity-score support under Cohorts III and IV. |

Figure 1 and Figure S1 are schematic diagrams with no numeric payload
and therefore have no corresponding worksheets.

---

## `generate_figures.R`

### Purpose

Rebuild every figure in the manuscript and supplement directly from
`Source Data.xlsx`, so a reviewer can verify any panel by sourcing a
single script. No access to the raw cohort microdata or ensemble
estimator is required.

### Inputs

- `Source Data.xlsx`, discovered automatically in either the script's
  own directory or a sibling `Figure and data/` folder, so the script
  runs correctly whether you source it from `Figure and data/` or from
  the `Revision_2/` parent.

### Outputs

Written to `generated_figures/` alongside the script (or in the
`Revision_2/` parent when sourced from the parent):

| PDF / PNG stem | Figure |
| --- | --- |
| `Figure_2_multimorbidity` | Fig. 2 |
| `Figure_3_subpopulations` | Fig. 3 |
| `Figure_4_aAR_regions`, `Figure_5_aRR_regions` | Figs. 4, 5 |
| `Figure_6_aAR_year`, `Figure_7_aRR_year` | Figs. 6, 7 |
| `Figure_S2_multimorbidity`, `Figure_S3_multimorbidity_other` | Figs. S2, S3 |
| `Figure_S4_aR0` … `Figure_S15_aRR_year` | Figs. S4–S15 |
| `Figure_S16_balance_cohort_I` … `Figure_S19_balance_cohort_IV` | Figs. S16–S19 |
| `Figure_S20_ps_support_IandII`, `Figure_S21_ps_support_IIIandIV` | Figs. S20, S21 |

Each figure is written as a vector PDF and a 300 dpi PNG with matching
page dimensions.

### Dependencies

R ≥ 4.2 and the CRAN packages
`readxl`, `dplyr`, `tidyr`, `ggplot2`, `stringr`, `maps`, `patchwork`,
`scales`, and `ggtext`. The script aborts with an installation hint if
any are missing.

### How it works

1. **Locate workbook and output folder.** Script resolves its own path,
   then searches for `Source Data.xlsx` in the two candidate locations
   above; creates `generated_figures/` if needed.
2. **Read a sheet (`read_figure_sheet`).** Uses `readxl::read_excel`
   with `col_names = FALSE` and locates the data header by scanning
   column A for the literal `"Cohort"`. This makes the reader robust to
   any number of preamble rows above the header, which is why the
   recent preamble cleanup (removal of `Source:` / `Selection:` lines)
   needed no R-side change.
3. **Shared transforms.** `metric_columns()` resolves estimate /
   CI-low / CI-high column names for a given metric (`aR0`, `aR1`,
   `aAR`, `aRR`) and `summarise_metric()` tidies them into long form.
4. **Figure builders.** Specialised helpers assemble each family of
   figures and share formatting, legends and colour scales:
   - `make_multimorbidity_figure()` for the multimorbidity stacked
     panels (Figs. 2, S2, S3).
   - `generate_radar_figure()` for the subpopulation radar (Fig. 3,
     and Figs. S4–S6 under different metrics).
   - `make_region_map()` for the U.S. census-region choropleths
     (Figs. 4, 5, and Figs. S7–S11).
   - `make_temporal_figure()` for the 2020–2024 year-by-year panels
     (Figs. 6, 7, and Figs. S12–S15).
   - `make_balance_figure()` and `make_ps_support_figure()` for the
     propensity diagnostics (Figs. S16–S21).
5. **Save (`save_figure`).** Writes both PDF and PNG using the figure's
   declared width / height so A4 constraints are respected.

### Running the script

From inside this folder:

```sh
Rscript generate_figures.R
```

From the `Revision_2/` parent:

```sh
Rscript "Figure and data/generate_figures.R"
```

For a quick sanity check that every sheet is still present and
well-formed without producing graphics:

```sh
GENERATE_FIGURES_DRY_RUN=1 Rscript generate_figures.R
```

### Changelog, briefly

- **R2b (current)**: a new main-text **Figure 2** (multimorbidity for
  four target comorbidities) was added; the former Figs. 2–6 shifted
  to Figs. 3–7; **Figure S2** was repopulated with seven
  neurological / psychiatric / diabetes groups and **Figure S3** with
  nine cardiovascular / cancer groups. The workbook was re-sorted into
  ascending order, the `Source:` / `Selection:` metadata lines were
  removed from every preamble, and the two-palette styling described
  above was applied.
- **R2a**: figure renumbering relative to R1 (former Figures S8–S11
  were promoted to main text, former Figures S12–S15 moved down, etc.).
  Historical comments inside `generate_figures.R` note the per-sheet
  mapping where it matters.

### Reproducibility notes

- All figures are rendered from the deposited workbook alone — the
  script never touches the raw cohort files, so results are stable
  under read-only distribution.
- Percentage metrics (`aR0`, `aR1`, `aAR`) in the workbook are stored
  on the percentage-point scale (already multiplied by 100), except
  the radar sheet (Fig. 3) which keeps the raw proportion scale used
  by the upstream `radar_charts.R`. `aRR` is always a dimensionless
  ratio.


