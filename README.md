# Comorbidity Exposure-Window Definitions and Multidimensional Disparities in Long COVID Risk
This Github page provides code and data for the manuscript: Chen Y, Chen Z, Yang G, et al. [Comorbidity Exposure Windows and Multidimensional Disparities in Long COVID Risk in a US National Cohort](https://www.medrxiv.org/content/10.64898/2026.07.11.26357794v1). medRxiv, 2026: 2026.07. 11.26357794.

This paper addresses two key scientific questions: (1) whether pre-infection comorbidities alone are sufficient to assess the risk of developing long COVID (LC); and (2) which comorbidities have risk estimates that are most sensitive to the definitions of comorbidity exposure windows. To address these questions, we employ ensemble double/debiased machine learning (EnDML) for binary exposures and outcomes to estimate the heterogeneous risks of developing LC associated with preexisting comorbidities.


# Software and R-package dependencies

The code has been tested using R version 4.4.0, with the required R packages and Java/H2O dependencies installed. The tested software versions include arrow (21.0.0.1), caret (6.0-94), data.table (1.16.0), dplyr (1.1.4), ggplot2 (4.0.3), ggtext (0.2.0), glmnet (4.1-9), Hmisc (5.1-3), maps (3.4.2), mice (3.18.0), patchwork (1.3.0), readxl (1.4.3), scales (1.4.0), stringr (1.5.1), tidyr (1.3.1), Java (23.0.1, 64-bit), and H2O (3.44.0.3). A properly configured Java environment is required to support the H2O machine-learning framework.

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
