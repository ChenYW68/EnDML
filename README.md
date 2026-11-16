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



