Code and Software Submission Checklist
======================================

This archive contains the R source code and a simulated data set for the
ensemble double/debiased machine-learning (DML) analysis of comorbidity-
associated Long COVID risk. 

Required content
----------------

* Source code: `R/dml.R` contains reusable cohort, DML, H2O ensemble, and
  diagnostic functions. `R/main_call.R` loads data and runs the analyses.
  Additional R scripts create figures, tables, sensitivity analyses, and
  diagnostic plots.
* Demonstration data: `data/Simulated_data_dml.parquet` contains 100,000
  simulated records with the same column structure required by the main
  analysis. `data/Simulated_data_dml_seed.parquet` contains the smaller seed
  data set used to generate it.
* Reproducible simulation: `R/simulate_100k.R` regenerates the demonstration
  data with a fixed random seed.

1. System requirements
----------------------

Operating system

The code has been tested on Windows 11 with R 4.4.0. It is expected to work
on macOS and Linux with R 4.4.x, subject to installation of the required R
packages and Java/H2O support. No non-standard hardware is required. A
64-bit operating system, at least 128 GB of RAM, and 4 or more CPU cores are
recommended for the demonstration. The production analysis is substantially
larger and benefits from additional memory and CPU resources.

Software and R-package dependencies

Install R 4.4.x and the following packages before running the code:

```r
install.packages(c(
  "arrow", "caret", "data.table", "dplyr", "ggplot2", "glmnet",
  "Hmisc", "mice", "ROSE"
))
install.packages("h2o", repos = c("https://h2o-release.s3.amazonaws.com/h2o/latest_stable_R.html", getOption("repos")))
```

The analysis uses H2O for stacked binary ensembles. H2O requires a supported
64-bit Java runtime; follow the H2O installation instructions for the
operating system. The code was tested with the R packages available in the
R 4.4.0 environment used for development, including `arrow`, `data.table`,
`h2o`, `glmnet`, `caret`, `Hmisc`, `ggplot2`, and `mice`.

2. Installation guide
---------------------

1. Install R 4.4.x and a supported Java runtime.
2. Unzip this archive and open a terminal in the `Revised_code` directory.
3. Install the R packages listed above.
4. Confirm that R can read the demonstration file:

```r
library(arrow)
dat <- read_parquet("data/Simulated_data_dml.parquet")
dim(dat)
```

On a normal desktop computer with a stable internet connection, R and the
packages typically install in 5–15 minutes. H2O may take additional time on
the first installation.

3. Demo
------

The default `main_call.R` configuration uses the simulated parquet file and
the complete-case, non-imputed testing scenario. Run:

```text
Rscript R/main_call.R
```

On Windows, use the full path to `Rscript.exe` if it is not on `PATH`, for
example:

```text
C:\R-4.4.2\bin\Rscript.exe R/main_call.R
```

The script runs the selected cohort and writes results below `results/`.
Depending on desktop memory, CPU count, Java startup time, and H2O model
settings, the 100,000-record demonstration generally requires approximately
30–120 minutes. The first run can take longer while H2O starts.

Expected output

For the default cohort setting, the main analysis writes parquet summaries,
cross-validation results, and (when enabled) diagnostic files under:

```text
results/Cohort_I/
```

The current code uses `before = 1` for Cohort IV but saves those outputs under
`results/Cohort_I/`, as required by the manuscript's file-naming convention.
`before = 4` represents Cohort I and saves under `results/Cohort_IV/`.

The output may include files named `overall.*.parquet`, `SUMMARY.*.parquet`,
`CV.*.parquet`, and diagnostic files under `diagnostics/`. Successful console
output reports each comorbidity and either a DML estimate or a readiness
decision.

4. Instructions for use
-----------------------

To analyze a data set with the same schema as the demonstration data:

1. Place the parquet file under `data/`, or provide its path through
   `DML_DATA_PATH`.
2. Confirm that the data contain the COVID-19 date, Long COVID date, target
   comorbidity dates, covariates, and subgroup variables expected by
   `main_call.R`.
3. Set the cohort code in `R/main_call.R`:

```r
before <- 1  # Cohort IV; outputs use the Cohort_I_* naming convention
```

The cohort codes are:

```text
1 = Cohort IV
2 = Cohort II
3 = Cohort III
4 = Cohort I
```

4. Run `Rscript R/main_call.R` from `Revised_code`.

The input path and sensitivity options can also be set without editing the
script. For example:

```powershell
$env:DML_DATA_PATH="./data/my_analysis_data.parquet"
$env:DML_NON_IMPUTED_EXC_DECEASED="TRUE"
$env:DML_IMPUTED_INC_DECEASED="FALSE"
& "C:\R-4.4.0\bin\Rscript.exe" "R/main_call.R"
```

The main script expects complete covariates for the non-imputed scenario. The
imputed-data branch is available in the code but requires an
`avg_data_imputed_all.RData` object with matching columns; that file is not
needed for the default simulated-data demo.

Optional reproduction instructions
-----------------------------------

Regenerate the 100,000-record demonstration data with:

```text
Rscript R/simulate_100k.R
```

The script uses a fixed seed (`20261016`) and writes
`data/Simulated_data_dml.parquet`. The remaining scripts can be run after the
main analysis to construct manuscript tables, figures, sensitivity summaries,
covariate-balance outputs, and propensity-score diagnostics. Their input and
output conventions are described in the header comments of each script.
