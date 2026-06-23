# Bayesian approximate factor model code

This repository contains simulation and empirical-analysis code for the paper

> Bayesian Estimation of the Eigenstructure in High-Dimensional Approximate Factor Models

The code implements a Bayesian approximate factor model based on a spiked-covariance/eigenstructure parameterization, together with simulation scripts and empirical analyses.

## Repository structure

```text
.
├── R/
│   ├── utils.R
│   └── pml_bai_ng_estimator.R
├── src/
│   ├── static_afm_sampler.cpp
│   ├── dynamic_afm_sampler.cpp
│   └── gsiw_baseline_sampler.cpp
├── simulations/
│   ├── simulation_covariance_recovery.R
│   └── summarize_simulation.R
└── real_data/
    ├── favar/
    │   └── favar_forecasting_analysis.R
    ├── sp500/
    │   └── sp500_sector_loading_analysis.R
    └── korean_macro/
        ├── README.md
        ├── ecos_variable_list.csv
        ├── 01_download_ecos_data.R
        ├── 02_preprocess_ecos_data.R
        └── korean_macro_dynamic_afm_analysis.R
```

## R package dependencies

The scripts install missing CRAN packages automatically where possible. Main dependencies are:

```r
Rcpp, RcppArmadillo, mvtnorm, POET, dplyr, tidyr, ggplot2, tibble, readr,
tidyverse, tidyquant, zoo, lubridate, dfms, BVAR, vars, frenchdata,
rvest, scales, viridis, stringr, purrr, forcats, rstiefel, GPArotation,
httr2, jsonlite
```

The C++ files require `Rcpp` and `RcppArmadillo`.

## Simulation

Run one simulation replication from the repository root:

```bash
Rscript simulations/simulation_covariance_recovery.R 1 1 1 outputs/simulation/case1
```

Arguments are:

```text
seed_id p_index n_index output_dir
```

where `p_index` selects `p` from `{300, 500}` and `n_index` selects `n` from `{30, 40, 50}`.

Summarize all saved simulation files:

```bash
Rscript simulations/summarize_simulation.R outputs/simulation/case1 outputs/simulation/case1_summary.csv
```

## Korean macro dynamic AFM analysis

The Korean macroeconomic variables are obtained from the Economic Statistics System (ECOS) of the Bank of Korea. ECOS statistical data may be used, reused, and redistributed for non-commercial purposes with proper attribution to the original source, subject to the ECOS terms of use.

This repository provides the retained ECOS variable list and scripts for downloading and preprocessing the data. It does not include an ECOS API key. To reproduce the Korean macroeconomic application, obtain an ECOS API key from the Bank of Korea and set it as the environment variable `ECOS_API_KEY`.

Run from the repository root:

```bash
Rscript real_data/korean_macro/01_download_ecos_data.R
Rscript real_data/korean_macro/02_preprocess_ecos_data.R
Rscript real_data/korean_macro/korean_macro_dynamic_afm_analysis.R
```

The download script writes raw data to `data/raw/`, the preprocessing script writes the processed panel to `data/processed/`, and the analysis script writes results to `outputs/dynamic_afm/`. These generated folders are excluded from version control.

## FAVAR analysis

Run from the repository root:

```bash
Rscript real_data/favar/favar_forecasting_analysis.R
```

This script uses FRED-MD as the macroeconomic panel and the ``12 Industry Portfolios'' excess returns from Kenneth French's Data Library as target variables. It includes the proposed AFM method when `src/static_afm_sampler.cpp` is available.

## S&P 500 sector loading analysis

Run from the repository root:

```bash
Rscript real_data/sp500/sp500_sector_loading_analysis.R
```

This script downloads the S&P 500 constituent table from Wikipedia and monthly stock prices through `tidyquant`, estimates the proposed AFM loadings, compares them with PCA loadings, and saves heatmaps and loading tables under `outputs/sp500/`.

## Included comparison baselines

The repository includes additional baseline implementations used by the simulation and empirical scripts:

- `src/gsiw_baseline_sampler.cpp` for the gSIW baseline function `new_algorithm_gSIW()`
- `R/pml_bai_ng_estimator.R` for the PML baseline function `em_mm_joint()`

The scripts source these files automatically when they are present.

## Privacy and data notes

This public version excludes private API keys. Generated raw data, processed data, and output files are not committed to the repository; they are created locally by the scripts when needed.

