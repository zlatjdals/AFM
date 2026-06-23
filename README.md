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
    └── dynamic_afm/
        └── korean_macro_dynamic_afm_analysis.R
```

## R package dependencies

The scripts install missing CRAN packages automatically where possible. Main dependencies are:

```r
Rcpp, RcppArmadillo, mvtnorm, POET, dplyr, tidyr, ggplot2, tibble, readr,
tidyverse, tidyquant, zoo, lubridate, dfms, BVAR, vars, frenchdata,
rvest, scales, viridis, stringr, purrr, forcats, rstiefel, GPArotation
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

The Korean macroeconomic data are not included in this repository. They can be accessed through the Bank of Korea ECOS Open API subject to the provider's terms of use and API-key registration.

To reproduce the Korean macro dynamic AFM analysis, place your local processed data object at the following path after creating the directory locally:

```text
data/processed/preprocessed_data.RData
```

Then run:

```bash
Rscript real_data/dynamic_afm/korean_macro_dynamic_afm_analysis.R
```

The script will stop with an explicit message if the processed data file is not found.

## FAVAR analysis

Run from the repository root:

```bash
Rscript real_data/favar/favar_forecasting_analysis.R
```

This script uses public data sources through the relevant R packages and includes the proposed AFM method when `src/static_afm_sampler.cpp` is available.

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

This public version excludes private API keys, raw Korean macroeconomic data, processed Korean macroeconomic data, and Korean-data-derived output files.
