# Korean macroeconomic ECOS application

This folder contains the code needed to reproduce the Korean macroeconomic dynamic AFM application in the paper.

## Data source

The Korean macroeconomic variables are obtained from the Economic Statistics System (ECOS) of the Bank of Korea. ECOS statistical data may be used, reused, and redistributed for non-commercial purposes with proper attribution to the original source, subject to the ECOS terms of use.

This repository does not include a private ECOS API key. To download the raw data, obtain an ECOS API key from the Bank of Korea and set it as an environment variable:

```r
Sys.setenv(ECOS_API_KEY = "your_key_here")
```

or add the following line to `~/.Renviron`:

```text
ECOS_API_KEY=your_key_here
```

## Files

- `ecos_variable_list.csv`: retained ECOS variables and transformations used in the paper.
- `01_download_ecos_data.R`: downloads raw ECOS series listed in `ecos_variable_list.csv`.
- `02_preprocess_ecos_data.R`: converts the raw data to the transformed monthly and quarterly panels.
- `korean_macro_dynamic_afm_analysis.R`: fits the dynamic AFM and produces the summary tables and figures.

## Reproduction order

Run the following commands from the repository root:

```bash
Rscript real_data/korean_macro/01_download_ecos_data.R
Rscript real_data/korean_macro/02_preprocess_ecos_data.R
Rscript real_data/korean_macro/korean_macro_dynamic_afm_analysis.R
```

The download script writes raw data to `data/raw/`, the preprocessing script writes the processed panel to `data/processed/`, and the analysis script writes results to `outputs/dynamic_afm/`.
