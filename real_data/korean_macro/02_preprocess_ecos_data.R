# Preprocess the ECOS data for the Korean macroeconomic dynamic AFM application.
#
# Run from the repository root after 01_download_ecos_data.R:
#   Rscript real_data/korean_macro/02_preprocess_ecos_data.R

library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(zoo)

variable_file <- file.path("real_data", "korean_macro", "ecos_variable_list.csv")
raw_file <- file.path("data", "raw", "ecos_korean_macro_raw.csv")
output_dir <- file.path("data", "processed")
output_file <- file.path(output_dir, "korean_macro_preprocessed.rds")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(variable_file)) stop("Variable list not found: ", variable_file)
if (!file.exists(raw_file)) {
  stop(
    "Raw ECOS file not found: ", raw_file, "\n",
    "Run real_data/korean_macro/01_download_ecos_data.R first."
  )
}

variables <- readr::read_csv(variable_file, show_col_types = FALSE)
raw_data <- readr::read_csv(raw_file, show_col_types = FALSE) |>
  mutate(date = as.Date(date))

first_non_missing <- function(x) {
  idx <- which(!is.na(x))
  if (length(idx) == 0) return(x[1])
  x[idx[1]]
}

last_non_missing <- function(x) {
  idx <- which(!is.na(x))
  if (length(idx) == 0) return(x[1])
  x[idx[length(idx)]]
}

impute_vector <- function(x) {
  x <- zoo::na.approx(x, na.rm = FALSE)
  x <- zoo::na.locf(x, na.rm = FALSE)
  x <- zoo::na.locf(x, fromLast = TRUE, na.rm = FALSE)
  x
}

apply_transform <- function(x, transform) {
  if (length(x) <= 1) return(numeric(0))
  if (transform == "logdiff") {
    if (any(x <= 0, na.rm = TRUE)) return(diff(x))
    return(diff(log(x)))
  }
  diff(x)
}

# Convert all series to monthly frequency. Daily series are aggregated by taking
# the last available observation within each month. Monthly series are kept as is.
monthly_long <- raw_data |>
  mutate(month_date = floor_date(date, "month")) |>
  arrange(series, date) |>
  group_by(short_code, series, month_date) |>
  summarise(
    value = if (first(cycle) == "D") last_non_missing(value) else first_non_missing(value),
    table_code = first_non_missing(table_code),
    item_code1 = first_non_missing(item_code1),
    item_code2 = first_non_missing(item_code2),
    unit_name = first_non_missing(unit_name),
    label_kr = first_non_missing(label_kr),
    transformation = first_non_missing(transformation),
    transform = first_non_missing(transform),
    .groups = "drop"
  ) |>
  rename(date = month_date) |>
  arrange(series, date)

panel_monthly <- monthly_long |>
  select(date, series, value) |>
  pivot_wider(names_from = series, values_from = value) |>
  arrange(date)

series_names <- setdiff(names(panel_monthly), "date")

# Keep the final retained variables reported in ecos_variable_list.csv.
series_names <- variables$series[variables$series %in% series_names]
panel_monthly <- panel_monthly |>
  select(date, all_of(series_names))

# Fill internal and boundary missing values series by series.
panel_monthly[series_names] <- lapply(panel_monthly[series_names], impute_vector)

# Apply the transformation specified in the variable list.
X_list <- lapply(series_names, function(s) {
  tr <- variables$transform[match(s, variables$series)]
  apply_transform(panel_monthly[[s]], tr)
})

X_mat_monthly <- do.call(cbind, X_list)
date_monthly <- panel_monthly$date[-1]
colnames(X_mat_monthly) <- series_names
rownames(X_mat_monthly) <- as.character(date_monthly)

X_scaled_monthly <- scale(X_mat_monthly)

# The empirical analysis uses quarterly observations selected from the monthly
# transformed panel: January, April, July, and October.
quarter_idx <- month(date_monthly) %in% c(1, 4, 7, 10)
X_mat_quarter_pick <- X_mat_monthly[quarter_idx, , drop = FALSE]
X_scaled_quarter_pick <- X_scaled_monthly[quarter_idx, , drop = FALSE]
date_quarter_pick <- date_monthly[quarter_idx]

trans_map <- variables |>
  filter(series %in% series_names) |>
  arrange(match(series, series_names)) |>
  transmute(
    series,
    short_code,
    short_label,
    group_code = sub("[0-9]+", "", short_code),
    group_name,
    label_kr,
    table_code,
    item_code1,
    item_name1,
    item_code2,
    item_name2,
    transform,
    unit_name
  )

preprocessed <- list(
  monthly_long = monthly_long,
  panel_monthly = panel_monthly,
  trans_map = trans_map,
  X_mat_monthly = X_mat_monthly,
  X_scaled_monthly = X_scaled_monthly,
  date_monthly = date_monthly,
  X_mat_quarter_pick = X_mat_quarter_pick,
  X_scaled_quarter_pick = X_scaled_quarter_pick,
  date_quarter_pick = date_quarter_pick,
  current_codes_raw = series_names
)

saveRDS(preprocessed, output_file)
message("Saved processed Korean macroeconomic data to: ", output_file)
message("Quarterly panel dimensions: n = ", nrow(X_scaled_quarter_pick), ", p = ", ncol(X_scaled_quarter_pick))
