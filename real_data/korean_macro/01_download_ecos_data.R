# Download the ECOS series used in the Korean macroeconomic application.
#
# Run from the repository root:
#   Rscript real_data/korean_macro/01_download_ecos_data.R
#
# Before running, set your ECOS API key as an environment variable:
#   Sys.setenv(ECOS_API_KEY = "your_key_here")
# or add ECOS_API_KEY=your_key_here to ~/.Renviron.

library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(httr2)
library(jsonlite)
library(lubridate)

variable_file <- file.path("real_data", "korean_macro", "ecos_variable_list.csv")
output_dir <- file.path("data", "raw")
output_file <- file.path(output_dir, "ecos_korean_macro_raw.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

api_key <- Sys.getenv("ECOS_API_KEY")
if (!nzchar(api_key)) {
  stop(
    "ECOS_API_KEY is not set. Obtain an ECOS API key from the Bank of Korea ",
    "and set it as an environment variable before running this script."
  )
}

if (!file.exists(variable_file)) {
  stop("Variable list not found: ", variable_file)
}

variables <- readr::read_csv(variable_file, show_col_types = FALSE)

format_start <- function(cycle) {
  if (cycle == "D") return("20100101")
  if (cycle == "M") return("201001")
  stop("Unsupported ECOS cycle: ", cycle)
}

format_end <- function(cycle) {
  # The paper uses quarterly observations through 2026Q1. For daily series,
  # the last available observation in January 2026 is enough to construct the
  # January monthly value used for 2026Q1.
  if (cycle == "D") return("20260131")
  if (cycle == "M") return("202601")
  stop("Unsupported ECOS cycle: ", cycle)
}

parse_ecos_date <- function(time, cycle) {
  if (cycle == "D") {
    return(as.Date(paste0(substr(time, 1, 4), "-", substr(time, 5, 6), "-", substr(time, 7, 8))))
  }
  if (cycle == "M") {
    return(as.Date(paste0(substr(time, 1, 4), "-", substr(time, 5, 6), "-01")))
  }
  as.Date(NA)
}

fetch_one_series <- function(row) {
  stat_code <- row$table_code
  cycle <- row$cycle
  start_time <- format_start(cycle)
  end_time <- format_end(cycle)

  path_parts <- c(
    "api", "StatisticSearch", api_key, "json", "kr", "1", "100000",
    stat_code, cycle, start_time, end_time, row$item_code1
  )

  if (!is.na(row$item_code2) && nzchar(row$item_code2)) {
    path_parts <- c(path_parts, row$item_code2)
  }

  url <- paste(c("https://ecos.bok.or.kr", path_parts), collapse = "/")

  message("Downloading ", row$short_code, " (", row$series, ")")
  Sys.sleep(0.35)

  resp <- request(url) |>
    req_user_agent("AFM replication package (academic use)") |>
    req_perform()

  txt <- resp_body_string(resp)
  dat <- jsonlite::fromJSON(txt, simplifyDataFrame = TRUE)

  if (is.null(dat$StatisticSearch$row)) {
    warning("No ECOS data returned for ", row$short_code, " (", row$series, ")")
    return(tibble())
  }

  out <- as_tibble(dat$StatisticSearch$row) |>
    transmute(
      short_code = row$short_code,
      series = row$series,
      table_code = row$table_code,
      item_code1 = row$item_code1,
      item_code2 = row$item_code2,
      cycle = row$cycle,
      date = parse_ecos_date(TIME, row$cycle),
      value = suppressWarnings(as.numeric(DATA_VALUE)),
      unit_name = row$unit_name,
      label_kr = row$label_kr,
      transformation = row$transformation,
      transform = row$transform
    ) |>
    filter(!is.na(date), !is.na(value)) |>
    arrange(date)

  out
}

raw_data <- purrr::map_dfr(seq_len(nrow(variables)), function(i) {
  fetch_one_series(variables[i, ])
})

if (nrow(raw_data) == 0) {
  stop("No ECOS data were downloaded. Check your API key and ECOS access.")
}

readr::write_csv(raw_data, output_file)
message("Saved raw ECOS data to: ", output_file)
