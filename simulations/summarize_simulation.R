# Summarize simulation outputs from simulations/simulation_case1_covariance_recovery.R
# Usage:
#   Rscript simulations/summarize_simulation_case1.R [input_dir] [output_csv]

rm(list = ls())
required <- c("dplyr", "tidyr")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}
source(file.path(getwd(), "R", "utils.R"))

args <- commandArgs(trailingOnly = TRUE)
base_dir <- ifelse(length(args) >= 1, args[1], file.path("outputs", "simulation", "case1"))
out_csv <- ifelse(length(args) >= 2, args[2], file.path("outputs", "simulation", "case1_summary.csv"))
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

files <- list.files(base_dir, pattern = "\\.Rdata$", full.names = TRUE)
if (length(files) == 0) stop("No .Rdata files found in: ", base_dir)

results <- list()
results_u <- list()

for (f in files) {
  load(f)  # mat_list, Sigma_u_list
  meta <- regmatches(basename(f), regexec("n([0-9]+)_p([0-9]+)_seed([0-9]+)\\.Rdata", basename(f)))[[1]]
  n_val <- if (length(meta) == 4) as.integer(meta[2]) else NA_integer_
  p_val <- if (length(meta) == 4) as.integer(meta[3]) else nrow(mat_list$True)
  seed_val <- if (length(meta) == 4) as.integer(meta[4]) else NA_integer_

  Sigtrue <- mat_list$True
  SigInvHalf <- inv_sqrt_from_cov(Sigtrue)
  for (mtd in setdiff(names(mat_list), "True")) {
    rn <- relative_norms(mat_list[[mtd]], Sigtrue, SigInvHalf)
    results[[length(results) + 1]] <- c(list(n = n_val, p = p_val, seed = seed_val, method = mtd), rn)
  }

  if (exists("Sigma_u_list") && !is.null(Sigma_u_list$True)) {
    Sigtrue_u <- Sigma_u_list$True
    for (mtd in setdiff(names(Sigma_u_list), "True")) {
      Delta_u <- Sigma_u_list[[mtd]] - Sigtrue_u
      results_u[[length(results_u) + 1]] <- list(
        n = n_val, p = p_val, seed = seed_val, method = mtd,
        frob = fro_norm(Delta_u), spec = spec_norm(Delta_u), max = max_norm(Delta_u)
      )
    }
  }
}

df <- dplyr::bind_rows(results)
summary_df <- df %>%
  dplyr::group_by(n, p, method) %>%
  dplyr::summarise(
    rel_spec_mean = mean(rel_spec, na.rm = TRUE),
    rel_spec_sd = sd(rel_spec, na.rm = TRUE),
    rel_frob_mean = mean(rel_frob, na.rm = TRUE),
    rel_frob_sd = sd(rel_frob, na.rm = TRUE),
    spec_mean = mean(spec, na.rm = TRUE),
    spec_sd = sd(spec, na.rm = TRUE),
    frob_mean = mean(frob, na.rm = TRUE),
    frob_sd = sd(frob, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_df, out_csv, row.names = FALSE)
print(summary_df, n = Inf)
message("Saved summary: ", out_csv)

if (length(results_u) > 0) {
  df_u <- dplyr::bind_rows(results_u)
  summary_df_u <- df_u %>%
    dplyr::group_by(n, p, method) %>%
    dplyr::summarise(
      frob_mean = mean(frob, na.rm = TRUE), frob_sd = sd(frob, na.rm = TRUE),
      spec_mean = mean(spec, na.rm = TRUE), spec_sd = sd(spec, na.rm = TRUE),
      max_mean = mean(max, na.rm = TRUE), max_sd = sd(max, na.rm = TRUE),
      .groups = "drop"
    )
  out_u <- sub("\\.csv$", "_idiosyncratic.csv", out_csv)
  write.csv(summary_df_u, out_u, row.names = FALSE)
  message("Saved idiosyncratic summary: ", out_u)
}
