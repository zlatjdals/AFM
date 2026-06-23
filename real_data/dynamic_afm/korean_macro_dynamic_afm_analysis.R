# Dynamic AFM empirical analysis for the Korean macroeconomic panel
# Run from the repository root:
#   Rscript real_data/dynamic_afm/korean_macro_dynamic_afm_analysis.R

# ============================================================
# 0) Packages
# ============================================================
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(lubridate)
library(zoo)
library(ggplot2)
library(forcats)
library(tibble)
library(readr)
library(Rcpp)
library(dfms)

# ============================================================
# 1) Load preprocessed data
# ============================================================
preproc_file <- file.path("data", "processed", "preprocessed_data.RData")
result_dir   <- file.path("outputs", "dynamic_afm")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(preproc_file)) {
  stop(
    "Processed Korean macroeconomic data not found: ", preproc_file, "\n",
    "This public repository does not include the Korean macroeconomic dataset or ",
    "the ECOS preprocessing script because the raw-data download requires a private API key. ",
    "Place your local preprocessed_data.RData file at data/processed/ before running this script."
  )
}

load(preproc_file)

# ------------------------------------------------------------
# Choose the dataset to analyze
#   - "quarter_pick": observations from Jan 1, Apr 1, Jul 1, and Oct 1
#   - "monthly": full monthly data
# ------------------------------------------------------------
use_data <- "quarter_pick"

if (use_data == "quarter_pick") {
  X_scaled <- X_scaled_quarter_pick
  X_mat    <- X_mat_quarter_pick
  date2    <- as.Date(date_quarter_pick)
} else if (use_data == "monthly") {
  X_scaled <- X_scaled_monthly
  X_mat    <- X_mat_monthly
  date2    <- as.Date(date_monthly)
} else {
  stop("use_data must be either 'quarter_pick' or 'monthly'.")
}

current_codes_raw <- colnames(X_scaled)

if (is.null(current_codes_raw) || length(current_codes_raw) == 0) {
  stop("X_scaled must have column names.")
}

# Clean trans_map
if (!"label_short" %in% names(trans_map)) {
  trans_map$label_short <- trans_map$label_kr
}
if (!"group_name" %in% names(trans_map)) {
  trans_map$group_name <- trans_map$group_code
}

trans_map <- trans_map %>%
  filter(series %in% current_codes_raw) %>%
  arrange(match(series, current_codes_raw))

# ============================================================
# 2) Build variable mapping for paper
# ============================================================
var_map <- trans_map %>%
  filter(series %in% current_codes_raw) %>%
  arrange(group_code, label_kr) %>%
  group_by(group_code) %>%
  mutate(
    within_group_id = row_number(),
    short_code = paste0(group_code, within_group_id),
    short_label = label_short
  ) %>%
  ungroup() %>%
  dplyr::select(
    series, short_code, short_label,
    group_code, group_name,
    label_kr, table_code,
    item_code1, item_name1,
    item_code2, item_name2,
    transform, unit_name
  )

# Convert raw series IDs to short codes
short_codes <- var_map$short_code[match(current_codes_raw, var_map$series)]

# Check unmatched variables
bad_idx <- which(is.na(short_codes))
bad_codes <- current_codes_raw[bad_idx]

if (length(bad_codes) > 0) {
  cat("\n[WARN] Some variables have no short_code and will be excluded:\n")
  print(bad_codes)
}

# Keep only matched variables
keep_idx <- which(!is.na(short_codes))

current_codes_raw <- current_codes_raw[keep_idx]
short_codes <- short_codes[keep_idx]

X_scaled <- X_scaled[, keep_idx, drop = FALSE]
X_mat    <- X_mat[, keep_idx, drop = FALSE]

var_map <- var_map %>%
  filter(series %in% current_codes_raw) %>%
  arrange(match(series, current_codes_raw))

# Assign final variable names
colnames(X_scaled) <- short_codes
colnames(X_mat)    <- short_codes

if (any(is.na(short_codes))) {
  stop("Some current_codes_raw entries are not matched to var_map.")
}

colnames(X_scaled) <- short_codes
colnames(X_mat)    <- short_codes

# Save variable names for the paper
var_names  <- short_codes
var_labels <- var_map$short_label[match(var_names, var_map$short_code)]
var_groups <- var_map$group_name[match(var_names, var_map$short_code)]

cat("\n=========================================\n")
cat("Data loading and variable mapping completed.\n")
cat("use_data =", use_data, "\n")
cat("n =", nrow(X_scaled), " p =", ncol(X_scaled), "\n")
cat("=========================================\n")

write.csv(
  var_map,
  file.path(result_dir, "variable_mapping_table.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 3) Fit dynamic factor model
# ============================================================
Rcpp::sourceCpp(file.path("src", "dynamic_afm_sampler.cpp"))

if (!exists("gibbs_factor_full_sigma_ar1_diagA_timed")) {
  stop("gibbs_factor_full_sigma_ar1_diagA_timed() was not found. Check sourceCpp().")
}

n <- nrow(X_scaled)
p <- ncol(X_scaled)
k <- dfms::ICr(X_scaled, max.r = 10)$r[2]


# PCA init
sv <- svd(X_scaled)
Gamma_init  <- sv$v[, 1:k, drop = FALSE]
Z_nk        <- sv$u[, 1:k, drop = FALSE] %*% diag(sv$d[1:k], k, k)
Z_init      <- t(Z_nk)

Lambda_init <- apply(Z_init, 1, var)
Lambda_init[Lambda_init <= 1e-8] <- 1
Sigma_init  <- diag(p)

# hyperparameters
a_vec   <- rep(2.0, k)
q_hyper <- 4.0
nu0     <- 2 * p + 1
S0      <- diag(p, p)
b0      <- 0.1
b1      <- 2

fit <- gibbs_factor_full_sigma_ar1_diagA_timed(
  X = X_scaled, k = k,
  iter = 20000, burnin = 10000, thin = 10,
  a_vec = a_vec, q_hyper = q_hyper,
  nu0 = nu0, S0 = S0,
  b0 = b0, b1 = b1,
  Gamma_init = Gamma_init,
  Lambda_init = Lambda_init,
  Z_init = Z_init,
  Sigma_init = Sigma_init,
  a_ar_init = 1,
  verbose = TRUE
)

# ============================================================
# 4) Posterior post-processing
# ============================================================
N_iter <- length(fit$Gamma)
if (N_iter == 0) stop("The fit object does not contain saved MCMC draws.")

get_G <- function(x) matrix(x, nrow = p, ncol = k)
get_S <- function(x) matrix(x, nrow = p, ncol = p)
get_L <- function(x) as.numeric(x)
get_a <- function(x) as.numeric(x)

# ------------------------------------------------------------
# 4-1) sign alignment
# ------------------------------------------------------------
Gamma_ref <- get_G(fit$Gamma[[1]])
fit_Gamma_aligned <- vector("list", N_iter)
sign_mat <- matrix(1, nrow = N_iter, ncol = k)

for (s in seq_len(N_iter)) {
  G_curr <- get_G(fit$Gamma[[s]])
  
  for (j in seq_len(k)) {
    if (sum(Gamma_ref[, j] * G_curr[, j]) < 0) {
      G_curr[, j] <- -G_curr[, j]
      sign_mat[s, j] <- -1
    }
  }
  fit_Gamma_aligned[[s]] <- G_curr
}

# ------------------------------------------------------------
# 4-2) posterior means
# ------------------------------------------------------------
Gamma_bar  <- Reduce("+", fit_Gamma_aligned) / N_iter
Lambda_bar <- Reduce("+", lapply(fit$Lambda, get_L)) / N_iter
Sigma_bar  <- Reduce("+", lapply(fit$Sigma, get_S)) / N_iter

if ("a" %in% names(fit)) {
  a_bar <- Reduce("+", lapply(fit$a, get_a)) / length(fit$a)
} else {
  a_bar <- rep(NA_real_, k)
}

eigS <- eigen(Sigma_bar, symmetric = TRUE, only.values = TRUE)$values
cat("\n========================================\n")
cat("Posterior mean summary\n")
cat("========================================\n")
cat("k =", k, " / saved draws =", N_iter, "\n")
cat("Sigma_u eigenvalue range: [", round(min(eigS), 4), ", ", round(max(eigS), 4), "]\n", sep = "")
cat("Posterior mean Lambda:\n")
print(round(Lambda_bar, 4))
if (all(!is.na(a_bar))) {
  cat("Posterior mean AR(1) coefficients a:\n")
  print(round(a_bar, 4))
}
cat("========================================\n")

# ============================================================
# 5) Loading summaries
# ============================================================
Gamma_array <- array(NA_real_, dim = c(N_iter, p, k))
for (s in seq_len(N_iter)) {
  Gamma_array[s, , ] <- fit_Gamma_aligned[[s]]
}

Gamma_mean  <- apply(Gamma_array, c(2, 3), mean)
Gamma_lower <- apply(Gamma_array, c(2, 3), quantile, probs = 0.025)
Gamma_upper <- apply(Gamma_array, c(2, 3), quantile, probs = 0.975)

df_gamma <- bind_rows(lapply(seq_len(k), function(j) {
  tibble(
    var_index   = seq_len(p),
    var_name    = var_names,
    short_label = var_labels,
    group       = var_groups,
    factor      = paste0("Factor ", j),
    mean        = Gamma_mean[, j],
    lower       = Gamma_lower[, j],
    upper       = Gamma_upper[, j],
    abs_mean    = abs(Gamma_mean[, j])
  )
}))

top_loading_plot_df <- df_gamma %>%
  group_by(factor) %>%
  slice_max(order_by = abs_mean, n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(display = paste0(var_name),
         display = fct_reorder(display, mean))

p_loading_top <- ggplot(top_loading_plot_df, aes(x = display, y = mean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_pointrange(aes(ymin = lower, ymax = upper), color = "darkblue", linewidth = 0.3) +
  coord_flip() +
  facet_wrap(~ factor, scales = "free_y", ncol = 2) +
  theme_minimal() +
  labs(
    title = "Top factor loadings with 95% credible intervals",
    x = NULL,
    y = "Loading"
  )

print(p_loading_top)

top_loading_table <- bind_rows(lapply(seq_len(k), function(j) {
  tibble(
    factor = paste0("Factor ", j),
    short_code = var_names,
    short_label = var_labels,
    group = var_groups,
    loading_mean = Gamma_mean[, j],
    loading_abs = abs(Gamma_mean[, j]),
    lower = Gamma_lower[, j],
    upper = Gamma_upper[, j]
  ) %>%
    arrange(desc(loading_abs)) %>%
    slice(1:min(10, n()))
}))

write.csv(
  top_loading_table,
  file.path(result_dir, "top_loading_table.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 6) Variance decomposition
# ============================================================
component_names <- c(paste0("Factor ", seq_len(k)), "Idiosyncratic")

share_total_arr <- array(
  NA_real_,
  dim = c(N_iter, p, k + 1),
  dimnames = list(NULL, var_names, component_names)
)

share_common_arr <- array(
  NA_real_,
  dim = c(N_iter, p, k),
  dimnames = list(NULL, var_names, paste0("Factor ", seq_len(k)))
)

communality_arr <- matrix(NA_real_, nrow = N_iter, ncol = p, dimnames = list(NULL, var_names))
common_var_arr  <- matrix(NA_real_, nrow = N_iter, ncol = p, dimnames = list(NULL, var_names))
idio_var_arr    <- matrix(NA_real_, nrow = N_iter, ncol = p, dimnames = list(NULL, var_names))
total_var_arr   <- matrix(NA_real_, nrow = N_iter, ncol = p, dimnames = list(NULL, var_names))

cat("\nComputing variance decomposition...\n")

for (s in seq_len(N_iter)) {
  G  <- fit_Gamma_aligned[[s]]
  Lm <- get_L(fit$Lambda[[s]])
  Su <- get_S(fit$Sigma[[s]])
  
  fac_contrib <- sweep(G^2, 2, Lm, `*`)
  common_var  <- rowSums(fac_contrib)
  idio_var    <- pmax(diag(Su), 1e-12)
  total_var   <- common_var + idio_var
  
  share_total_arr[s, , 1:k] <- sweep(fac_contrib, 1, total_var, "/")
  share_total_arr[s, , k + 1] <- idio_var / total_var
  
  share_common_arr[s, , ] <- sweep(fac_contrib, 1, pmax(common_var, 1e-12), "/")
  
  communality_arr[s, ] <- common_var / total_var
  common_var_arr[s, ]  <- common_var
  idio_var_arr[s, ]    <- idio_var
  total_var_arr[s, ]   <- total_var
}

cat("Variance decomposition completed.\n")

share_total_mean  <- apply(share_total_arr, c(2, 3), mean)
share_total_lower <- apply(share_total_arr, c(2, 3), quantile, probs = 0.025)
share_total_upper <- apply(share_total_arr, c(2, 3), quantile, probs = 0.975)
share_common_mean <- apply(share_common_arr, c(2, 3), mean)

communality_mean  <- colMeans(communality_arr)
communality_lower <- apply(communality_arr, 2, quantile, probs = 0.025)
communality_upper <- apply(communality_arr, 2, quantile, probs = 0.975)

common_var_mean <- colMeans(common_var_arr)
idio_var_mean   <- colMeans(idio_var_arr)
total_var_mean  <- colMeans(total_var_arr)

df_comm <- tibble(
  var_name = var_names,
  short_label = var_labels,
  group = var_groups,
  communality = communality_mean,
  uniqueness = 1 - communality_mean,
  common_var = common_var_mean,
  idio_var = idio_var_mean,
  total_var = total_var_mean,
  lower = communality_lower,
  upper = communality_upper
) %>%
  arrange(desc(communality))

cat("\nTop 15 variables by communality\n")
print(
  df_comm %>%
    slice(1:min(15, n())) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
)

factor_tables <- lapply(seq_len(k), function(j) {
  tibble(
    factor = paste0("Factor ", j),
    var_name = var_names,
    short_label = var_labels,
    group = var_groups,
    share_total = share_total_mean[, j],
    share_common = share_common_mean[, j],
    communality = communality_mean
  ) %>%
    arrange(desc(share_total))
})

for (j in seq_len(k)) {
  cat("\n========================================\n")
  cat(sprintf("Factor %d: top variables by total variance share\n", j))
  cat("========================================\n")
  print(
    factor_tables[[j]] %>%
      slice(1:min(10, n())) %>%
      mutate(across(where(is.numeric), ~ round(.x, 3)))
  )
}

# ------------------------------------------------------------
# 6-1) Communality plot
# ------------------------------------------------------------
top_n_comm <- min(30, p)

df_comm_plot <- df_comm %>%
  slice(1:top_n_comm) %>%
  mutate(display = paste0(var_name),
         display = fct_reorder(display, communality))

p_comm <- ggplot(df_comm_plot, aes(x = display, y = communality)) +
  geom_pointrange(aes(ymin = lower, ymax = upper), color = "darkblue") +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "Top variables by communality",
    x = NULL,
    y = "Communality"
  )

print(p_comm)

# ------------------------------------------------------------
# 6-2) Total variance decomposition (top communality)
# ------------------------------------------------------------
top_n_stack <- min(15, p)
sel_vars <- df_comm %>% slice(1:top_n_stack) %>% pull(var_name)

df_vd_total <- map_dfr(sel_vars, function(vn) {
  i <- match(vn, var_names)
  tibble(
    var_name = vn,
    short_label = var_labels[i],
    display = paste0(vn),
    component = component_names,
    share = share_total_mean[i, ],
    lower = share_total_lower[i, ],
    upper = share_total_upper[i, ]
  )
}) %>%
  mutate(
    display = factor(display, levels = rev(unique(display))),
    component = factor(component, levels = component_names)
  )

p_vd_total <- ggplot(df_vd_total, aes(x = display, y = share, fill = component)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal() +
  labs(
    title = "Total variance decomposition for top-communality variables",
    x = NULL,
    y = "Share of total variance",
    fill = NULL
  )

print(p_vd_total)

# ------------------------------------------------------------
# 6-3) Group-level heatmap
# ------------------------------------------------------------
df_heat <- map_dfr(seq_len(k), function(j) {
  tibble(
    var_name = var_names,
    group = var_groups,
    factor = paste0("Factor ", j),
    share_common = share_common_mean[, j],
    share_total = share_total_mean[, j],
    communality = communality_mean
  )
}) %>%
  group_by(group, factor) %>%
  summarise(
    mean_common_share = mean(share_common, na.rm = TRUE),
    mean_total_share = mean(share_total, na.rm = TRUE),
    mean_communality = mean(communality, na.rm = TRUE),
    n_var = n(),
    .groups = "drop"
  )

p_heat <- ggplot(df_heat, aes(x = factor, y = fct_rev(group), fill = mean_total_share)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    title = "Group-level variance decomposition heatmap",
    x = NULL,
    y = NULL,
    fill = "Mean share"
  )

print(p_heat)

# ------------------------------------------------------------
# 6-4) Variable-level heatmap
# ------------------------------------------------------------
df_heat_var <- map_dfr(seq_len(k), function(j) {
  tibble(
    var_name = var_names,
    short_label = var_labels,
    factor = paste0("Factor ", j),
    share_total = share_total_mean[, j],
    communality = communality_mean
  )
})

top_vars_for_heat <- df_comm %>% slice(1:min(40, p)) %>% pull(var_name)

p_heat_var <- df_heat_var %>%
  filter(var_name %in% top_vars_for_heat) %>%
  mutate(display = paste0(var_name)) %>%
  group_by(var_name, display) %>%
  mutate(comm_rank = first(match(var_name, top_vars_for_heat))) %>%
  ungroup() %>%
  mutate(display = factor(display, levels = rev(unique(display[order(comm_rank)])))) %>%
  ggplot(aes(x = factor, y = display, fill = share_total)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(
    title = "Variable-level variance decomposition heatmap",
    x = NULL,
    y = NULL,
    fill = "Share"
  )

print(p_heat_var)

# ============================================================
# 7) Reconstruct latent factor paths draw-by-draw
# ============================================================
kalman_smoother_diagA <- function(X, Gamma, Lambda, Sigma, a_vec) {
  X <- as.matrix(X)
  Tn <- nrow(X)
  p  <- ncol(X)
  k  <- ncol(Gamma)
  
  A <- diag(a_vec, k, k)
  
  q_vec <- Lambda * pmax(1 - a_vec^2, 1e-8)
  Q <- diag(q_vec, k, k)
  
  m0 <- rep(0, k)
  C0 <- diag(pmax(Lambda, 1e-8), k, k)
  
  m_pred <- matrix(NA_real_, Tn, k)
  m_filt <- matrix(NA_real_, Tn, k)
  C_pred <- vector("list", Tn)
  C_filt <- vector("list", Tn)
  
  for (t in seq_len(Tn)) {
    if (t == 1) {
      a_t <- m0
      R_t <- C0
    } else {
      a_t <- A %*% m_filt[t - 1, ]
      R_t <- A %*% C_filt[[t - 1]] %*% t(A) + Q
    }
    
    y_t <- X[t, ]
    F_t <- Gamma %*% R_t %*% t(Gamma) + Sigma
    F_t <- (F_t + t(F_t)) / 2
    
    K_t <- R_t %*% t(Gamma) %*% solve(F_t)
    v_t <- y_t - as.numeric(Gamma %*% a_t)
    
    m_t <- as.numeric(a_t + K_t %*% v_t)
    C_t <- R_t - K_t %*% Gamma %*% R_t
    C_t <- (C_t + t(C_t)) / 2
    
    m_pred[t, ] <- a_t
    m_filt[t, ] <- m_t
    C_pred[[t]] <- R_t
    C_filt[[t]] <- C_t
  }
  
  Z_smooth <- matrix(NA_real_, Tn, k)
  Z_smooth[Tn, ] <- m_filt[Tn, ]
  
  if (Tn >= 2) {
    for (t in (Tn - 1):1) {
      J_t <- C_filt[[t]] %*% t(A) %*% solve(C_pred[[t + 1]])
      Z_smooth[t, ] <- m_filt[t, ] + J_t %*% (Z_smooth[t + 1, ] - m_pred[t + 1, ])
    }
  }
  
  t(Z_smooth)
}

cat("\nStarting draw-by-draw latent factor reconstruction...\n")

Tn <- nrow(X_scaled)
Z_smooth_array <- array(NA_real_, dim = c(N_iter, k, Tn))

for (s in seq_len(N_iter)) {
  Gs <- fit_Gamma_aligned[[s]]
  Ls <- get_L(fit$Lambda[[s]])
  Ss <- get_S(fit$Sigma[[s]])
  as <- if ("a" %in% names(fit)) get_a(fit$a[[s]]) else rep(0, k)
  
  Zs <- kalman_smoother_diagA(
    X = X_scaled,
    Gamma = Gs,
    Lambda = Ls,
    Sigma = Ss,
    a_vec = as
  )
  
  Z_smooth_array[s, , ] <- Zs
}

cat("Draw-by-draw latent factor reconstruction completed.\n")

Z_mean  <- apply(Z_smooth_array, c(2, 3), mean)
Z_lower <- apply(Z_smooth_array, c(2, 3), quantile, probs = 0.025)
Z_upper <- apply(Z_smooth_array, c(2, 3), quantile, probs = 0.975)

factor_dates <- date2

df_factor_path <- bind_rows(lapply(seq_len(k), function(j) {
  tibble(
    date = factor_dates,
    factor = paste0("Factor ", j),
    mean = Z_mean[j, ],
    lower = Z_lower[j, ],
    upper = Z_upper[j, ]
  )
}))

p_factor_path <- ggplot(df_factor_path, aes(x = date, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey70") +
  geom_line(linewidth = 0.5, color = "darkblue") +
  facet_wrap(~ factor, scales = "free_y", ncol = 1) +
  theme_minimal() +
  labs(
    title = "Reconstructed latent factor paths",
    x = NULL,
    y = "Factor value"
  )

print(p_factor_path)

write.csv(
  df_factor_path,
  file.path(result_dir, "factor_path_reconstructed.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 8) Factor persistence summary
# ============================================================
if ("a" %in% names(fit)) {
  a_mat <- do.call(rbind, lapply(fit$a, get_a))
  a_df <- tibble(
    factor = paste0("Factor ", seq_len(k)),
    mean  = colMeans(a_mat),
    lower = apply(a_mat, 2, quantile, probs = 0.025),
    upper = apply(a_mat, 2, quantile, probs = 0.975)
  )
  
  p_a <- ggplot(a_df, aes(x = factor, y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_pointrange(aes(ymin = lower, ymax = upper), color = "darkblue") +
    theme_minimal() +
    labs(
      title = "Posterior summary of AR(1) persistence parameters",
      x = NULL,
      y = "AR coefficient"
    )
  
  print(p_a)
  
  write.csv(
    a_df,
    file.path(result_dir, "ar_persistence_summary.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

# ============================================================
# 9) Final summary tables
# ============================================================
vd_total_table <- tibble(
  var_name = var_names,
  short_label = var_labels,
  group = var_groups,
  communality = communality_mean,
  uniqueness = 1 - communality_mean,
  common_var = common_var_mean,
  idio_var = idio_var_mean,
  total_var = total_var_mean
)

for (j in seq_len(k)) {
  vd_total_table[[paste0("Factor", j, "_share_total")]]  <- share_total_mean[, j]
  vd_total_table[[paste0("Factor", j, "_share_common")]] <- share_common_mean[, j]
}

vd_total_table$Idiosyncratic_share <- share_total_mean[, k + 1]

vd_total_table <- vd_total_table %>%
  arrange(desc(communality))

write.csv(
  vd_total_table,
  file.path(result_dir, "variance_decomposition_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

paper_factor_table <- bind_rows(lapply(seq_len(k), function(j) {
  tibble(
    factor = paste0("Factor ", j),
    short_code = var_names,
    short_label = var_labels,
    group = var_groups,
    loading_mean = Gamma_mean[, j],
    total_share = share_total_mean[, j],
    common_share = share_common_mean[, j],
    communality = communality_mean
  ) %>%
    arrange(desc(abs(loading_mean))) %>%
    slice(1:min(8, n()))
}))

write.csv(
  paper_factor_table,
  file.path(result_dir, "paper_factor_interpretation_table.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

group_summary_table <- df_comm %>%
  group_by(group) %>%
  summarise(
    n_var = n(),
    mean_communality = mean(communality, na.rm = TRUE),
    median_communality = median(communality, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_communality))

write.csv(
  group_summary_table,
  file.path(result_dir, "group_summary_table.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# 10) Save figures
# ============================================================
ggsave(
  file.path(result_dir, "plot_loading_top.png"),
  p_loading_top, width = 12, height = 8, dpi = 300, bg = "white"
)

ggsave(
  file.path(result_dir, "plot_communality.png"),
  p_comm, width = 10, height = 8, dpi = 300, bg = "white"
)

ggsave(
  file.path(result_dir, "plot_variance_decomposition_top.png"),
  p_vd_total, width = 11, height = 8, dpi = 300, bg = "white"
)

ggsave(
  file.path(result_dir, "plot_group_heatmap.png"),
  p_heat, width = 9, height = 6, dpi = 300, bg = "white"
)

ggsave(
  file.path(result_dir, "plot_variable_heatmap.png"),
  p_heat_var, width = 10, height = 10, dpi = 300, bg = "white"
)

ggsave(
  file.path(result_dir, "plot_factor_paths.png"),
  p_factor_path, width = 10, height = 9, dpi = 300, bg = "white"
)

if (exists("p_a")) {
  ggsave(
    file.path(result_dir, "plot_ar_persistence.png"),
    p_a, width = 7, height = 4, dpi = 300, bg = "white"
  )
}

# ============================================================
# 11) Final message
# ============================================================
cat("\n========================================\n")
cat("Paper-ready results were generated.\n")
cat("Generated files:\n")
cat("- variable_mapping_table.csv\n")
cat("- variance_decomposition_summary.csv\n")
cat("- top_loading_table.csv\n")
cat("- paper_factor_interpretation_table.csv\n")
cat("- group_summary_table.csv\n")
cat("- factor_path_reconstructed.csv\n")
if ("a" %in% names(fit)) cat("- ar_persistence_summary.csv\n")
cat("and the main PNG figure files.\n")
cat("========================================\n")