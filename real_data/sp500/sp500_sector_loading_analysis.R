# S&P 500 sector loading analysis for the Bayesian approximate factor model
# Run from the repository root:
#   Rscript real_data/sp500/sp500_sector_loading_analysis.R

rm(list = ls())

output_dir <- file.path("outputs", "sp500")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c("tidyverse", "tidyquant", "lubridate", "reshape2", "viridis",
              "dfms", "rvest", "Rcpp", "scales")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

mean_Gamma_sqrtLambda <- function(draws, align = c("none","sign"), p){
  align <- match.arg(align)
  S <- length(draws$Gamma)
  G_ref <- diag(1, nrow = p, ncol = length(draws$Gamma[[1]]) / p)
  
  p2 <- nrow(G_ref)
  k  <- ncol(G_ref)
  acc <- matrix(0, p2, k)
  
  align_sign <- function(G, G0){
    s <- sign(colSums(G * G0))
    s[s == 0] <- 1
    G %*% diag(s, k)
  }
  
  for (s in seq_len(S)) {
    G <- matrix(draws$Gamma[[s]], nrow = p)
    L <- draws$Lambda[[s]]
    Ld <- if (is.null(dim(L))) diag(sqrt(as.vector(L)), k) else diag(sqrt(diag(L)), k)
    if (align == "sign") G <- align_sign(G, G_ref)
    acc <- acc + (G %*% Ld)
  }
  acc / S
}

# 1. Retrieve the S&P 500 constituent list from Wikipedia
get_sp500_structure <- function() {
  url <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies"
  
  sp500_raw <- read_html(url) %>%
    html_node("table#constituents") %>%
    html_table()
  
  sp500_clean <- sp500_raw %>%
    dplyr::select(
      symbol = Symbol,
      sector = `GICS Sector`,
      sub_industry = `GICS Sub-Industry`
    ) %>%
    mutate(symbol = str_replace(symbol, "\\.", "-"))
  
  return(sp500_clean)
}

all_sp500 <- get_sp500_structure()
target_sectors <- unique(all_sp500$sector)

selected_stocks <- all_sp500 %>%
  filter(sector %in% target_sectors) %>%
  arrange(sector, symbol)

tickers <- selected_stocks$symbol
sector_map <- selected_stocks %>%
  dplyr::select(symbol, sector)

print(sector_map)

# Download price data
stock_data <- tq_get(
  tickers,
  get  = "stock.prices",
  from = "2015-01-01",
  to   = "2023-12-31"
)

# Monthly log returns
returns_df <- stock_data %>%
  group_by(symbol) %>%
  tq_transmute(
    select     = adjusted,
    mutate_fun = periodReturn,
    period     = "monthly",
    type       = "log"
  ) %>%
  ungroup() %>%
  pivot_wider(names_from = symbol, values_from = monthly.returns) %>%
  dplyr::select(where(~ !any(is.na(.))))

X_ret <- returns_df %>%
  dplyr::select(-date) %>%
  as.matrix()

colnames(X_ret) <- returns_df %>% dplyr::select(-date) %>% colnames()

# Demean only
X_ret <- scale(X_ret, center = TRUE, scale = FALSE)

# Bai & Ng
bn_result <- dfms::ICr(X_ret, max.r = 20)
cat("Suggested number of factors (IC2) for Returns:", bn_result$r[2], "\n")
k_factors <- bn_result$r[2]

# =========================
# AFM
# =========================
Rcpp::sourceCpp(file.path("src", "static_afm_sampler.cpp"))

n <- nrow(X_ret)
p <- ncol(X_ret)

S   <- crossprod(X_ret) / n
evS <- eigen(S)
sp  <- evS$values[1:k_factors]

iter   <- 20000
burnin <- floor(iter / 2)
thin   <- max(1, iter / 2000)

nu <- 2 * p + 1
S0 <- diag(p)

res <- gibbs_factor_full_sigma_timed(
  scale(X_ret),
  k_factors,
  iter,
  burnin,
  thin,
  a_vec      = rep(2, k_factors),
  q_hyper    = 4,
  nu0        = nu,
  S0         = S0,
  b0         = 0.5,
  b1         = 2,
  evS$vectors[, 1:k_factors],
  sp,
  Sigma_init = diag(1, p),
  Z_init     = t(evS$vectors[, 1:k_factors]) %*% t(X_ret),
  verbose    = TRUE
)

# AFM loadings
afm_loadings_mat <- mean_Gamma_sqrtLambda(res, align = "sign", p = p)
colnames(afm_loadings_mat) <- paste0("V", 1:k_factors)

afm_loadings_df <- as.data.frame(afm_loadings_mat)
afm_loadings_df$symbol <- colnames(X_ret)

# =========================
# PCA
# =========================
pca_fit <- prcomp(X_ret, center = FALSE, scale. = FALSE)

# Include sqrt(eigenvalue), not only the rotation, for comparison
pca_loadings_mat <- pca_fit$rotation[, 1:k_factors, drop = FALSE] %*%
  diag(pca_fit$sdev[1:k_factors], nrow = k_factors)

colnames(pca_loadings_mat) <- paste0("V", 1:k_factors)

pca_loadings_df <- as.data.frame(pca_loadings_mat)
pca_loadings_df$symbol <- rownames(pca_loadings_mat)

# =========================
# long format
# =========================
plot_data_afm <- afm_loadings_df %>%
  pivot_longer(cols = starts_with("V"), names_to = "Factor", values_to = "Loading") %>%
  left_join(sector_map, by = "symbol") %>%
  mutate(Model = "AFM")

plot_data_pca <- pca_loadings_df %>%
  pivot_longer(cols = starts_with("V"), names_to = "Factor", values_to = "Loading") %>%
  left_join(sector_map, by = "symbol") %>%
  mutate(Model = "PCA")

# =========================
# Common ordering criterion
# =========================
# Use one reference ordering (AFM V1) for all stocks
order_df <- plot_data_afm %>%
  filter(Factor == "V1") %>%
  group_by(sector) %>%
  arrange(sector, desc(Loading), .by_group = TRUE) %>%
  ungroup()

symbol_levels <- order_df$symbol

plot_data_afm$symbol <- factor(plot_data_afm$symbol, levels = symbol_levels)
plot_data_pca$symbol <- factor(plot_data_pca$symbol, levels = symbol_levels)

# =========================
# Use a common color scale
# =========================
afm_min <- min(plot_data_afm$Loading, na.rm = TRUE)
afm_max <- max(plot_data_afm$Loading, na.rm = TRUE)

pca_min <- min(plot_data_pca$Loading, na.rm = TRUE)
pca_max <- max(plot_data_pca$Loading, na.rm = TRUE)

global_min <- min(afm_min, pca_min)
shared_max <- min(afm_max, pca_max)
global_max <- max(afm_max, pca_max)

# =========================
# plotting function
# =========================
plot_loading_heatmap_piecewise <- function(plot_data, global_min, shared_max, global_max, plot_title = "") {
  ggplot(plot_data, aes(x = Factor, y = symbol, fill = Loading)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = c(
        "#D73027",  # low: red
        "white",    # zero
        "#4575B4",  # shared positive range
        "#08306B"   # extra dark blue for AFM-only upper range
      ),
      values = scales::rescale(
        c(global_min, 0, shared_max, global_max),
        from = c(global_min, global_max)
      ),
      limits = c(global_min, global_max),
      oob = scales::squish,
      name = "Loading"
    ) +
    facet_grid(sector ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = plot_title,
      x = "Factors",
      y = "Stocks"
    ) +
    theme_minimal() +
    theme(
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank(),
      strip.text.y = element_text(angle = 0, hjust = 0, size = 10)
    )
}

# =========================
# Generate the two figures
# =========================
p_afm <- plot_loading_heatmap_piecewise(
  plot_data_afm, global_min, shared_max, global_max
)

p_pca <- plot_loading_heatmap_piecewise(
  plot_data_pca, global_min, shared_max, global_max
)

print(p_afm)
print(p_pca)

ggsave(file.path(output_dir, "sp500_afm_loading_heatmap.png"), p_afm, width = 8, height = 10, dpi = 300)
ggsave(file.path(output_dir, "sp500_pca_loading_heatmap.png"), p_pca, width = 8, height = 10, dpi = 300)
write.csv(sector_map, file.path(output_dir, "sp500_selected_sector_map.csv"), row.names = FALSE)
write.csv(afm_loadings_df, file.path(output_dir, "sp500_afm_loadings.csv"), row.names = FALSE)
write.csv(pca_loadings_df, file.path(output_dir, "sp500_pca_loadings.csv"), row.names = FALSE)
