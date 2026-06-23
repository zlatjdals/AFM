# Simulation case 1 for Bayesian approximate factor model
# Usage:
#   Rscript simulations/simulation_case1_covariance_recovery.R <seed_id> <p_index> <n_index> [output_dir]
# Example:
#   Rscript simulations/simulation_case1_covariance_recovery.R 1 1 1 outputs/simulation/case1

rm(list = ls())

required <- c("Rcpp", "mvtnorm", "POET")
for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

source(file.path(getwd(), "R", "utils.R"))
Rcpp::sourceCpp(file.path(getwd(), "src", "static_afm_sampler.cpp"))

optional_gsiw <- file.path(getwd(), "src", "gsiw_baseline_sampler.cpp")
if (file.exists(optional_gsiw)) Rcpp::sourceCpp(optional_gsiw)

optional_pml <- file.path(getwd(), "R", "pml_bai_ng_estimator.R")
if (file.exists(optional_pml)) source(optional_pml)

args <- commandArgs(trailingOnly = TRUE)
seed_id <- ifelse(length(args) >= 1, as.integer(args[1]), 1L)
p_index <- ifelse(length(args) >= 2, as.integer(args[2]), 1L)
n_index <- ifelse(length(args) >= 3, as.integer(args[3]), 1L)
out_dir <- ifelse(length(args) >= 4, args[4], file.path("outputs", "simulation", "case1"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1000 + seed_id)
p <- c(300, 500)[p_index]
n <- c(30, 40, 50)[n_index]

c1 <- 0.2; c2 <- 0.5; c3 <- 1
lambda <- p / n / c(c1, c2, c3)
k <- length(lambda)

generate_AFM_obs <- function(n, p, k, lambda) {
  B <- matrix(NA_real_, nrow = p, ncol = k)
  for (i in seq_len(k)) {
    B[, i] <- mvtnorm::rmvnorm(1, mean = rep(0, p), sigma = diag(p))
    B[, i] <- sqrt(lambda[i]) / norm(B[, i], "2") * B[, i]
  }
  F_mat <- t(mvtnorm::rmvnorm(n, mean = rep(0, k), sigma = diag(k)))
  Q_mat <- qr.Q(qr(matrix(rnorm(p^2), p, p)))
  Sigma_u <- Q_mat %*% diag(rgamma(p, 100, 100)) %*% t(Q_mat)
  U_noise <- t(mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = Sigma_u))
  Y <- B %*% F_mat + U_noise
  Sigma <- B %*% t(B) + Sigma_u
  list(Y = Y, Sigma = Sigma, Sigma_u = Sigma_u)
}

AFM_data <- generate_AFM_obs(n, p, k, lambda)
x <- t(AFM_data$Y)
S <- crossprod(x) / n
eigen_S <- eigen(S)
spiked_eigenvalue <- eigen_S$values[seq_len(k)]
Sigma_true <- AFM_data$Sigma
Sigma_u_true <- AFM_data$Sigma_u

## Proposed AFM model
iter <- as.integer(Sys.getenv("AFM_SIM_ITER", unset = "20000"))
burnin <- floor(iter / 2)
thin <- max(1L, floor(iter / 2000))
nu <- 2 * p + 1
S0 <- diag(p)

result_AFM <- gibbs_factor_full_sigma_timed(
  x, k, iter, burnin, thin,
  a_vec = rep(2, k), q_hyper = 4, nu0 = nu, S0 = S0, b0 = 0.5, b1 = 2,
  eigen_S$vectors[, seq_len(k)], spiked_eigenvalue,
  Sigma_init = diag(1, p), Z_init = diag(1, k, n), verbose = TRUE
)
Sigma_u_afm <- matrix(Reduce("+", result_AFM$Sigma) / length(result_AFM$Sigma), nrow = p)
Sigma_afm <- mean_spike_cov(result_AFM, p) + Sigma_u_afm
rm(result_AFM)

mat_list <- list(AFM = Sigma_afm, Sample = S, True = Sigma_true)
Sigma_u_list <- list(AFM = Sigma_u_afm, True = Sigma_u_true)

## Optional gSIW baseline. Requires src/gsiw_baseline_sampler.cpp.
if (exists("new_algorithm_gSIW")) {
  result_gSIW <- new_algorithm_gSIW(p, n, x, h = 4, a_vec = rep(2, p),
                                  Gamma0 = eigen_S$vectors, iter, burnin, thin,
                                  modified = TRUE, k)
  mat_list$gSIW <- mean_spike_cov(result_gSIW, p)
  rm(result_gSIW)
} else {
  message("Skipping gSIW baseline: new_algorithm_gSIW() was not found.")
}

## SPOET baseline
shrinkage_lambda <- eigen_S$values[seq_len(k)] -
  p / (n * p - n * k - p * k) * sum(eigen_S$values[(k + 1):p])
poet <- POET::POET(t(x) - colMeans(x), k)
mat_list$SPOET <- poet$SigmaU + eigen_S$vectors[, seq_len(k)] %*%
  diag(shrinkage_lambda, k) %*% t(eigen_S$vectors[, seq_len(k)])
Sigma_u_list$SPOET <- poet$SigmaU

## Optional PML baseline. Requires R/pml_bai_ng_estimator.R.
if (exists("em_mm_joint")) {
  pml_result <- em_mm_joint(t(x), k,
                            Lambda_init = t(mvtnorm::rmvnorm(k, mean = rep(0, p), sigma = diag(p))),
                            Sigma_init = diag(p), tol = 1e-7, max_iter = 500,
                            verbose = FALSE)
  mat_list$PML <- pml_result$Lambda %*% t(pml_result$Lambda) + pml_result$Sigma
  Sigma_u_list$PML <- pml_result$Sigma
} else {
  message("Skipping PML baseline: em_mm_joint() was not found.")
}

out_file <- file.path(out_dir, sprintf("n%d_p%d_seed%d.Rdata", n, p, 1000 + seed_id))
save(mat_list, Sigma_u_list, file = out_file)
message("Saved: ", out_file)
