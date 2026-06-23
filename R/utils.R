# Utility functions for AFM simulation and empirical analyses

symmetrize <- function(A) (A + t(A)) / 2

spec_norm <- function(M) {
  eig <- eigen(symmetrize(M), symmetric = TRUE, only.values = TRUE)$values
  max(abs(eig))
}

fro_norm <- function(M) sqrt(sum(M * M))

max_norm <- function(M) max(abs(M))

safe_eig <- function(A) {
  A <- symmetrize(A)
  ev <- eigen(A)
  ev$values[ev$values < 0] <- pmax(ev$values[ev$values < 0], 0)
  ev
}

safe_solve <- function(A) {
  A <- symmetrize(A)
  ok <- TRUE
  R <- tryCatch(chol(A + diag(1e-8, nrow(A))), error = function(e) { ok <<- FALSE; NULL })
  if (ok) return(chol2inv(R))
  ev <- safe_eig(A)
  V <- ev$vectors
  d <- pmax(ev$values, 1e-12)
  V %*% diag(1 / d, nrow(A)) %*% t(V)
}

inv_sqrt_from_cov <- function(Sigma, eps = 1e-8) {
  Ssym <- symmetrize(Sigma)
  es <- eigen(Ssym, symmetric = TRUE)
  d <- pmax(es$values, eps)
  V <- es$vectors
  V %*% diag(1 / sqrt(d), nrow(Sigma)) %*% t(V)
}

relative_norms <- function(Sighat, Sigtrue, SigInvHalf = NULL) {
  p <- nrow(Sigtrue)
  if (is.null(SigInvHalf)) SigInvHalf <- inv_sqrt_from_cov(Sigtrue)
  Delta <- Sighat - Sigtrue
  A <- SigInvHalf %*% Delta %*% SigInvHalf
  list(
    rel_spec = spec_norm(A),
    rel_frob = fro_norm(A) / sqrt(p),
    spec = spec_norm(Delta),
    frob = fro_norm(Delta)
  )
}

mean_spike_cov <- function(draws, p, weights = NULL) {
  S <- length(draws$Gamma)
  if (is.null(weights)) weights <- rep(1 / S, S)
  weights <- weights / sum(weights)
  G1 <- matrix(draws$Gamma[[1]], nrow = p)
  k <- ncol(G1)
  acc <- matrix(0, p, p)
  for (s in seq_len(S)) {
    G <- matrix(draws$Gamma[[s]], nrow = p)
    L <- draws$Lambda[[s]]
    Ld <- if (is.null(dim(L))) diag(as.vector(L), nrow = k, ncol = k) else diag(diag(L), nrow = k, ncol = k)
    acc <- acc + weights[s] * G %*% Ld %*% t(G)
  }
  acc
}

mean_Gamma_sqrtLambda <- function(draws, align = c("none", "sign"), p) {
  align <- match.arg(align)
  S <- length(draws$Gamma)
  G_ref <- matrix(draws$Gamma[[1]], nrow = p)
  k <- ncol(G_ref)
  acc <- matrix(0, nrow = p, ncol = k)
  align_sign <- function(G, G0) {
    s <- sign(colSums(G * G0))
    s[s == 0] <- 1
    G %*% diag(s, nrow = ncol(G))
  }
  for (s in seq_len(S)) {
    G <- matrix(draws$Gamma[[s]], nrow = p)
    L <- draws$Lambda[[s]]
    Ld <- if (is.null(dim(L))) diag(sqrt(as.vector(L)), nrow = k, ncol = k) else diag(sqrt(diag(L)), nrow = k, ncol = k)
    if (align == "sign") G <- align_sign(G, G_ref)
    acc <- acc + G %*% Ld
  }
  acc / S
}

gls_scores <- function(L, Sigma_u, x, identity = FALSE) {
  if (!identity) {
    Su_inv <- safe_solve(Sigma_u)
    A <- t(L) %*% Su_inv %*% L
    safe_solve(A) %*% t(L) %*% Su_inv %*% t(x)
  } else {
    A <- t(L) %*% L
    safe_solve(A) %*% t(L) %*% t(x)
  }
}

repo_path <- function(...) {
  root <- Sys.getenv("AFM_REPO_ROOT", unset = getwd())
  file.path(root, ...)
}
