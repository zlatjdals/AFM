# ============================================================
# R reimplementation of Bai–Liao (2012) PML simulation + SFM init + EM+MM
# - Ports MATLAB code you provided line-by-line.
# - Orientation kept: Y is N x T (rows = cross-section, cols = time).
# - Requires only base R.
# ============================================================

# ---------- helper funcs: exact MATLAB ports ----------

likelihoodTrue <- function(Sy, P, Sigmau, Lambda, lambda) {
  A <- Lambda %*% t(Lambda) + Sigmau
  A <- (A + t(A)) / 2
  
  dec <- determinant(A, logarithm = TRUE)
  if (!is.finite(dec$modulus)) return(NA_real_)
  
  Ainv <- tryCatch(solve(A), error = function(e) matrix(NA_real_, nrow(A), ncol(A)))
  if (anyNA(Ainv)) return(NA_real_)
  
  f1 <- as.numeric(dec$modulus)
  f2 <- sum(Sy * Ainv)
  B  <- Sigmau * P
  f3 <- sum(abs(B))
  
  (f1 + f2 + lambda * f3) / ncol(Sy)
}

likelihoodlambda <- function(Sy, Sigmau, Lambda) {
  # f = log|ΛΛ' + Σu|/N + tr(Sy (ΛΛ'+Σu)^{-1})/N
  A  <- Lambda %*% t(Lambda) + Sigmau
  f1 <- log(abs(det(A)))
  f2 <- sum(Sy * solve(A))
  f  <- (f1 + f2) / ncol(Sy)
  as.numeric(f)
}

Pmatrix <- function(S, gamma) {
  N <- ncol(S)
  P <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:i) {
      if (j < i) {
        if (abs(S[i, j]) > 1e-9) P[i, j] <- abs(S[i, j])^(-gamma) else P[i, j] <- 1e9
      } else {
        P[i, j] <- 0
      }
      P[j, i] <- P[i, j]
    }
  }
  P
}

soft <- function(A, B) {
  # elementwise soft-threshold: sign(A) * max(|A| - B, 0)
  stopifnot(all(dim(A) == dim(B)))
  H <- sign(A)
  D <- abs(A) - B
  D[D < 0] <- 0
  H * D
}

# ---------- PCA init (as in your MATLAB block) ----------
pca_init <- function(Y, r) {
  # Y: N x T
  Tn <- ncol(Y)
  ev <- eigen(t(Y) %*% Y, symmetric = TRUE)
  ord <- order(ev$values)                    # ascending
  F   <- matrix(0, Tn, r)
  for (i in 1:r) {
    F[, i] <- sqrt(Tn) * ev$vectors[, ord[Tn - i + 1]]
  }
  Lambda <- Y %*% F / Tn                     # N x r
  uhat   <- Y - Lambda %*% t(F)              # N x T
  SuPCA  <- uhat %*% t(uhat) / Tn
  list(F = F, Lambda = Lambda, Su = SuPCA)
}

# ---------- SFM (diagonal ML) initial value, per your MATLAB loop ----------
sfm_diagonal_ml <- function(Y, r, tol = 1e-7, max_iter = 4000) {
  N  <- nrow(Y); Tn <- ncol(Y)
  Sy <- Y %*% t(Y) / Tn
  
  # PCA as initial value
  pca     <- pca_init(Y, r)
  Lambda  <- pca$Lambda
  Su      <- pca$Su
  Sigma1  <- diag(diag(Su))
  Sigmaold <- diag(N) * 100
  Lambda0 <- matrix(10, N, r)
  
  kk <- 1
  while ((likelihoodlambda(Sy, Sigmaold, Lambda0) - likelihoodlambda(Sy, Sigma1, Lambda) > tol) &&
         kk < max_iter) {
    Sigmaold <- Sigma1
    Lambda0  <- Lambda
    A   <- solve(Lambda0 %*% t(Lambda0) + Sigmaold)
    C   <- Sy %*% A %*% Lambda0
    Eff <- diag(r) - t(Lambda0) %*% A %*% Lambda0 + t(Lambda0) %*% A %*% C
    Lambda <- C %*% solve(Eff)                 # N x r
    M      <- Sy - Lambda %*% t(Lambda0) %*% A %*% Sy
    Sigma1 <- diag(diag(M))
    kk <- kk + 1
  }
  
  Ybar   <- matrix(rowMeans(Y), nrow = N, ncol = Tn)
  Factor <- t(solve(t(Lambda) %*% solve(Sigma1) %*% Lambda) %*%
                t(Lambda) %*% solve(Sigma1) %*% (Y - Ybar)) # T x r
  list(Lambda = Lambda, Sigma = Sigma1, Factor = Factor)
}

# ---------- EM + Majorize-Minimize joint estimation (your main paper loop) ----------
em_mm_joint <- function(Y, r, lambda = 0.08, gamma = 5, tstep = 0.1,
                        Lambda_init, Sigma_init,
                        tol = 1e-7, max_iter = 5000, verbose = FALSE) {
  N  <- nrow(Y); Tn <- ncol(Y)
  Sy <- Y %*% t(Y) / Tn
  
  Sigmaold <- Sigma_init
  Lambda0  <- Lambda_init
  
  # warm update (one step before while-loop), exactly as in MATLAB
  A   <- solve(Lambda0 %*% t(Lambda0) + Sigmaold)
  C   <- Sy %*% A %*% Lambda0
  Eff <- diag(r) - t(Lambda0) %*% A %*% Lambda0 + t(Lambda0) %*% A %*% C
  Lambda <- C %*% solve(Eff)
  Su     <- Sy - C %*% t(Lambda) - Lambda %*% t(C) + Lambda %*% Eff %*% t(Lambda)
  KML    <- Sigmaold - tstep * (solve(Sigmaold) - solve(Sigmaold) %*% Su %*% solve(Sigmaold))
  P      <- Pmatrix(Su, gamma)
  Bmat   <- lambda * tstep * P
  Sigma1 <- soft(KML, Bmat)
  Sigma1 <- (Sigma1 + t(Sigma1)) / 2
  diag(Sigma1) <- pmax(diag(Sigma1), 1e-6)
  
  kk <- 1
  while (abs( likelihoodTrue(Sy, P, Sigmaold, Lambda0, lambda) -
          likelihoodTrue(Sy, P, Sigma1,  Lambda,  lambda) ) > tol &&
         kk < max_iter) {
    Sigmaold <- Sigma1
    Lambda0  <- Lambda
    A   <- solve(Lambda0 %*% t(Lambda0) + Sigmaold)
    C   <- Sy %*% A %*% Lambda0
    Eff <- diag(r) - t(Lambda0) %*% A %*% Lambda0 + t(Lambda0) %*% A %*% C
    Lambda <- C %*% solve(Eff)
    Su     <- Sy - C %*% t(Lambda) - Lambda %*% t(C) + Lambda %*% Eff %*% t(Lambda)
    
    KML    <- Sigmaold - tstep * (solve(Sigmaold) - solve(Sigmaold) %*% Su %*% solve(Sigmaold))
    P      <- Pmatrix(Su, gamma)
    Bmat   <- lambda * tstep * P
    Sigma1 <- soft(KML, Bmat)
    
    kk <- kk + 1
    if (verbose && kk %% 50 == 0) message(sprintf("iter %d", kk))
  }
  message(paste0("error : ", likelihoodTrue(Sy, P, Sigmaold, Lambda0, lambda) -
                    likelihoodTrue(Sy, P, Sigma1,  Lambda,  lambda)))
  # keep "old" ones after stopping criterion (as in MATLAB)
  Lambda_hat <- Lambda0
  Sigma_hat  <- Sigmaold
  
  Ybar <- matrix(rowMeans(Y), nrow = N, ncol = Tn)
  Factor <- t(solve(t(Lambda_hat) %*% solve(Sigma_hat) %*% Lambda_hat) %*%
                t(Lambda_hat) %*% solve(Sigma_hat) %*% (Y - Ybar))  # T x r
  
  list(Lambda = Lambda_hat, Sigma = Sigma_hat, Factor = Factor, iters = kk)
}

# ---------------- demo (smaller reps) ----------------
# res <- run_simulation(r = 2, Tn = 100, N = 150, nrep = 10,
#                       lambda = 0.08, gamma = 5, tstep = 0.1, seed = 123)
# print(res)

# After fitting, the model-implied covariance for demeaned Y is:
# Sigma_y_hat <- Lambda_hat %*% t(Lambda_hat) + Sigma_hat
