# FAVAR empirical analysis
# Run from the repository root:
#   Rscript real_data/favar/favar_forecasting_analysis.R

## ============================== CLEAN START ==============================
rm(list = ls())

## ============================== CONFIG ==================================
cfg <- list(
  start_year  = 2018,      # start year
  start_month = 1,         # start month (1--12)
  n_months    = 48,       # number of months to use; NA uses all available months
  r_max       = 12,        # maximum number of factors considered by Bai--Ng
  use_wpc     = FALSE,     # TRUE uses SPOET-based WPC factors; FALSE uses PC factors
  k_factors   = NULL,      # NULL selects k by Bai--Ng; otherwise use a fixed k
  var_h       = 12,        # forecast horizon H
  var_lagmax  = 12,        # maximum lag for VAR lag selection
  bvar_draw   = 8000L,
  bvar_burn   = 2000L,
  bvar_thin   = 5L,
  # Optional external code paths; they are sourced automatically when available
  path_cpp_static_afm = file.path("src", "static_afm_sampler.cpp"),
  # Optional baseline files. Put them in these paths if you want to reproduce
  # the gSIW and PML comparisons. They are not required for the proposed AFM method.
  path_cpp_gsiw     = file.path("src", "gsiw_baseline_sampler.cpp"),
  path_r_pml       = file.path("R", "pml_bai_ng_estimator.R")
)

## ============================== PACKAGES =================================
need <- c(
  "tidyverse","zoo","lubridate","dfms","POET","BVAR","vars",
  "frenchdata","purrr","mvtnorm","GPArotation","rstiefel"
)
for(p in need) if(!requireNamespace(p, quietly = TRUE)) install.packages(p)
invisible(lapply(need, library, character.only = TRUE))

## ============================== UTILS ====================================
symmetrize <- function(A) (A + t(A))/2
safe_eig <- function(A){ A <- symmetrize(A); ev <- eigen(A); ev$values[ev$values < 0] <- pmax(ev$values[ev$values < 0], 0); ev }
sym_sqrtm <- function(A, inverse = FALSE) {
  ev <- safe_eig(A); vals <- pmax(ev$values, 1e-12)
  if (inverse) ev$vectors %*% diag(1/sqrt(vals),nrow(A)) %*% t(ev$vectors)
  else         ev$vectors %*% diag(   sqrt(vals),nrow(A)) %*% t(ev$vectors)
}
safe_solve <- function(A){
  A <- symmetrize(A); ok <- TRUE
  R <- tryCatch(chol(A + diag(1e-8, nrow(A))), error = function(e){ok <<- FALSE; NULL})
  if (ok) return(chol2inv(R))
  ev <- safe_eig(A); V <- ev$vectors; d <- pmax(ev$values, 1e-12)
  V %*% diag(1/d,nrow(A)) %*% t(V)
}

qnm <- function(q){
  s <- as.character(q); s2 <- paste0(100*q, "%")
  if (exists("bands", inherits = TRUE)) {
    if (s  %in% bands) return(s)
    if (s2 %in% bands) return(s2)
  }
  s2
}

gls_scores <- function(L, Sigma_u, x, identity = FALSE){  # x: n x p, L: p x k
  if(!identity){
    Su_inv <- safe_solve(Sigma_u)
    A <- t(L) %*% Su_inv %*% L
    safe_solve(A) %*% t(L) %*% Su_inv %*% t(x)   # k x n
  }else{
    A <- t(L) %*%  L
    safe_solve(A) %*% t(L) %*%  t(x)   # k x n
  }
}
## ============================== 1) DATA LOADER ===========================
load_panels <- function(cfg){
  # FRED-MD data included in the BVAR package
  data("fred_md", package = "BVAR")
  fred_md <- as.data.frame(t(na.omit(t(fred_md))))
  md <- BVAR::fred_transform(fred_md, type = "fred_md") %>% as.data.frame()
  
  n_raw <- nrow(fred_md)
  dates_raw <- as.yearmon(seq(as.Date("1959-01-01"),by="month",length.out=n_raw))
  start_pt <- as.numeric(rownames(md)[1]) -1
  md <- cbind(date = dates_raw[start_pt:length(dates_raw)], md)  # adjust dates after transformation trimming
  
  # French 12 Industry (Value-Weighted, Monthly) + RF
  ind12 <- frenchdata::download_french_data("12 Industry Portfolios")
  ind12_m <- ind12$subsets |>
    dplyr::filter(grepl("Average Value Weighted Returns -- Monthly", name)) |>
    tidyr::unnest(data) |>
    dplyr::mutate(date = lubridate::ymd(paste(date, "01")))|>
    dplyr::mutate(date = zoo::as.yearmon(date)) |>
    as.data.frame()
  
  ff3 <- frenchdata::download_french_data("Fama/French 3 Factors")
  rf_m <- ff3$subsets |>
    dplyr::mutate(n = purrr::map_int(data, nrow)) |>
    arrange(desc(n)) |>
    slice(1) |>
    tidyr::unnest(data) |>
    dplyr::mutate(date = lubridate::ymd(paste(date, "01"))) |>
    dplyr::mutate(date = zoo::as.yearmon(date) ) |>
    dplyr::mutate(RF = RF/100) |>   # convert RF to decimal units
    dplyr::select(date, RF)
  
  ind_names <- setdiff(colnames(ind12_m), c("date","name"))
  ind_ex <- ind12_m %>%
    left_join(rf_m, by = "date") %>%
    # dplyr::mutate(across(all_of(ind_names), ~log(1 + .x/100 )- log(1 + RF))) %>%           # log gross return
    dplyr::mutate(across(all_of(ind_names), ~.x/100  -  RF)) %>% # excess return
    dplyr::mutate(across(all_of(ind_names), ~ .x  , .names = "{.col}_ex"))%>% # log excess
    dplyr::select(date, ends_with("_ex"))
  
  both <- inner_join(md, ind_ex, by = "date") %>% tidyr::drop_na()
  # Slice the target period
  start_ym <- zoo::as.yearmon(sprintf("%d-%02d", cfg$start_year, cfg$start_month))
  both <- dplyr::filter(both, date >= start_ym)
  if (!is.null(cfg$n_months) && is.finite(cfg$n_months)) {
    both <- head(both, cfg$n_months + cfg$var_h)
  }
  
  dates <- both$date
  Y <- both %>% dplyr::select(ends_with("_ex"))              # targets: industry log excess returns
  X <- both %>% dplyr::select(-date, -all_of(colnames(Y)))   # factor panel: FRED-MD
  list(X = as.matrix(X)[1:cfg$n_months,], Y = as.matrix(Y) , dates = dates)
}

## ============================== 2) FACTORING (X -> F) ===================
make_F_from_X <- function(X, dates,ic_idx, k = NA){
  # Select the number of factors by Bai--Ng IC2
  if(is.na(k)){
    r_star <- (dfms::ICr(X, max.r = cfg$r_max))$r[ic_idx]
  }else{
    r_star <- k 
  }
  
  
  if (is.null(r_star) || r_star < 1) r_star <- 6
  
  eigen_S <- eigen(t(X)%*%X/nrow(X))
  spiked_eigval <- eigen_S$values[1:r_star,drop=FALSE]
  spiked_eigvec <- eigen_S$vectors[,1:r_star,drop=FALSE]
  
  L <- spiked_eigvec %*% diag(sqrt(spiked_eigval),r_star)
  
  F_use <- diag(1/sqrt(spiked_eigval),r_star) %*% t(spiked_eigvec) %*% t(X)
  colnames(F_use) <- paste0("F", seq_len(ncol(F_use)))
  
  list(F = F_use, r = r_star)
}

## ============================== 3) AFM on Y =============================
# Optional external functions/files are used when available and skipped otherwise.
maybe_source <- function(path, type = c("R","cpp")){
  type <- match.arg(type)
  if (isTRUE(file.exists(path))) {
    if (type == "cpp") {
      if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
      Rcpp::sourceCpp(path); message("sourceCpp: ", path)
    } else {
      source(path); message("source: ", path)
    }
    TRUE
  } else { message("skip (not found): ", path); FALSE }
}

mean_Gamma_sqrtLambda <- function(draws, align = c("none","sign"), p){
  align <- match.arg(align)
  S <- length(draws$Gamma); G_ref <-diag( 1 , nrow = p ,ncol = length(draws$Gamma[[1]])/p)
  #G_ref <- matrix(draws$Gamma[[1]], nrow = p)
  p2 <- nrow(G_ref); k <- ncol(G_ref); acc <- matrix(0, p2, k)
  align_sign <- function(G, G0){ s <- sign(colSums(G*G0)); s[s==0] <- 1; G %*% diag(s, k) }
  for (s in seq_len(S)) {
    G <- matrix(draws$Gamma[[s]], nrow = p)
    L <- draws$Lambda[[s]]
    Ld <- if (is.null(dim(L))) diag(sqrt(as.vector(L)), k) else diag(sqrt(diag(L)), k)
    if (align == "sign") G <- align_sign(G, G_ref)
    acc <- acc + (G %*% Ld)
  }
  acc / S
}

run_afm_methods_on_X <- function(X, k, cfg){
  n <- nrow(X); p <- ncol(X); x <- X; S <- crossprod(x)/n
  evS <- eigen(S); sp <- evS$values[1:k]
  
  out <- list(); Sigma_u <- list(); Factor_gls <- list()
  
  ## 1) Proposed static AFM sampler
  have1 <- maybe_source(cfg$path_cpp_static_afm, "cpp")
  if (have1 && exists("gibbs_factor_full_sigma_timed")) {
    iter <- 20000; burnin <- floor(iter/2); thin <- max(1, iter/10000)
    nu <- 2*p + 1; S0 <- diag(p,p)
    #nu <- p + 1 + sqrt(n*p); S0 <- diag(sqrt(n*p)/2,p)
    res <- gibbs_factor_full_sigma_timed(
      x, k, iter, burnin, thin, a_vec = rep(2,k), q_hyper = 4,
      nu0 = nu, S0 = S0, b0 = 0.1, b1 = 2,
      evS$vectors[,1:k], sp, Sigma_init = diag(1,p),Z_init = t(evS$vectors[,1:k])%*%t(x),
      verbose = TRUE
    )
    Sigma_u_afm <- matrix(Reduce("+", res$Sigma) / length(res$Sigma), nrow = p)
    L_afm <- mean_Gamma_sqrtLambda(res, align = "sign", p = p)
    out$AFM   <- L_afm %*% t(L_afm) + Sigma_u_afm
    Sigma_u$AFM <- Sigma_u_afm
    Factor_gls$AFM <-  gls_scores(L_afm, Sigma_u_afm, x,identity = FALSE)
    rm(res)
  }
  
  ## 3) SPOET (always available)
  shrink_lambda <- evS$values[1:k] - p / (n * p - n * k - p * k) * sum(evS$values[(k + 1):p])
  poet <- POET::POET(t(x) - colMeans(x), k)
  out$SPOET <- poet$SigmaU + evS$vectors[,1:k] %*% diag(shrink_lambda,k) %*% t(evS$vectors[,1:k])
  Sigma_u$SPOET <- poet$SigmaU
  L_poet <- evS$vectors[,1:k] %*% diag(sqrt(shrink_lambda),k)
  Factor_gls$SPOET <- gls_scores(L_poet, poet$SigmaU, x,identity = FALSE)
  
  ## 4) PML (only when the optional file is available)
  have4 <- maybe_source(cfg$path_r_pml, "R")
  if (have4 && exists("em_mm_joint")) {
    set.seed(42)
    pml <- em_mm_joint(t(x), k, Lambda_init = t(mvtnorm::rmvnorm(k, mean = rep(0,p), sigma = diag(1,p))),
                       Sigma_init = diag(1,p), tol = 1e-9, max_iter = 5000, verbose = FALSE)
    out$PML <- pml$Lambda %*% t(pml$Lambda) + pml$Sigma
    Sigma_u$PML <- pml$Sigma
    Factor_gls$PML <- gls_scores(pml$Lambda, pml$Sigma, x,identity = FALSE)
  }
  
  ## 5) Sample
  out$Sample <- S
  
  list(Sigma = out, Sigma_u = Sigma_u, Factor_gls = Factor_gls)
}

## ============================== 4) VAR/BVAR =============================
safe_bvar_predict <- function(bfit, H, qs = c(0.05,0.16,0.50,0.84,0.95)){
  # Call predict() from the BVAR namespace directly because bvartools may mask predict.bvar.
  f <- getFromNamespace("predict.bvar", "BVAR")
  pred <- f(bfit, BVAR::bv_fcast(as.integer(H)), conf_bands = qs, n_thin = 1L)
  sb <- summary(pred)
  qu <- if (is.list(sb) && !is.null(sb$quants)) sb$quants else sb
  if (is.null(dimnames(qu)[[1]]) || length(dimnames(qu)[[1]]) != length(qs))
    dimnames(qu)[[1]] <- paste0(round(qs*100), "%")
  assign("bands", dimnames(qu)[[1]], inherits = TRUE)
  qu
}

run_var_pipeline <- function(Y, Factor, dates, cfg, lag, pred_num , scale = FALSE){
  # Z <- cbind(Y[1:cfg$n_months,], Factor); colnames(Z) <- c(colnames(Y), colnames(Factor)); rownames(Z) <- as.character(dates[1:cfg$n_months])
  # Z: [Y , Factor]  or just [Y] if Factor is NULL/empty
  
  if (!is.na(Factor[1,1])) {
    F_train <- Factor[1:cfg$n_months, , drop = FALSE]
    muF  <- colMeans(F_train)
    sdF  <- apply(F_train, 2, sd)
    sdF[sdF == 0 | is.na(sdF)] <- 1
    Factor <- sweep(Factor, 2, muF, "-")
    Factor <- sweep(Factor, 2, sdF, "/")
  }
  
  if(scale){
    Y<- scale(Y)
    Factor <- scale(Factor)
  }
  
  if (is.na(Factor[1,1]) ) {
    Z <- Y[1:cfg$n_months, , drop = FALSE]
  } else {
    Z <- cbind(Y[1:cfg$n_months, , drop = FALSE],
               Factor[1:cfg$n_months, , drop = FALSE])
  }
  
  colnames(Z) <- colnames(Z)  # no-op, keeps names
  rownames(Z) <- as.character(dates[1:cfg$n_months])
  
  if(is.na(pred_num)){
    H <- cfg$var_h
  }else{
    H <- pred_num
  }
  n_train <- nrow(Z)
  
  if(is.na(lag)){
    sel <- vars::VARselect(Z[,], lag.max = cfg$var_lagmax, type = "const")
    lag <- as.integer(sel$selection["AIC(n)"]); if (is.na(lag) || lag < 2) lag <- 4
  }
  
  
  set.seed(42)
  bfit <- BVAR::bvar(
    data   = Z[,],
    lags   = lag,
    n_draw = cfg$bvar_draw,
    n_burn = cfg$bvar_burn,
    n_thin = cfg$bvar_thin,
    priors = BVAR::bv_priors()
  )
  
  qq <- safe_bvar_predict(bfit, H , qs = c(0.005,0.025,0.5,0.975,0.995))
  dimnames(qq) <- list(dimnames(qq)[[1]],NULL,colnames(Z))
  oos_true <- Y[(cfg$n_months+1):(cfg$n_months+H), colnames(Y), drop = FALSE] 
  
  # if(scale){
  #   oos_true <- oos_true * attr(Y,"scaled:scale")  + attr(Y,"scaled:center")
  #   for(i in 1:ncol(Y)){
  #     qq[,,i] <- qq[,,i] * attr(Y,"scaled:scale")[i] + attr(Y,"scaled:center")[i]
  #   }
  # }
  
  # qq: [bands x (H+1) x variables]; keep only the final H forecast intervals
  t_idx <- tail(seq_len(dim(qq)[2]), H)
  med   <- aperm(qq[qnm(0.5), t_idx, , drop = FALSE], c(2,3,1))[,,1]
  
  if( H == 1){
    rmse_oos <- abs(oos_true - med[1:ncol(Y)])
  }else {
    rmse_oos <- sqrt(colMeans((oos_true - med[, 1:ncol(Y), drop = FALSE])^2))
    #rmse_oos <- (oos_true - med[, 1:ncol(Y), drop = FALSE])^2
  }
  
  list(bfit = bfit, quants = qq, med = med, rmse = rmse_oos, Z = Z, n_train = n_train , Y = Y)
}

## ============================== 5) MAIN =================================
pan <- load_panels(cfg)
# X <- scale(pan$X,center = TRUE,scale = TRUE); Y <- scale(pan$Y,center=FALSE,scale =TRUE); dates <- pan$dates
X <- scale(pan$X,center = TRUE,scale = TRUE); Y <- scale(pan$Y,center=TRUE,scale =TRUE); dates <- pan$dates
eigen(t(X)%*%X/nrow(X))$values


fx <- make_F_from_X(X, dates[1:nrow(X)],ic_idx=2)
F_use <- fx$F; k_x <- fx$r 
# cat(sprintf("Selected r (X factors) = %d (use_wpc = %s)\n", k_x, cfg$use_wpc))

set.seed(1234)
afm <- run_afm_methods_on_X(X, k = k_x, cfg = cfg)
cat("AFM methods computed:", paste(names(afm$Sigma), collapse = ", "), "\n")

# Add PCA factors (from step 2) and a VAR-only baseline
afm$Factor_gls$PCA    <- F_use
afm$Factor_gls$VARonly <- NA

var_out_list <- list()
lag <- 4
pred_num <- 1

set.seed(1234)
for(i in 1:length(afm$Factor_gls)){ 
  Factor <- t(afm$Factor_gls[[i]])
  if ( !is.na(Factor[1,1]) ) {
    colnames(Factor) <- paste0("F", seq_len(ncol(Factor)))
  }
  var_out <- run_var_pipeline(Y, Factor, dates, cfg, lag, pred_num , scale = FALSE)
  var_out_list <- append(var_out_list,list(var_out))
  paste0("OOS RMSE method ",names(afm$Factor_gls)[i]," :\n"); print(var_out$rmse)
}
names(var_out_list)<-names(afm$Factor_gls)

## ============================== 6) QUICK PLOTS ===========================
suppressPackageStartupMessages(require(ggplot2))
var_out <- var_out_list$AFM
if(is.na(pred_num)){
  H <- cfg$var_h
}else{
  H <- pred_num
}
n_train <- var_out$n_train; Z <- var_out$Z; qq <- var_out$quants
scale_Y <- attr(Y,"scaled:scale")


vars_show <- colnames(Y)[1:6]
idx_hist  <- max(1, n_train - 500):(n_train+1)
idx_test  <- (n_train+1):(n_train+H)
split_date <- as.Date(zoo::as.Date(dates[n_train]))

med_lab <- qnm(0.5); l95 <- qnm(0.025); u95 <- qnm(0.975)
l99 <- qnm(0.005); u99 <- qnm(0.995)

hist_df <- purrr::map_dfr(vars_show, function(v){
  tibble(date = as.Date(zoo::as.Date(dates[idx_hist])),
         actual = as.numeric(Y[idx_hist, v]* scale_Y[v]),
         var = v, piece = "history")
})
fc_df <- purrr::map_dfr(vars_show, function(v){
  tibble(date = as.Date(zoo::as.Date(dates[idx_test])),
         actual = as.numeric(Y[idx_test, v] * scale_Y[v]),
         med    = as.numeric(aperm(qq[med_lab, idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]),
         lo95   = as.numeric(aperm(qq[l95,    idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]),
         hi95   = as.numeric(aperm(qq[u95,    idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]),
         lo99   = if("0.5%" %in% dimnames(qq)[[1]]) as.numeric(aperm(qq[l99, idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]) else NA_real_,
         hi99   = if("99.5%"%in% dimnames(qq)[[1]]) as.numeric(aperm(qq[u99, idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]) else NA_real_,
         var = v, piece = "forecast")
})
ggplot() +
  geom_ribbon(data = subset(fc_df, !is.na(lo95)), aes(date, ymin = lo95, ymax = hi95), alpha = 0.12) +
  geom_ribbon(data = subset(fc_df, !is.na(lo99)), aes(date, ymin = lo99, ymax = hi99), alpha = 0.20) +
  geom_line  (data = hist_df, aes(date, y = actual)) +
  geom_line  (data = fc_df,   aes(date, y = actual)) +
  geom_line  (data = fc_df,   aes(date, y = med), linetype = "dashed") +
  geom_vline (xintercept = split_date, linetype = "dotted") +
  facet_wrap(~var, ncol = 1, scales = "free_y") +
  labs(title = NULL, x = NULL, y = NULL)

var_out_list$AFM$rmse ; var_out_list$SPOET$rmse ; var_out_list$PML$rmse ; var_out_list$PCA$rmse ; var_out_list$VARonly$rmse
sqrt(mean(var_out_list$AFM$rmse^2 )) ; sqrt( mean(var_out_list$SPOET$rmse^2)) ; sqrt(mean(var_out_list$PML$rmse^2)) ; sqrt(mean(var_out_list$PCA$rmse^2)) ; sqrt(mean(var_out_list$VARonly$rmse^2))
sqrt(mean(var_out_list$AFM$rmse^2 * attr(Y,"scaled:scale")^2)) ;  sqrt(mean(var_out_list$SPOET$rmse^2* attr(Y,"scaled:scale")^2)) ;sqrt( mean(var_out_list$PML$rmse^2* attr(Y,"scaled:scale")^2)) ; sqrt(mean(var_out_list$PCA$rmse^2* attr(Y,"scaled:scale")^2)) ; sqrt(mean(var_out_list$VARonly$rmse^2* attr(Y,"scaled:scale")^2))

### plotting
colnames(Y) <- c("Nondurables","Durables","Manufacturing","Energy","Chemicals","Business Equipment","Telecom","Utilities","Shops","Health","Money","Other") 
dimnames(qq) <- list(dimnames(qq)[[1]],NULL,c(colnames(Y),colnames(Z)[13:(12+k_x)]) )
names(scale_Y) <- colnames(Y) 
plot_forcast<-function(start, end){
  vars_show <- colnames(Y)[start: end]
  idx_hist  <- max(1, n_train - 500):(n_train+1)
  idx_test  <- (n_train+1):(n_train+H)
  split_date <- as.Date(zoo::as.Date(dates[n_train]))
  
  med_lab <- qnm(0.5); l95 <- qnm(0.025); u95 <- qnm(0.975)
  l99 <- qnm(0.005); u99 <- qnm(0.995)
  
  hist_df <- purrr::map_dfr(vars_show, function(v){
    tibble(date = as.Date(zoo::as.Date(dates[idx_hist])),
           actual = as.numeric(Y[idx_hist, v]* scale_Y[v]),
           var = v, piece = "history")
  })
  fc_df <- purrr::map_dfr(vars_show, function(v){
    tibble(date = as.Date(zoo::as.Date(dates[idx_test])),
           actual = as.numeric(Y[idx_test, v] * scale_Y[v]),
           med    = as.numeric(aperm(qq[med_lab, idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]),
           lo95   = as.numeric(aperm(qq[l95,    idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]),
           hi95   = as.numeric(aperm(qq[u95,    idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]),
           lo99   = if("0.5%" %in% dimnames(qq)[[1]]) as.numeric(aperm(qq[l99, idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]) else NA_real_,
           hi99   = if("99.5%"%in% dimnames(qq)[[1]]) as.numeric(aperm(qq[u99, idx_test - n_train, v, drop = FALSE]* scale_Y[v], c(2,1,3))[,1,1]) else NA_real_,
           var = v, piece = "forecast")
  })
  h1<- ggplot() +
    geom_ribbon(data = subset(fc_df, !is.na(lo95)), aes(date, ymin = lo95, ymax = hi95), alpha = 0.12) +
    geom_ribbon(data = subset(fc_df, !is.na(lo99)), aes(date, ymin = lo99, ymax = hi99), alpha = 0.20) +
    geom_line  (data = hist_df, aes(date, y = actual)) +
    geom_line  (data = fc_df,   aes(date, y = actual)) +
    geom_line  (data = fc_df,   aes(date, y = med), linetype = "dashed") +
    geom_vline (xintercept = split_date, linetype = "dotted") +
    facet_wrap(~var, ncol = 1, scales = "free_y") +
    labs(title = NULL, x = NULL, y = NULL)
  return(h1)
}

library(patchwork)
h1<- plot_forcast(1,4)
h2<- plot_forcast(5,8)
h3<- plot_forcast(9,12)
h1+h2+h3 