#include <RcppArmadillo.h>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =================================================================================
// 0. Utilities
// =================================================================================

// Inverse-Gamma
double rinv_gamma(double alpha, double beta){
  return 1.0 / R::rgamma(alpha, 1.0 / beta);
}

// Givens Rotation for Columns (Q updated)
inline void givens_rotate_cols(arma::mat& Q, int u, int v, double c, double s, double eps1){
  arma::vec cu = Q.col(u);
  arma::vec cv = Q.col(v);
  Q.col(u) =  eps1 * ( c * cu - s * cv );
  Q.col(v) =  s * cu + c * cv;
}

inline double logsumexp2(double a, double b){
  double m = std::max(a,b);
  return m + std::log(std::exp(a-m) + std::exp(b-m));
}

inline double wrap_angle(double th){
  th = std::fmod(th + M_PI, 2.0*M_PI);
  if (th < 0) th += 2.0*M_PI;
  return th - M_PI; // (-pi, pi]
}

// Efficient update: rotate symmetric matrix M by 2D Givens on (u,v)
inline void rotate_sym_matrix(arma::mat& M, int u, int v, double c, double s){
  int p = M.n_rows;
  
  arma::vec Mu = M.col(u);
  arma::vec Mv = M.col(v);
  
  double m_uu = M(u,u);
  double m_vv = M(v,v);
  double m_uv = M(u,v);
  
  double m_uu_new = c*c*m_uu - 2.0*c*s*m_uv + s*s*m_vv;
  double m_vv_new = s*s*m_uu + 2.0*c*s*m_uv + c*c*m_vv;
  double m_uv_new = c*s*(m_uu - m_vv) + (c*c - s*s)*m_uv;
  
  for(int kk=0; kk<p; ++kk){
    if(kk == u || kk == v) continue;
    double muk = Mu(kk);
    double mvk = Mv(kk);
    
    M(kk,u) = M(u,kk) = c * muk - s * mvk;
    M(kk,v) = M(v,kk) = s * muk + c * mvk;
  }
  
  M(u,u) = m_uu_new;
  M(v,v) = m_vv_new;
  M(u,v) = M(v,u) = m_uv_new;
}

// Random Pairs generator
inline std::vector<std::pair<int,int>> make_random_pairs(int p){
  arma::uvec perm = arma::randperm(p);
  std::vector<std::pair<int,int>> pairs;
  for(int i=0; i+1<p; i+=2)
    pairs.emplace_back( (int)perm[i], (int)perm[i+1] );
  if(p % 2 == 1){
    int u = (int)perm[p-1];
    int v = (int)perm[(int)std::floor(R::runif(0.0, (double)(p-1)))];
    if(v==u) v = (v+1) % p;
    pairs.emplace_back(u,v);
  }
  return pairs;
}

// Slice Sampler for Circular Domain (Theta) – anchored at 0 (no-rotation)
inline double slice_sample_theta(const std::function<double(double)>& logp,
                                 double w = 0.1, int m = 50,
                                 int max_shrink = 1000)
{
  const double theta0 = 0.0;
  const double logp0  = logp(theta0);
  const double logy   = logp0 + std::log(R::runif(0.0, 1.0));
  
  double U  = R::runif(0.0,  1.0);
  double L  = theta0 - w * U;
  double Rr = L + w;
  
  int J = static_cast<int>(std::floor(m * R::runif(0.0, 1.0)));
  int K = (m - 1) - J;
  while (J-- > 0 && logp(L) > logy)  L -= w;
  while (K-- > 0 && logp(Rr) > logy) Rr += w;
  
  for (int it = 0; it < max_shrink; ++it) {
    double th = R::runif(L, Rr);
    double lp = logp(th);
    if (lp >= logy) return th;
    if (th < theta0) L = th; else Rr = th;
    if (Rr - L < 1e-9) break;
  }
  return theta0;
}

// Slice Sampler for Interval Domain
inline double slice_sample_interval(const std::function<double(double)>& logp,
                                    double x0, double lb, double ub, double w = 0.1)
{
  if(x0 < lb) x0 = lb + 1e-6;
  if(x0 > ub) x0 = ub - 1e-6;
  
  double logp0 = logp(x0);
  if(!std::isfinite(logp0)){
    x0 = (lb + ub)/2.0;
    logp0 = logp(x0);
  }
  
  double logy = logp0 + std::log(R::runif(0.0, 1.0));
  double U = R::runif(0.0, 1.0);
  double L = x0 - w * U;
  double Rr = L + w;
  
  L  = std::max(L,  lb);
  Rr = std::min(Rr, ub);
  
  for(int i=0; i<250; ++i){
    double x_prop  = R::runif(L, Rr);
    double lp_prop = logp(x_prop);
    if(lp_prop >= logy) return x_prop;
    if(x_prop < x0) L = x_prop; else Rr = x_prop;
    if((Rr - L) < 1e-10) break;
  }
  return x0;
}

// =================================================================================
// 1. MAIN GIBBS FUNCTION (Factor-wise AR(1): A = diag(a_j))
// =================================================================================

// [[Rcpp::export]]
Rcpp::List gibbs_factor_full_sigma_ar1_diagA_timed(
    const arma::mat &X, int k,
    int iter, int burnin, int thin,
    arma::vec &a_vec,               // IG shape hyper for Lambda (existing)
    double q_hyper,                 // IG scale hyper for Lambda (h)
    double nu0, const arma::mat &S0,
    double b0, double b1,            // Eigenvalue constraints
    const arma::mat &Gamma_init,
    const arma::vec &Lambda_init,
    const arma::mat &Z_init,
    const arma::mat &Sigma_init,
    arma::vec a_ar_init,             // NEW: initial a_j (length k); if empty -> 0.5
    double sigma_a_prior = 1.0,       // NEW: prior sd for each a_j
    bool sample_a = true,             // NEW
    bool verbose = false
){
  using clk = std::chrono::steady_clock;
  auto wall_start = clk::now();
  
  // Timings
  double t_Z = 0.0;
  double t_a = 0.0;
  double t_Lambda = 0.0;
  double t_Gamma = 0.0;
  double t_Sigma_Setup = 0.0;
  double t_Sigma_Vec = 0.0;
  double t_Sigma_Val = 0.0;
  
  int n = X.n_rows; // T
  int p = X.n_cols;
  arma::mat Y = X.t(); // p x n
  
  // Initialize a_vec for Lambda prior if needed (your original logic)
  arma::mat S_emp = Y * Y.t();
  arma::vec eigen_S;
  arma::eig_sym(eigen_S, S_emp / n);
  eigen_S = arma::reverse(eigen_S);
  if (p > n && (a_vec.n_elem == 0 || a_vec(0) == 0)) {
    double mean_non_sp = arma::mean(eigen_S.subvec(k, n-1));
    a_vec = n * (mean_non_sp) / (eigen_S.subvec(0, k-1) - mean_non_sp) / 2 + 2;
  }
  
  // States
  arma::mat Gamma = Gamma_init;
  arma::vec Lambda = Lambda_init;
  arma::mat Z = Z_init;
  
  // a_ar init
  arma::vec a_ar(k, fill::zeros);
  if(a_ar_init.n_elem == (unsigned)k) {
    a_ar = a_ar_init;
  } else {
    a_ar.fill(0.5);
  }
  // clip to (-1,1)
  for(int j=0; j<k; ++j){
    if(a_ar(j) <= -0.999) a_ar(j) = -0.5;
    if(a_ar(j) >=  0.999) a_ar(j) =  0.5;
  }
  
  // Sigma eigendecomp
  arma::vec d_u;
  arma::mat Gamma_u;
  eig_sym(d_u, Gamma_u, Sigma_init);
  arma::vec inv_d_u = 1.0 / d_u;
  arma::mat Sigma = Sigma_init;
  
  // Gamma null space
  arma::mat Gamma_perp = null(Gamma.t());
  
  // Outputs
  std::vector<arma::mat> out_Gamma, out_Sigma;
  std::vector<arma::vec> out_Lambda;
  std::vector<arma::vec> out_a;
  
  // Workspace
  arma::mat SInvY(p, n);        // Sigma^{-1} Y
  arma::mat GtSInvY(k, n);      // Gamma' Sigma^{-1} Y
  
  // For FFBS storage (k small)
  arma::mat m_filt(k, n, fill::zeros);
  arma::mat a_pred(k, n, fill::zeros);
  arma::cube C_filt(k, k, n, fill::zeros);
  arma::cube R_pred(k, k, n, fill::zeros);
  
  Rcpp::RNGScope scope;
  
  // Helper to build A and Q quickly
  auto build_A = [&](const arma::vec& a)->arma::mat{
    return diagmat(a);
  };
  auto build_Q = [&](const arma::vec& a, const arma::vec& lam)->arma::mat{
    arma::vec om = 1.0 - square(a);
    // numerical guard
    for(uword j=0; j<om.n_elem; ++j) if(om(j) < 1e-8) om(j) = 1e-8;
    return diagmat(lam % om);
  };
  
  for(int it=0; it < iter; ++it){
    
    // -------------------------------------------------------
    // 1. Update Z: AR(1) via FFBS with A=diag(a_j), Q=diag(lambda_j(1-a_j^2))
    // -------------------------------------------------------
    auto t0 = clk::now();
    
    // Sigma^{-1} Y using eigendecomp: SInvY = V D^{-1} V' Y
    {
      arma::mat tmp = Gamma_u.t() * Y;  // p x n
      tmp.each_col() %= inv_d_u;        // scale by 1/d
      SInvY = Gamma_u * tmp;            // Sigma^{-1} Y
    }
    
    // Compute GtSigInvG = Gamma' Sigma^{-1} Gamma
    arma::vec inv_sqrt_d = sqrt(inv_d_u);
    arma::mat W = Gamma_u.t() * Gamma;      // p x k
    W.each_col() %= inv_sqrt_d;             // D^{-1/2} * (V' Gamma)
    arma::mat GtSigInvG = W.t() * W;        // k x k
    
    // Precompute Gamma' Sigma^{-1} Y
    GtSInvY = Gamma.t() * SInvY;            // k x n
    
    arma::mat A = build_A(a_ar);
    arma::mat Q = build_Q(a_ar, Lambda);
    arma::mat LambdaMat = diagmat(Lambda);
    
    // Forward filter (information form)
    for(int t=0; t<n; ++t){
      arma::vec at(k, fill::zeros);
      arma::mat Rt(k,k, fill::zeros);
      
      if(t == 0){
        at.zeros();
        Rt = LambdaMat; // Z1 ~ N(0, Lambda)
      } else {
        at = A * m_filt.col(t-1);
        // Rt = A C A' + Q; A diagonal -> row/col scaling
        arma::mat Cp = C_filt.slice(t-1);
        Rt = Cp;
        Rt.each_row() %= a_ar.t(); // left-multiply by diag(a)
        Rt.each_col() %= a_ar;     // right-multiply by diag(a)
        Rt += Q;
      }
      
      arma::mat Rt_sym = symmatu(Rt);
      arma::mat Rt_inv = inv_sympd(Rt_sym);
      
      arma::mat Ct = inv_sympd(symmatu(Rt_inv + GtSigInvG));
      arma::vec bt = Rt_inv * at + GtSInvY.col(t);
      arma::vec mt = Ct * bt;
      
      a_pred.col(t) = at;
      R_pred.slice(t) = Rt_sym;
      m_filt.col(t) = mt;
      C_filt.slice(t) = symmatu(Ct);
    }
    
    // Backward sampling
    {
      arma::mat Ct = C_filt.slice(n-1);
      arma::mat L;
      bool ok = chol(L, Ct + 1e-10 * eye(k,k), "lower");
      if(!ok) stop("Cholesky failed in FFBS (final state).");
      Z.col(n-1) = m_filt.col(n-1) + L * randn<vec>(k);
    }
    
    for(int t=n-2; t>=0; --t){
      arma::mat Ct = C_filt.slice(t);
      arma::mat Rn = R_pred.slice(t+1);
      
      arma::mat Rn_inv = inv_sympd(Rn);
      
      // J_t = C_t A' R_{t+1}^{-1}; A' = A (diag)
      arma::mat J = Ct * A.t() * Rn_inv;
      
      arma::vec mean = m_filt.col(t) + J * (Z.col(t+1) - a_pred.col(t+1));
      arma::mat Cov  = symmatu(Ct - J * Rn * J.t());
      
      arma::mat L;
      bool ok = chol(L, Cov + 1e-10 * eye(k,k), "lower");
      if(!ok) {
        ok = chol(L, Cov + 1e-6 * eye(k,k), "lower");
        if(!ok) stop("Cholesky failed in FFBS (backward).");
      }
      Z.col(t) = mean + L * randn<vec>(k);
    }
    
    t_Z += std::chrono::duration<double>(clk::now() - t0).count();
    
    // -------------------------------------------------------
    // 1.5 Update a_j by slice sampling: a_j ~ N(0,sigma^2) I(|a_j|<1)
    //     Prior of Z given (a_j,lambda_j) factorizes across j (state eq diagonal)
    // -------------------------------------------------------
    if(sample_a){
      t0 = clk::now();
      
      const double inv_sig2 = 1.0 / (sigma_a_prior * sigma_a_prior);
      
      for(int j=0; j<k; ++j){
        // precompute sufficient stats for factor j
        double s00 = 0.0, s11 = 0.0, s10 = 0.0;
        for(int t=1; t<n; ++t){
          double z0 = Z(j, t-1);
          double z1 = Z(j, t);
          s00 += z0*z0;
          s11 += z1*z1;
          s10 += z1*z0;
        }
        
        double lamj = Lambda(j);
        
        auto logp_aj = [&](double a)->double{
          if(std::abs(a) >= 1.0) return -std::numeric_limits<double>::infinity();
          double om = 1.0 - a*a;
          if(om <= 0.0) return -std::numeric_limits<double>::infinity();
          
          double quad = s11 - 2.0*a*s10 + a*a*s00;
          // likelihood (ignoring constants and lambda-only terms)
          double ll = -0.5 * (double)(n-1) * std::log(om)
            -0.5 * quad / (lamj * om);
          double lp = -0.5 * a*a * inv_sig2;
          return ll + lp;
        };
        
        a_ar(j) = slice_sample_interval(logp_aj, a_ar(j), -0.999, 0.999, 0.05);
      }
      
      t_a += std::chrono::duration<double>(clk::now() - t0).count();
    }
    
    // -------------------------------------------------------
    // 2. Update Lambda (factor-wise AR(1) sufficient stats)
    // -------------------------------------------------------
    t0 = clk::now();
    
    for(int j=0; j<k; j++){
      double aj = a_ar(j);
      double om = 1.0 - aj*aj;
      if(om < 1e-8) om = 1e-8;
      
      double sum_innov = 0.0;
      for(int t=1; t<n; ++t){
        double diff = Z(j,t) - aj * Z(j,t-1);
        sum_innov += diff*diff;
      }
      double Sj = Z(j,0)*Z(j,0) + sum_innov / om;
      
      double alpha_post = a_vec(j) - 1.0 + 0.5 * n;
      double beta_post  = 0.5 * q_hyper + 0.5 * Sj;
      Lambda(j) = rinv_gamma(alpha_post, beta_post);
    }
    
    t_Lambda += std::chrono::duration<double>(clk::now() - t0).count();
    
    // -------------------------------------------------------
    // 3. Update Gamma using the existing logic with the updated Z
    // -------------------------------------------------------
    t0 = clk::now();
    
    arma::mat S_Z = Z * Z.t();
    arma::vec eigval_Sz; arma::mat eigvec_Sz;
    eig_sym(eigval_Sz, eigvec_Sz, S_Z);
    eigval_Sz = arma::reverse(eigval_Sz);
    eigvec_Sz = arma::fliplr(eigvec_Sz);
    
    arma::mat F = SInvY * Z.t(); // Sigma^{-1} Y Z'
    arma::mat mod_F = join_horiz(F * eigvec_Sz, arma::zeros<arma::mat>(p, p - k));
    arma::mat mod_Q = join_horiz(Gamma * eigvec_Sz, Gamma_perp);
    arma::vec diag_D_k = arma::join_cols(eigval_Sz , arma::zeros<arma::vec>(p - k));
    
    arma::mat T_A = Gamma_u.t() * mod_Q;
    T_A.each_col() %= inv_d_u;
    arma::mat AQ = Gamma_u * T_A;
    arma::mat QtF  = mod_Q.t() * mod_F;
    
    auto pairs_g = make_random_pairs(p);
    for(const auto& pr : pairs_g){
      int u = std::min(pr.first, pr.second);
      int v = std::max(pr.first, pr.second);
      
      if(u >= k && v >= k) {
        double th = R::runif(-M_PI, M_PI);
        double eps1  = (R::runif(0.0, 1.0) < 0.5) ? -1.0 : 1.0;
        double c = std::cos(th), s = std::sin(th);
        givens_rotate_cols(mod_Q, u, v, c, s, eps1);
        givens_rotate_cols(AQ, u, v, c, s, eps1);
        
        arma::uvec idx = { (arma::uword)u, (arma::uword)v };
        arma::mat R(2,2); R << eps1*c << -eps1*s << arma::endr
                            <<     s  <<      c  << arma::endr;
        QtF.rows(idx) = R.t() * QtF.rows(idx);
        continue;
      }
      
      const arma::vec& g1 = mod_Q.col(u);
      const arma::vec& g2 = mod_Q.col(v);
      double t11 = arma::dot(g1, AQ.col(u));
      double t22 = arma::dot(g2, AQ.col(v));
      double t12 = arma::dot(g1, AQ.col(v));
      double Buu = QtF(u,u), Bvu = QtF(v,u), Buv = QtF(u,v), Bvv = QtF(v,v);
      
      auto logp_cached = [&](double th, double eps1)->double{
        double c = std::cos(th), s = std::sin(th);
        double c2=c*c, s2=s*s, sc=s*c;
        double dL  = eps1*(Buu*c - Bvu*s) + (Buv*s + Bvv*c);
        double dK  = ( c2*t11 - 2.0*sc*t12 + s2*t22 ) * diag_D_k(u)
          + ( s2*t11 + 2.0*sc*t12 + c2*t22 ) * diag_D_k(v);
        return dL - 0.5*dK;
      };
      auto logp_theta = [&](double th)->double{
        th = wrap_angle(th);
        return logsumexp2(logp_cached(th, +1.0), logp_cached(th, -1.0));
      };
      
      double theta = slice_sample_theta(logp_theta, 0.01, 50, 1000);
      theta = wrap_angle(theta);
      
      double lp_p = logp_cached(theta,+1.0), lp_m = logp_cached(theta,-1.0);
      double den = std::exp(lp_p) + std::exp(lp_m);
      double prp = (den > 0) ? (std::exp(lp_p)/den) : 0.5;
      double eps1 = (R::runif(0,1) < prp) ? 1.0 : -1.0;
      
      double c = std::cos(theta), s = std::sin(theta);
      givens_rotate_cols(mod_Q, u, v, c, s, eps1);
      givens_rotate_cols(AQ,   u, v, c, s, eps1);
      
      arma::uvec idx = { (arma::uword)u, (arma::uword)v };
      arma::mat R(2,2); R << eps1*c << -eps1*s << arma::endr
                          <<     s  <<      c  << arma::endr;
      QtF.rows(idx) = R.t() * QtF.rows(idx);
    }
    
    Gamma = mod_Q.cols(0, k-1) * eigvec_Sz.t();
    Gamma_perp = mod_Q.cols(k, p-1);
    
    t_Gamma += std::chrono::duration<double>(clk::now() - t0).count();
    
    // -------------------------------------------------------
    // 4. Update Sigma_u using the existing eigen-decomposition Gibbs step
    // -------------------------------------------------------
    
    // [4-0] Setup
    t0 = clk::now();
    arma::mat E = Y - Gamma * Z;
    arma::mat Spost = S0 + E * E.t();
    double nu_post = nu0 + n;
    t_Sigma_Setup += std::chrono::duration<double>(clk::now() - t0).count();
    
    // [4-A] Eigenvectors update
    t0 = clk::now();
    arma::mat M = Gamma_u.t() * Spost * Gamma_u;
    
    auto pairs_s = make_random_pairs(p);
    for(const auto& pr : pairs_s){
      int u = pr.first;
      int v = pr.second;
      if(u == v) continue;
      
      auto logp_V = [&](double th) -> double {
        double c = std::cos(th), s = std::sin(th);
        double m_uu = M(u,u), m_vv = M(v,v), m_uv = M(u,v);
        double m_uu_new = c*c*m_uu - 2.0*c*s*m_uv + s*s*m_vv;
        double m_vv_new = s*s*m_uu + 2.0*c*s*m_uv + c*c*m_vv;
        return -0.5 * ( m_uu_new * inv_d_u(u) + m_vv_new * inv_d_u(v) );
      };
      
      double theta = slice_sample_theta(logp_V, 0.1, 20);
      if(std::abs(theta) > 1e-9){
        double c = std::cos(theta), s = std::sin(theta);
        givens_rotate_cols(Gamma_u, u, v, c, s, 1.0);
        rotate_sym_matrix(M, u, v, c, s);
      }
    }
    t_Sigma_Vec += std::chrono::duration<double>(clk::now() - t0).count();
    
    // [4-B] Eigenvalues update with Jacobian
    t0 = clk::now();
    for(int ii=0; ii<p; ++ii){
      auto logp_lam = [&](double lam) -> double {
        if(lam < b0 || lam > b1) return -std::numeric_limits<double>::infinity();
        double val = -0.5 * (nu_post + p + 1.0) * std::log(lam) - 0.5 * M(ii,ii) / lam;
        double jac = 0.0;
        for(int j=0; j<p; ++j){
          if(j == ii) continue;
          jac += std::log(std::abs(lam - d_u(j)));
        }
        return val + jac;
      };
      double w = 0.5 * std::max(d_u(ii), 1e-8);
      d_u(ii) = slice_sample_interval(logp_lam, d_u(ii), b0, b1, w);
    }
    inv_d_u = 1.0 / d_u;
    t_Sigma_Val += std::chrono::duration<double>(clk::now() - t0).count();
    
    Sigma = Gamma_u * diagmat(d_u) * Gamma_u.t();
    
    // -------------------------------------------------------
    // Save draws
    // -------------------------------------------------------
    if(it >= burnin && ((it - burnin) % thin == 0)){
      out_Gamma.push_back(Gamma);
      out_Lambda.push_back(Lambda);
      out_Sigma.push_back(Sigma);
      out_a.push_back(a_ar);
    }
    
    if(verbose && (it % 1000 == 0)){
      Rcpp::Rcout << "Iter: " << it << " completed." << std::endl;
    }
    Rcpp::checkUserInterrupt();
  }
  
  // Timing report
  double total_core = t_Z + t_a + t_Lambda + t_Gamma + t_Sigma_Setup + t_Sigma_Vec + t_Sigma_Val;
  auto wall_sec = std::chrono::duration<double>(clk::now() - wall_start).count();
  
  Rcpp::Rcout << "\n============================================\n";
  Rcpp::Rcout << "           Timing Report (Seconds)          \n";
  Rcpp::Rcout << "============================================\n";
  Rcpp::Rcout << "Total Wall Time  : " << wall_sec << " s\n";
  Rcpp::Rcout << "Total Core Time  : " << total_core << " s (overhead: " << wall_sec - total_core << " s)\n\n";
  
  auto pct = [&](double t){ return (total_core > 0) ? 100.0 * t / total_core : 0.0; };
  
  Rcpp::Rcout << "1. Z step (FFBS) : " << t_Z << " (" << pct(t_Z) << "%)\n";
  Rcpp::Rcout << "1.5 a step       : " << t_a << " (" << pct(t_a) << "%)\n";
  Rcpp::Rcout << "2. Lambda step   : " << t_Lambda << " (" << pct(t_Lambda) << "%)\n";
  Rcpp::Rcout << "3. Gamma step    : " << t_Gamma << " (" << pct(t_Gamma) << "%)\n";
  Rcpp::Rcout << "4. Sigma (Total) : " << (t_Sigma_Setup + t_Sigma_Vec + t_Sigma_Val)
              << " (" << pct(t_Sigma_Setup + t_Sigma_Vec + t_Sigma_Val) << "%)\n";
  Rcpp::Rcout << "   - Setup (E,S) : " << t_Sigma_Setup << "\n";
  Rcpp::Rcout << "   - Vectors (V) : " << t_Sigma_Vec << "\n";
  Rcpp::Rcout << "   - Values (d)  : " << t_Sigma_Val << "\n";
  Rcpp::Rcout << "============================================\n";
  
  return Rcpp::List::create(
    Named("Gamma")  = out_Gamma,
    Named("Lambda") = out_Lambda,
    Named("Sigma")  = out_Sigma,
    Named("a")      = out_a
  );
}
