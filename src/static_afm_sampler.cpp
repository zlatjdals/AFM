#include <RcppArmadillo.h>
#include <chrono>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// =================================================================================
// 0. Utility functions (slice sampler, math, and rotations)
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
  return th - M_PI;                 // now th ∈ (-π, π]
}

// Givens Rotation for Symmetric Matrix M update (M = R^T M R)
// Efficient O(p) update
inline void rotate_sym_matrix(arma::mat& M, int u, int v, double c, double s){
  int p = M.n_rows;
  
  arma::vec Mu = M.col(u); // copy
  arma::vec Mv = M.col(v); // copy
  
  double m_uu = M(u,u);
  double m_vv = M(v,v);
  double m_uv = M(u,v);
  
  double m_uu_new = c*c*m_uu - 2.0*c*s*m_uv + s*s*m_vv;
  double m_vv_new = s*s*m_uu + 2.0*c*s*m_uv + c*c*m_vv;
  double m_uv_new = c*s*(m_uu - m_vv) + (c*c - s*s)*m_uv;
  
  for(int k=0; k<p; ++k){
    if(k == u || k == v) continue;
    double muk = Mu(k);
    double mvk = Mv(k);
    
    M(k,u) = M(u,k) = c * muk - s * mvk;
    M(k,v) = M(v,k) = s * muk + c * mvk;
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

// Slice Sampler for Circular Domain (Theta)
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

// Slice Sampler for Interval Domain (Eigenvalues)
inline double slice_sample_interval(const std::function<double(double)>& logp,
                                    double x0, double lb, double ub, double w = 0.5)
{
  if(x0 < lb) x0 = lb + 1e-6;
  if(x0 > ub) x0 = ub - 1e-6;
  
  double logp0 = logp(x0);
  if(!std::isfinite(logp0)){ x0 = (lb + ub)/2.0; logp0 = logp(x0); }
  
  double logy = logp0 + std::log(R::runif(0.0, 1.0));
  double U = R::runif(0.0, 1.0);
  double L = x0 - w * U;
  double Rr = L + w;
  
  L = std::max(L, lb);
  Rr = std::min(Rr, ub);
  
  for(int i=0; i<100; ++i){
    double x_prop = R::runif(L, Rr);
    double lp_prop = logp(x_prop);
    if(lp_prop >= logy) return x_prop;
    if(x_prop < x0) L = x_prop; else Rr = x_prop;
    if((Rr - L) < 1e-9) break;
  }
  return x0;
}

// =================================================================================
// 1. MAIN GIBBS FUNCTION (With Timing)
// =================================================================================

// [[Rcpp::export]]
Rcpp::List gibbs_factor_full_sigma_timed(
    const arma::mat &X, int k,
    int iter, int burnin, int thin,
    arma::vec &a_vec, double q_hyper,
    double nu0, const arma::mat &S0,
    double b0, double b1, // Eigenvalue constraints
    const arma::mat &Gamma_init,
    const arma::vec &Lambda_init,
    const arma::mat &Z_init,
    const arma::mat &Sigma_init,
    bool verbose = false
){
  // Set up timers
  using clk = std::chrono::steady_clock;
  auto wall_start = clk::now(); // start wall-clock timer
  
  // Accumulated timing variables (seconds)
  double t_Z = 0.0;
  double t_Lambda = 0.0;
  double t_Gamma = 0.0;
  double t_Sigma_Setup = 0.0; // E and Spost calculations
  double t_Sigma_Vec = 0.0;   // eigenvector rotations (Step A)
  double t_Sigma_Val = 0.0;   // eigenvalue sampling (Step B)
  
  // Basic setup
  int n = X.n_rows, p = X.n_cols;
  arma::mat Y = X.t(); 
  
  // Initialize a_vec
  arma::mat S_emp = Y * Y.t();
  arma::vec eigen_S;
  arma::eig_sym(eigen_S, S_emp / n);
  eigen_S = arma::reverse(eigen_S);
  if (p > n && (a_vec.n_elem == 0 || a_vec(0) == 0)) {
    double mean_non_sp = arma::mean(eigen_S.subvec(k,n-1));
    a_vec = n * (mean_non_sp) / (eigen_S.subvec(0, k-1) - mean_non_sp) / 2 + 2;
  }
  
  // Initial values
  arma::mat Gamma = Gamma_init;
  arma::vec Lambda = Lambda_init;
  arma::mat Z = Z_init;
  
  // Sigma State (Eigen-decomposition)
  arma::vec d_u; 
  arma::mat Gamma_u; 
  eig_sym(d_u, Gamma_u, Sigma_init); 
  arma::vec inv_d_u = 1.0 / d_u;
  arma::mat Sigma = Sigma_init; 
  
  // Gamma Null Space
  arma::mat Gamma_perp = null(Gamma.t());
  
  // Storage
  std::vector<arma::mat> out_Gamma, out_Sigma;
  std::vector<arma::vec> out_Lambda;
  
  // Workspace
  arma::mat SInvY(p, n), Vmat(k, k), Ck(k, k), MU(k, n), S_Z(k, k);
  
  Rcpp::RNGScope scope;
  
  for(int it=0; it < iter; ++it){
    
    // -------------------------------------------------------
    // 1. Update Z (latent factors)
    // -------------------------------------------------------
    auto t0 = clk::now(); // start timer
    
    arma::mat tmp = Gamma_u.t() * Y; 
    tmp.each_col() %= inv_d_u;       
    SInvY = Gamma_u * tmp;           
    
    arma::vec inv_sqrt_d = sqrt(inv_d_u);
    arma::mat W = Gamma_u.t() * Gamma; 
    W.each_col() %= inv_sqrt_d;
    arma::mat GtSigInvG = W.t() * W;
    
    arma::mat Lminv = diagmat(1.0 / Lambda);
    Vmat = inv_sympd(Lminv + GtSigInvG);
    if(!chol(Ck, Vmat, "lower")) stop("Cholesky failed in Z-step");
    
    MU = Vmat * (Gamma.t() * SInvY);
    for(int i=0; i<n; i++) Z.col(i) = MU.col(i) + Ck * randn<vec>(k);
    
    t_Z += std::chrono::duration<double>(clk::now() - t0).count(); // stop timer and accumulate
    
    // -------------------------------------------------------
    // 2. Update Lambda
    // -------------------------------------------------------
    t0 = clk::now();
    
    for(int j=0; j<k; j++){
      double alpha_post = a_vec(j) - 1.0 + 0.5 * n;
      double beta_post  = 0.5 * q_hyper + 0.5 * accu(square(Z.row(j)));
      Lambda(j) = rinv_gamma(alpha_post, beta_post);
    }
    
    t_Lambda += std::chrono::duration<double>(clk::now() - t0).count();
    
    // -------------------------------------------------------
    // 3. Update Gamma (factor loadings)
    // -------------------------------------------------------
    t0 = clk::now();
    
    // Keep the basic Gamma update structure
    S_Z = Z * Z.t();
    arma::vec eigval_Sz; arma::mat eigvec_Sz;
    eig_sym(eigval_Sz, eigvec_Sz, S_Z);
    eigval_Sz = arma::reverse(eigval_Sz); 
    eigvec_Sz = arma::fliplr(eigvec_Sz);
    
    arma::mat F = SInvY * Z.t();
    arma::mat mod_F = join_horiz(F * eigvec_Sz, arma::zeros<arma::mat>(p, p - k));
    arma::mat mod_Q = join_horiz(Gamma * eigvec_Sz, Gamma_perp);
    arma::vec diag_D_k = arma::join_cols(eigval_Sz , arma::zeros<arma::vec>(p - k));
    
    arma::mat T_A = Gamma_u.t() * mod_Q;
    T_A.each_col() %= inv_d_u;
    arma::mat AQ = Gamma_u * T_A;
    arma::mat QtF  = mod_Q.t() * mod_F;   // p×p
    
    // Detailed rotation logic follows the same structure as before
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
      // Loading rotation logic; extend here if needed
      // Log-likelihood contribution using only scalar summaries
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
      double eps1 = (R::runif(0,1) < std::exp(lp_p)/(std::exp(lp_p)+std::exp(lp_m))) ? 1.0 : -1.0;
      
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
    // 4. Update Sigma_u (eigen-decomposition Gibbs step)
    // -------------------------------------------------------
    
    // [4-0] Setup (E, Spost)
    t0 = clk::now();
    arma::mat E = Y - Gamma * Z;
    arma::mat Spost = S0 + E * E.t(); 
    double nu_post = nu0 + n;
    t_Sigma_Setup += std::chrono::duration<double>(clk::now() - t0).count();
    
    // [4-A] Eigenvectors (Gamma_u) Update
    t0 = clk::now();
    arma::mat M = Gamma_u.t() * Spost * Gamma_u; // V^T S V
    
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
        givens_rotate_cols(Gamma_u, u, v, c, s, 1);
        rotate_sym_matrix(M, u, v, c, s);
      }
    }
    t_Sigma_Vec += std::chrono::duration<double>(clk::now() - t0).count();
    
    // [4-B] Eigenvalues (d_u) Update with Jacobian
    t0 = clk::now();
    for(int k=0; k<p; ++k){
      auto logp_lam = [&](double lam) -> double {
        if(lam < b0 || lam > b1) return -std::numeric_limits<double>::infinity();
        double val = -0.5 * (nu_post + p + 1.0) * std::log(lam) - 0.5 * M(k,k) / lam;
        double jac = 0.0;
        for(int j=0; j<p; ++j){
          if(j == k) continue;
          jac += std::log(std::abs(lam - d_u(j)));
        }
        return val + jac;
      };
      d_u(k) = slice_sample_interval(logp_lam, d_u(k), b0, b1, 0.5 * d_u(k));
    }
    inv_d_u = 1.0 / d_u;
    t_Sigma_Val += std::chrono::duration<double>(clk::now() - t0).count();
    
    // Final reconstruction for output
    Sigma = Gamma_u * diagmat(d_u) * Gamma_u.t();
    
    // -------------------------------------------------------
    // Save draws and print progress
    // -------------------------------------------------------
    if(it >= burnin && ((it - burnin) % thin == 0)){
      out_Gamma.push_back(Gamma);
      out_Lambda.push_back(Lambda);
      out_Sigma.push_back(Sigma);
    }
    
    if(verbose && (it % 1000 == 0)){
      Rcpp::Rcout << "Iter: " << it << " completed." << std::endl;
    }
    Rcpp::checkUserInterrupt();
  }
  
  // === Final timing report ===
  double total_core = t_Z + t_Lambda + t_Gamma + t_Sigma_Setup + t_Sigma_Vec + t_Sigma_Val;
  auto wall_sec = std::chrono::duration<double>(clk::now() - wall_start).count();
  
  Rcpp::Rcout << "\n============================================\n";
  Rcpp::Rcout << "           Timing Report (Seconds)          \n";
  Rcpp::Rcout << "============================================\n";
  Rcpp::Rcout << "Total Wall Time  : " << wall_sec << " s\n";
  Rcpp::Rcout << "Total Core Time  : " << total_core << " s (overhead: " << wall_sec - total_core << " s)\n\n";
  
  auto pct = [&](double t){ return (total_core > 0) ? 100.0 * t / total_core : 0.0; };
  
  Rcpp::Rcout << "1. Z step        : " << t_Z << " (" << pct(t_Z) << "%)\n";
  Rcpp::Rcout << "2. Lambda step   : " << t_Lambda << " (" << pct(t_Lambda) << "%)\n";
  Rcpp::Rcout << "3. Gamma step    : " << t_Gamma << " (" << pct(t_Gamma) << "%)\n";
  Rcpp::Rcout << "4. Sigma (Total) : " << (t_Sigma_Setup + t_Sigma_Vec + t_Sigma_Val) 
              << " (" << pct(t_Sigma_Setup + t_Sigma_Vec + t_Sigma_Val) << "%)\n";
  Rcpp::Rcout << "   - Setup (E,S) : " << t_Sigma_Setup << "\n";
  Rcpp::Rcout << "   - Vectors (V) : " << t_Sigma_Vec << "\n";
  Rcpp::Rcout << "   - Values (d)  : " << t_Sigma_Val << "\n";
  Rcpp::Rcout << "============================================\n";
  
  return Rcpp::List::create(
    Named("Gamma") = out_Gamma,
    Named("Lambda") = out_Lambda,
    Named("Sigma") = out_Sigma
  );
}