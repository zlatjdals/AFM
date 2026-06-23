#include <RcppArmadillo.h>
#include <chrono>
using namespace Rcpp;
using namespace arma;
using namespace std;

using std::sort;

// [[Rcpp::depends(RcppArmadillo)]]
arma::mat rbing(arma::vec Lam, int t, bool ratio_print) {
  // Sort eigenvalues in descending order and add zero
  Lam.insert_rows(Lam.n_rows, 1);  // Add a zero at the end
  
  const int qa = Lam.n_elem;
  
  // Initialize matrices and vectors
  arma::vec X(qa, fill::zeros);  // Output matrix
  arma::vec sigacginv = 1 + 2 * Lam;
  arma::vec SigACG = sqrt(1 / sigacginv);
  arma::vec y2(qa), yp(qa), y(qa);
  
  const double qa2 = qa / 2.0;
  const double tmp = -qa2 * log(qa) + 0.5 * (qa - 1);
  double lratio;
  
  // Rejection sampling loop
  for (int i = 0; i < 1;) {
    // Generate normal samples and normalize them
    for (int j = 0; j < qa; ++j) {
      yp[j] = R::rnorm(0, SigACG[j]);
    }
    y = yp / sqrt(sum(square(yp)));  // Normalize
    y2 = square(y);
    
    // Compute the log acceptance ratio
    lratio = -dot(y2, Lam) + tmp + qa2 * log(dot(y2, sigacginv));
    if(ratio_print == TRUE){
      Rcpp::Rcout << "Log Acceptance rate (rbing) : " << lratio << std::endl ;
    }
    
    // Accept or reject the sample
    if (log(R::runif(0, 1)) < lratio) {
      X = y;
      i++ ;
      
    }
    Rcpp::checkUserInterrupt();
  }
  
  if (t % 1000 == 0) {
    Rcpp::Rcout << "Log Acceptance Rate (Bingham) : " << lratio << std::endl;
  }
  
  return X;
}

// [[Rcpp::export]]
List new_algorithm_FM(int p, int n, const arma::mat& x, double& h, arma::vec a_vec, 
                      const arma::mat& Gamma0, int iter, int burnin, int thin, 
                      bool modified, int k, double b0, double b1, double sigma2_0, arma::vec Lambda_0, bool ratio_print) {
  auto start_time = std::chrono::steady_clock::now();
  auto while_time_Lambda = std::chrono::duration<double>::zero();
  auto while_time_Gamma = std::chrono::duration<double>::zero();
  auto while_time_bing = std::chrono::duration<double>::zero();
  
  arma::mat S = x.t() * x;
  arma::vec eigen_S;
  arma::eig_sym(eigen_S, S / n);
  // Reverse the order of eigenvalues to be in descending order
  eigen_S = arma::reverse(eigen_S);
  
  if (modified && p>n) {
    double mean_non_sp = arma::mean(eigen_S.subvec(k,n-1));
    
    a_vec.subvec(0, k-1) = n * (mean_non_sp+h/n) / (eigen_S.subvec(0, k-1) - mean_non_sp) / 2 + 2;
    
  }
  
  arma::mat eigen_diag = diagmat(eigen_S);
  arma::mat eigen_partial_sqrt = diagmat(sqrt(eigen_S.subvec(0,n-1)));
  arma::vec r_vec = a_vec + n / 2.0;
  
  
  arma::mat H_0 = h * eye<mat>(p, p) + S;
  
  
  //std::vector<arma::mat> mod_Gamma_list;
  std::vector<arma::mat> Gamma_list;
  std::vector<arma::vec> Lambda_list;
  std::vector<double> sigma2_list;
  
  arma::vec eigval_H_0;
  arma::mat eigvec_H_0;
  arma::eig_sym(eigval_H_0, eigvec_H_0, H_0);
  
  // Reverse the order of eigenvalues and eigenvectors to be in descending order
  eigval_H_0 = arma::reverse(eigval_H_0);
  eigvec_H_0 = arma::fliplr(eigvec_H_0);
  
  arma::mat mod_Gamma_t = eigvec_H_0.t() * Gamma0;
  
  double sigma2_new = sigma2_0 ;
  double sigma2_mh = sigma2_0 ;
  
  arma::vec Lambda_new = Lambda_0 ;
  arma::vec Lambda_mh = Lambda_0 ;
  
  
  for (int t = 0; t < iter; ++t) {
    /////// Sampling Lambda
    arma::vec c = arma::sum(mod_Gamma_t.each_col() % eigval_H_0 % mod_Gamma_t, 0).t() / 2;
  
    auto while_Lambda_start = std::chrono::steady_clock::now();
    
    for (int i = 0; i < k; ++i) {
      Lambda_new( i ) = 1 / R::rgamma( r_vec( i ) - 1, 1 / c( i ) );
    }
    
    double log_accept_ratio = ( n / 2.0 ) * log( prod( Lambda_new ) / prod( Lambda_new + sigma2_new ) ) ;
    
    while ( log( R::runif( 0 , 1 ) ) > log_accept_ratio ) {
      for (int i = 0; i < k; ++i) {
        Lambda_new( i ) = 1 / R::rgamma( r_vec( i ) - 1, 1 / c( i ) );
      }
      log_accept_ratio = ( n / 2.0 ) * log( prod( Lambda_new ) / prod( Lambda_new + sigma2_new ) ) ;
      
      if(ratio_print == TRUE){
        Rcpp::Rcout << "Log Acceptance rate (lambda) : " << log_accept_ratio << std::endl ;
      }
      Rcpp::checkUserInterrupt();
    }
    auto while_Lambda_end = std::chrono::steady_clock::now();
    while_time_Lambda += (while_Lambda_end - while_Lambda_start);
    
    if (t % 1000 == 0) {
      Rcpp::Rcout << "Log Acceptance rate (lambda) : " << log_accept_ratio << std::endl ;
    }
    

    
    /////// Sampling Gamma
    
    arma::uvec idx_perm = arma::randperm(k, k) ;
    arma::vec B = - ( n / 2.0 ) * ( 1 / (Lambda_new + sigma2_new) - 1 / sigma2_new ) ;
    
    arma::mat N , Z , Y ;
    int idx ;
    arma::mat Gamma_temp ;
    
    arma::vec eigval_NtWN;
    arma::mat eigvec_NtWN;
    
    arma::vec neg_eigval_NtWN ;           // Vector with signs reversed
    arma::vec sorted_neg_eigval ; // Sorted vector
    arma::vec eigval_minus_NtWN ;
    
    
    for (int i = 0; i < k ;++i ) {
      idx = idx_perm(i) ;
      
      // orthogonal basis for the null space of Gamma[,-idx]

      if (idx > 0 & idx < k - 1) {
        // Combine columns before and after idx
        Gamma_temp = join_horiz(
          mod_Gamma_t.cols(0, idx - 1),
          mod_Gamma_t.cols(idx + 1, mod_Gamma_t.n_cols - 1)
        );
      } else if (idx == k - 1) {
        Gamma_temp = mod_Gamma_t.cols(0, idx - 1) ;
      } else {
        // If idx is 0, keep only the columns after idx
        Gamma_temp = mod_Gamma_t.cols(1, mod_Gamma_t.n_cols - 1);
      }
      
      N = null(Gamma_temp.t()) ;
      
      
      auto while_Gamma_start = std::chrono::steady_clock::now();
      
      arma::mat scaled_N = B(idx) * N.t() * eigen_diag * N ; // 5s/90s
      
      arma::eig_sym( eigval_NtWN , eigvec_NtWN , scaled_N ) ; 
      
      
      auto while_Gamma_end = std::chrono::steady_clock::now();
      while_time_Gamma += (while_Gamma_end - while_Gamma_start);

      sorted_neg_eigval = -eigval_NtWN; // Sorted vector
  
      eigval_minus_NtWN = sorted_neg_eigval.subvec(0, eigval_NtWN.n_elem - 2) - sorted_neg_eigval(eigval_NtWN.n_elem - 1);
      
      
      auto while_bing_start = std::chrono::steady_clock::now();
      
      Y = rbing(eigval_minus_NtWN , t , ratio_print) ; // 9s/90s
      auto while_bing_end = std::chrono::steady_clock::now();
      
      while_time_bing += (while_bing_end - while_bing_start);

      
      mod_Gamma_t.col(idx) = N * eigvec_NtWN * Y ;
      
    }
    

    
    if (t >= burnin && t % thin == 0) {
      arma::uvec ord = sort_index(Lambda_new, "descend");   // Alternatively, sort by s_j if desired
      Lambda_list.push_back(Lambda_new(ord));
      Gamma_list.push_back( eigvec_H_0 * mod_Gamma_t.cols(ord) ) ;
    }
    
    /////// Sampling sigma^2
    
    arma::mat GtWG = mod_Gamma_t.t() * eigen_diag * mod_Gamma_t ;
    double trace_GtWG = trace(GtWG);
    
    double IW_beta =  1 / (( sum( eigen_S ) - trace_GtWG  ) * n / 2.0 );
    
    sigma2_mh = 1 / R::rgamma( n * ( p - k ) / 2.0 , IW_beta ) ; 
    
    
    while(sigma2_mh < b0 || sigma2_mh > b1){
      sigma2_mh = 1 / R::rgamma( n * ( p - k ) / 2.0 , IW_beta ) ;
      Rcpp::checkUserInterrupt();
    }
    
    log_accept_ratio = ( n / 2.0 ) * log( prod( Lambda_new + sigma2_new ) / prod( Lambda_new + sigma2_mh ) ) 
      + ( n / 2.0 ) * trace( diagmat( 1 / ( Lambda_new + sigma2_new  ) -1 / ( Lambda_new + sigma2_mh  ) ) * GtWG ) ;
    
    if(log( R::runif( 0 , 1 ) ) < log_accept_ratio){
      sigma2_new = sigma2_mh ; 
    }    

    
    if (t % 1000 == 0) {
      Rcpp::Rcout << "Log Acceptance Rate (sigma2) : " << log_accept_ratio << std::endl;
    }
    
    if (t >= burnin && t % thin == 0) {
      sigma2_list.push_back( sigma2_new ) ;
    }
    
    if (t % 1000 == 0) {
      Rcpp::Rcout << "Iteration : " << t << std::endl;
      auto end_time = std::chrono::steady_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
      Rcpp::Rcout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
      
    }
    Rcpp::checkUserInterrupt(); 
  }
  
  Rcpp::Rcout << "Iteration : " << iter << std::endl;
  auto end_time = std::chrono::steady_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
  Rcpp::Rcout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
  
  Rcpp::Rcout << "Time spent in Lambda while : " << while_time_Lambda.count() << " seconds" << std::endl;
  Rcpp::Rcout << "Time spent in Gamma while : " << while_time_Gamma.count() << " seconds" << std::endl;
  Rcpp::Rcout << "Time spent in Bingham while : " << while_time_bing.count() << " seconds" << std::endl;
  
  return List::create(Named("Lambda") = Lambda_list, 
                      Named("Gamma") = Gamma_list,
                      Named("sigma2") = sigma2_list);
}

// [[Rcpp::export]]
List mat_result(arma::vec matr){
  arma::mat matrix(2,2,fill::zeros) ;
  matrix.col(0) = matr.subvec(0,1) ;
  matrix.col(1) = matr.subvec(2,3) ;
  
  std::vector<arma::mat> matrix_list ;
  for (int i = 0; i < 5 ; ++i ) {
    matrix_list.push_back( matrix ) ;
  }
  
  return List::create(Named("mat") = matrix_list);
}

// [[Rcpp::export]]
List new_algorithm_gSIW(int p, int n, const arma::mat& x, double& h, arma::vec a_vec, 
                       const arma::mat& Gamma0, int iter, int burnin, int thin, 
                       bool modified, int k) {
  auto start_time = std::chrono::steady_clock::now();
  
  arma::mat S = x.t() * x;
  if (modified && p>n) {
    arma::vec eigen_S;
    arma::eig_sym(eigen_S, S / n);
    // Reverse the order of eigenvalues to be in descending order
    eigen_S = arma::reverse(eigen_S);
    double mean_non_sp = arma::mean(eigen_S.subvec(k,n-1));
    
    a_vec.subvec(0, k-1) = n * (mean_non_sp+h/n) / (eigen_S.subvec(0, k-1) - mean_non_sp  )/2 + 2;
    a_vec.subvec(k, n-1).fill(p/2);
    a_vec.subvec(n, p-1).fill(p*2);
  }
  arma::vec r_vec = a_vec + n / 2;
  
  arma::mat H_0 = h * diagmat(arma::ones<arma::vec>(p)) + S;
  
  //std::vector<arma::mat> mod_Gamma_list;
  std::vector<arma::mat> Gamma_list;
  std::vector<arma::vec> Lambda_list;
  
  arma::vec eigval_H_0;
  arma::mat eigvec_H_0;
  arma::eig_sym(eigval_H_0, eigvec_H_0, H_0);
  
  // Reverse the order of eigenvalues and eigenvectors to be in descending order
  eigval_H_0 = arma::reverse(eigval_H_0);
  eigvec_H_0 = arma::fliplr(eigvec_H_0);
  
  arma::mat mod_Gamma_t = eigvec_H_0.t() * Gamma0;
  
  for (int t = 0; t < iter; ++t) {
    arma::vec c = arma::sum(mod_Gamma_t.each_col() % eigval_H_0 % mod_Gamma_t, 0).t() / 2;
    
    arma::vec Lambda_new(p);
    
    for (int i = 0; i < p; ++i) {
      Lambda_new(i) = 1 / R::rgamma(r_vec(i) - 1, 1 / c(i));
    }
    
    
    arma::uvec var_idx = arma::regspace<arma::uvec>(0, p-1);
    while (!var_idx.is_empty()) {
      arma::uvec two_idx;
      if (var_idx.n_elem == 1) {
        arma::uvec remaining_idx = arma::regspace<arma::uvec>(0, p - 1);
        remaining_idx.shed_rows(var_idx); // Remove var_idx from the remaining indices
        two_idx = join_vert(var_idx, remaining_idx(arma::randi<arma::uvec>(1, arma::distr_param(0, remaining_idx.n_elem - 1))));
      } else {
        two_idx = arma::randperm(var_idx.n_elem, 2);
      }
      two_idx = sort(two_idx);
      
      arma::mat sub_mod_Gamma = mod_Gamma_t.rows(two_idx);
      arma::mat sub_mat = sub_mod_Gamma.each_row() / Lambda_new.t() * sub_mod_Gamma.t();
      
      arma::vec eigval_sub_mat;
      arma::mat eigvec_sub_mat;
      arma::eig_sym(eigval_sub_mat, eigvec_sub_mat, sub_mat);
      
      
      double diff_s = eigval_sub_mat(1) - eigval_sub_mat(0);
      double diff_h = eigval_H_0(two_idx(1)) - eigval_H_0(two_idx(0));
      
      double c0 = -0.5 * std::abs(diff_s * diff_h);
      
      
      double omega = std::atan(eigvec_sub_mat(1, 1)/ eigvec_sub_mat(0, 1)); // Eigenvalues are returned in ascending order
      
      
      double alpha0 = R::rbeta(0.5, 0.5);
      while (R::runif(0, 1) > std::exp(c0 * alpha0)) {
        alpha0 = R::rbeta(0.5, 0.5);
      }
      
      
      
      double theta = std::acos(2 * alpha0 - 1) / 2 - omega;
      if (theta < -M_PI / 2) {
        theta += M_PI;
      } else if (theta > M_PI / 2) {
        theta -= M_PI;
      }
      
      arma::vec random_signs = 2 * arma::randi<arma::vec>(2, distr_param(0, 1)) - 1;
      arma::mat rot = {{std::cos(theta), -std::sin(theta)}, {std::sin(theta), std::cos(theta)}};
      mod_Gamma_t.rows(two_idx) = arma::diagmat(random_signs) * rot * mod_Gamma_t.rows(two_idx);
      
      if(var_idx.n_elem == 1){
        break;
      }
      var_idx.shed_rows(two_idx);
    }
    
    
    if (t >= burnin && t % thin == 0) {
      arma::uvec ord = sort_index(Lambda_new, "descend");   // Alternatively, sort by s_j if desired
      Lambda_list.push_back(Lambda_new(ord));
      Gamma_list.push_back( eigvec_H_0 * mod_Gamma_t.cols(ord) ) ;
    }
    
    if (t % 1000 == 0) {
      Rcpp::Rcout << "Iteration : " << t << std::endl;
      auto end_time = std::chrono::steady_clock::now();
      auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
      Rcpp::Rcout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
      
    }
    Rcpp::checkUserInterrupt();
  }
  
  Rcpp::Rcout << "Iteration : " << iter << std::endl;
  auto end_time = std::chrono::steady_clock::now();
  auto elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
  Rcpp::Rcout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
  
  
  return List::create(Named("Lambda") = Lambda_list, 
                      Named("Gamma") = Gamma_list);
}
