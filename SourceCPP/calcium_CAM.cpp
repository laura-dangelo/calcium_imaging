#include <RcppArmadillo.h>
#include <math.h>
#include <time.h>
#include <mvnorm.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>
 
using namespace Rcpp;


// log-likelihood
/*
 * y = (y_1, ... , y_n)  vec length = n
 * cc = (c_0, ..., c_n-1, c_n)  vec length = n + 1
 * AA = (A_0 = 0, A_1, ... , A_n)
 */
double loglik(const arma::vec& y, const arma::vec& cc, const arma::vec& AA, 
              double & b, double & gamma, double & sigma2, double & tau2)
{
  int n = y.n_elem;
  arma::vec llik(n);
  
  for(int k = 0; k < n; k++) { 
    llik(k) = R::dnorm(y(k), b + gamma * cc(k) + AA(k+1), std::sqrt(sigma2 + tau2), true) ;
  }
  return arma::accu(llik); 
}

// prior on gamma
/*
 * gamma ~ Beta(hyp_gamma1, hyp_gamma2)
 */
double logprior_gamma(double & gamma, double & hyp_gamma1, double & hyp_gamma2) 
{
  double out;
  out = R::dbeta( gamma, hyp_gamma1, hyp_gamma2, true);
  return(out);
}

// log-posterior (MH step on gamma)
double logpost_gamma(const arma::vec& y, const arma::vec& cc, const arma::vec& AA, 
                     double & b, double & gamma, double & sigma2, double & tau2, 
                     double & hyp_gamma1, double & hyp_gamma2)
{
  double out;
  out = loglik(y, cc, AA, b, gamma, sigma2, tau2) + logprior_gamma(gamma, hyp_gamma1, hyp_gamma2) ;
  return(out);
}

// prior on A for non-zero values
/*
 * A ~ G0 = Gamma(hyp_A1, hyp_A2) 
 * the mean is hyp_A1/hyp_A2
 */
double logprior_A(double & A, double & hyp_A1, double & hyp_A2) 
{
  double out;
  out = R::dgamma( A, hyp_A1, 1/hyp_A2, true );
  return(out);
}

// log-posterior on A
// input of the function are the y_i s.t. A_i = A_k for k = 1,..., K (fixed).
double logpost_A(const arma::vec& y, const arma::vec& cc, double & A, 
                 double & b, double & gamma, double & sigma2, double & tau2, 
                 double & hyp_A1, double & hyp_A2)
{
  double out ;
  arma::vec AAA(y.n_elem + 1) ; AAA.fill(A) ;
  out = loglik(y, cc, AAA, b, gamma, sigma2, tau2) + logprior_A(A, hyp_A1, hyp_A2) ;
  return(out) ;
}

// sampling from prior on A (spike and slab)
double sample_mix(double & p, double & hyp_A1, double & hyp_A2)
{
  double out = 0 ;
  int k = Rcpp::rbinom( 1, 1, p )[0] ;
  if( k == 1 ) {
    out = R::rgamma(hyp_A1, 1/hyp_A2) ;
  }
  return(out) ;
}



/* 
 * functions for slice sampler
 */

// compute coefficients xi_D = (xi_1, xi_2, ...)
arma::vec fun_xi_D(double & kappa_D, int & max_xiK)
{
  arma::vec out(max_xiK);
  for(int k = 0; k < max_xiK; k++)
  {
    out(k) = (1-kappa_D) * pow(kappa_D, k) ;
  }
  return(out) ;
}

// compute coefficients xi_O = (xi_1, xi_2, ...)
arma::vec fun_xi_O(double & kappa_O, int & max_xiL)
{
  arma::vec out(max_xiL);
  for(int l = 0; l < max_xiL; l++)
  {
    out(l) = (1-kappa_O) * pow(kappa_O, l) ;
  }
  return(out) ;
}

// compute stick-breaking weights starting from the vector of beta r.v. Beta(a_k,b_k)
arma::vec stick_breaking(arma::vec beta_var)
{
  int len = beta_var.n_elem ;
  arma::vec out(len) ;
  
  out(0) = beta_var(0) ;
  for(int k = 1; k < len; k++)
  {
    out(k) = beta_var(k) * arma::prod( 1 - beta_var.head(k) ) ;
  }
  return(out) ;
}



/*
 * Slice sampler
 */
Rcpp::List slice_sampler(const arma::vec& y, const arma::vec& g, 
                          arma::vec clusterD, arma::vec clusterO, 
                          const arma::vec& cc,
                          arma::vec A, double b, double gamma, 
                          double p,
                          double sigma2, double tau2,
                          double & alpha, double & beta, 
                          double & hyp_A1, double & hyp_A2,  
                          double & kappa_D, double & kappa_O, // 0.5 per far sì che sia il valor medio
                          const arma::vec& xi_D, const arma::vec& xi_O, 
                          double & eps_A,
                          int & check)
{
  /*
   * clusterD (lunghezza = J) cluster allocation on the distributions: takes value in {1,2,3,...}
   * clusterO (lunghezza = T) cluster allocation on the observations: takes value in {0,1,2,3,...}, where 0 is a non-spike
   * 
   * xi_O = (xi^O_1, xi^O_2, ...)
   * xi_D = (xi^D_1, xi^D_2, ...)
   */
 
  int T = y.n_elem ;
  arma::vec clusterD_long(T) ;
  for(int t = 0; t < T; t++) { clusterD_long(t) = clusterD( g(t)-1 ) ; }
  
  arma::vec unique_g = arma::unique(g) ;
  int J = unique_g.n_elem ;
  
  arma::vec u_O(T) ;
  arma::vec u_D(J) ;
  
  double oldA ;
  double newA ;
  double ratio ;
 
  // step 1: sample latent uniform on the distributions
  for(int j = 0; j < J; j++)
  {
    u_D(j) = R::runif( 0, xi_D(clusterD(j)-1) ) ;
  }
  arma::vec maxK_j(J) ;
  maxK_j = 1 + arma::floor( (log(u_D) - log(1 - kappa_D)) / log(kappa_D) ) ;
  int maxK = std::max(clusterD.max(), maxK_j.max()) ; // upper bound on the number of distributions

  
  // step 2: sample latent uniform on the observations
  for(int t = 0; t < T; t++)
  {
    u_O(t) = R::runif( 0, xi_O(clusterO(t)) ) ;
  } 
  arma::vec maxL_t(T) ;
  maxL_t = 1 + arma::floor( (log(u_O) - log(1 - kappa_O)) / log(kappa_O) ) ;
  int maxL = std::max(clusterO.max(), maxL_t.max()) ; // upper bound on the number of clusters on the observations
 
  
  // step 3: sample the stick-breaking weights on the distributions
  arma::vec v_k(maxK) ;
  arma::vec pi_k(maxK) ;
  for(int k = 1; k < maxK + 1; k++)
  {
    arma::uvec ind_k = find(clusterD == k) ;
    arma::uvec ind_mk = find(clusterD > k) ;
    
    int a_k = 1 + ind_k.n_elem ;
    double b_k = alpha + ind_mk.n_elem ;
    v_k(k - 1) = R::rbeta(a_k, b_k) ;
  }
  pi_k = stick_breaking( v_k ) ;

  
  // step 4: sample the stick-breaking weights on the observations
  // for each k in {1,...,maxK} we have maxL weights: matrix(maxL, maxK)
  arma::vec v_lk(maxL) ;
  arma::mat omega_lk(maxL, maxK) ;
  for(int k = 1; k < maxK + 1; k++)
  {
    arma::uvec ind_clusterD_k = find(clusterD_long == k) ; // indici delle y_t t.c. clusterD_t = k
    arma::vec subcluster = clusterO.elem( ind_clusterD_k ) ; // labels dei cluster sulle osservazioni per le y ~ G*k
      
    for(int l = 0; l < maxL; l++)
    {
      int a_lk = 1 ;
      double b_lk = alpha ;
      
      if( subcluster.n_elem > 0 )
      {
        arma::uvec ind_lk = find(subcluster == l) ; 
        arma::uvec ind_mlk = find(subcluster > l) ;
        
        a_lk = 1 + ind_lk.n_elem ;
        b_lk = beta + ind_mlk.n_elem ;
      }
      
      v_lk(l) = R::rbeta(a_lk, b_lk) ;
    }
     omega_lk.col(k-1) = stick_breaking( v_lk ) ;
  }
  
  
  // step 5: sample the distributional cluster indicator
  NumericVector probK(maxK) ;
  IntegerVector clusterD_id =  Rcpp::seq(1, maxK);
  
  
  for(int j = 0; j < J; j++)
  {
    arma::uvec ind_t = find(g == (j+1)) ;
    arma::vec mixcomp(ind_t.n_elem) ;
    for(int k = 0; k < maxK; k++)
    {
      probK[k] = -999 ;
      if( u_D(j) < xi_D(k) )
      {
        arma::vec omega_col = omega_lk.col(k) ;
        for(int t = 0; t < ind_t.n_elem; t++)
        {
          int cl = clusterO( ind_t(t) ) ;
          mixcomp(t) = log( omega_col(cl) ) ;
        }
        probK[k] = log( pi_k(k) ) - log( xi_D(k) ) + arma::accu(mixcomp) ;
      }
    }
    probK =  probK - max(probK)  ;
    probK =  exp(probK)  ;
    clusterD(j) = Rcpp::sample(clusterD_id, 1, false, probK)[0] ;
  }
  

  // step 6: sample the observational cluster indicator
  NumericVector probL(maxL) ;
  IntegerVector clusterO_id = Rcpp::seq(0, maxL-1) ;
  for(int t = 0; t < T; t++)
  {
    for(int l = 0; l < maxL; l++)
    {
      probL[l] = 0 ;
      if( u_O(t) < xi_O(l) )
      {
        probL[l] = omega_lk(l, clusterD(g(t)-1)-1) / xi_O(l) * R::dnorm(y(t) - b - gamma * cc(t), A(l), std::sqrt(sigma2 + tau2), false) ;
      }
    }
    
    clusterO(t) = Rcpp::sample(clusterO_id, 1, false, probL)[0] ;
    if( A(clusterO(t)) == 0 ) { clusterO(t) = 0 ; }
  }

  
  // step 7: sample the cluster parameters
  for(int l = 1; l < maxL; l++)
  {
    arma::uvec ind_l = find( clusterO == l ) ;
    if( ind_l.n_elem > 0 )
    {
      arma::vec sub_y = y(ind_l) ;
      // MH step su A(l)
      oldA = A(l) ;
      newA = oldA ;
      
      newA = oldA + R::runif(-eps_A, eps_A) ;
      ratio = exp( logpost_A(sub_y, cc(ind_l), 
                             newA, b, 
                             gamma, sigma2, tau2,
                             hyp_A1, hyp_A2) -
                    logpost_A(sub_y, cc(ind_l), 
                              oldA, b, 
                              gamma, sigma2, tau2,
                              hyp_A1, hyp_A2) ) ;
      if(R::runif(0, 1) < ratio) oldA = newA ;
      A(l) = oldA ; 
    } else {
      // sample from the prior
      A(l) = sample_mix(p, hyp_A1, hyp_A2) ;
    }
  }
 
      
  return Rcpp::List::create(Rcpp::Named("clusterO") = clusterO,
                            Rcpp::Named("clusterD") = clusterD,
                            Rcpp::Named("A") = A);
}


// Gibbs sampler: main function
// [[Rcpp::export]]
Rcpp::List calcium_gibbs(int Nrep, 
                         arma::vec y, arma::vec g,
                         arma::vec cal,
                         arma::vec clO, arma::vec clD, 
                         arma::vec A_start,
                         double b_start, double gamma_start, 
                         double sigma2_start, 
                         double tau2_start, // varianza su equazione del calcio
                         double p_start,
                         double alpha, // concentration param DP on distributions
                         double beta, // concentration param DP on observations
                         int max_xiK, int max_xiL, // max length of deterministic sequences
                         double kappa_D, double kappa_O, // parameters of the deterministic sequences
                         double c0, double varC0, // media e varianza c0 (calcio istante 0)
                         double hyp_A1, double hyp_A2, // shape e rate Gamma su A
                         double hyp_b1, double hyp_b2, // media e varianza Normale su b
                         double hyp_sigma21, double hyp_sigma22, // shape e rate Gamma su 1/sigma2 (precision)
                         double hyp_tau21, double hyp_tau22, // shape e rate Gamma su 1/tau2 (precision state equation)
                         double hyp_gamma1, double hyp_gamma2, // shape1 e shape2 Beta su gamma (parametro AR)
                         double hyp_p1, double hyp_p2, // shape1 e shape2 Beta su p (prob di spike)
                         double eps_gamma, // MH step size
                         double eps_A) 
{
  // allocate output matrices
  int n = y.n_elem ;
  arma::vec unique_g = arma::unique(g) ;
  int J = unique_g.n_elem ;
  
  arma::mat out_c(n+1, Nrep) ; // 0 1 ... n
  arma::mat out_A = arma::zeros(A_start.n_elem, Nrep) ;
  arma::vec out_b(Nrep) ;
  arma::vec out_gamma = arma::zeros(Nrep) ;
  arma::vec out_sigma2 = arma::zeros(Nrep) ;
  arma::vec out_tau2 = arma::zeros(Nrep) ;
  arma::mat clusterO = arma::zeros(n, Nrep)  ;
  arma::mat clusterD = arma::zeros(J, Nrep)  ;
  arma::vec out_p(Nrep) ;
  
  // initialize the chains
  out_c(0,0) = c0 ;
  out_b(0) = b_start ;
  out_gamma(0) = gamma_start ;
  out_sigma2(0) = sigma2_start ;
  out_tau2(0) = tau2_start ;
  out_p(0) = p_start ;
  
  // kalman filter quantities
  out_c.col(0) = cal;
  arma::vec filter_mean(n+1) ; // 0 1 ... n
  arma::vec filter_var(n+1) ; // 0 1 ... n
  arma::vec R(n+1) ; arma::vec a(n+1) ; // 0 1 ... n
  double back_mean; double back_var ;
  
  // MH quantities
  double oldgamma ; double newgamma ;
  double ratio ;
  
  // CAM quantities
  double n_clus ;
  arma::vec AA = arma::zeros(n+1) ;
  Rcpp::List out_slice ;
  
  clusterO.col(0) = clO ;
  clusterD.col(0) = clD ;
  out_A.col(0) = A_start; 
  
  arma::vec xi_D(max_xiK) ;
  arma::vec xi_O(max_xiL) ;
  xi_D = fun_xi_D(kappa_D, max_xiK) ;
  xi_O = fun_xi_O(kappa_O, max_xiL) ;

  
  bool display_progress = true ;
  Progress p(Nrep, display_progress) ;
  
  //// Main loop ////
  for(int i = 0; i < Nrep -1 ; i++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
    int check = 0 ;
    
    AA(0) = 0 ;
    for(int j = 1; j < n+1; j++) { AA(j) =  out_A(clusterO(j-1,i), i) ; }
    
    /*
     * Sampling calcium level
     * Kalman filter + backward sampling
     */
    a(0) = 0 ; // a0
    R(0) = varC0 ; // R0
    filter_mean(0) = 0 ;
    filter_var(0) = varC0 ;
    
    for(int j = 1; j < n +1 ; j++)
    {
      a(j) = out_gamma(i) * filter_mean(j-1) + AA(j);
      R(j) = out_gamma(i) * out_gamma(i) * filter_var(j-1) + out_tau2(i) ;
      
      filter_mean(j) = a(j) + R(j) / (R(j) + out_sigma2(i) ) * (y(j-1) - out_b(i) - a(j)) ;
      filter_var(j) = out_sigma2(i) * R(j) / (R(j) + out_sigma2(i)) ;
    }
    out_c(n, i+1) = R::rnorm(filter_mean(n), filter_var(n)) ;
    
    for(int j = n-1; j > -1; j--)
    {
      back_mean = filter_mean(j) + out_gamma(i) * filter_var(j) / R(j+1) * (out_c(j+1, i+1) - a(j+1)) ;
      back_var = filter_var(j) - pow(out_gamma(i) * filter_var(j), 2) / R(j+1) ;
      
      out_c(j, i+1) = R::rnorm(back_mean, back_var) ;
    }
    
    
    /*
     * Sampling b
     */
    arma::vec z(n) ; 
    for(int j = 0; j < n; j++) { z(j) = y(j) - out_c(j+1, i+1) ; } 
    out_b(i+1) = R::rnorm( (hyp_b2 * hyp_b1 + 1/out_sigma2(i) * arma::accu(z)) / (hyp_b2 + n / out_sigma2(i)) , 
          std::sqrt( 1/ (hyp_b2 + n / out_sigma2(i)) ) ) ;
    
    
    /*
     * Sampling sigma2
     */
    arma::vec sq(n) ; 
    for(int j = 0; j < n; j++) { sq(j) = (z(j) - out_b(i+1)) * (z(j) - out_b(i+1)) ; }
    out_sigma2(i+1) = 1/ R::rgamma(hyp_sigma21 + n/2, 1/(hyp_sigma22 + 0.5 * arma::accu(sq)) ) ;
    
    
    /*
     * Sampling tau2
     */
    arma::vec sq2(n) ; 
    for(int j = 0; j < n ; j++) { sq2(j) = (out_c(j+1,i+1) - out_gamma(i) * out_c(j,i+1) - AA(j+1)) * 
       (out_c(j+1,i+1) - out_gamma(i) * out_c(j,i+1) - AA(j+1)) ; }
    out_tau2(i+1) = 1/ R::rgamma(hyp_tau21 + n/2, 1/(hyp_tau22 + 0.5 * arma::accu(sq2)) ) ;
    
    
    /*
     * Sampling gamma
     * MH step with uniform random walk
     */
    oldgamma = out_gamma(i) ;
    newgamma = oldgamma + R::runif(-eps_gamma, eps_gamma) ;
    ratio = exp( logpost_gamma(y, out_c.col(i+1), 
                         AA, out_b(i+1), 
                         newgamma, out_sigma2(i+1), out_tau2(i+1),
                         hyp_gamma1, hyp_gamma2) -
                  logpost_gamma(y, out_c.col(i+1), 
                          AA, out_b(i+1), 
                          oldgamma, out_sigma2(i+1), out_tau2(i+1),
                          hyp_gamma1, hyp_gamma2) ) ;
    if(R::runif(0, 1) < ratio) oldgamma = newgamma ;
    out_gamma(i+1) = oldgamma ;
    
    
    /*
     * Sampling of clusters and cluster parameters
     * Slice sampler
     */
    out_slice = slice_sampler(y, g, 
                              clusterD.col(i), clusterO.col(i), 
                              out_c.col(i+1),
                              out_A.col(i), out_b(i+1), out_gamma(i+1), 
                              out_p(i),
                              out_sigma2(i+1), out_tau2(i+1),
                              alpha, beta, 
                              hyp_A1, hyp_A2,  
                              kappa_D, kappa_O, 
                              xi_D, xi_O, 
                              eps_A,
                              check) ;

    arma::vec out_slice_clusO = out_slice["clusterO"] ;
    clusterO.col(i+1) = out_slice_clusO ;
    
    arma::vec out_slice_clusD = out_slice["clusterD"] ;
    clusterD.col(i+1) = out_slice_clusD ;
    
    arma::vec out_slice_A = out_slice["A"] ;
    out_A.col(i+1) = out_slice_A ;
    

    /*
    * Sampling p (spike and slab proportion)
    */
    arma::vec line(n) ; line = clusterO.col(i+1) ;
    double n0 = std::count(line.begin(), line.end(), 0) ;
      
    out_p(i+1) = R::rbeta(hyp_p2 + n - n0, hyp_p1 + n0) ;
    // if(out_p(i+1) > 0.4) { check = 1 ; }
    
    
    //// END Gibbs sampler ////
    if(check == 1) { 
      Rcout << "Stop at iter. " << i << "\n" ;
      i = Nrep - 2 ; 
    }
  }
  return Rcpp::List::create(Rcpp::Named("calcium") = out_c,
                            Rcpp::Named("A") = out_A,
                            Rcpp::Named("b") = out_b,
                            Rcpp::Named("gamma") = out_gamma,
                            Rcpp::Named("sigma2") = out_sigma2,
                            Rcpp::Named("clusterO") = clusterO,
                            Rcpp::Named("clusterD") = clusterD,
                            Rcpp::Named("p") = out_p,
                            Rcpp::Named("tau2") = out_tau2);
}



