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
 * la media è hyp_A1/hyp_A2
 */
double logprior_A(double & A, double & hyp_A1, double & hyp_A2) 
{
  double out;
  out = R::dgamma( A, hyp_A1, 1/hyp_A2, true );
  return(out);
}

// log-posterior on A
// la funzione prende in input solo le y_i t.c. A_i = A_k per un certo k = 1,..., K (fissato).
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
 * funzioni per slice sampler
 */

// calcola coefficienti xi_D = (xi_1, xi_2, ...)
arma::vec xi_D(double & kappa_D, int & max_xiK)
{
  arma::vec out(max_xiK);
  for(int k = 0; k < max_xiK; k++)
  {
    out(k) = (1-kappa_D) * pow(kappa_D, k) ;
  }
  return(out) ;
}

// calcola coefficienti xi_O = (xi_1, xi_2, ...)
arma::vec xi_O(double & kappa_O, int & max_xiL)
{
  arma::vec out(max_xiL);
  for(int l = 0; l < max_xiL; l++)
  {
    out(l) = (1-kappa_O) * pow(kappa_O, l) ;
  }
  return(out) ;
}

// calcola pesi costruiti con stick breaking process a partire dal vettore di variabili Beta(a_k,b_k)
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
                          double & alpha, 
                          double & hyp_A1, double & hyp_A2,  
                          double & kappa_D, double & kappa_O,
                          const arma::vec& xi_D, const arma::vec& xi_O, 
                          int & check)
{
  /*
   * clusterD (lunghezza = J) allocazione del cluster sulle distribuzioni assume valori in {1,2,3,...}
   * clusterO (lunghezza = T) allocazione del cluster delle osservazioni assume valori in {0,1,2,3,...}, dove 0 indica lo spike
   */
  
  int T = y.n_elem ;
  arma::vec clusterD_long(T) ;
  for(int t = 0; t < T; t++) { clusterD_long = clusterD(g(t)) ; }
  
  arma::vec unique_g = arma::unique(g) ;
  int J = unique_g.n_elem ;
  
  arma::vec u_O(T) ;
  arma::vec u_D(J) ;
  
  // step 1: sample latent uniform on the distribution
  for(int j = 0; j < J; j++)
  {
    u_D(j) = R::runif( 0, xi_D(clusterD(j)-1) ) ;
  }
  arma::vec maxK_j(J) ;
  maxK_j = 1 + arma::floor( (log(u_D) - log(1 - kappa_D)) / log(kappa_D) ) ;
  int maxK = std::max(clusterD.max(), maxK_j.max()) ; // upper bound su numero di distribuzioni
  
  // step 2: sample latent uniform on the observations
  for(int t = 0; t < T; t++)
  {
    u_O(t) = R::runif( 0, xi_O(clusterO(t)) ) ;
  } 
  arma::vec maxL_t(T) ;
  maxL_t = 1 + arma::floor( (log(u_O) - log(1 - kappa_O)) / log(kappa_O) ) ;
  int maxL = std::max(clusterO.max(), maxL_t.max()) ; // upper bound su numero di parametri per ogni distribuzione
  
  // step 3: sample the stick-breaking weights on the distribution
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
      arma::uvec ind_lk = find(subcluster == l) ; //potrebbe creare problemi su label vuote
      arma::uvec ind_mlk = find(subcluster > l) ;
      
      int a_lk = 1 + ind_lk.n_elem ;
      double b_lk = alpha + ind_mlk.n_elem ;
      v_lk(l) = R::rbeta(a_lk, b_lk) ;
    }
     omega_lk.col(k-1) = stick_breaking( v_lk ) ;
  }
  
  return Rcpp::List::create(Rcpp::Named("cluster") = clusterO,
                            Rcpp::Named("A") = A);
}



// Gibbs sampler: funzione principale
// [[Rcpp::export]]
Rcpp::List calcium_gibbs(int Nrep, 
                         arma::vec y, arma::vec g,
                         arma::vec cal,
                         arma::vec cl, 
                         arma::vec A_start,
                         double b_start, double gamma_start, 
                         double sigma2_start, 
                         double tau2_start, // varianza su equazione del calcio
                         double p_start,
                         double c0, double varC0, // media e varianza c0 (calcio istante 0)
                         double alpha, // concentration param DP
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
  arma::mat out_c(n+1, Nrep) ; // 0 1 ... n
  arma::mat out_A = arma::zeros(A_start.n_elem, Nrep) ;
  arma::vec out_b(Nrep) ;
  arma::vec out_gamma = arma::zeros(Nrep) ;
  arma::vec out_sigma2 = arma::zeros(Nrep) ;
  arma::vec out_tau2 = arma::zeros(Nrep) ;
  arma::mat cluster = arma::zeros(n, Nrep)  ;
  arma::vec out_p(Nrep) ;
  
  // initialize the chain
  out_c(0,0) = c0 ;
  out_b(0) = b_start ;
  out_gamma(0) = gamma_start ;
  out_sigma2(0) = sigma2_start ;
  out_tau2(0) = tau2_start ;
  out_p(0) = p_start ;
  
  // other quantities
  arma::vec filter_mean(n+1) ; // 0 1 ... n
  arma::vec filter_var(n+1) ; // 0 1 ... n
  arma::vec R(n+1) ; arma::vec a(n+1) ; // 0 1 ... n
  double back_mean; double back_var ;
  
  double oldgamma ; double newgamma ;
  double ratio ;
  
  double n_clus ;
  arma::vec AA = arma::zeros(n+1) ;
  Rcpp::List polyaurn ;
  double oldA ; double newA ;
  
  cluster.col(0) = cl ;
  out_A.col(0) = A_start; 
  out_c.col(0) = cal;
  
  bool display_progress = true ;
  Progress p(Nrep, display_progress) ;
  
  for(int i = 0; i < Nrep -1 ; i++)
  {
    if( Progress::check_abort() ) { return -1.0 ; }
    p.increment();
    int check = 0 ;
    
    AA(0) = 0 ;
    for(int j = 1; j < n+1; j++) { AA(j) =  out_A(cluster(j-1,i), i) ; }
    
    // sampling calcium: kalman filter
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
    
    //  out_c.col(i+1) = out_c.col(i) ;
    
    
    
    // sampling b
    
    arma::vec z(n) ; 
    for(int j = 0; j < n; j++) { z(j) = y(j) - out_c(j+1, i+1) ; } 
    
    out_b(i+1) = R::rnorm( (hyp_b2 * hyp_b1 + 1/out_sigma2(i) * arma::accu(z)) / (hyp_b2 + n / out_sigma2(i)) , 
          std::sqrt( 1/ (hyp_b2 + n / out_sigma2(i)) ) ) ;
    //    out_b(i+1) = out_b(i) ;
    
    // sampling sigma2 
    
    arma::vec sq(n) ; 
    for(int j = 0; j < n; j++) { sq(j) = (z(j) - out_b(i+1)) * (z(j) - out_b(i+1)) ; }
    
    out_sigma2(i+1) = 1/ R::rgamma(hyp_sigma21 + n/2, 1/(hyp_sigma22 + 0.5 * arma::accu(sq)) ) ;
    
    //out_sigma2(i+1) = out_sigma2(i) ;
    
    arma::vec sq2(n) ; 
    for(int j = 0; j < n ; j++) { sq2(j) = (out_c(j+1,i+1) - out_gamma(i) * out_c(j,i+1) - AA(j+1)) * 
       (out_c(j+1,i+1) - out_gamma(i) * out_c(j,i+1) - AA(j+1)) ; }
    out_tau2(i+1) = 1/ R::rgamma(hyp_tau21 + n/2, 1/(hyp_tau22 + 0.5 * arma::accu(sq2)) ) ;
    //   out_tau2(i+1) = out_tau2(i) ;
    
    
    
    //MH per gamma: random walk
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
    
    
    // Sampling of clusters and cluster parameters
    // Polya-Urn
/*    polyaurn = polya_urn(y, cluster.col(i), out_c.col(i+1),
                          out_A.col(i), out_b(i+1), out_gamma(i+1), 
                          out_p(i),
                          out_sigma2(i+1), out_tau2(i+1),
                          alpha,
                          hyp_A1, hyp_A2, check) ; 
    arma::vec out_polya_clus = polyaurn["cluster"] ;
    cluster.col(i+1) = out_polya_clus ;
 */
    cluster.col(i+1) = cluster.col(i) ;
    
    // sampling of parameters A1,...,Ak
    arma::vec out_polya_A =out_A.col(i) ;
//    arma::vec out_polya_A = polyaurn["A"] ;
    out_A.col(i+1) = out_polya_A ;
    
    arma::vec line(n) ; line = cluster.col(i+1) ;
    n_clus = line.max() ;
    
    for(int l = 1; l < n_clus + 1; l++)
    {
      // MH step su A(l)
      oldA = out_A(l, i+1) ;
      newA = oldA ;
      
      arma::uvec idk = find( line == l ) ;
      arma::vec sub_y = y(idk) ;
      arma::vec col_c = out_c.col(i+1) ;
      
      newA = oldA + R::runif(-eps_A, eps_A) ;
      ratio = exp( logpost_A(sub_y, col_c(idk), 
                             newA, out_b(i+1), 
                             out_gamma(i+1), out_sigma2(i+1), out_tau2(i+1),
                             hyp_A1, hyp_A2) -
                    logpost_A(sub_y, col_c(idk), 
                              oldA, out_b(i+1), 
                              out_gamma(i+1), out_sigma2(i+1), out_tau2(i+1),
                              hyp_A1, hyp_A2) ) ;
      if(R::runif(0, 1) < ratio) oldA = newA ;
      out_A(l, i+1) = oldA ;
    }
    
    // Update p
    double n0 = std::count(line.begin(), line.end(), 0) ;
    out_p(i+1) = R::rbeta(hyp_p2 + n - n0, hyp_p1 + n0) ;
    //   out_p(i+1) = out_p(i) ;
    if(out_p(i+1) > 0.1) { check = 1 ; }
    
    
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
                            Rcpp::Named("cluster") = cluster,
                            Rcpp::Named("p") = out_p,
                            Rcpp::Named("tau2") = out_tau2
  );
}



