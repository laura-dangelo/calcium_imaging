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

// prior on A
double logprior_A(double & A, double & hyp_A1, double & hyp_A2) 
{
  double out;
  out = R::dgamma( A, hyp_A1, hyp_A2, true );
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


// log-posterior (MH step on A)
/*
 * A ~ G0 = Gamma(hyp_A1, hyp_A2)
 * la funzione prende in input solo le y_i t.c. A_i = A_k per un certo k = 1,..., K (fissato).
 */
double logpost_A(const arma::vec& y, const arma::vec& cc, double & A, 
                     double & b, double & gamma, double & sigma2, double & tau2, 
                     double & hyp_A1, double & hyp_A2)
{
  double out;
  arma::vec AAA(y.n_elem + 1) ; AAA.fill(A) ;
  out = loglik(y, cc, AAA, b, gamma, sigma2, tau2) + logprior_A(A, hyp_A1, hyp_A2) ;
  return(out);
}



// trova differenza tra due vettori
double missing_value(const arma::vec& x, const arma::vec& y)
{
  double out = 0;
  int n = x.n_elem ;
  for(int i = 0; i < n; i++)
  {
    double lookfor = x(i) ;
    if(std::count(y.begin(), y.end(), lookfor) == 0) { out = lookfor ; }
  }
  return(out) ;
}


// Polya-Urn
Rcpp::List polya_urn_nonc(const arma::vec& y, arma::vec cluster, const arma::vec& cc,
                    arma::vec A, double b, double gamma, 
                    double p,
                    double sigma2, double tau2,
                    double & alpha, int & m,
                    double & hyp_A1, double & hyp_A2,  
                    int & check)
{
  int n = y.n_elem ;
  
  int n_clus = 0 ; 
  int n_j ;
  int max_clus = 0 ;
  arma::vec A_tmp(A.n_elem) ;
  
  int old_c ;
  double old_A ;
  
  arma::vec non0 ; // potrebbero essere trasformati in vectors of integers
  arma::vec unique_non0 ;
  arma::vec labels ; 
  double diff ; arma::uvec ids ;
  double denom = 0 ;
    
  for(int j = 0; j < n; j++)
  {
    old_c = cluster(j) ;
    old_A = A(old_c) ;
    A_tmp = A ;
    A.fill(0) ;
    
    cluster(j) = (-99) ; // remove the j-th element
    
    non0 = cluster.elem( find( cluster > 0 ) ); // subvector of observations assigned to some cluster
    n_clus = 0 ;
    
    if(non0.n_elem > 0) 
    {
      unique_non0 = arma::unique( non0 ) ; // unique labels of existing clusters (-j)
      max_clus = unique_non0.max() ; // max label of the cluster
      A_tmp(max_clus + 1) = 0 ; // if the removed cluster was a singleton with the highest label
      n_clus = unique_non0.n_elem ; // number of clusters 
      
      if(max_clus > n_clus) 
      {
        // it means that there is a "hole" in the label sequence
        // if cluster(j) was a singleton
        
        labels = arma::linspace(1, max_clus, max_clus) ;
        diff = missing_value(labels, unique_non0); // find missing label
        
        for(int k = diff + 1 ; k < max_clus +1; k++)
        {
          // relabeling
          ids = find( cluster == k ) ; 
          cluster.elem(ids).fill(k - 1) ;
          
          A_tmp(k-1) = A_tmp(k) ; // sort cluster parameters
        }
        
        // assign the singleton the last label
        cluster(j) = max_clus ;
        A_tmp(max_clus) = old_A ;
        
        // sample (m-1) temporary parameters
        for(int l = 0; l < m-1; l++) { A_tmp(max_clus + 1 + l) = R::rgamma(hyp_A1, hyp_A2) ; }
      }
      
      if(max_clus == n_clus)
      {
        // sample m temporary parameters
        for(int l = 0; l < m; l++) { A_tmp(max_clus + 1 + l) = R::rgamma(hyp_A1, hyp_A2) ; }
      }
    }
    // end if( non0.elem>0 )
    
    double pr0 = log(p) + R::dnorm(y(j), b + gamma * cc(j), std::sqrt(sigma2 + tau2), true) ;
    
    NumericVector prob(n_clus + 1 + m); // cluster assignment probabilities: [ pr(0), pr(cl 1), ..., pr(cl K), pr(cl K+1), ..., pr(cl K+m) ]
    prob[0] = pr0 ;
    
    if(n_clus > 0)
    {
      for(int k = 1; k < n_clus + 1; k++)
      {
        double nj = std::count(non0.begin(), non0.end(), k) ;
        prob[k] = log(nj) - log(non0.n_elem - 1 + alpha) + log(1-p) + R::dnorm(y(j), b + gamma * cc(j) + A_tmp(k), std::sqrt(sigma2 + tau2), true) ;
      }
    }
    
    for(int k = n_clus + 1; k < n_clus + m + 1; k++)
    {
      prob[k] = log(1-p) + log(alpha/m)  - log(non0.n_elem - 1 + alpha) + R::dnorm(y(j), b + gamma * cc(j) + A_tmp(k), std::sqrt(sigma2 + tau2), true) ;
    }
    
    double max_p = max(prob) ;
    prob = exp( prob - max_p ) ;
    
    IntegerVector clusters_id =  Rcpp::seq(0, n_clus + m);
    cluster(j) = Rcpp::sample(clusters_id, 1, false, prob)[0] ;

    for(int k = 0; k < n_clus +1 ; k++) { A(k) = A_tmp(k); }
    if(cluster(j) > n_clus) 
    {
      A(n_clus + 1) = A_tmp(cluster(j)) ;
      cluster(j) = n_clus + 1 ;
    }
  } 
  return Rcpp::List::create(Rcpp::Named("cluster") = cluster,
                            Rcpp::Named("A") = A);
}



// Gibbs sampler: funzione principale
// [[Rcpp::export]]
Rcpp::List calcium_gibbs_nonc(int Nrep, arma::vec y,
                               arma::vec cal,
                               arma::vec cl, 
                               arma::vec A_start,
                               double b_start, double gamma_start, 
                               double lambda_start, double p_start,
                               double c0, double varC0,
                               double tau2,
                               double alpha, int m,
                               double hyp_A1, double hyp_A2,
                               double hyp_b1, double hyp_b2,
                               double hyp_lambda1, double hyp_lambda2,
                               double hyp_gamma1, double hyp_gamma2,
                               double hyp_p1, double hyp_p2,
                               double eps_gamma,
                               double eps_A)
{
  // allocate output matrices
  int n = y.n_elem ;
  arma::mat out_c(n+1, Nrep) ; // 0 1 ... n
  arma::mat out_A = arma::zeros(A_start.n_elem, Nrep) ;
  arma::vec out_b(Nrep) ;
  arma::vec out_gamma = arma::zeros(Nrep) ;
  arma::vec out_lambda = arma::zeros(Nrep) ;
  arma::mat cluster = arma::zeros(n, Nrep)  ;
  arma::vec out_p(Nrep) ;
  
  // initialize the chain
  out_c(0,0) = c0 ;
  out_b(0) = b_start ;
  out_gamma(0) = gamma_start ;
  out_lambda(0) = lambda_start ;
  out_p(0) = p_start ;
  
  // other quantities
  arma::vec filter_mean(n+1) ; // 0 1 ... n
  arma::vec filter_var(n+1) ; // 0 1 ... n
  arma::vec R(n+1) ; arma::vec a(n+1) ; // 0 1 ... n
  double back_mean; double back_var ;
  
  double sigma2 = 1/lambda_start ;
  
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
      R(j) = pow(out_gamma(i), 2) * filter_var(j-1) + tau2 ;
      
      filter_mean(j) = a(j) + R(j) / (R(j) + sigma2) * (y(j-1) - out_b(i) - a(j)) ;
      filter_var(j) = sigma2 * R(j) / (R(j) + sigma2) ;
    }
    out_c(n, i+1) = R::rnorm(filter_mean(n), filter_var(n)) ;
    
    for(int j = n-1; j > -1; j--)
    {
      back_mean = filter_mean(j) + out_gamma(i) * filter_var(j) / R(j+1) * (out_c(j+1, i+1) - a(j+1)) ;
      back_var = filter_var(j) - pow(out_gamma(i) * filter_var(j), 2) / R(j+1) ;
      
      out_c(j, i+1) = R::rnorm(back_mean, back_var) ;
    }
    //out_c.col(i+1) = out_c.col(i) ;

    
    // sampling b
    arma::vec z(n) ; 
    for(int j = 0; j < n; j++) { z(j) = y(j) - out_c(j+1, i+1) ; } 
    
    out_b(i+1) = R::rnorm( (hyp_b2 * hyp_b1 + out_lambda(i) * arma::accu(z)) / (hyp_b2 + n * out_lambda(i)) , 
          std::sqrt( 1/ (hyp_b2 + n * out_lambda(i)) ) ) ;
    //out_b(i+1) = out_b(i) ;
    
    // sampling lambda (precision)
    arma::vec sq(n) ; 
    for(int j = 0; j < n; j++) { sq(j) = pow(z(j) - out_b(i+1), 2) ; }
    
    out_lambda(i+1) = R::rgamma(hyp_lambda1 + n/2, 1/(hyp_lambda2 + 0.5 * arma::accu(sq)) ) ;
    //out_lambda(i+1) = out_lambda(i) ;
    sigma2 = 1/out_lambda(i+1) ;
    
    if( !out_b.is_finite() ) { check = 1 ; }
    if( !out_lambda.is_finite() ) { check = 1 ; }
  
    
    //MH per gamma: random walk
    oldgamma = out_gamma(i) ;
    newgamma = oldgamma + R::runif(-eps_gamma, eps_gamma) ;
    ratio = exp( logpost_gamma(y, out_c.col(i+1), 
                         AA, out_b(i+1), 
                         newgamma, sigma2, tau2,
                         hyp_gamma1, hyp_gamma2) -
                  logpost_gamma(y, out_c.col(i+1), 
                          AA, out_b(i+1), 
                          oldgamma, sigma2, tau2,
                          hyp_gamma1, hyp_gamma2) ) ;
    
    if(R::runif(0, 1) < ratio) oldgamma = newgamma ;
    out_gamma(i+1) = oldgamma ;
    
    
    // Sampling of clusters and cluster parameters
    // Polya-Urn
    polyaurn = polya_urn_nonc(y, cluster.col(i), out_c.col(i+1),
                    out_A.col(i), out_b(i+1), out_gamma(i+1), 
                    out_p(i),
                    sigma2, tau2,
                    alpha, m,
                    hyp_A1, hyp_A2, check) ; 
    arma::vec out_polya_clus = polyaurn["cluster"] ;
    cluster.col(i+1) = out_polya_clus ;
    //cluster.col(i+1) = cluster.col(i) ;
    
    // sampling of parameters A1,...,Ak
    arma::vec out_polya_A = polyaurn["A"] ;
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
                            out_gamma(i+1), sigma2, tau2,
                            hyp_A1, hyp_A2) -
                    logpost_A(sub_y, col_c(idk), 
                            oldA, out_b(i+1), 
                            out_gamma(i+1), sigma2, tau2,
                            hyp_A1, hyp_A2) ) ;
      if(R::runif(0, 1) < ratio) oldA = newA ;
      out_A(l, i+1) = oldA ;
    }
    
    
    // Update p
    double n0 = std::count(line.begin(), line.end(), 0) ;
    out_p(i+1) = R::rbeta(hyp_p1 + n0, hyp_p2 + n - n0) ;
    
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
                            Rcpp::Named("lambda") = out_lambda,
                            Rcpp::Named("cluster") = cluster,
                            Rcpp::Named("p") = out_p);
}

































