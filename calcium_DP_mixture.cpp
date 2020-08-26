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
 */v
double logprior_gamma(double & gamma, double & hyp_gamma1, double & hyp_gamma2) 
{
  double out;
  out = R::dbeta( gamma, hyp_gamma1, hyp_gamma2, true);
  return(out);
}

// log-posterior
double logpost(const arma::vec& y, const arma::vec& cc, const arma::vec& AA, 
               double & b, double & gamma, double & sigma2, double & tau2, 
               double & hyp_gamma1, double & hyp_gamma2)
{
  double out;
  out = loglik(y, cc, AA, b, gamma, sigma2, tau2) + logprior_gamma(gamma, hyp_gamma1, hyp_gamma2);
  return(out);
}

// generazione da normale troncata tra 0 e Inf
double gen_truncnorm(double mean, double sd)
{
  double out = -1;
  while(out < 0) { out = R::rnorm(mean, sd) ; }
  return(out) ;
}

// trova differenza tra due vettori
double missing_value(const arma::vec& x, const arma::vec& y)
{
  double out = 0;
  int n = x.n_elem ;
  for(int i = 0; i < n; i++)
  {
    double lookfor = x(i) ;
    if(std::count(y.begin(), y.end(), lookfor) == 0) { out = lookfor ;}
  }
  return(out) ;
}


// Polya-Urn
arma::vec polya_urn(const arma::vec& y, arma::vec cluster, const arma::vec& cc,
                    arma::vec A, double b, double gamma, 
                    double p,
                    double sigma2, double tau2,
                    double & alpha, 
                    double & hyp_A1, double & hyp_A2,  
                    int & check)
{
  int n = y.n_elem ;
  int n_clus = 0 ; 
  double max_clus = 0 ;
    
  for(int j = 0; j < n; j++)
  {
    cluster(j) = (-99) ; // remove the j-th element
    n_clus = 0 ; 
    
    arma::vec non0 = cluster.elem( find( cluster > 0 ) ); // subvector of observations assigned to some cluster

    if(non0.n_elem > 0) 
    {
      arma::vec unique_non0 = arma::unique( non0 ) ; // unique labels of existing clusters
      max_clus = non0.max() ; // max label of the cluster
      A(max_clus + 1) = 0 ; // if the removed cluster was a singleton with the highest label
      n_clus = unique_non0.n_elem ; // number of clusters 
      
      if(max_clus > n_clus) // it means that there is a "hole" in the label sequence
      {
        arma::vec labels = arma::linspace(1, max_clus, max_clus) ;
        double diff = missing_value(labels, unique_non0); // find missing label
        
        for(int k = diff + 1 ; k < max_clus +1; k++)
        {
          // rename labels 
          arma::uvec ids = find( cluster == k ) ; 
          arma::vec elem_k = cluster.elem(ids) ;
          cluster.elem(ids).fill(k - 1) ;
          A(k-1) = A(k) ;
        }
        
        A(max_clus) = 0 ;
      }
    }
    
    NumericVector prob(n_clus + 2); // cluster assignment probabilities: [ pr(0), pr(cl 1), pr(cl 2), ..., pr(cl K), pr(cl K+1) ]

    double z =  y(j) - b - gamma * cc(j) ;
    
    // prob di non-spike
    double pr0 = log(p) + R::dnorm(z, 0, std::sqrt(sigma2 + tau2), true) ;
    prob[0] = pr0 ;
  
    // prob di spike e assegno a un cluster esistente
    if(n_clus > 0)
    {
      for(int k = 1; k < n_clus + 1; k++)
      {
        double nj = std::count(non0.begin(), non0.end(), k) ;
        prob[k] = log(nj) - log(non0.n_elem - 1 + alpha)  + log(1-p) + 
                    R::dnorm(y(j), b + gamma * cc(j) + A(k), std::sqrt(sigma2 + tau2), true) ;
      }
    }
    
    // prob di spike e assegno a un nuovo cluster
    double pr_new = log(1-p) + log(alpha) - log(non0.n_elem - 1 + alpha) + 
                      R::dnorm( z, hyp_A1, std::sqrt(sigma2 + tau2 + hyp_A2), true ) ;
    prob[n_clus + 1] = pr_new ; 
    
    // normalizzo e esponenziale
    double max_p = max(prob) ;
    prob = exp( prob - max_p ) ;
    
    IntegerVector clusters_id =  seq(0, n_clus + 1);
    
    if( exp(max_p) == 0 ) 
    { 
      prob.fill(0) ;
      prob[0] = 1 ;
      if(check == 0) { Rcout << "Probabilità nulle \n" ; }
      check = 1 ;
    }
    
    cluster(j) = Rcpp::sample(clusters_id, 1, false, prob)[0] ;
    
    if(cluster(j) == n_clus + 1)
    {
      double precA = 1/(sigma2 + tau2) + 1/hyp_A2 ;
      double meanA = 1/precA * ( (y(j) - b - gamma * cc(j))/(sigma2 + tau2) + hyp_A1/hyp_A1 ) ;
      
      A(n_clus + 1) = gen_truncnorm( meanA , std::sqrt(1/precA) ) ;
    }
  } 
  return(cluster) ;
}



// Gibbs sampler: funzione principale
// [[Rcpp::export]]
Rcpp::List calcium_gibbs(int Nrep, arma::vec y,
                          arma::vec cal,
                          arma::vec cl, 
                          arma::vec A_start,
                          double b_start, double gamma_start, 
                          double lambda_start, double p_start,
                          double c0, double varC0,
                          double tau2,
                          double alpha, 
                          double hyp_A1, double hyp_A2,
                          double hyp_b1, double hyp_b2,
                          double hyp_lambda1, double hyp_lambda2,
                          double hyp_gamma1, double hyp_gamma2,
                          double hyp_p1, double hyp_p2,
                          double eps_gamma)
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
    
  //  out_c.col(i+1) = out_c.col(i) ;
    
    
    
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
   // out_lambda(i+1) = out_lambda(i) ;
    sigma2 = 1/out_lambda(i+1) ;
    
    if( !out_b.is_finite() ) { check = 1 ; }
    if( !out_lambda.is_finite() ) { check = 1 ; }
    

    
    
    //MH per gamma: random walk
    oldgamma = out_gamma(i) ;
    newgamma = oldgamma + R::runif(-eps_gamma, eps_gamma) ;
    ratio = exp( logpost(y, out_c.col(i+1), 
                         AA, out_b(i+1), 
                         newgamma, sigma2, tau2,
                         hyp_gamma1, hyp_gamma2) -
                  logpost(y, out_c.col(i+1), 
                          AA, out_b(i+1), 
                          oldgamma, sigma2, tau2,
                          hyp_gamma1, hyp_gamma2) ) ;
    
    if(R::runif(0, 1) < ratio) oldgamma = newgamma ;
    out_gamma(i+1) = oldgamma ;
    /*
    
    // Sampling of clusters and cluster parameters
    // Polya-Urn
    cluster.col(i+1) = polya_urn(y, cluster.col(i), out_c.col(i+1),
                out_A.col(i), out_b(i+1), out_gamma(i+1), 
                out_p(i),
                sigma2, tau2,
                alpha, 
                hyp_A1, hyp_A2, check) ; 
    //cluster.col(i+1) = cluster.col(i) ;
*/
    // sampling of parameters A1,...,Ak
    arma::vec line(n); line = cluster.col(i+1) ;
    arma::vec non0 = line.elem( find( line > 0 ) );
    arma::vec tmp2 = arma::unique( non0 ) ; // unique labels 
    
    n_clus = tmp2.n_elem ; // number of clusters (n.unique \{0}) ...
    double nj = 0 ;
    double ssum = 0;
    double meanA = 0; double precA = 0;
    
    /*
    arma::vec lincomb(n) ;
    for(int l = 0; l < n; l++) { lincomb(l) = y(l) - out_b(i+1) - out_gamma(i+1) * out_c(l, i+1) ; }
    
    for(int k = 1; k < n_clus + 1; k++) 
    { 
      nj = std::count(non0.begin(), non0.end(), k) ; 
      arma::uvec idk = find( line == k ) ;
      ssum = arma::accu( lincomb(idk) ) ;
      
      precA = nj/(sigma2 + tau2) + 1/hyp_A2 ;
      meanA = 1/precA * ( ssum/(sigma2 + tau2) + hyp_A1/hyp_A2 ) ;
      
      out_A(k, i+1) = gen_truncnorm( meanA, std::sqrt(1/precA) ) ;
    }
    
    */
    // Update p
    double n0 = std::count(line.begin(), line.end(), 0) ;
    out_p(i+1) = R::rbeta(hyp_p1 + n0, hyp_p2 + n - n0) ;
    
    if(out_p(i+1) < 0.5) { 
      Rcout << "Stop at iter. " << i << "\n" ;
      i = Nrep - 2 ; 
    }
    
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

































