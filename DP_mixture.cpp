#include <RcppArmadillo.h>
#include <math.h>
#include <time.h>
#include <mvnorm.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]


using namespace Rcpp;

// log-likelihood
double loglik(arma::vec y, arma::vec cc, arma::vec A, 
              double b, double gamma, double lambda)
{
  int n = y.n_elem;
  arma::vec llik(n);
  
  arma::vec mean;
  for(int j = 0; j < n; j++) { mean(j) = b + gamma * cc(j) + A(j); }
  
  for(int k = 0; k < n; k++) { 
    llik(k) = 0.5 * (std::log(lambda) - std::log(2 * M_PI) - lambda * pow(y[k] - mean[k], 2) );
  }
  return arma::accu(llik);
}

// prior on gamma
double logprior_gamma(double gamma, double hyp_gamma1, double hyp_gamma2) 
{
  double out;
  out = R::dbeta( gamma, hyp_gamma1, hyp_gamma2, true);
  return(out);
}

// log-posterior
double logpost(arma::vec y, arma::vec cc, arma::vec A, 
               double b, double gamma, double lambda,
               double hyp_gamma1, double hyp_gamma2)
{
  double out;
  out = loglik(y, cc, A, b, gamma, lambda) + logprior_gamma(gamma, hyp_gamma1, hyp_gamma2);
  return(out);
}


// marginale per Polya Urn
double marginal(double y, double cc, double b, double gamma, 
                double sigma2, double psi2)
{
  double out;
  double mean1 = b + gamma * cc;
  double sd1 = std::sqrt(psi2 + sigma2);
  double mean2 = psi2 / (sigma2 + psi2) * (y - b - gamma * cc);
  double sd2 = std::sqrt(sigma2 * psi2 / (sigma2 + psi2));
  
  out = 2 * R::dnorm(y, mean1, sd1, false) * R::pnorm(0, mean2, sd2, false, false);
  return(out);
}


// generazione da normale troncata tra 0 e Inf
// [[Rcpp::export]]
double gen_truncnorm(double mean, double sd)
{
  double out = -1;
  while(out < 0) { out = R::rnorm(mean, sd) ; }
  return(out) ;
}


Rcpp::NumericVector arma_setdiff(arma::vec x, arma::vec y)
{
  x = arma::unique(x);
  y = arma::unique(y);
  for (size_t j = 0; j < y.n_elem; j++) 
  {
    arma::uvec q1 = arma::find(x == y[j]);
    if (!q1.empty()) { x.shed_row(q1(0)); }
  }
  Rcpp::NumericVector x2 = Rcpp::wrap(x);
  x2.attr("dim") = R_NilValue;
  return x2;
}




Rcpp::List calcium_gibbs(int Nrep, arma::vec y,  
                         double gamma_start, double lambda_start,
                         double c0, double varC0,
                         double tau2,
                         double p,
                         double alpha, double psi2,
                         double hyp_A1, double hyp_A2,
                         double hyp_b1, double hyp_b2,
                         double hyp_gamma1, double hyp_gamma2,
                         double hyp_lambda1, double hyp_lambda2,
                         double eps_gamma)
{
  double b_start = median(y);
  
  // allocate output matrices
  int n = y.n_elem ;
  arma::mat out_c(Nrep, n+1) ;
  arma::mat out_A = arma::zeros(Nrep, n) ;
  arma::vec out_b(Nrep) ;
  arma::vec out_gamma(Nrep) ;
  arma::vec out_lambda(Nrep) ;
  arma::mat cluster = arma::zeros(Nrep, n) ;
  
  // initialize the chain
  out_c(0,0) = c0 ;
  out_b(0) = b_start ;
  out_gamma(0) = gamma_start ;
  out_lambda(0) = lambda_start ;
  
  // other quantities
  arma::vec filter_mean(n) ;
  arma::vec filter_var(n) ;
  arma::vec R; arma::vec a ;
  double back_mean; double back_var ;
  double sigma2 ;
  
  arma::vec z; arma::vec sq ; 
  
  double oldgamma; double newgamma ;
  double ratio ;
  
  arma::vec AA = arma::zeros(n) ;
  arma::vec clus_tmp = arma::zeros(n); 
  arma::vec A_tmp = arma::zeros(n) ;

    
  for(int i = 0; i < Nrep; i++)
  {
    sigma2 = 1/out_lambda(i) ;
    for(int j = 0; j < n; j++) { AA(j) =  out_A(i, cluster(i,j)) ; }
    
    // sampling c
    R(0) = tau2 + pow(out_gamma(i), 2) * varC0 ;
    filter_mean(0) = (sigma2 * AA(0) + R(0) * (y(0) - out_b(i))) / (sigma2 + R(0)) ;
    filter_var(0) = sigma2 * R(0) /  (sigma2 + R(0)) ;
    
    for(int j = 1; j < n; j++)
    {
      R(j) = pow(out_gamma(i), 2) * filter_var(j-1) + tau2 ;
      a(j) = out_gamma(i) * filter_mean(j-1) + AA(j) ;
      
      filter_mean(j) = a(j) + R(j) / (R(j) + sigma2) * (y(j) - out_b(i) - a(j)) ;
      filter_var(j) = R(j) - pow(R(j), 2) / (R(j) + sigma2) ;
    }
    out_c(i+1, n) = R::rnorm(filter_mean(n-1), filter_var(n-1)) ;
      
    for(int j = n-1; j > 0; j--)
    {
      back_mean = filter_mean(j) + out_gamma(i) * filter_var(j) / R(j+1) * (out_c(i+1, j+2) - a(j+1)) ;
      back_var = filter_var(j) - pow(out_gamma(i) * filter_var(j), 2) / R(j+1) ;
        
      out_c(i+1, j+1) = R::rnorm(back_mean, back_var) ;
    }
    back_mean = out_gamma(i) * varC0 / (varC0 + tau2) * out_c(i+1, 1) ;
    back_var = varC0 - pow(out_gamma(i) * varC0, 2) / R(1) ;
    out_c(i+1, 0) = R::rnorm(back_mean, back_var) ;
      
    // sampling lambda (precision)
    for(int j = 0; j < n; j++) { z(j) = y(j) - out_gamma(i) * out_c(i+1, j) - AA(j) ; }
    for(int j = 0; j < n; j++) { sq(j) = pow(z(j) - mean(z), 2) ; } 
    out_lambda(i+1) = R::rgamma(hyp_lambda1 + n/2, 
                                hyp_lambda2 + 0.5 * accu(sq) + 0.5 * n * hyp_b2 / (n + hyp_b2) * pow(mean(z) - hyp_b1, 2) ) ;
    // sampling b
    out_b(i+1) = R::rnorm((n * mean(z) + hyp_b2 * hyp_b1) / (n + hyp_b2), 1/ std::sqrt((n + hyp_b2) * out_lambda(i+1)) ) ;
      
    //MH per gamma: random walk
    oldgamma = out_gamma(i) ;
    newgamma = oldgamma + R::runif(-eps_gamma, eps_gamma) ;
    ratio = exp( logpost(y, out_c.row(i+1), 
                         AA, out_b(i+1), 
                         newgamma, out_lambda(i+1),
                         hyp_gamma1, hyp_gamma2) -
                 logpost(y, out_c.row(i+1), 
                         AA, out_b(i+1), 
                         oldgamma, out_lambda(i+1),
                         hyp_gamma1, hyp_gamma2) ) ;
    
    if(R::runif(0, 1) < ratio) oldgamma = newgamma ;
    out_gamma(i+1) = oldgamma ;
    
    
    // Sampling of cluster and cluster parameters
    cluster.row(i+1) = cluster.row(i) ; 
    
    // Polya-Urn
    for(int j = 0; j < n; j++)
    {
      clus_tmp = cluster.row(i+1) ;
      A_tmp = out_A.row(i) ;
      clus_tmp(j) = -99 ; // remove the j-th element
      
      double max_clus = clus_tmp.max() ; // max label of the cluster
      arma::vec tmp = arma::unique(clus_tmp) ; // unique labels
      double n_clus = tmp.n_elem -2; // number of clusters (n.unique \{0,-99})
      if(max_clus > n_clus) // it means that there is a "hole" in the label sequence
      {
        arma::vec v = arma::linspace(1, max_clus, max_clus) ;
        arma::vec diff = arma_setdiff(tmp, v) ;
        for(int k = min(diff) +1 ; k <= max_clus + 1; k++)
        {
          arma::uvec ids = find(clus_tmp == k) ;
          clus_tmp.elem(ids).fill(k - 1) ;
          A_tmp(k-1) = A_tmp(k) ;
        }
      }
      
      double pr0 = p * R::dnorm(y(j), out_b(i+1) + out_gamma(i+1) * out_c(i+1, j), std::sqrt(sigma2 + tau2), false) ;
        
      double pr_new = (1-p) * alpha * marginal(y(j), out_c(i+1, j), out_b(i+1), out_gamma(i+1), sigma2 + tau2, psi2) ;
      
      arma::vec nj(n_clus);
      arma::vec pr_clus(n_clus) ;
      for(int k = 1; k < n_clus + 1; k++)
      {
        nj(k-1) = std::count(clus_tmp.begin(), clus_tmp.end(), k) ;
        pr_clus(k-1) =  nj(k-1) * (1-p) * R::dnorm(y(j), out_b(i+1) + out_gamma(i+1) * out_c(i+1, j) + A_tmp(k), std::sqrt(sigma2 + tau2), false) ;
      }
      
      arma::vec prob(n_clus + 2);
      prob(0) = pr0 ;
      for(int k = 0; k < n_clus; k++)
      {
        prob(k+1) = pr_clus(k) ; 
      }
      prob(n_clus + 1) = pr_new ; 
      
      arma::vec clusters_id = arma::linspace(0, n_clus + 1, n_clus + 2);
      clus_tmp(j) = RcppArmadillo::sample( clusters_id, 1, false, prob )[0] ;
      
      if(clus_tmp(j) == n_clus + 1)
      {
        A_tmp(n_clus + 1) = gen_truncnorm(psi2 / (psi2 + sigma2 + tau2) * (y(j) - out_b(i+1) - out_gamma(i+1) * out_c(i+1, j)), 
                                            std::sqrt((sigma2 + tau2) * psi2 / (psi2 + sigma2 + tau2) )) ;
      }
      cluster.row(i+1) = clus_tmp ;
      out_A.row(i+1) = A_tmp ;
    }
    
    // sampling of parameters A1,...,Ak
    arma::vec tmp = arma::unique(cluster.row(i+1)) ; // unique labels
    double n_clus = tmp.n_elem -1 ; // number of clusters (n.unique \{0})
    arma::vec nj(n_clus) ;
    arma::vec ssum(n_clus) ;
    for(int k = 1; k < n_clus + 1; k++) 
    { 
      nj(k-1) = std::count(cluster.row(i+1).begin(), cluster.row(i+1).end(), k) ; 
      arma::vec lincomb ;
      lincomb = Rcpp::sapply
    }

    ///
    // manca sum 
    // sampling A
    
  }
  return Rcpp::List::create(Rcpp::Named("beta") = out_c,
                            Rcpp::Named("acceptance_rate") = out_A) ;
}














