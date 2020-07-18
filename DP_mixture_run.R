library(Rcpp)
library(RcppDist)
sourceCpp('DP_mixture.cpp')

calcium_gibbs(Nrep = 10, y = c(1,2,3), 
              gamma_start = 0.3, lambda_start = 10, 
              c0 = 0, varC0 = 1, tau2 = 0.0001, p = 0.995, 
              alpha = 1, psi2 = 1, 
              hyp_A1 = 1, hyp_A2 = 1, hyp_b1 = 1, hyp_b2 = 1, 
              hyp_gamma1 = 1, hyp_gamma2 = 2, hyp_lambda1 = 1, hyp_lambda2 = 1, eps_gamma = 0.3)


out2 <- sapply(1:10000, function(x) gen_truncnorm(5,2))

out1 <- rtruncnorm(10000, a = 0, b = Inf, 5, 2)

hist(out1)
hist(out2, add=T, border=2)
