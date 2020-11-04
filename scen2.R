library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)

## mixture ##
sim_data <- function(n, sigma2, tau2, time_spike, b, gamma, prob, par)
{
  c = rep(0,n)
  s = rep(0,n)
  A = rep(0,n)
  s[time_spike] = 1
  
  k <- sample(1:length(prob), sum(s), prob, replace = TRUE)
  A[time_spike] <- par[k]
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i] + rnorm(1, 0, sqrt(tau2))
  }
  return(list("y" = b + c + rnorm(n, 0, sqrt(sigma2)), "c" = c, "s" = s, "A" = A, "k" = k))
}

times_spike <- function(n, ns, p1, p2)
{
  times = rbinom(n, 1, p1)
  for(i in 1:n)
  {
    if(times[i] == 1)
    {
      times[(i+1):(i+ns)] = rbinom(ns, 1, p2)
    }
  }
  times 
}


#--------# #--------# #--------# #--------# #--------#
#--------# #--------# Scenario 1 #--------# #--------#
# distribuzioni ben separate, atomi ben distinti

set.seed(1234)
sigma2 = 0.004
tau2 = 0.0001

n1 = 5000
n2 = 5000
n3 = 5000
n4 = 5000
n5 = 5000
n6 = 5000

gamma = 0.7
b = 0


spp = c( times_spike(n1, 9, 0.016, 0.11)[1:n1], # clD1
         times_spike(n2, 7, 0.01, 0.14)[1:n2], ## clD2
         times_spike(n3, 6, 0.018, 0.13)[1:n3], ### cld3
         times_spike(n4, 9, 0.016, 0.11)[1:n4], # cld1
         times_spike(n5, 9, 0.016, 0.11)[1:n5], ## cld1
         times_spike(n6, 7, 0.01, 0.14)[1:n6]) ### clD2
spp = which(spp == 1)
length(spp)

sum(spp<n1)
sum((spp>n1)&(spp<n1+n2))
sum((spp>n1+n2)&(spp<n1+n2+n3))
sum((spp>n1+n2+n3)&(spp<n1+n2+n3+n4))
sum((spp>n1+n2+n3+n4)&(spp<n1+n2+n3+n4+n5))
sum((spp>n1+n2+n3+n4+n5)&(spp<n1+n2+n3+n4+n5+n6))

prob1 = rep(1, 3)/3
par1 = c(0.4, 0.6, 0.9)

prob2 = rep(0.25, 4)
par2 = c(0.4, 0.75, 1.15, 1.5)

prob3 = c(0.2, 0.4, 0.4)
par3 = c(0.75, 1.5, 1.75)

### parametri: 
unip <- c(0.4, 0.6, 0.75, 0.9, 1.15, 1.5, 1.75)
### n clus: 8 = (7 spike + 0)

### n distr: 3

group1 <- sim_data(n = n1, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp<=n1],
                   gamma = gamma, b = b,
                   prob = prob1, par = par1)

group2 <- sim_data(n = n2, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>n1)&(spp<=(n1 + n2))] - n1,
                   gamma = gamma, b = b,
                   prob = prob2, par = par2)

group3 <- sim_data(n = n3, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2))&(spp<=(n1 + n2 + n3))] - (n1 + n2),
                   gamma = gamma, b = b,
                   prob = prob3, par = par3)

group4 <- sim_data(n = n4, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2 + n3))&(spp<=(n1 + n2 + n3 + n4))] - (n1 + n2 + n3),
                   gamma = gamma, b = b,
                   prob = prob1, par = par1)

group5 <- sim_data(n = n5, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2 + n3 + n4))&(spp<=(n1 + n2 + n3 + n4 + n5))] - (n1 + n2 + n3 + n4),
                   gamma = gamma, b = b,
                   prob = prob1, par = par1)

group6 <- sim_data(n = n5, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp>(n1 + n2 + n3 + n4 + n5)] - (n1 + n2 + n3 + n4 + n5),
                   gamma = gamma, b = b,
                   prob = prob2, par = par2)

y = c(group1$y, group2$y, group3$y, group4$y, group5$y, group6$y)
# save(y, file="y_scen2.Rdata")
# save(spp, file="spp_scen2.Rdata")
g = c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4), rep(5,n5), rep(6,n6))
plot(y, type = "l")


A_start = rep(0,100)
# A_start[2:7] = c(0.5, 0.8, 1.1, 1.4, 1.75, 1.9)
# AA = c(group1$A, group2$A, group3$A, group4$A)
# 
# cluster = rep(0, length(y))
# cluster[AA == A_start[2]] = 1
# cluster[AA == A_start[3]] = 2
# cluster[AA == A_start[4]] = 3
# cluster[AA == A_start[5]] = 4
# cluster[AA == A_start[6]] = 5
# cluster[AA == A_start[7]] = 6


clus = kmeans(y[y>0.5], centers = 5)

A_start = rep(0,100)
A_start[2:6] = c(clus$centers)
cluster = numeric(length(y))
cluster[y>0.5] = clus$cluster

sum(length(spp))/length(y)

n = length(y)
J = length(unique(g))
nrep = 2000
burnin = 1:50


sourceCpp('./SourceCPP/calcium_nested_gMFM2.cpp')
set.seed(1234)
run_gMFM = calcium_gibbs(Nrep = nrep, 
                    y = y, g = g,
                    cal = c(0,y),
                    clO = cluster, clD = c(1,2,3,4,5,6), 
                    A_start = A_start,
                    b_start = 0,
                    gamma_start = 0.7, 
                    sigma2_start = 0.004, 
                    tau2_start = 0.0001, 
                    p_start = 0.01, 
                    alpha_start = 0.4, beta_start = 0.6,
                    maxK_start = 5,
                    maxL_start = 7,
                    c0 = 0, varC0 = 0.1, 
                    hyp_A1 = 5, hyp_A2 = 7, 
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                    hyp_tau21 = 1000, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    hyp_alpha1 = 3, hyp_alpha2 = 3,
                    hyp_beta1 = 3, hyp_beta2 = 6,
                    hyp_maxK1 = 2, hyp_maxK2 = 4, hyp_maxK3 = 3,
                    hyp_maxL1 = 2, hyp_maxL2 = 4, hyp_maxL3 = 3,
                    eps_alpha = 0.5, eps_beta = 0.7,
                    eps_gamma = 0.005,
                    eps_A = 0.02,
                    eps_maxK = 4, eps_maxL = 5)

barplot(table(apply(run_gMFM$clusterD[,-burnin], 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
barplot(table(apply(run_gMFM$clusterO[,-burnin], 2, function(x) length(unique(x)) )))
# save(run_gMFM, file="scen2_res_gMFM.Rdata")




sourceCpp('./SourceCPP/calcium_CAM.cpp')
run_DP = calcium_gibbs(Nrep = nrep, 
                    y = y, g = g,
                    cal = c(0,y),
                    clO = cluster, clD = c(1,2,3,4,5,6), 
                    A_start = A_start,
                    b_start = 0,
                    gamma_start = 0.8, 
                    sigma2_start = 0.004, 
                    tau2_start = 0.0001, 
                    p_start = 0.01, 
                    alpha = 0.5, beta = 1,
                    max_xiK = 700, max_xiL = 500,
                    kappa_D = 0.5, kappa_O = 0.5, 
                    c0 = 0, varC0 = 0.1, 
                    hyp_A1 = 5, hyp_A2 = 7, 
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 100, hyp_sigma22 = 1, 
                    hyp_tau21 = 100, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 99,
                    eps_gamma = 0.008,
                    eps_A = 0.02)
save(run_DP, file="scen2_res_DP.Rdata")
barplot(table(apply(run_DP$clusterD, 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
barplot(table(apply(run_DP$clusterO[,-burnin], 2, function(x) length(unique(x)) )))



