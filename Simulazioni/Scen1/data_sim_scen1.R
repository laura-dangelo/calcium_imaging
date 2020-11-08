###################################################
#----------# #----------# #---------- # #----------#
#----------#         Scenario 1         #----------#
# J = 6
# cluster su distr = 4
# parametri ben distinti
# distribuzioni ben distinte

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


#----------# parametri
sigma2 = 0.01
tau2 = 0.0003
n1 = n2 = n3 = n4 = n5 = n6 = 5000
gamma = 0.8
b = 0

# prob di spike "distanti"
pp1 = 0.018
pp2 = 0.010
pp3 = 0.014
pp4 = 0.012

# per quanti istanti successivi ho spikes
m1 = 5
m2 = 8
m3 = 7
m4 = 10

# probabilità di spike per m istanti successivi a uno spike
pm1 = 0.17
pm2 = 0.13
pm3 = 0.13
pm4 = 0.14

set.seed(1234)
spp = c( times_spike(n1, m1, pp1, pm1)[1:n1], # clD1
         times_spike(n2, m2, pp2, pm2)[1:n2], ## clD2
         times_spike(n3, m3, pp3, pm3)[1:n3], ### cld3
         times_spike(n4, m4, pp4, pm4)[1:n4], #### clD4
         times_spike(n5, m1, pp1, pm1)[1:n5], # clD1
         times_spike(n6, m2, pp2, pm2)[1:n6]) ## clD2

spp = which(spp == 1)
length(spp)

sum(spp<n1)
sum((spp>n1)&(spp<n1+n2))
sum((spp>n1+n2)&(spp<n1+n2+n3))
sum((spp>n1+n2+n3)&(spp<n1+n2+n3+n4))
sum((spp>n1+n2+n3+n4)&(spp<n1+n2+n3+n4+n5))
sum((spp>n1+n2+n3+n4+n5)&(spp<n1+n2+n3+n4+n5+n6))

prob1 = c(0.25, 0.25, 0.2, 0.15, 0.15)
par1 = c(0.5, 1.2, 1.7, 2, 2.3)

prob2 = rep(0.25, 4)
par2 = c(0.9, 1.2, 1.7, 2.6)

prob3 = c(0.4, 0.2, 0.4)
par3 = c(0.5, 1.2, 2)

prob4 = rep(1, 3)/3
par4 = c(0.9, 1.7, 2)


### parametri: 
unip <- sort(unique(c(par1,par2,par3,par4)))
unip
length(unip)

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
                   prob = prob4, par = par4)

group5 <- sim_data(n = n5, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2 + n3 + n4))&(spp<=(n1 + n2 + n3 + n4 + n5))] - (n1 + n2 + n3 + n4),
                   gamma = gamma, b = b,
                   prob = prob1, par = par1)

group6 <- sim_data(n = n6, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp>(n1 + n2 + n3 + n4 + n5)] - (n1 + n2 + n3 + n4 + n5),
                   gamma = gamma, b = b,
                   prob = prob2, par = par2)



y = c(group1$y, group2$y, group3$y, group4$y, group5$y, group6$y)
g = c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4), rep(5,n5), rep(6,n6))
plot(y, type = "l")

# save(y, file="y_scen1.Rdata")
# save(spp, file="spp_scen1.Rdata")




clus = kmeans(y[(diff(y,1))>0.5], centers = 5)
A_start = rep(0,50)
A_start[2:6] = c(clus$centers)
cluster = numeric(length(y))
cluster[ (diff(y,1))>0.5] = clus$cluster

n = length(y)
J = length(unique(g))
nrep = 3000
burnin = 1:1000
set.seed(1234)

sourceCpp('./SourceCPP/calcium_nested_gMFM2.cpp')
run_gMFM = calcium_gibbs(Nrep = nrep, 
                         y = y, g = g,
                         cal = c(0,y),
                         clO = cluster, clD = c(1,2,3,4,5,6), 
                         A_start = A_start,
                         b_start = 0,
                         gamma_start = 0.6, 
                         sigma2_start = 0.004, 
                         tau2_start = 0.0001, 
                         p_start = 0.001, 
                         alpha_start = 1, beta_start = 1,
                         maxK_start = 7,
                         maxL_start = 20,
                         c0 = 0, varC0 = 0.1, 
                         hyp_A1 = 4, hyp_A2 = 5, 
                         hyp_b1 = 0, hyp_b2 = 1, 
                         hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                         hyp_tau21 = 1000, hyp_tau22 = 1, 
                         hyp_gamma1 = 1, hyp_gamma2 = 1,
                         hyp_p1 = 1, hyp_p2 = 999,
                         hyp_alpha1 = 3, hyp_alpha2 = 3,
                         hyp_beta1 = 3, hyp_beta2 = 3,
                         hyp_maxK1 = 4, hyp_maxK2 = 4, hyp_maxK3 = 3,
                         hyp_maxL1 = 4, hyp_maxL2 = 4, hyp_maxL3 = 3,
                         eps_alpha = 0.7, eps_beta = 0.7,
                         eps_gamma = 0.005,
                         eps_A = 0.002,
                         eps_maxK = 7, eps_maxL = 10)

# save(run_gMFM, file="scen1_res_gMFM.Rdata")
barplot(table(apply(run_gMFM$clusterD, 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
barplot(table(apply(run_gMFM$clusterO[,-burnin], 2, function(x) length(unique(x)) )))



sourceCpp('./SourceCPP/calcium_CAM.cpp')
run_DP = calcium_gibbs(Nrep = nrep, 
                       y = y, g = g,
                       cal = c(0,y),
                       clO = cluster, clD = c(1,2,3,4,5,6), 
                       A_start = A_start,
                       b_start = 0,
                       gamma_start = 0.6, 
                       sigma2_start = 0.004, 
                       tau2_start = 0.0001, 
                       p_start = 0.01, 
                       alpha = 1, beta = 1,
                       max_xiK = 700, max_xiL = 500,
                       kappa_D = 0.5, kappa_O = 0.5, 
                       c0 = 0, varC0 = 0.1, 
                       hyp_A1 = 4, hyp_A2 = 5, 
                       hyp_b1 = 0, hyp_b2 = 1, 
                       hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                       hyp_tau21 = 1000, hyp_tau22 = 1, 
                       hyp_gamma1 = 1, hyp_gamma2 = 1,
                       hyp_p1 = 1, hyp_p2 = 999,
                       eps_gamma = 0.008,
                       eps_A = 0.02)
# save(run_DP, file="scen1_res_DP.Rdata")
barplot(table(apply(run_DP$clusterD, 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
barplot(table(apply(run_DP$clusterO[,-burnin], 2, function(x) length(unique(x)) )))


