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


#--------# #--------# #--------# #--------# #--------#
#--------# #--------# Scenario 1 #--------# #--------#
# distribuzioni ben separate, atomi ben distinti


sigma2 = 0.004
tau2 = 0.0001
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
pm4 = 0.09

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

prob1 = rep(1, 3)/3
par1 = c(0.5, 0.8, 1.1)

prob2 = rep(0.25, 4)
par2 = c(0.8, 1.1, 1.4, 1.75)

prob3 = c(0.2, 0.4, 0.4)
par3 = c(0.5, 1.4, 1.9)

### parametri: 
unip <- c(0.5, 0.8, 1.1, 1.4, 1.75, 1.9)
### n clus: 9 = (8 spike + 0)

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

group4 <- sim_data(n = n4, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp>(n1 + n2 + n3 )] - (n1 + n2 + n3 ),
                   gamma = gamma, b = b,
                   prob = prob3, par = par3)

y = c(group1$y, group2$y, group3$y, group4$y)
# save(y, file="y_scen1.Rdata")
# save(spp, file="spp_scen1.Rdata")
g = c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4))
plot(y, type = "l")
