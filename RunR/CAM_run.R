library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)
# install.packages("RcppProgress")
sourceCpp('./SourceCPP/calcium_CAM.cpp')


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



spp = c(380,385,
        500,553,
        985,
        1000, 1001, 1003, 1004, 1005, 1007, 
        1012,
        1021, 1022, 1024,
        1100, 1820, 1825,
        2800, 
        2990, 2991, 2993, 2995, 2999, 3001, 3002,
        3101,3104,
        3500,3502,3503,
        4000, 4004,
        4700, 4701, 4703,
        5500,
        6000,6003,
        6233, 
        6250,
        7100,7120,
        7700, 7703,
        8700,
        8800, 8802, 8803, 8804,
        9200, 9300, 9340)

set.seed(1234)
sigma2 = 0.004
tau2 = 0.00002
n1 = 3000
n2 = 1500
n3 = 3000
n4 = 2500
gamma = 0.6
b = 0
group1 <- sim_data(n = n1, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp<=n1],
                 gamma = gamma, b = b,
                 prob = c(0.13, 0.5, 0.43), par = c(0.41, 0.5, 0.75))

group2 <- sim_data(n = n2, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>n1)&(spp<=(n1 + n2))] - n1,
                   gamma = gamma, b = b,
                   prob = c(0.3, 0.43), par = c(0.6, 0.8))

group3 <- sim_data(n = n3, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2))&(spp<=(n1 + n2 + n3))] - (n1 + n2),
                   gamma = gamma, b = b,
                   prob = c(0.13, 0.5, 0.43), par = c(0.41, 0.5, 0.75))

group4 <- sim_data(n = n4, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp>(n1 + n2 + n3)] - (n1 + n2 + n3),
                   gamma = gamma, b = b,
                   prob = c(0.42, 0.3, 0.2), par = c(0.41, 0.75, 0.8))


y = c(group1$y, group2$y, group3$y, group4$y)
g = c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4))
plot(y, type = "l")

A_start = rep(0,50)
A_start[2:6] = c(0.41,0.5,0.6,0.75,0.8)
AA = c(group1$A, group2$A, group3$A, group4$A)

cluster = rep(0, length(y))
cluster[AA == A_start[2]] = 1
cluster[AA == A_start[3]] = 2
cluster[AA == A_start[4]] = 3
cluster[AA == A_start[5]] = 4
cluster[AA == A_start[6]] = 5

1-sum(cluster == 0)/length(y)

n = length(y)
nrep = 500
set.seed(1234)

run = calcium_gibbs(Nrep = nrep, 
                    y = y, g = g,
                    cal = c(0,y),
                    clO = cluster, clD = c(1,2,1,3), 
                    A_start = A_start,
                    b_start = 0,
                    gamma_start = 0.6, 
                    sigma2_start = 0.002, 
                    tau2_start = 0.001, 
                    p_start = 0.001, 
                    alpha = 1, beta = 1,
                    max_xiK = 100, max_xiL = 100,
                    kappa_D = 0.5, kappa_O = 0.5, 
                    c0 = 0, varC0 = 0.1, 
                    hyp_A1 = 10, hyp_A2 = 10, 
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                    hyp_tau21 = 1000, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    eps_gamma = 0.01,
                    eps_A = 0.002)



# burnin = 1:500
# run$calcium = run$calcium[,-burnin]
# run$cluster = run$cluster[,-burnin]
# run$b = run$b[-burnin]
# run$gamma = run$gamma[-burnin]
# run$sigma2 = run$sigma2[-burnin]
# run$tau2 = run$tau2[-burnin]
# run$A = run$A[,-burnin]
# run$p = run$p[-burnin]
# save(run, file = "res_sim_1709_par4_low.Rdata")

burnin = 1
plot(1:length(run$p[-burnin]), run$p[-burnin], type = "l")
lines(1:length(run$p[-burnin]), cumsum(run$p[-burnin])/1:length(run$p[-burnin]), col =2)

plot(1:length(run$sigma2[-burnin]), run$sigma2[-burnin], type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(run$sigma2[-burnin]), cumsum(run$sigma2[-burnin])/1:length(run$sigma2[-burnin]), col =2)

plot(1:length(run$tau[-burnin]), run$tau[-burnin], type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(run$tau[-burnin]), cumsum(run$tau[-burnin])/1:length(run$tau[-burnin]), col =2)

plot(1:length(run$b[-burnin]), run$b[-burnin], type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run$b[-burnin]), cumsum(run$b[-burnin])/1:length(run$b[-burnin]), col =2)

plot(1:length(run$gamma[-burnin]), run$gamma[-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run$gamma[-burnin]), cumsum(run$gamma[-burnin])/1:length(run$gamma[-burnin]), col =2)


mean(run$sigma2[-burnin])
mean(run$tau2[-burnin])
mean(run$b[-burnin])
mean(run$gamma[-burnin])

apply(run$A, 2, function(x) max(which(x>0)) )
apply(run$A, 2, function(x) length(which(x>0)) )
maxx = max(apply(run$A, 2, function(x) max(which(x>0)) ))
t(run$A[1:maxx,100:120]) 

run$clusterO[377:387,1:10]
run$clusterO[,150]
colSums(run$clusterO>0)
run$clusterD[,1:10]

