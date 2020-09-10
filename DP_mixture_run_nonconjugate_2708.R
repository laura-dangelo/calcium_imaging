library(Rcpp)
library(RcppDist)
# install.packages("RcppProgress")
sourceCpp('calcium_DP_mixture_nonconjugate_2708.cpp')

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


set.seed(123)
data <- sim_data(n = 6000, sigma2 = 0.002, tau2 = 0.0005, time_spike = c(380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500),
                 gamma = 0.9, b = 0,
                 prob = c(0.23, 0.44, 0.33), par = c(2.6, 1.1, 0.6))


y = data$y
plot(y, type = "l")


clus = data$s
clus[clus>0] = data$k
data$k
c(2.6, 1, 0.5)[data$k]

1-length(data$k)/length(data$y)


A_start = rep(0,50)
#A_start[2:4] = c(2.6, 1.1, 0.6)
n = length(y)

plot(function(x) dgamma(x, 10,5),0,4)

nrep = 500
set.seed(123)
run = calcium_gibbs(Nrep = nrep, y = y, 
                    cal = c(0,y),
                    cl = rep(0,n), 
                    A_start = A_start,
                    b_start = 0,
                    gamma_start = 0.8, 
                    sigma2_start = 0.002, 
                    tau2_start = 0.001, 
                    p_start = 0.99, 
                    c0 = 0, varC0 = 0.1, 
                    alpha = 1, m = 5,
                    hyp_A1 = 10, hyp_A2 = 5, 
                    kappa = 0,
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                    hyp_tau21 = 1000, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 99, hyp_p2 = 1,
                    eps_gamma = 0.005,
                    eps_A = 0.002)


str(run) 
burnin = 1:300
which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# (380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500)

plot(1:length(run$p), run$p, type = "l")
lines(1:length(run$p), cumsum(run$p)/1:length(run$p), col =2)

plot(1:length(run$sigma2), run$sigma2, type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(run$sigma2), cumsum(run$sigma2)/1:length(run$sigma2), col =2)

plot(1:length(run$tau), run$tau, type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(run$tau), cumsum(run$tau)/1:length(run$tau), col =2)

plot(1:length(run$b), run$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run$b), cumsum(run$b)/1:length(run$b), col =2)

plot(1:length(run$gamma), run$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run$gamma), cumsum(run$gamma)/1:length(run$gamma), col =2)




iter = 1:length(run$gamma)
obs = 1:500
image(1:(length(iter)), 1:(length(obs)), t(run$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)

obs = 501:1000
image(1:(length(iter)), 1:(length(obs)), t(run$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)


burnin = 1:300

mean(run$sigma2[-burnin])
mean(run$tau2[-burnin])
mean(run$b[-burnin])
mean(run$gamma[-burnin])

plot(0:n, c(0,data$c), type = "l")
lines(0:n, rowMeans(run$calcium[,-burnin]), col = "turquoise")



plot(1:n, y, type = "l")
AA = matrix(0,nrep,n)
for(i in 1:nrep)
{
  AA[i, t(run$clus)[i,] >0] = run$A[run$clus[run$clus[,i] >0,i]+1,i]
}

calcium = matrix(0, nrep, n+1)
for(i in 1:nrep)
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i, j-1] * mean(run$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y, type = "l")
lines(1:n, mean(run$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# (380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500)

hist(apply(run$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


### number of clusters
# head(t(run$cluster[375:383,-burnin]))
# barplot( table(unlist(apply(run$cluster[,-burnin], 2, function(x) unique(x[x>0])))), main="number of cluster")
plot(1:nrep, apply(run$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
minA = min( which(apply(run$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA -1
run$A = run$A[1:(minA +1),]
str(run$A)

out_A = t(run$A)
out_A[400:500,1:(minA)]










set.seed(1234)
data2 <- sim_data(n = 6000, sigma2 = 0.002, tau2 = 0.0003, time_spike = c(380,385,
                                                                          1000, 1002, 1005, 1006, 1007, 1008, 1011, 1012, 1014,
                                                                          1300, 
                                                                          2000, 
                                                                          2990, 2991, 2993, 2995, 2999, 3001, 3002,
                                                                          4000, 4004,
                                                                          4700, 4701, 4703,
                                                                          5500),
                 gamma = 0.9, b = 0,
                 prob = c(0.23, 0.44, 0.33), par = c(0.35, 0.8, 0.6))


y2 = data2$y
plot(y2, type = "l")


clus = data2$s
clus[clus>0] = data2$k
data2$k
c(2.6, 1, 0.5)[data2$k]

1-length(data2$k)/length(data2$y)


A_start = rep(0,50)
#A_start[2:4] = c(2.6, 1.1, 0.6)
n = length(y2)

plot(function(x) dgamma(x, 10,5),0,4)
# media (hyp_A1 + 2*kappa) / hyp_A2
# varianza (hyp_A1 + 2*kappa) / hyp_A2^2


hyp_A1 = 5
hyp_A2 = 10

hyp_A1 / hyp_A2
hyp_A1 / hyp_A2^2


nrep = 500
set.seed(123)
run2 = calcium_gibbs(Nrep = nrep, y = y2, 
                    cal = c(0,y2),
                    cl = rep(0,n), 
                    A_start = A_start,
                    b_start = 0,
                    gamma_start = 0.8, 
                    sigma2_start = 0.002, 
                    tau2_start = 0.001, 
                    p_start = 0.99, 
                    c0 = 0, varC0 = 0.1, 
                    alpha = 1, m = 5,
                    hyp_A1 = hyp_A1, hyp_A2 = hyp_A2, 
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                    hyp_tau21 = 1000, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 99, hyp_p2 = 1,
                    eps_gamma = 0.005,
                    eps_A = 0.002)


str(run2) 
burnin = 1:300

plot(1:length(run2$p), run2$p, type = "l")
lines(1:length(run2$p), cumsum(run2$p)/1:length(run2$p), col =2)

plot(1:length(run2$sigma2), run2$sigma2, type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(run2$sigma2), cumsum(run2$sigma2)/1:length(run2$sigma2), col =2)

plot(1:length(run2$tau), run2$tau, type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(run2$tau), cumsum(run2$tau)/1:length(run2$tau), col =2)

plot(1:length(run2$b), run2$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run2$b), cumsum(run2$b)/1:length(run2$b), col =2)

plot(1:length(run2$gamma), run2$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run2$gamma), cumsum(run2$gamma)/1:length(run2$gamma), col =2)




iter = 1:length(run2$gamma)
obs = 1:500
image(1:(length(iter)), 1:(length(obs)), t(run2$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)

obs = 501:1000
image(1:(length(iter)), 1:(length(obs)), t(run2$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)


burnin = 1:300

mean(run2$sigma2[-burnin])
mean(run2$tau2[-burnin])
mean(run2$b[-burnin])
mean(run2$gamma[-burnin])

plot(0:n, c(0,data2$c), type = "l")
lines(0:n, rowMeans(run2$calcium[,-burnin]), col = "turquoise")



plot(1:n, y2, type = "l")
AA = matrix(0,nrep,n)
for(i in 1:nrep)
{
  AA[i, t(run2$clus)[i,] >0] = run2$A[run2$clus[run2$clus[,i] >0,i]+1,i]
}

calcium = matrix(0, nrep, n+1)
for(i in 1:nrep)
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i, j-1] * mean(run2$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y2, type = "l")
lines(1:n, mean(run2$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(t(run2$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
which( apply(t(run2$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
which( apply(t(run2$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# c(380,385,
#   1000, 1002, 1005, 1006, 1007, 1008, 1011, 1012, 1014,
#   1300, 
#   2000, 
#   2990, 2991, 2993, 2995, 2999, 3001, 3002,
#   4000, 4004,
#   4700, 4701, 4703,
#   5500)


hist(apply(run2$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


### number of clusters
# head(t(run2$cluster[375:383,-burnin]))
# barplot( table(unlist(apply(run2$cluster[,-burnin], 2, function(x) unique(x[x>0])))), main="number of cluster")
plot(1:nrep, apply(run2$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
minA = min( which(apply(run2$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA -1
run2$A = run2$A[1:(minA +1),]
str(run2$A)

out_A = t(run2$A)
out_A[400:500,1:(minA)]

