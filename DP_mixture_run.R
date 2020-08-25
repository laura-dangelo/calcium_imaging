library(Rcpp)
library(RcppDist)
# install.packages("RcppProgress")
sourceCpp('calcium_DP_mixture.cpp')

## mixture ##
sim_data <- function(n, lambda, time_spike, b, gamma, prob, par)
{
  c = rep(0,n)
  s = rep(0,n)
  A = rep(0,n)
  s[time_spike] = 1
  
  k <- sample(1:length(prob), sum(s), prob, replace = TRUE)
  A[time_spike] <- par[k]
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i] + rnorm(1,0,0.01)
  }
  return(list("y" = b + c + rnorm(n, 0, 1/sqrt(lambda)), "c" = c, "s" = s, "A" = A, "k" = k))
}


set.seed(123)
data <- sim_data(n = 6000, lambda = 500, time_spike = c(380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500),
                 gamma = 0.9, b = 0,
                 prob = c(0.23, 0.44, 0.33), par = c(2.6, 1.1, 0.5))
y = data$y
plot(y, type = "l")


clus = data$s
clus[clus>0] = data$k
data$k
c(2.6, 1, 0.5)[data$k]

1-length(data$k)/length(data$y)


A_start = rep(0,50)
n = length(y)
A_start[2:4] = c(2.6, 1, 0.5)
clus = clus[1:n]

plot(function(x) dnorm(x, 1.8, 0.4), 0, 3)

nrep = 500
run = calcium_gibbs(Nrep = nrep, y = y, 
                      cal = #c(0, data$c ), 
                    rep(0,n+1),
                      cl = #clus, #
                    rep(0,n), 
                      A_start = A_start,
                      b_start = 0,
                      gamma_start = 0.9, lambda_start = 500, 
                      p_start = 0.9985, 
                      c0 = 0, varC0 = 0.4, 
                      tau2 = 0.0001,
                      alpha = 1, 
                      hyp_A1 = 2.2, hyp_A2 = 0.4^2, 
                      hyp_b1 = 0, hyp_b2 = 1, 
                      hyp_lambda1 = 50, hyp_lambda2 = 1, 
                      hyp_gamma1 = 1, hyp_gamma2 = 1,
                      hyp_p1 = 999, hyp_p2 = 100,
                      eps_gamma = 0.008)

while( sum(run$gamma == 0) > (nrep/2) )
{
  run = NULL
  run = calcium_gibbs(Nrep = nrep, y = y, 
                      cal = #c(0, data$c ), 
                        rep(0,n+1),
                      cl = #clus, #
                        rep(0,n), 
                      A_start = A_start,
                      b_start = 0,
                      gamma_start = 0.9, lambda_start = 500, 
                      p_start = 0.9985, 
                      c0 = 0, varC0 = 0.4, 
                      tau2 = 0.0001,
                      alpha = 1, 
                      hyp_A1 = 2.2, hyp_A2 = 0.4^2, 
                      hyp_b1 = 0, hyp_b2 = 1, 
                      hyp_lambda1 = 50, hyp_lambda2 = 1, 
                      hyp_gamma1 = 1, hyp_gamma2 = 1,
                      hyp_p1 = 999, hyp_p2 = 100,
                      eps_gamma = 0.008)
}




str(run) #1202


plot(1:length(run$p), run$p, type = "l")
lines(1:length(run$p), cumsum(run$p)/1:length(run$p), col =2)

plot(1:length(run$lambda), run$lambda, type = "l", xlab = "iterazioni", ylab = "lambda")
lines(1:length(run$lambda), cumsum(run$lambda)/1:length(run$lambda), col =2)

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

mean(run$lambda[-burnin])
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
head(t(run$cluster[375:383,-burnin]))
barplot( table(unlist(apply(run$cluster[,-burnin], 2, function(x) unique(x[x>0])))), main="number of cluster")
plot(1:nrep, apply(run$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
minA = min( which(apply(run$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA
run$A = run$A[1:(minA +1),]
str(run$A)

out_A = t(run$A)
out_A[900:1000,2:7]

