library(Rcpp)
library(RcppDist)
# install.packages("RcppProgress")
sourceCpp('calcium_DP_mixture_nonconjugate.cpp')

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

plot(diff(y,1), type = "l")
which(diff(y,1) > 0.5)+1

clus = data$s
clus[clus>0] = data$k
data$k
c(2.6, 1, 0.5)[data$k]

1-length(data$k)/length(data$y)


A_start = rep(0,50)
n = length(y)
A_start[2:4] = c(2.6, 1, 0.5)
clus = clus[1:n]

plot(function(x) dgamma(x, 200, 150), xlim=c(0,3) )
plot(function(x) dnorm(x, 1.5, sqrt(150/100^2)), xlim=c(0,3), add=TRUE, col=2 )


nrep = 1000
run_nonc = calcium_gibbs_nonc(Nrep = nrep, y = y, 
                            cal = #rep(0,n+1),
                              c(0, data$c ),
                            cl = #rep(0,n),
                              clus, 
                            A_start = A_start,
                            b_start = 0,
                            gamma_start = 0.9, lambda_start = 500, 
                            p_start = 0.9985, 
                            c0 = 0, varC0 = 0.4, 
                            tau2 = 0.0001, 
                            alpha = 1, m = 20,
                            hyp_A1 = 200, hyp_A2 = 1/150,  # funziona se parte dai punti giusti
                            hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_lambda1 = 50, hyp_lambda2 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 1,
                            hyp_p1 = 999, hyp_p2 = 100,
                            eps_gamma = 0.01,
                            eps_A = 0.2)

# while( sum(run_nonc$gamma == 0) > (nrep/2) )
# {
#   run_nonc = NULL
#   run_nonc = calcium_gibbs_run_nonc(Nrep = nrep, y = y, 
#                               cal = #rep(0,n+1),#
#                                 c(0, data$c ),
#                               cl = #rep(0,n),#
#                                 clus, 
#                               A_start = A_start,
#                               b_start = 0,
#                               gamma_start = 0.6, lambda_start = 200, 
#                               p_start = 0.996, 
#                               c0 = 0, varC0 = 0.4, 
#                               tau2 = 0.001, 
#                               alpha = 1, m = 3,
#                               hyp_A1 = 100, hyp_A2 = 1/50, 
#                               hyp_b1 = 0, hyp_b2 = 1, 
#                               hyp_lambda1 = 50, hyp_lambda2 = 1, 
#                               hyp_gamma1 = 1, hyp_gamma2 = 1,
#                               hyp_p1 = 999, hyp_p2 = 100,
#                               eps_gamma = 0.012,
#                               eps_A = 0.2)
# }



str(run_nonc) 

plot(1:length(run_nonc$p), run_nonc$p, type = "l")
lines(1:length(run_nonc$p), cumsum(run_nonc$p)/1:length(run_nonc$p), col =2)

plot(1:length(run_nonc$lambda), run_nonc$lambda, type = "l", xlab = "iterazioni", ylab = "lambda")
lines(1:length(run_nonc$lambda), cumsum(run_nonc$lambda)/1:length(run_nonc$lambda), col =2)

plot(1:length(run_nonc$b), run_nonc$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run_nonc$b), cumsum(run_nonc$b)/1:length(run_nonc$b), col =2)

plot(1:length(run_nonc$gamma), run_nonc$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run_nonc$gamma), cumsum(run_nonc$gamma)/1:length(run_nonc$gamma), col =2)




iter = 1:length(run_nonc$gamma)

obs = 1:500
image(1:(length(iter)), 1:(length(obs)), t(run_nonc$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)

obs = 501:1000
image(1:(length(iter)), 1:(length(obs)), t(run_nonc$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)


burnin = 1:500

mean(run_nonc$lambda[-burnin])
mean(run_nonc$b[-burnin])
mean(run_nonc$gamma[-burnin])

plot(0:n, c(0,data$c), type = "l")
lines(0:n, rowMeans(run_nonc$calcium[,-burnin]), col = "turquoise")



plot(1:n, y, type = "l")
AA = matrix(0,nrep,n)
for(i in 1:nrep)
{
  AA[i, t(run_nonc$clus)[i,] >0] = run_nonc$A[run_nonc$clus[run_nonc$clus[,i] >0,i]+1,i]
}

calcium = matrix(0, nrep, n+1)
for(i in 1:nrep)
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i, j-1] * mean(run_nonc$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y, type = "l")
lines(1:n, mean(run_nonc$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(t(run_nonc$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
which( apply(t(run_nonc$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# (380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500)

hist(apply(run_nonc$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


### number of clusters
head(t(run_nonc$cluster[375:383,-burnin]))
barplot( table(unlist(apply(run_nonc$cluster[,-burnin], 2, function(x) unique(x[x>0])))), main="number of cluster")
plot(1:nrep, apply(run_nonc$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
minA = min( which(apply(run_nonc$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA
run_nonc$A = run_nonc$A[1:(minA +1),]
str(run_nonc$A)

out_A_nonc = t(run_nonc$A)
out_A_nonc[900:1000,2:7]

