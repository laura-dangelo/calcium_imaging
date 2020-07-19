library(Rcpp)
library(RcppDist)
sourceCpp('DP_mixture.cpp')



## mixture ##
sim_data <- function(n, lambda, time_spike, b, gamma, prob, par)
{
  c = rep(0,n)
  s = rep(0,n)
  A = rep(0,n)
  s[time_spike] = 1
  
  k <- sample(1:length(prob), sum(s), prob, replace = TRUE)
  A[time_spike] <- rnorm(sum(s), par[k], 0.5)
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i]
  }
  return(list("y" = b + c + rnorm(n, 0, 1/sqrt(lambda)), "c" = c, "s" = s, "A" = A, "k" = k))
}

data <- sim_data(n = 500, lambda = 10, time_spike = c(50,52, 140, 180, 250, 350, 420, 421, 460),
                 gamma = 0.8, b = 0,
                 prob = c(0.4, 0.6), par = c(4, 10))
y = data$y
#plot(y, type = "l")
clus = data$s
clus[clus>0] = data$k

debug = calcium_gibbs_debug(Nrep = 3, y = data$y, 
                    clus = t(clus),
                    gamma_start = 0.3, lambda_start = 13, 
                    c0 = 0, varC0 = 1, tau2 = 0.0001, p = 0.995, 
                    alpha = 1, psi2 = 1, 
                    hyp_A1 = 1, hyp_A2 = 1, hyp_b1 = 1, hyp_b2 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 2, hyp_lambda1 = 10, hyp_lambda2 = 1, eps_gamma = 0.2)


str(debug)
debug$calcium[2,1:100]
debug$A[2,1:10]
debug$cluster[1,]
c(debug$AA)

n = length(y)
burnin = 1:400
plot(0:n, colMeans(debug$calcium[-c(1:400),]), type = "l")
lines(0:n, c(0,data$c), col = 2)

mean(debug$lambda[-burnin])
plot(1:length(debug$lambda[-burnin]), debug$lambda[-burnin], type = "l")
lines(1:length(debug$lambda[-burnin]), cumsum(debug$lambda[-burnin])/1:length(debug$lambda[-burnin]), col =2)

mean(debug$b[-burnin])
plot(1:length(debug$b[-burnin]), debug$b[-burnin], type = "l")
lines(1:length(debug$b[-burnin]), cumsum(debug$b[-burnin])/1:length(debug$b[-burnin]), col =2)

mean(debug$gamma[-burnin])
plot(1:length(debug$gamma[-burnin]), debug$gamma[-burnin], type = "l")
lines(1:length(debug$gamma[-burnin]), cumsum(debug$gamma[-burnin])/1:length(debug$gamma[-burnin]), col =2)



nrep=1000
burnin = 1:500
plot(1:n, y, type = "l")
AA = matrix(0,nrep,n)
for(i in 1:nrep)
{
  AA[i, debug$clus[i,] >0] = debug$A[i, debug$clus[i, debug$clus[i,] >0]+1]
}

calcium = matrix(0, nrep, n+1)
for(i in 1:nrep)
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i,j-1] * mean(debug$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y, type = "l")
lines(1:n, mean(debug$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(debug$clus[-burnin,], 2, function(x) mean(x != 0))>0.95), lty = 3, col = "salmon")


which(apply(debug$A, 2, function(x) sum(!is.na(x))>0))
debug$A = prova$A[,1:20]

debug$A[500:600,1:4]
debug$A[is.na(debug$A)] = 0
maxA = apply(debug$A, 1, max)
mean(maxA)

max2A = t(apply(debug$A, 1, sort))
colMeans(max2A)













data <- sim_data(n = 200, lambda = 10, time_spike = c(10,25,50,52, 140, 180),
                 gamma = 0.8, b = 0,
                 prob = c(0.4, 0.6), par = c(4, 10))
y = data$y
#plot(y, type = "l")
clus = data$s
clus[clus>0] = data$k


pu <- polya_urn(j=180, y = y, cc = data$c, clus = clus, A = c(0,4,10,0,0,0), sigma2 = 1,
          tau2 = 1, b = 0, gamma = 0.8, lambda = 10, p = 0.995, alpha = 1, psi2 = 1)
pu$A
clus[1:200]
c(pu$cluster)











