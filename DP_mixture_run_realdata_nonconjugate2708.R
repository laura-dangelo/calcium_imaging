library(Rcpp)
library(RcppDist)
sourceCpp('calcium_DP_mixture_nonconjugate_2708.cpp')

data <- read.csv("data.csv", header = FALSE)


y_real = c(data$V1[35000:45000])
length(y_real)

rm(list = ("data"))
plot(1:length(y_real), y_real, type = "l")
y_real = y_real[4000:10000]
str(y_real)

n = length(y_real)
nrep = 1500


out = list()
out$calcium = matrix(c(0,y_real),length(y_real)+1, 1)
out$cluster = matrix(0, length(y_real), 1)
out$A = matrix(0, 100, 1)
out$AA = matrix(NA, length(y_real), 1)
out$b = 0
out$gamma = 0.8
out$sigma2 = 0.001
out$tau2 = 0.001
out$p = 0.99

plot(function(x) dgamma(x, 9, 4), 0, 5)
debug = calcium_gibbs(Nrep = nrep, y = y_real, 
                      cal = c(out$calcium[,length(out$b)]),
                      cl = c(out$cluster[,length(out$b)]), 
                      A_start = c(out$A[,length(out$b)]),
                      b_start = c(out$b[length(out$b)]),
                      gamma_start = c(out$gamma[length(out$b)]), 
                      sigma2_start = c(out$sigma2[length(out$b)]), 
                      tau2_start = c(out$tau2[length(out$b)]), 
                      p_start = c(out$p[length(out$b)]), 
                      c0 = 0, varC0 = 0.1, 
                      alpha = 1, m = 5,
                      hyp_A1 = 9, hyp_A2 = 4, 
                      hyp_b1 = 0, hyp_b2 = 1, 
                      hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                      hyp_tau21 = 1000, hyp_tau22 = 1, 
                      hyp_gamma1 = 1, hyp_gamma2 = 1,
                      hyp_p1 = 99, hyp_p2 = 1,
                      eps_gamma = 0.005,
                      eps_A = 0.002)



out$calcium = cbind(out$calcium, debug$calcium)
out$cluster = cbind(out$cluster, debug$cluster)
out$b = c(out$b, debug$b)
out$gamma = c(out$gamma, debug$gamma)
out$sigma2 = c(out$sigma2, debug$sigma2)
out$tau2 = c(out$tau2, debug$tau2)
out$A = cbind(out$A, debug$A)
out$p = c(out$p, debug$p)




str(out) 


plot(1:length(out$p), out$p, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$p), cumsum(out$p)/1:length(out$p), col =2)

plot(1:length(out$sigma2), out$sigma2, type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(out$sigma2), cumsum(out$sigma2)/1:length(out$sigma2), col =2)

plot(1:length(out$tau2), out$tau2, type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(out$tau2), cumsum(out$tau2)/1:length(out$tau2), col =2)

plot(1:length(out$b), out$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(out$b), cumsum(out$b)/1:length(out$b), col =2)

plot(1:length(out$gamma), out$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$gamma), cumsum(out$gamma)/1:length(out$gamma), col =2)



burnin = 1:700
mean(out$sigma2[-burnin])
mean(out$tau2[-burnin])
mean(out$b[-burnin])
mean(out$gamma[-burnin])

plot(1:n, y_real, col="turquoise", type = "l")
lines(0:n, rowMeans(out$calcium[,-burnin]))


obs = 1:500
iter = 1:length(out$gamma)
image(1:(length(iter)), 1:(length(obs)), t(out$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = FALSE)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,max(obs), by = 1))


obs = 501:1000
image(1:(length(iter)), 1:(length(obs)), t(out$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = FALSE)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)

obs = 1001:1500
image(1:(length(iter)), 1:(length(obs)), t(out$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = FALSE)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)

obs = 1501:2000
image(1:(length(iter)), 1:(length(obs)), t(out$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = FALSE)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)




plot(1:n, y_real, type = "l")
AA = matrix(0,(length(out$b)-max(burnin)),n)
for(i in 1:(length(out$b)-max(burnin)))
{
  ii = i + (max(burnin)) 
  AA[i, t(out$clus)[ii,] >0] = out$A[out$clus[out$clus[,ii] >0,ii]+1,ii]
}

calcium = matrix(0, nrow(AA), n+1)
for(i in 1:nrow(AA))
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i, j-1] * mean(out$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y_real, type = "l", col = "turquoise3")
lines(1:n, mean(out$b[-burnin]) + colMeans(calcium)[2:(n+1)], lwd = 1.5)

abline(v = which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)

hist(apply(out$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


### number of clusters
plot(1:length(out$b), apply(out$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
minA = min( which(apply(out$A[-1,-burnin], 1, function(x) sum(x == 0)) == length(out$b)-(max(burnin))) )
minA -1
out$A = out$A[1:(minA +1),]
str(out$A)

out_A = t(out$A)
out_A[1100:1200,1:(minA)]

