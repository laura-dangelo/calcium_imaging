library(Rcpp)
library(RcppDist)
sourceCpp('calcium_DP_mixture.cpp')

data <- read.csv("data.csv", header = FALSE)
str(data)
head(data)
# plot(1:nrow(data), data$V1, type = "l")
# plot(1:length(data$V1[35000:45000]), data$V1[35000:45000], type = "l")
y = c(data$V1[35000:45000])
length(y)
str(y)
rm(list = ("data"))
plot(1:length(y), y, type = "l")
y = y[4000:10000]
plot(1:length(y), y, type = "l")
#abline(v = which(y >0.5), col = 2)
plot(diff(y,1), type = "l")

pr = which(diff(y,1)>0.22) + 1
plot(1:length(y), y, type = "l")
abline(v = pr, col = 2)

clus = rep(0, length(y))
clus[pr] = 1
mean(y[pr])


out = list()
out$calcium = matrix(0,length(y)+1, 1)
out$cluster = matrix(clus, length(y), 1)
out$A = matrix(0, length(y), 1)
out$A[2,1] = mean(y[pr])
out$AA = matrix(NA, length(y), 1)
out$b = 0
out$gamma = 0.4
out$lambda = 40
out$p = 0.997

nrep = 20
debug = calcium_gibbs_debug(Nrep = nrep, y = y, 
                            cal = out$calcium[,ncol(out$cluster)],
                            cl = out$cluster[,ncol(out$cluster)], A_start = out$A[,ncol(out$cluster)],
                            b_start = out$b[length(out$b)],
                            gamma_start = out$gamma[length(out$b)], lambda_start = out$lambda[length(out$b)], 
                            p_start = out$p[length(out$b)],
                            c0 = 0, varC0 = 0.4, tau2 = 0.00001, 
                            alpha = 1, psi2 = 1, 
                            hyp_A1 = 1, hyp_A2 = 1, 
                            hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 2, 
                            hyp_lambda1 = 10, hyp_lambda2 = 1, 
                            hyp_p1 = 99, hyp_p2 = 1,
                            eps_gamma = 0.2)


str(debug) #1202

out$calcium = cbind(out$calcium, debug$calcium)
out$cluster = cbind(out$cluster, debug$cluster)
out$b = c(out$b, debug$b)
out$gamma = c(out$gamma, debug$gamma)
out$lambda = c(out$lambda, debug$lambda)
out$A = cbind(out$A, debug$A)

n = length(y)

plot(1:length(out$lambda), out$lambda, type = "l", xlab = "iterazioni", ylab = "lambda")
lines(1:length(out$lambda), cumsum(out$lambda)/1:length(out$lambda), col =2)

plot(1:length(out$b), out$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(out$b), cumsum(out$b)/1:length(out$b), col =2)

plot(1:length(out$gamma), out$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$gamma), cumsum(out$gamma)/1:length(out$gamma), col =2)


burnin = 1:800
mean(out$lambda[-burnin])
mean(out$b[-burnin])
mean(out$gamma[-burnin])

plot(0:n, rowMeans(out$calcium[,-burnin]), type = "l")


obs = 1:500
iter = 800:1800
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



burnin = 1:600
plot(1:n, y, type = "l")
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

plot(1:n, y, type = "l")
lines(1:n, mean(out$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.5), lty = 3, col = "salmon")
which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)


hist(apply(out$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")



minA = min( which(  apply(out$A[-1,-burnin], 1, function(x) sum(x == 0)) == length(out$b) - (max(burnin))) )
minA
out$A = out$A[1:(minA +1),]
str(out$A)

out_A = t(out$A)
str(out_A)

sortA = t(apply(out_A, 1, sort, decreasing = TRUE))
str(sortA)
sortA[1:20,1:minA]
round(colMeans(sortA),2)


cluster = t(out$cluster[,-burnin])
str(cluster)
cluster[, which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))<0.8)] = 0
table(apply(cluster,1, function(x) length(unique(x))))

spike_times =  which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
est_spike = matrix(NA, nrep - max(burnin), length(spike_times) )
for(i in 1:nrow(est_spike))
{
  est_spike[i,] = out_A[-burnin,][i, cluster[i, spike_times ]+1 ]
}

colMeans(est_spike)
