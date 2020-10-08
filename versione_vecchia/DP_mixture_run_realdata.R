library(Rcpp)
library(RcppDist)
sourceCpp('calcium_DP_mixture_nonconjugate.cpp')

data <- read.csv("data.csv", header = FALSE)
str(data)
head(data)
# plot(1:nrow(data), data$V1, type = "l")
# plot(1:length(data$V1[35000:45000]), data$V1[35000:45000], type = "l")
y_real = c(data$V1[35000:45000])
length(y_real)
str(y_real)
rm(list = ("data"))
plot(1:length(y_real), y_real, type = "l")
y_real = y_real[4000:10000]

n = length(y_real)
plot(1:length(y_real), y_real, type = "l")
#plot(diff(y_real,1), type = "l")

pr = which(diff(y_real,1)>0.2) + 1
abline(v=pr, col = "salmon", lty = 3)

clus = rep(0, length(y_real))
clus[pr] = 1
mean(y_real[pr])

out = list()
out$calcium = matrix(0,length(y_real)+1, 1)
out$cluster = matrix(clus, length(y_real), 1)
out$A = matrix(0, 100, 1)
out$A[2,1] = mean(y_real[pr])+0.5
out$AA = matrix(NA, length(y_real), 1)
out$b = 0
out$gamma = 0.5
out$lambda = 100
out$p = 0.999

plot(function(x) dgamma(x, 30, 15), xlim=c(-0.2,3))

nrep = 500
while(length(out$gamma)==1)
{
  debug = calcium_gibbs_debug(Nrep = nrep, y = y_real, 
                              cal = out$calcium[,ncol(out$cluster)],
                              cl = out$cluster[,ncol(out$cluster)], 
                              A_start = out$A[,ncol(out$cluster)],
                              b_start = out$b[length(out$b)],
                              gamma_start = out$gamma[length(out$b)], 
                              lambda_start = out$lambda[length(out$b)], 
                              p_start = out$p[length(out$b)],
                              c0 = 0, varC0 = 0.4, 
                              tau2 = 0.01, 
                              alpha = 1, m = 5,
                              hyp_A1 = 30, hyp_A2 = 1/15, 
                              hyp_b1 = 0, hyp_b2 = 1, 
                              hyp_lambda1 = 100, hyp_lambda2 = 1, 
                              hyp_gamma1 = 1, hyp_gamma2 = 2,
                              hyp_p1 = 999, hyp_p2 = 1,
                              eps_gamma = 0.012,
                              eps_A = 0.05)
  
  if(sum(debug$gamma==0) < nrep/2)
  {
    out$calcium = cbind(out$calcium, debug$calcium)
    out$cluster = cbind(out$cluster, debug$cluster)
    out$b = c(out$b, debug$b)
    out$gamma = c(out$gamma, debug$gamma)
    out$lambda = c(out$lambda, debug$lambda)
    out$A = cbind(out$A, debug$A)
    out$p = c(out$p, debug$p)
  }
}



str(out) 


plot(1:length(out$p), out$p, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$p), cumsum(out$p)/1:length(out$p), col =2)

plot(1:length(out$lambda), out$lambda, type = "l", xlab = "iterazioni", ylab = "lambda")
lines(1:length(out$lambda), cumsum(out$lambda)/1:length(out$lambda), col =2)

plot(1:length(out$b), out$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(out$b), cumsum(out$b)/1:length(out$b), col =2)

plot(1:length(out$gamma), out$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$gamma), cumsum(out$gamma)/1:length(out$gamma), col =2)


burnin = 1:200
mean(out$lambda[-burnin])
mean(out$b[-burnin])
mean(out$gamma[-burnin])

plot(0:n, rowMeans(out$calcium[,-burnin]), type = "l")
lines(1:n, y_real, col="turquoise", lty = 2)

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

plot(1:n, y_real, type = "l")
lines(1:n, mean(out$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)


hist(apply(out$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


burnin = 1:300

### number of clusters
head(t(out$cluster[1040:1060,-burnin]))
plot(1:(nrep+1), apply(out$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
# minA = min( which( apply(out$A[-1,-burnin], 1, function(x) sum(x == 0) )> (200)  ) )
minA =4
out$A = out$A[1:(minA +1),]
str(debug$A)

out_A = t(debug$A)
out_A[1:200,2:3]
