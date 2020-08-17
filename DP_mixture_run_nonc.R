library(Rcpp)
library(RcppDist)
# install.packages("RcppProgress")
sourceCpp('calcium_DP_mixture_nonconjugate.cpp')

# data <- read.csv("data.csv", header = FALSE)
# str(data)
# head(data)
# # plot(1:nrow(data), data$V1, type = "l")
# # plot(1:length(data$V1[35000:45000]), data$V1[35000:45000], type = "l")
# y = c(data$V1[35000:45000])
# length(y)
# str(y)
# rm(list = ("data"))
# plot(1:length(y), y, type = "l")
# 
# #abline(v = which(y >0.5), col = 2)
# pr = c(1829, 2716 , 4145, 4150 , 5048, 5081 , 5091, 6860 , 6877 , 9580, 9581, 9582, 9583, 9590)
# abline(v = c(1829, 2716 , 4145, 4150 , 5048, 5081 , 5091, 6860 , 6877 , 9580, 9581, 9582, 9583, 9590), col = 2)

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

data <- sim_data(n = 500, lambda = 10, time_spike = c(50,52, 140, 180, 250, 350, 420, 421, 460),
                 gamma = 0.8, b = 0,
                 prob = c(0.23, 0.44), par = c(4, 10))
set.seed(123)
data <- sim_data(n = 3000, lambda = 200, time_spike = c(380, 550, 820, 821, 1060,1600, 1850,1852, 2662, 2802,2803,2904),
                 gamma = 0.6, b = 0,
                 prob = c(0.23, 0.44, 0.33), par = c(2.6, 1, 0.5))
y = data$y
plot(y, type = "l")


clus = data$s
clus[clus>0] = data$k
data$k
c(2.6, 1, 0.5)[data$k]

1-length(data$k)/length(data$y)


A_start = rep(0,50)


plot(function(x) dnorm(x, 2.5, 0.3), xlim=c(-0.2,3))
plot(function(x) dgamma(x, 80, 40), xlim=c(-0.2,3))


y = y[1:1000]
n = length(y)
A_start[2:4] = c(2.6, 1, 0.5)
clus = clus[1:n]

nrep = 500

debug = calcium_gibbs_debug(Nrep = nrep, y = y, 
                            cal = c(0, data$c[1:1000]),
                            cl = clus, 
                            A_start = A_start,
                            b_start = 0,
                            gamma_start = 0.6, lambda_start = 200, 
                            p_start = 0.996, 
                            c0 = 0, varC0 = 0.4, 
                            tau2 = 0.001, 
                            alpha = 1, m = 3,
                            hyp_A1 = 80, hyp_A2 = 1/40, 
                            hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_lambda1 = 50, hyp_lambda2 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 1,
                            hyp_p1 = 999, hyp_p2 = 1,
                            eps_gamma = 0.012)


str(debug) #1202


plot(1:length(debug$p), debug$p, type = "l")
lines(1:length(debug$p), cumsum(debug$p)/1:length(debug$p), col =2)

plot(1:length(debug$lambda), debug$lambda, type = "l", xlab = "iterazioni", ylab = "lambda")
lines(1:length(debug$lambda), cumsum(debug$lambda)/1:length(debug$lambda), col =2)

plot(1:length(debug$b), debug$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(debug$b), cumsum(debug$b)/1:length(debug$b), col =2)

plot(1:length(debug$gamma), debug$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(debug$gamma), cumsum(debug$gamma)/1:length(debug$gamma), col =2)




iter = 1:length(debug$gamma)

obs = 1:500
image(1:(length(iter)), 1:(length(obs)), t(debug$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)

obs = 501:1000
image(1:(length(iter)), 1:(length(obs)), t(debug$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)


burnin = 1:200

mean(debug$lambda[-burnin])
mean(debug$b[-burnin])
mean(debug$gamma[-burnin])

plot(0:n, c(0,data$c), type = "l")
lines(0:n, rowMeans(debug$calcium[,-burnin]), col = "turquoise")



plot(1:n, y, type = "l")
AA = matrix(0,nrep,n)
for(i in 1:nrep)
{
  AA[i, t(debug$clus)[i,] >0] = debug$A[debug$clus[debug$clus[,i] >0,i]+1,i]
}

calcium = matrix(0, nrep, n+1)
for(i in 1:nrep)
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i, j-1] * mean(debug$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y, type = "l")
lines(1:n, mean(debug$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(t(debug$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
which( apply(t(debug$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)

hist(apply(debug$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


### number of clusters
#head(t(debug$cluster[375:383,-burnin]))
barplot( table(unlist(apply(debug$cluster[,-burnin], 2, function(x) unique(x[x>0])))), horiz = T, main="number of cluster")


### cluster parameter
minA = min( which(apply(debug$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA
debug$A = debug$A[1:(minA +1),]
str(debug$A)

out_A = t(debug$A)
head(out_A)

sortA = t(apply(out_A[-burnin,], 1, sort, decreasing = TRUE))
apply(sortA, 2, function(x) mean(x>0))
round(colMeans(sortA),2)

colMeans(AA[-burnin,which( apply(t(debug$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)])

