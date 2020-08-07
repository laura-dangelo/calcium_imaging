library(Rcpp)
library(RcppDist)
# install.packages("RcppProgress")
sourceCpp('calcium_DP_mixture.cpp')

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
  A[time_spike] <- rnorm(sum(s), par[k], 0.2)

  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i] + rnorm(1,0,0.01)
  }
  return(list("y" = b + c + rnorm(n, 0, 1/sqrt(lambda)), "c" = c, "s" = s, "A" = A, "k" = k))
}

data <- sim_data(n = 500, lambda = 10, time_spike = c(50,52, 140, 180, 250, 350, 420, 421, 460),
                 gamma = 0.8, b = 0,
                 prob = c(0.23, 0.44), par = c(4, 10))
# data <- sim_data(n = 200, lambda = 10, time_spike = c(15,50,52,66,82,130,168),
#                  gamma = 0.3, b = 0,
#                  prob = c(0.4, 0.6), par = c(4, 10))
y = data$y
plot(y, type = "l")
clus = data$s
clus[clus>0] = data$k
data$k

1-length(data$k)/length(data$y)


A_start = rep(0,50)
# A_start[2] = 4
# A_start[3] = 10


n = length(y)
nrep = 1000
debug = calcium_gibbs_debug(Nrep = nrep, y = y, 
                            cal = rep(0,n+1),
                            cl = rep(0,n), A_start = A_start,
                            b_start = -1,
                            gamma_start = 0.5, lambda_start = 5, 
                            p_start = 0.999, #1-length(data$k)/length(data$y), 
                            c0 = 0, varC0 = 0.4, tau2 = 0.01, 
                            alpha = 1, psi2 = 1, 
                            hyp_A1 = 1, hyp_A2 = 3, hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 2, hyp_lambda1 = 1, hyp_lambda2 = 1, 
                            hyp_p1 = 99, hyp_p2 = 1,
                            eps_gamma = 0.012)


str(debug) #1202





plot(1:length(debug$lambda), debug$lambda, type = "l", xlab = "iterazioni", ylab = "lambda")
lines(1:length(debug$lambda), cumsum(debug$lambda)/1:length(debug$lambda), col =2)

plot(1:length(debug$b), debug$b, type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(debug$b), cumsum(debug$b)/1:length(debug$b), col =2)

plot(1:length(debug$gamma), debug$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(debug$gamma), cumsum(debug$gamma)/1:length(debug$gamma), col =2)

plot(1:length(debug$p), debug$p, type = "l")
lines(1:length(debug$p), cumsum(debug$p)/1:length(debug$p), col =2)


obs = 1:250
iter = 1:length(debug$gamma)
image(1:(length(iter)), 1:(length(obs)), t(debug$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,max(obs), by = 1))

obs = 251:500
image(1:(length(iter)), 1:(length(obs)), t(debug$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)




burnin = 1:600

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


burnin = 1:600
minA = min( which(apply(debug$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA
debug$A = debug$A[1:(minA +1),]
str(debug$A)

out_A = t(debug$A)
str(out_A)

sortA = t(apply(out_A, 1, sort, decreasing = TRUE))
apply(sortA, 2, function(x) mean(x==0))
round(colMeans(sortA),2)
