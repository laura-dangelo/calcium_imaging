data <- read.csv("data.csv", header = FALSE)
str(data)
head(data)
plot(1:nrow(data), data$V1, type = "l")
plot(1:length(data$V1[22000:28000]), data$V1[22000:28000], type = "l")
y = data$V1[22000:28000]


sim_data <- function(n, lambda, time_spike, A, gamma)
{
  c = rep(0,n)
  s = rep(0,n)
  s[time_spike] = 1
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A * s[i]
  }
  return(list("y" = c + rnorm(n, 0, 1/sqrt(lambda)), "c" = c, "s" = s))
}

data <- sim_data(n = 1000, lambda = 10, time_spike = c(100, 101, 103, 200, 300, 600), A = 5, gamma = 0.8)
y = data$y 
plot(y, type = "l")



loglik <- function(y, cc, c_0, s, A, gamma, lambda)
{
  n = length(y)
  sum(dnorm(y, mean = gamma * c(c_0, cc[1:(n-1)]) + A * s, sd = 1/sqrt(lambda), log = TRUE))
}

# ### verosimiglianza per A
# A.seq = seq(0.1, 10, length.out = 500)
# lik.seq = sapply(A.seq, function(x) loglik(y = y, cc = data$c, c_0 = 0,
#                                            s = data$s,
#                                            A = x, gamma = 0.5, lambda = 10))
# plot(A.seq, lik.seq, type = "l")

logprior_A <- function(A, a, b) dgamma(A, a, b, log = TRUE)
logprior_gamma <- function(gamma, c, d) dbeta(gamma, c, d, log = TRUE)
logprior_lambda <- function(lambda, e, f) dgamma(lambda, e, f, log = TRUE)

logpost <- function(y, cc, c_0, s, A, gamma, lambda, 
                    a = 1, b = 1, c = 1, d = 1, e = 1, f = 1)
{
  loglik(y, cc, c_0, s, A, gamma, lambda) + logprior_A(A, a, b) + 
    logprior_gamma(gamma, c, d) + logprior_lambda(lambda, e,f)
}

### posteriori per A
# A.seq = seq(0.1, 10, length.out = 500)
# post.seq = sapply(A.seq, function(x) logpost(y = y, cc = data$c, c_0 = 0,
#                                            s = data$s,
#                                            A = x, gamma = 0.3, lambda = 10,
#                                            a = 3, b = 1, c = 5.5, d = 2,
#                                            e = 10, f = 1))
# plot(A.seq, post.seq, type = "l")
# abline(v = A.seq[which.max(post.seq)], lty = 2)

### posteriori per gamma
# gamma.seq = seq(0.2, 0.998, length.out = 500)
# post.seq = sapply(gamma.seq, function(x) logpost(y = y, cc = data$c, c_0 = 0,
#                                            s = data$s,
#                                            A = 5, gamma = x, lambda = 10,
#                                            a = 3, b = 1, c = 5.5, d = 2,
#                                            e = 10, f = 1))
# plot(gamma.seq, post.seq, type = "l")
# abline(v = gamma.seq[which.max(post.seq)], lty = 3)

#---------------------# #---------------------# #---------------------# #---------------------# 
#---------------------# #---------------------# #---------------------# #---------------------# 

gibbs_calcium <- function(nrep, y,
                          gamma_start = 0.5, A_start = 5,
                          tau2 = 0.0001, lambda_start = 10, c0 = 0,
                          p = 0.005, 
                          a = 3, b = 1, c = 5.5, d = 2,
                          e = 10, f = 1,
                          eps_gamma, eps_A)
{
  n = length(y)
  # creo matrici di output
  out_c = matrix(NA, nrep, n)
  out_s = matrix(NA, nrep, n)
  out_p1 = matrix(NA, nrep, n)
  out_p0 = matrix(NA, nrep, n)
  out_A = rep(NA, nrep) 
  out_gamma = rep(NA, nrep) 
  out_lambda = rep(NA, nrep)
  
  filter_mean = numeric(n)
  filter_var = numeric(n)
  
  # inizializzazione dei parametri
  out_c[1,] = 0
  out_s[1,] = 0
  out_A[1] = A_start
  out_gamma[1] = gamma_start
  out_lambda[1] = lambda_start

  
  for(i in 1:(nrep-1))
  {
    sigma2 = 1/out_lambda[i]
    
    # sampling di s
    out_p1[i+1,] = exp(-.5 / sigma2 * (y - out_gamma[i] * c(c0, out_c[i,1:(n-1)]) - out_A[i])^2 ) * p
    out_p0[i+1,] = exp(-.5 / sigma2 * (y - out_gamma[i] * c(c0, out_c[i,1:(n-1)]))^2 ) * (1-p)
    out_s[i+1,] = apply(cbind(out_p0[i+1,], out_p1[i+1,]), 1, function(x) sample(c(0,1), 1, prob = x))
    
    # sampling di c
    filter_mean[1] = (sigma2 * (out_gamma[i] * c0 + out_A[i] * out_s[i,1]) + tau2 * y[1]) / (tau2 + sigma2)
    filter_var[1] = tau2 * sigma2 / (tau2 + sigma2)
    
    for(j in 2:n)
    {
      filter_mean[j] = (sigma2 * (out_gamma[i] * filter_mean[j-1] + out_A[i] * out_s[i+1,j]) + 
                          (filter_var[j-1] + tau2) * y[j]) / (filter_var[j-1] + tau2 + sigma2)
      filter_var[j] = ((filter_var[j-1] + tau2) * sigma2) / (filter_var[j-1] + tau2 + sigma2)
    }
    out_c[i+1,n] = rnorm(1, filter_mean[n], filter_var[n])
    
    for(j in (n-1):1)
    {
      back_mean = filter_mean[j] + out_gamma[i] * filter_var[j]/(filter_var[j] + tau2) * 
        (out_c[i+1,j+1] - out_gamma[i] * filter_mean[j] - out_A[i] * out_s[i+1,j+1])
      back_var = filter_var[j] - out_gamma[i]^2 * filter_var[j]^2 / (filter_var[j] + tau2)
      out_c[i+1,j] = rnorm(1, back_mean, back_var)
    }
    
    # sampling di lambda (precision)
    mean_y = out_gamma[i] * c(c0, out_c[i+1,1:(n-1)]) + out_A[i] * out_s[i+1,]
    ssum = sum(apply(cbind(y,mean_y), 1, function(x) (x[1] - x[2])^2))
    out_lambda[i+1] = rgamma(1, e + n/2, f + .5 * ssum)
    
    # # MH per gamma: random walk
    oldgamma = out_gamma[i]
    newgamma = oldgamma + runif(1, -eps_gamma, eps_gamma)
    alpha = exp( logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                         A = out_A[i], gamma = newgamma, lambda = out_lambda[i+1],
                         a = a, b = b, c = c, d = d, e = e, f = f) -
                  logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                          A = out_A[i], gamma = oldgamma, lambda = out_lambda[i+1],
                          a = a, b = b, c = c, d = d, e = e, f = f) )
    if(runif(1) < alpha) oldgamma = newgamma
    out_gamma[i+1] = oldgamma
    
    # MH per A: uso RW normale
    oldA = out_A[i]
    newA = rnorm(1, oldA, eps_A)
    alpha = exp( logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                         A = newA, gamma = out_gamma[i+1], lambda = out_lambda[i+1],
                         a = a, b = b, c = c, d = d, e = e, f = f) -
                   logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                           A = oldA, gamma = out_gamma[i+1], lambda = out_lambda[i+1],
                           a = a, b = b, c = c, d = d, e = e, f = f) )
    if(runif(1) < alpha) oldA = newA
    out_A[i+1] = oldA
  }
  return(list(c = out_c, s = out_s, lambda = out_lambda, A = out_A, gamma = out_gamma))
}
## magari il kalman filter lo posso fare solo per un tot di iterazioni - burnin ?

nrep = 1500
prova <- gibbs_calcium(nrep = nrep, y = y, 
                       lambda_start = 10,
                       gamma_start = 0.5, A_start = 5,
                       eps_gamma = 0.07, 
                       eps_A = 0.2)
str(prova)

burnin = 1:700
plot(1:nrep, prova$lambda, type = "l", main = "lambda")
lines(1:nrep, cumsum(prova$lambda)/1:nrep, col = 2)
abline(h=10, col = "blue")

plot(1:nrep, prova$gamma, type = "l", main = "gamma")
lines(1:nrep, cumsum(prova$gamma)/1:nrep, col = 2)
abline(h=0.8, col = "blue")

plot(1:nrep, prova$A, type = "l", main = "A")
lines(1:nrep, cumsum(prova$A)/1:nrep, col = 2)
abline(h=5, col = "blue")

plot(1:1000, colMeans(prova$c[-burnin,]), type = "l")
lines(1:1000, data$c, col = "blue", lty = 3)

plot(1:1000, y, type = "l")
lines(1:1000, colMeans(prova$c[-burnin,]), col = 2)
abline(v = which(colMeans(prova$s[-burnin,])>0.5), col = "salmon", lty = 2)













