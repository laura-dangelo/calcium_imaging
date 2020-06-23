data <- read.csv("data.csv", header = FALSE)
str(data)
head(data)
plot(1:nrow(data), data$V1, type = "l")
plot(1:length(data$V1[22000:28000]), data$V1[22000:28000], type = "l")
y = data$V1[22000:28000]

## mixture ##
sim_data <- function(n, lambda, time_spike, A, b, gamma, mix_K, )
{
  c = rep(0,n)
  s = rep(0,n)
  s[time_spike] = 1
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A * s[i]
  }
  return(list("y" = b + c + rnorm(n, 0, 1/sqrt(lambda)), "c" = c, "s" = s))
}

data <- sim_data(n = 1000, lambda = 10, time_spike = c(100, 101, 103, 700, 300, 600), 
                 A = 5, gamma = 0.8, b = 2)
y = data$y 
plot(y, type = "l", main = "")



loglik <- function(y, cc, c_0, s, A, b, gamma, lambda)
{
  n = length(y)
  sum(dnorm(y, mean = b + gamma * c(c_0, cc[1:(n-1)]) + A * s, sd = 1/sqrt(lambda), log = TRUE))
}

# ### verosimiglianza per A
# A.seq = seq(0.1, 10, length.out = 500)
# lik.seq = sapply(A.seq, function(x) loglik(y = y, cc = data$c, c_0 = 0,
#                                            s = data$s, b = 2,
#                                            A = x, gamma = 0.5, lambda = 10))
# plot(A.seq, lik.seq, type = "l")
# abline(v = A.seq[which.max(lik.seq)])

logprior_A <- function(A, hyp_A1, hyp_A2) dgamma(A, hyp_A1, hyp_A2, log = TRUE)
logprior_gamma <- function(gamma, hyp_gamma1, hyp_gamma2) dbeta(gamma, hyp_gamma1, hyp_gamma2, log = TRUE)
# logprior_lambda <- function(lambda, hyp_lambda1, hyp_lambda2) dgamma(lambda, hyp_lambda1, hyp_lambda2, log = TRUE)
# logprior_b <- function(b, hyp_b1, hyp_b2) dnorm(b, mean = hyp_b1, sd = hyp_b2, log = TRUE)

logprior_b_lambda <- function(b, lambda, hyp_b1, hyp_b2, hyp_lambda1, hyp_lambda2)
{
  dnorm(b, mean = hyp_b1, sd = 1 / sqrt(hyp_b2 * lambda), log = TRUE) + dgamma(lambda, hyp_lambda1, hyp_lambda2, log = TRUE)
}

logpost <- function(y, cc, c_0, s, A, b, gamma, lambda, 
                    hyp_A1, hyp_A2, hyp_gamma1, hyp_gamma2, 
                    hyp_lambda1, hyp_lambda2, hyp_b1, hyp_b2)
{
  loglik(y, cc, c_0, s, A, b, gamma, lambda) + logprior_A(A, hyp_A1, hyp_A2) + 
    logprior_gamma(gamma, hyp_gamma1, hyp_gamma2) + 
    logprior_b_lambda(b, lambda, hyp_b1, hyp_b2, hyp_lambda1, hyp_lambda2)
}

### posteriori per A
# A.seq = seq(0.1, 10, length.out = 500)
# post.seq = sapply(A.seq, function(x) logpost(y = y, cc = data$c, c_0 = 0,
#                                            s = data$s,
#                                            A = x, gamma = 0.8, 
#                                            b = 2, lambda = 10,
#                                            hyp_A1 = 3, hyp_A2 = 1, hyp_gamma1 = 5, hyp_gamma2 = 2, 
#                                            hyp_lambda1 = 10, hyp_lambda2 = 1, hyp_b1 = 1, hyp_b2 = 1))
# plot(A.seq, post.seq, type = "l")
# abline(v = A.seq[which.max(post.seq)], lty = 2)

### posteriori per gamma
# gamma.seq = seq(0.2, 0.998, length.out = 500)
# post.seq = sapply(gamma.seq, function(x) logpost(y = y, cc = data$c, c_0 = 0,
#                                            s = data$s,
#                                            A = 5, gamma = x, b = 2, lambda = 10,
#                                            hyp_A1 = 3, hyp_A2 = 1, hyp_gamma1 = 5, hyp_gamma2 = 2, 
#                                            hyp_lambda1 = 10, hyp_lambda2 = 1, hyp_b1 = 1, hyp_b2 = 1))
# plot(gamma.seq, post.seq, type = "l")
# abline(v = gamma.seq[which.max(post.seq)], lty = 3)

#---------------------# #---------------------# #---------------------# #---------------------# 
#---------------------# #---------------------# #---------------------# #---------------------# 

gibbs_calcium <- function(nrep, y,
                          gamma_start = 0.5, A_start = 5, b_start = 2,
                          tau2 = 0.0001, lambda_start = 10, c0 = 0,
                          p = 0.005, 
                          hyp_A1 = 3, hyp_A2 = 1, hyp_gamma1 = 5, hyp_gamma2 = 2,
                          hyp_lambda1 = 10, hyp_lambda2 = 1, hyp_b1 = 1, hyp_b2 = 1,
                          eps_gamma, eps_A)
{
  n = length(y)
  # creo matrici di output
  out_c = matrix(NA, nrep, n)
  out_s = matrix(NA, nrep, n)
  out_p1 = matrix(NA, nrep, n)
  out_p0 = matrix(NA, nrep, n)
  out_A = rep(NA, nrep) 
  out_b = rep(NA, nrep)
  out_gamma = rep(NA, nrep) 
  out_lambda = rep(NA, nrep)
  
  filter_mean = numeric(n)
  filter_var = numeric(n)
  
  # inizializzazione dei parametri
  out_c[1,] = 0
  out_s[1,] = 0
  out_A[1] = A_start
  out_b[1] = b_start
  out_gamma[1] = gamma_start
  out_lambda[1] = lambda_start

  
  for(i in 1:(nrep-1))
  {
    sigma2 = 1/out_lambda[i]
    
    # sampling di s
    out_p1[i+1,] = exp(-.5 / sigma2 * (y - out_b[i] - out_gamma[i] * c(c0, out_c[i,1:(n-1)]) - out_A[i])^2 ) * p
    out_p0[i+1,] = exp(-.5 / sigma2 * (y - out_b[i] - out_gamma[i] * c(c0, out_c[i,1:(n-1)]))^2 ) * (1-p)
    out_s[i+1,] = apply(cbind(out_p0[i+1,], out_p1[i+1,]), 1, function(x) sample(c(0,1), 1, prob = x))
    
    # sampling di c
    filter_mean[1] = (sigma2 * (out_b[i] + out_gamma[i] * c0 + out_A[i] * out_s[i,1]) + tau2 * y[1]) / (tau2 + sigma2)
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
    z = y - out_gamma[i] * c(c0, out_c[i+1,1:(n-1)]) - out_A[i] * out_s[i+1,]
    ssum = sum((z - mean(z))^2)
    out_lambda[i+1] = rgamma(1, hyp_lambda1 + n/2, 
                             hyp_lambda2 + .5 * ssum + .5 * n * hyp_b2 / (n + hyp_b2) * (mean(z) - hyp_b1)^2 )

    # sampling di b
    out_b[i+1] = rnorm(1, mean = (n * mean(z) + hyp_b2 * hyp_b1) / (n + hyp_b2),
                       sd = 1/ sqrt((n + hyp_b2) * out_lambda[i+1]) )
    #out_b[i+1] = out_b[i]
    
    # # MH per gamma: random walk
    oldgamma = out_gamma[i]
    newgamma = oldgamma + runif(1, -eps_gamma, eps_gamma)
    alpha = exp( logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                         A = out_A[i], gamma = newgamma, 
                         b = out_b[i+1], lambda = out_lambda[i+1],
                         hyp_A1 = hyp_A1, hyp_A2 = hyp_A2, 
                         hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2, 
                         hyp_lambda1 = hyp_lambda1, hyp_lambda2 = hyp_lambda2, 
                         hyp_b1 = hyp_b1, hyp_b2 = hyp_b2) -
                  logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                          A = out_A[i], gamma = oldgamma, 
                          b = out_b[i+1], lambda = out_lambda[i+1],
                          hyp_A1 = hyp_A1, hyp_A2 = hyp_A2, 
                          hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2, 
                          hyp_lambda1 = hyp_lambda1, hyp_lambda2 = hyp_lambda2, 
                          hyp_b1 = hyp_b1, hyp_b2 = hyp_b2) )
    if(runif(1) < alpha) oldgamma = newgamma
    out_gamma[i+1] = oldgamma
    
    # MH per A: uso RW normale
    oldA = out_A[i]
    newA = rnorm(1, oldA, eps_A)
    alpha = exp( logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                           A = newA, gamma = out_gamma[i+1], 
                           b = out_b[i+1], lambda = out_lambda[i+1],
                           hyp_A1 = hyp_A1, hyp_A2 = hyp_A2, 
                           hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2, 
                           hyp_lambda1 = hyp_lambda1, hyp_lambda2 = hyp_lambda2, 
                           hyp_b1 = hyp_b1, hyp_b2 = hyp_b2) -
                   logpost(y = y, cc = out_c[i+1,], c_0 = c0, s = out_s[i+1,],
                           A = oldA, gamma = out_gamma[i+1], 
                           b = out_b[i+1], lambda = out_lambda[i+1],
                           hyp_A1 = hyp_A1, hyp_A2 = hyp_A2, 
                           hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2, 
                           hyp_lambda1 = hyp_lambda1, hyp_lambda2 = hyp_lambda2, 
                           hyp_b1 = hyp_b1, hyp_b2 = hyp_b2) )
    if(runif(1) < alpha) oldA = newA
    out_A[i+1] = oldA
  }
  return(list(c = out_c, s = out_s, lambda = out_lambda, A = out_A, gamma = out_gamma, b = out_b))
}
## magari il kalman filter lo posso fare solo per un tot di iterazioni - burnin ?

nrep = 5000
start <- Sys.time()
prova <- gibbs_calcium(nrep = nrep, y = y, 
                       lambda_start = 10, b_start = 0,
                       gamma_start = 0.5, A_start = 5,
                       eps_gamma = 0.05, 
                       eps_A = 0.17)
end <- Sys.time()
end - start
str(prova)

burnin = 1:1500
plot(1:nrep, prova$lambda, type = "l", main = "lambda")
lines(1:nrep, cumsum(prova$lambda)/1:nrep, col = 2)
#abline(h=10, col = "blue")

plot(1:nrep, prova$b, type = "l", main = "b")
lines(1:nrep, cumsum(prova$b)/1:nrep, col = 2)
abline(h=2, col = "blue")

plot(1:nrep, prova$gamma, type = "l", main = "gamma")
lines(1:nrep, cumsum(prova$gamma)/1:nrep, col = 2)
abline(h=0.8, col = "blue")

plot(1:nrep, prova$A, type = "l", main = "A")
lines(1:nrep, cumsum(prova$A)/1:nrep, col = 2)
abline(h=5, col = "blue")

n = length(y)
plot(1:n, colMeans(prova$c[-burnin,]), type = "l")
lines(1:n, data$c, col = "blue", lty = 2)

plot(1:n, y, type = "l")
lines(1:n, mean(prova$b[-burnin]) + colMeans(prova$c[-burnin,]), col = 2)
abline(v = which(colMeans(prova$s[-burnin,])>0.5), col = "salmon", lty = 2)

which(colMeans(prova$s[-burnin,])>0.5)











