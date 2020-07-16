

library(truncnorm)

# data <- read.csv("data.csv", header = FALSE)
# str(data)
# head(data)
# plot(1:nrow(data), data$V1, type = "l")
# plot(1:length(data$V1[22000:28000]), data$V1[22000:28000], type = "l")
# y = data$V1[22000:28000]

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

# data <- sim_data(n = 1000, lambda = 8, time_spike = c(100, 101, 150, 300, 400,402, 500, 700, 720, 850), 
#                  gamma = 0.8, b = 2,
#                  prob = c(0.1, 0.7, 0.2), par = c(4, 10, 6))

data <- sim_data(n = 500, lambda = 10, time_spike = c(50,52, 140, 180, 250, 350, 420, 421, 460),
                 gamma = 0.8, b = 0,
                 prob = c(0.4, 0.6), par = c(4, 10))

y = data$y 
plot(y, type = "l", main = "")

data$k
data$A[data$A>0]


loglik <- function(y, cc, A, b, gamma, lambda)
{
  n = length(y)
  sum(dnorm(y, mean = b + gamma * cc + A , sd = 1/sqrt(lambda), log = TRUE))
}


logprior_gamma <- function(gamma, hyp_gamma1, hyp_gamma2) dbeta(gamma, hyp_gamma1, hyp_gamma2, log = TRUE)


logpost <- function(y, cc, s, A, b, gamma, lambda, 
                    hyp_gamma1, hyp_gamma2
                    )
{
  loglik(y = y, cc = cc, A = A, b = b, gamma = gamma, lambda = lambda) +
    logprior_gamma(gamma, hyp_gamma1, hyp_gamma2)
}



#---------------------# #---------------------# #---------------------# #---------------------# 
#---------------------# #---------------------# #---------------------# #---------------------# 

marg <- function(y, b, c, gamma, sigma2, psi2)
{
  2 * dnorm(y, mean = b + gamma * c, sd = sqrt(psi2 + sigma2)) * 
    (1 - pnorm(0, mean = psi2 / (sigma2 + psi2) * (y-b-gamma*c), sd = sqrt(sigma2 * psi2 / (sigma2 + psi2)) ) )
}

# 
# gamma_start = 0.8; b_start = 0
# tau2 = 0.0001
# lambda_start = 10
# c0 = 0; p = 0.995; alpha = 1
# hyp_A1 = 3; hyp_A2 = 1; hyp_gamma1 = 5; hyp_gamma2 = 2
# hyp_lambda1 = 10; hyp_lambda2 = 1; hyp_b1 = 1; hyp_b2 = 1
# psi2 = 2
# eps_gamma = 0.3
# C0 = 1
# 
# nrep=50
gibbs_calcium <- function(nrep, y,
                          # chain starting points
                          b_start = median(y),
                          gamma_start = 0.8,  lambda_start = 10,  
                          # (fized) variance of the underlying calcium
                          tau2 = 0.0001,
                          # calcium initial value and variance
                          c0 = 0, C0 = 1,
                          # prior probability of a spike
                          p = 0.995, 
                          # DP concentration parameter, variance of G0
                          alpha = 1, psi2 = 2,
                          # hyperparameters
                          hyp_A1 = 3, hyp_A2 = 1, hyp_gamma1 = 5, hyp_gamma2 = 2,
                          hyp_lambda1 = 10, hyp_lambda2 = 1, hyp_b1 = 1, hyp_b2 = 1,
                          # MH step
                          eps_gamma)
{
  n = length(y)
  # creo matrici di output
  out_c = matrix(NA, nrep, n+1)
  out_p1 = matrix(NA, nrep, n)
  out_p0 = matrix(NA, nrep, n)
  out_A = matrix(NA, nrep, n) # contiene i valori A1,...,Ak
  out_b = rep(NA, nrep)
  out_gamma = rep(NA, nrep) 
  out_lambda = rep(NA, nrep)
  
  filter_mean = numeric(n)
  filter_var = numeric(n)
  R = numeric(n)
  a = numeric(n)
  
  cluster = matrix(0, nrep, n) 
  
  # inizializzazione dei parametri
  out_c[1,] = 0
  out_c[,1] = c0
  out_A[,1] = 0
  out_b[1] = b_start
  out_gamma[1] = gamma_start
  out_lambda[1] = lambda_start
  
  cluster[1,] = 0 # vettore di lunghezza n con il cluster Zj

  for(i in 1:(nrep-1))
  {
    sigma2 = 1/out_lambda[i]

    # sampling di c
    R[1] = tau2 + out_gamma[i]^2 * C0
    filter_mean[1] = (sigma2 * out_A[i, cluster[i,1]+1] + R[1] * (y[1] - out_b[i])) / (sigma2 + R[1])
    filter_var[1] = sigma2 * R[1] /  (sigma2 + R[1])
    
    for(j in 2:n)
    {
      R[j] = out_gamma[i]^2 * filter_var[j-1] + tau2
      a[j] = out_gamma[i] * filter_mean[j-1] + out_A[i, cluster[i,j]+1]
      
      filter_mean[j] = a[j] + R[j]/(R[j] + sigma2) * (y[j] - out_b[i] - a[j])
      filter_var[j] = R[j] - R[j]^2/(R[j] + sigma2)
    }
    out_c[i+1, n+1] = rnorm(1, filter_mean[n], filter_var[n])
    
    for(j in (n-1):1)
    {
      back_mean = filter_mean[j] + out_gamma[i] * filter_var[j] / (R[j+1]) * (out_c[i+1,j+2] - a[j+1])
      back_var = filter_var[j] - out_gamma[i]^2 * filter_var[j]^2 / (R[j+1])

      out_c[i+1,j+1] = rnorm(1, back_mean, back_var)
    }
    
    back_mean = out_gamma[i] * C0/(C0 + tau2) * (out_c[i+1,2])
    back_var = C0 - out_gamma[i]^2 * C0^2 / (R[1])
    out_c[i+1, 1] = rnorm(1, back_mean, back_var)
    
    
    # sampling di lambda (precision)
    z = y - out_gamma[i] * out_c[i+1,1:n] - out_A[i, cluster[i,]+1] 
    ssum = sum((z - mean(z))^2)
    out_lambda[i+1] = rgamma(1, hyp_lambda1 + n/2, 
                             hyp_lambda2 + .5 * ssum + .5 * n * hyp_b2 / (n + hyp_b2) * (mean(z) - hyp_b1)^2 )

    # sampling di b
    out_b[i+1] = rnorm(1, mean = (n * mean(z) + hyp_b2 * hyp_b1) / (n + hyp_b2),
                       sd = 1/ sqrt((n + hyp_b2) * out_lambda[i+1]) )
    
    # # MH per gamma: random walk
    oldgamma = out_gamma[i] 
    newgamma = oldgamma + runif(1, -eps_gamma, eps_gamma)
    ratio = exp( logpost(y = y, cc = out_c[i+1,1:n], 
                         A = out_A[i, cluster[i,]+1], gamma = newgamma,
                         b = out_b[i+1], lambda = out_lambda[i+1],
                         hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2) -
                  logpost(y = y, cc = out_c[i+1,1:n],
                          A = out_A[i, cluster[i,]+1], gamma = oldgamma,
                          b = out_b[i+1], lambda = out_lambda[i+1],
                          hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2) )
    if(runif(1) < ratio) oldgamma = newgamma
    out_gamma[i+1] = oldgamma
    
    cluster[i+1,] = cluster[i,]
    out_A[i+1,] = out_A[i,]

    for(j in 1:n)
    {
      clus_tmp = cluster[i+1,]
      A_tmp = out_A[i+1,]
      clus_tmp[j] = NA

      if(length(unique(clus_tmp))>2)
      {
        clus_tmp[!is.na(clus_tmp) & clus_tmp>0] = as.numeric( factor( clus_tmp[!is.na(clus_tmp) & clus_tmp>0],
                                                                      labels = 1:(length(unique(clus_tmp[clus_tmp>0]))-1) ) )

        out_A[i+1, 2:n] = NA
        out_A[i+1, 2:(length(unique(clus_tmp[!is.na(clus_tmp) & clus_tmp>0]))+1)] = A_tmp[2:max(sort(unique(clus_tmp[!is.na(clus_tmp) & clus_tmp>0]))+1) ]

      }

      pr0 = p * dnorm(y[j], mean = out_b[i+1] + out_gamma[i+1] * out_c[i+1,j], sd = sqrt(sigma2 + tau2))

      pr_new = (1-p) * alpha * marg(y = y[j], b = out_b[i+1], c = out_c[i+1,j], gamma = out_gamma[i+1],
                                    sigma2 = (sigma2 + tau2), psi2 = psi2)

      nj = sapply(sort(unique(clus_tmp[!is.na(clus_tmp) & clus_tmp>0])), function(x) sum(clus_tmp[-j][clus_tmp[-j]>0] == x))
      pr_clus = 0
      if(length(nj)>0)
      {
        pr_clus = nj  * (1-p) * sapply(sort(unique(clus_tmp[-j][clus_tmp[-j]>0])),
                                          function(x) dnorm(y[j], mean = out_b[i+1] + out_gamma[i+1] * out_c[i+1,j] + out_A[i+1,x+1],
                                           sd = sqrt(sigma2 + tau2)) ) 
      }
      clus_tmp[j] = sample( -1:length(pr_clus), 1, prob = c(pr_new, pr0, pr_clus) )

      if(clus_tmp[j] == -1)
      {
        clus_tmp[j] = max(clus_tmp) +1
        out_A[i+1, clus_tmp[j]+1] = rtruncnorm(1, mean = psi2 / (psi2 + sigma2 + tau2) * (y[j] - out_b[i+1] - out_gamma[i+1] * out_c[i+1, j]),
                                     sd = sqrt((sigma2 + tau2) * psi2 / (psi2 + sigma2 + tau2) ))
      }

      cluster[i+1, ] = clus_tmp
    }
    
    # sampling di A
    out_A[i+1,] = NA
    out_A[i+1, 1] = 0
    nj = sapply(sort(unique(cluster[i+1, cluster[i+1,]>0])), function(x) sum(cluster[i+1, cluster[i+1,]>0] == x))
    sumj = sapply(sort(unique(cluster[i+1, cluster[i+1,]>0])), 
                  function(x) sum(y[cluster[i+1,] == x] - out_b[i+1] - out_gamma[i+1] * out_c[i+1,1:n][cluster[i+1,] == x]) )
    out_A[i+1, 2:(length(unique(cluster[i+1, cluster[i+1,]>0]))+1)] = rtruncnorm( length(unique(cluster[i+1, cluster[i+1,]>0])), 
                                                                           a = 0, b = Inf, mean = psi2 * sumj / (nj * psi2 + sigma2 + tau2), 
                                                                           sd = sqrt((sigma2 + tau2) * psi2 / (nj*psi2 + sigma2 + tau2) ) )
    
    
  }
  return(list(c = out_c, lambda = out_lambda, A = out_A, gamma = out_gamma, b = out_b, clus = cluster))
}




nrep = 1500
start <- Sys.time()
prova <- gibbs_calcium(nrep = nrep, y = y, 
                       C0 = 0.5,
                       p = 0.995,
                       lambda_start = 5, b_start = 1,
                       gamma_start = 0.5, 
                       eps_gamma = 0.03, 
                       alpha = 1,
                       psi2 = 2, tau2 = 0.00001)
end <- Sys.time()
end - start
str(prova)

n = length(y)
plot(1:nrep, prova$lambda, type = "l", main = "lambda")
lines(1:nrep, cumsum(prova$lambda)/1:nrep, col = 2)
abline(h=10, col = "blue")

plot(1:nrep, prova$b, type = "l", main = "b")
lines(1:nrep, cumsum(prova$b)/1:nrep, col = 2)
abline(h=0, col = "blue")

plot(1:nrep, prova$gamma, type = "l", main = "gamma")
lines(1:nrep, cumsum(prova$gamma)/1:nrep, col = 2)
abline(h=0.8, col = "blue")



burnin = 1:500
plot(1:n, y, type = "l")
AA = matrix(0,nrep,n)
for(i in 1:nrep)
{
  AA[i, prova$clus[i,] >0] = prova$A[i, prova$clus[i, prova$clus[i,] >0]+1]
}

calcium = matrix(0, nrep, n+1)
for(i in 1:nrep)
{
  for(j in 2:(n+1))
  {
    calcium[i,j] = calcium[i,j-1] * mean(prova$gamma[-burnin]) + AA[i,j-1]
  }
}

plot(1:n, y, type = "l")
lines(1:n, mean(prova$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)

abline(v = which( apply(prova$clus[-burnin,], 2, function(x) mean(x != 0))>0.95), lty = 3, col = "salmon")

hist(apply(prova$clus[-burnin,], 2, function(x) mean(x != 0)), main = "Distr. of spike probabilities")
which( apply(prova$clus[-burnin,], 2, function(x) mean(x != 0))>0.95 )
#(50,52, 140, 180, 250, 350, 420, 421, 460)

which(apply(prova$A, 2, function(x) sum(!is.na(x))>0))
prova$A = prova$A[,1:20]
  

prova$A[500:540,]
prova$A[is.na(prova$A)] = 0
maxA = apply(prova$A, 1, max)
mean(maxA)

max2A = t(apply(prova$A, 1, sort))
colMeans(max2A)


