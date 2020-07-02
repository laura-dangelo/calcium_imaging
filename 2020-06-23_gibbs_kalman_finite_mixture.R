

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
  A[time_spike] <- rnorm(sum(s), par[k], 0.8)
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i]
  }
  return(list("y" = b + c + rnorm(n, 0, 1/sqrt(lambda)), "c" = c, "s" = s, "A" = A, "k" = k))
}

data <- sim_data(n = 1000, lambda = 10, time_spike = c(100, 150, 400, 500, 700, 300, 600,850), 
                 gamma = 0.8, b = 2,
                 prob = c(0.1, 0.7, 0.2), par = c(0.6, 1, 2))

data <- sim_data(n = 500, lambda = 10, time_spike = c(50,140,180,250,350,420,460), 
                 gamma = 0.8, b = 0,
                 prob = c(0.4, 0.6), par = c(4, 10))

y = data$y 
plot(y, type = "l", main = "")

data$k
data$A[data$A>0]


loglik <- function(y, cc, s, A, b, gamma, lambda)
{
  n = length(y)
  sum(dnorm(y, mean = b + gamma * cc + A * s, sd = 1/sqrt(lambda), log = TRUE))
}


#logprior_A <- function(A, hyp_A1, hyp_A2) dgamma(A, hyp_A1, hyp_A2, log = TRUE)
# modificare!!!

logprior_gamma <- function(gamma, hyp_gamma1, hyp_gamma2) dbeta(gamma, hyp_gamma1, hyp_gamma2, log = TRUE)
# logprior_b_lambda <- function(b, lambda, hyp_b1, hyp_b2, hyp_lambda1, hyp_lambda2)
# {
#   dnorm(b, mean = hyp_b1, sd = 1 / sqrt(hyp_b2 * lambda), log = TRUE) + dgamma(lambda, hyp_lambda1, hyp_lambda2, log = TRUE)
# }

logpost <- function(y, cc, s, A, b, gamma, lambda, 
                    # hyp_A1, hyp_A2, 
                    hyp_gamma1, hyp_gamma2
                    # hyp_lambda1, hyp_lambda2, hyp_b1, hyp_b2
                    )
{
  loglik(y = y, cc = cc, s = s, A = A, b = b, gamma = gamma, lambda = lambda) + #logprior_A(A, hyp_A1, hyp_A2) + 
    logprior_gamma(gamma, hyp_gamma1, hyp_gamma2)# + 
   # logprior_b_lambda(b, lambda, hyp_b1, hyp_b2, hyp_lambda1, hyp_lambda2)
}



#---------------------# #---------------------# #---------------------# #---------------------# 
#---------------------# #---------------------# #---------------------# #---------------------# 

marg <- function(y, b, c, gamma, s, sigma2, psi2)
{
  2 * dnorm(y, mean = b + gamma * c, sd = sqrt(psi2 + sigma2)) * 
    (1 - pnorm(0, mean = psi2 / (sigma2 + psi2) * (y-b-gamma*c), sd = sqrt(sigma2 * psi2 / (sigma2 + psi2)) ) )
}


gamma_start = 0.8; b_start = 2
A_start = 5; tau2 = 0.0001
lambda_start = 10
c0 = 0; p = 0.005; alpha = 1
hyp_A1 = 3; hyp_A2 = 1; hyp_gamma1 = 5; hyp_gamma2 = 2
hyp_lambda1 = 10; hyp_lambda2 = 1; hyp_b1 = 1; hyp_b2 = 1
psi2 = 2

nrep=50
gibbs_calcium <- function(nrep, y,
                          gamma_start = 0.8,  b_start = 2, A_start = 5,
                          tau2 = 0.0001, lambda_start = 10, c0 = 0,
                          p = 0.005, 
                          alpha = 1, psi2 = 2,
                          hyp_A1 = 3, hyp_A2 = 1, hyp_gamma1 = 5, hyp_gamma2 = 2,
                          hyp_lambda1 = 10, hyp_lambda2 = 1, hyp_b1 = 1, hyp_b2 = 1,
                          eps_gamma, eps_A)
{
  n = length(y)
  # creo matrici di output
  out_c = matrix(NA, nrep, n+1) # calcio (c_0, c_1, ..., c_n)
  out_s = matrix(NA, nrep, n) # spikes (s_1, s_2, ..., s_n)
  out_p1 = matrix(NA, nrep, n) # probabilita' di osservare uno spike (pi_1, ..., pi_n)
  out_p0 = matrix(NA, nrep, n) # (1-pi_1, ..., 1-pi_n) -- a meno di costanti
  
  out_A = matrix(NA, nrep, n) # contiene i valori (A1,...,Ak)
  out_b = rep(NA, nrep) # posteriori di b
  out_gamma = rep(NA, nrep) # posteriori di gamma
  out_lambda = rep(NA, nrep) # posteriori di lambda
  
  filter_mean = numeric(n) # vettori di supporto per il kalman filter
  filter_var = numeric(n)
  
  cluster = matrix(0, nrep, n) # cluster per gli y_t tali che s_t = 1: (Z_1,...,Z_n)
  
  # inizializzazione dei parametri
  out_c[1,] = 0
  out_c[,1] = c0
  out_s[1,] = 0
  out_A[1,1] = A_start
  out_b[1] = b_start
  out_gamma[1] = gamma_start
  out_lambda[1] = lambda_start
  
  cluster[1,] = 1 # inizializzo tutti nello stesso cluster
  
  AA = rep(A_start, n) # contiene un vettore di lunghezza n che vale A_Zj anche quando s=0

  for(i in 1:(nrep-1))
  {
    sigma2 = 1/out_lambda[i]

    # sampling di c
    filter_mean[1] = (sigma2 * (out_b[i] + out_gamma[i] * out_c[i,1] + AA[1]*out_s[i,1] ) + tau2 * y[1]) / (tau2 + sigma2)
    filter_var[1] = tau2 * sigma2 / (tau2 + sigma2)
    
    for(j in 2:n)
    {
      filter_mean[j] = (sigma2 * (out_gamma[i] * filter_mean[j-1] + AA[j] * out_s[i,j] ) + 
                          (filter_var[j-1] + tau2) * y[j]) / (filter_var[j-1] + tau2 + sigma2)
      filter_var[j] = ((filter_var[j-1] + tau2) * sigma2) / (filter_var[j-1] + tau2 + sigma2)
    }
    out_c[i+1,n+1] = rnorm(1, filter_mean[n], filter_var[n])
    
    for(j in (n-1):1)
    {
      back_mean = filter_mean[j] + out_gamma[i] * filter_var[j]/(filter_var[j] + tau2) * 
        (out_c[i+1,j+2] - out_gamma[i] * filter_mean[j] - AA[j] * out_s[i,j])
      back_var = filter_var[j] - out_gamma[i]^2 * filter_var[j]^2 / (filter_var[j] + tau2)
      out_c[i+1,j+1] = rnorm(1, back_mean, back_var)
    }
    
    # sampling di lambda (precision)
    z = y - out_gamma[i] * out_c[i+1,1:n] - AA * out_s[i,]
    ssum = sum((z - mean(z))^2)
    out_lambda[i+1] = rgamma(1, hyp_lambda1 + n/2, 
                             hyp_lambda2 + .5 * ssum + .5 * n * hyp_b2 / (n + hyp_b2) * (mean(z) - hyp_b1)^2 )

    # sampling di b
    out_b[i+1] = rnorm(1, mean = (n * mean(z) + hyp_b2 * hyp_b1) / (n + hyp_b2),
                       sd = 1/ sqrt((n + hyp_b2) * out_lambda[i+1]) )
    
    # # MH per gamma: random walk
    oldgamma = out_gamma[i] 
    newgamma = oldgamma + runif(1, -eps_gamma, eps_gamma)
    ratio = exp( logpost(y = y, cc = out_c[i+1,1:n], s = out_s[i,],
                         A = AA, gamma = newgamma,
                         b = out_b[i+1], lambda = out_lambda[i+1],
                         hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2) -
                  logpost(y = y, cc = out_c[i+1,1:n], s = out_s[i,],
                          A = AA, gamma = oldgamma,
                          b = out_b[i+1], lambda = out_lambda[i+1],
                          hyp_gamma1 = hyp_gamma1, hyp_gamma2 = hyp_gamma2) )
    if(runif(1) < ratio) oldgamma = newgamma
    out_gamma[i+1] = oldgamma
    
    # sampling di s
    out_p1[i+1,] = exp(-.5 / sigma2 * (y - out_b[i] - out_gamma[i] * out_c[i,1:n] - AA)^2 ) * p
    out_p0[i+1,] = exp(-.5 / sigma2 * (y - out_b[i] - out_gamma[i] * out_c[i,1:n])^2 ) * (1-p)
    out_s[i+1,] = apply(cbind(out_p0[i+1,], out_p1[i+1,]), 1, function(x) sample(c(0,1), 1, prob = x))
    
    cluster[i+1,] = cluster[i,]
    cluster[i+1,][out_s[i+1,] == 0] = 0
    out_A[i+1,] = out_A[i,]
    
    # "Polya Urn" finite mixture
    for(j in 1:sum(out_s[i+1,]>0))
    {
      # sampling del cluster
      clus_tmp = cluster[i+1,out_s[i+1,]>0]
      clus_tmp[j] = NA

      out_A[i+1, 1:(length(unique(clus_tmp[!is.na(clus_tmp) & clus_tmp>0])))] = out_A[i+1, sort(unique(clus_tmp[!is.na(clus_tmp) & clus_tmp>0])) ]
      
      
      clus_tmp[!is.na(clus_tmp) & clus_tmp>0] = as.numeric( factor( clus_tmp[!is.na(clus_tmp) & clus_tmp>0], 
                                                                    labels = 1:(length(unique(clus_tmp[clus_tmp>0]))-1) ) )
      
      nj = sapply(sort(unique(clus_tmp[!is.na(clus_tmp) & clus_tmp>0])), function(x) sum(clus_tmp[-j][clus_tmp[-j]>0] == x))
      prob_c = sapply(sort(unique(clus_tmp[-j][clus_tmp[-j]>0])), 
                      function(x) dnorm(y[out_s[i+1,]>0][j], mean = out_b[i+1] + out_gamma[i+1] * out_c[i+1,j] + out_A[i+1,x], 
                                        sd = sqrt(sigma2)) ) * nj / (n-1+alpha)
    
      prob_new = alpha / (n-1+alpha) * marg(y = y[out_s[i+1,]>0][j], b = out_b[i+1], c = out_c[i+1,j], gamma = out_gamma[i+1], 
                                      s = out_s[i+1,j], sigma2 = sigma2, psi2 = 1)
      
      clus_tmp[j] = sample( 1:(max(unique(clus_tmp[-j][clus_tmp[-j]>0]))+1), 1, prob = c(prob_c, prob_new) )
      
      cluster[i+1, out_s[i+1,]>0] = clus_tmp
      if(clus_tmp[j] == max(unique(clus_tmp[-j][clus_tmp[-j]>0]))+1) 
        {out_A[i+1,clus_tmp[j]] = rtruncnorm( 1,
                                                a = 0, b = Inf, 
                                                mean = psi2 * (y[out_s[i+1,]>0][j] - out_b[i+1] - out_gamma[i+1] * out_c[i+1, j]) / (psi2 + sigma2),
                                                sd = sqrt(sigma2 * psi2 / (psi2 + sigma2) ) ) }
      
    }
    
    # sampling di A
    out_A[i+1,] = NA
    nj = sapply(sort(unique(cluster[i+1,out_s[i+1,]>0])), function(x) sum(cluster[i+1,out_s[i+1,]>0] == x))
    sumj = sapply(sort(unique(cluster[i+1,out_s[i+1,]>0])), 
                  function(x) sum(y[cluster[i+1,] == x] - out_b[i+1] - out_gamma[i+1] * out_c[i+1, cluster[i+1,] == x]) )
    out_A[i+1, 1:length(unique(cluster[i+1,out_s[i+1,]>0]))] = rtruncnorm( length(unique(cluster[i+1,out_s[i+1,]>0])), 
                                                                           a = 0, b = Inf, mean = psi2 * sumj / (nj * psi2 + sigma2), 
                                                                           sd = sqrt(sigma2 * psi2 / (nj*psi2 + sigma2) ) )

    AA[out_s[i+1,]==1] = out_A[i+1, cluster[i+1, out_s[i+1,]==1]]
    for(j in which(out_s[i+1,]==0))
    {
      clus_tmp = c(NA, cluster[i+1, out_s[i+1,]==1])

      nj = sapply(sort(unique(clus_tmp[-1])), function(x) sum(clus_tmp[-1] == x))
      prob_c = sapply(sort(unique(clus_tmp[-1])),
                      function(x) dnorm(y[j], mean = out_b[i+1] + out_gamma[i+1] * out_c[i+1,j] + out_A[i+1,x], 
                                        sd = sqrt(sigma2)) ) * nj / (n-1+alpha)

      # prob_new = alpha / (n-1+alpha) * marg(y = y[j], b = out_b[i+1], c = out_c[i+1,j], gamma = out_gamma[i+1],
      #                                 s = out_s[i+1,j], sigma2 = sigma2, psi2 = 1)

      clus_tmp[1] = sample(1:(max(unique(clus_tmp[-1][clus_tmp[-1]>0]))), 1, prob = c(prob_c) )
      AA[j] = out_A[i+1, clus_tmp[1]]
      
      # sampling di A
      # if(clus_tmp[1] == max(unique(clus_tmp[-1])+1) )
      # AA[j] = rtruncnorm( 1, a = 0, b = Inf,
      #                      mean = psi2 * (y[j] - out_b[i+1] - out_gamma[i+1] * out_c[i+1, j]) / (psi2 + sigma2),
      #                       sd = sqrt(sigma2 * psi2 / (psi2 + sigma2) ) ) 
    }
    
    
    
  }
  return(list(c = out_c, s = out_s, lambda = out_lambda, A = out_A, gamma = out_gamma, b = out_b, clus = cluster))
}




nrep = 2000
start <- Sys.time()
prova <- gibbs_calcium(nrep = nrep, y = y, 
                       lambda_start = 10, b_start = 0,
                       gamma_start = 0.8, A_start = 3,
                       eps_gamma = 0.05, 
                       eps_A = 0.17)
end <- Sys.time()
end - start
str(prova)

n = length(y)
burnin = 1:500
plot(1:nrep, prova$lambda, type = "l", main = "lambda")
lines(1:nrep, cumsum(prova$lambda)/1:nrep, col = 2)
abline(h=10, col = "blue")

plot(1:nrep, prova$b, type = "l", main = "b")
lines(1:nrep, cumsum(prova$b)/1:nrep, col = 2)
abline(h=0, col = "blue")

plot(1:nrep, prova$gamma, type = "l", main = "gamma")
lines(1:nrep, cumsum(prova$gamma)/1:nrep, col = 2)
abline(h=0.8, col = "blue")


apply(prova$A, 2, function(x) sum(!is.na(x)))
plot(1:nrep, prova$A[,1], type = "l", main = "A")
lines(1:nrep, cumsum(prova$A[,1])/1:nrep, col = 2)

plot(1:nrep, prova$A[,2], type = "l", main = "A")
lines(1:nrep, cumsum(prova$A[,2])/1:nrep, col = 2)

plot(1:n, y, type = "l")
AA = matrix(0,nrep,n)
AA[prova$clus>0] = prova$A[prova$clus[prova$clus>0]]
lines(1:n, mean(prova$gamma[-burnin]) * colMeans(prova$c[-burnin,1:n]) + mean(prova$b[-burnin]) + 
        (colMeans(prova$s[-burnin,])>0.6) * colMeans(AA[-burnin,]),
     col = 2)

str(AA)

n = length(y)
plot(1:n, y, type = "l")
lines(1:n, colMeans(prova$c[-burnin,1:n]) + mean(prova$b[-burnin]), col = "blue", lty = 2)
abline(v = which(colMeans(prova$s[-burnin,])>0.6), lty = 3, col = 2)


which(colMeans(prova$s[-burnin,])>0.6)
#c(50,140,180,250,350,420,460)










