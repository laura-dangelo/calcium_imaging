library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)

# sourceCpp('./SourceCPP/calcium_CAM.cpp')
sourceCpp('./SourceCPP/calcium_CAM_gMFM.cpp')

## mixture ##
sim_data <- function(n, sigma2, tau2, time_spike, b, gamma, prob, par)
{
  c = rep(0,n)
  s = rep(0,n)
  A = rep(0,n)
  s[time_spike] = 1
  
  k <- sample(1:length(prob), sum(s), prob, replace = TRUE)
  A[time_spike] <- par[k]
  
  for(i in 2:n)
  {
    c[i] = gamma * c[i-1] + A[i] * s[i] + rnorm(1, 0, sqrt(tau2))
  }
  return(list("y" = b + c + rnorm(n, 0, sqrt(sigma2)), "c" = c, "s" = s, "A" = A, "k" = k))
}


times_spike <- function(n, ns, p1, p2)
{
  times = rbinom(n, 1, p1)
  for(i in 1:n)
  {
    if(times[i] == 1)
    {
      times[(i+1):(i+ns)] = rbinom(ns, 1, p2)
    }
  }
  times 
}


set.seed(1234)
sigma2 = 0.004
tau2 = 0.00001
n1 = 4000
n2 = 2500
n3 = 5000
n4 = 3500
n5 = 2000
gamma = 0.6
b = 0


spp = c( times_spike(n1, 9, 0.007, 0.09)[1:n1], # clD1
         times_spike(n2, 6, 0.01, 0.08)[1:n2], ## clD2
         times_spike(n3, 5, 0.01, 0.1)[1:n3], ### cld3
         times_spike(n4, 9, 0.007, 0.09)[1:n4], # clD1
         times_spike(n5, 5, 0.01, 0.1)[1:n5]) ### clD3
spp = which(spp == 1)
length(spp)


prob1 = rep(0.2, 5)
par1 = c(0.5, 0.9, 1.5, 1.9, 2.5)

prob2 = rep(0.25, 4)
par2 = c(0.5, 0.9, 1.5, 2.5)

prob3 = c(0.33, 0.33, 0.34)
par3 = c(0.5, 1.5, 1.9)

group1 <- sim_data(n = n1, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp<=n1],
                 gamma = gamma, b = b,
                 prob = prob1, par = par1)

group2 <- sim_data(n = n2, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>n1)&(spp<=(n1 + n2))] - n1,
                   gamma = gamma, b = b,
                   prob = prob2, par = par2)

group3 <- sim_data(n = n3, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2))&(spp<=(n1 + n2 + n3))] - (n1 + n2),
                   gamma = gamma, b = b,
                   prob = prob3, par = par3)

group4 <- sim_data(n = n4, sigma2 = sigma2, tau2 = tau2, time_spike = spp[(spp>(n1 + n2 + n3))&(spp<=(n1 + n2 + n3 + n4))] - (n1 + n2 + n3),
                   gamma = gamma, b = b,
                   prob = prob1, par = par1)

group5 <- sim_data(n = n5, sigma2 = sigma2, tau2 = tau2, time_spike = spp[spp>(n1 + n2 + n3 + n4)] - (n1 + n2 + n3 + n4),
                   gamma = gamma, b = b,
                   prob = prob3, par = par3)



y = c(group1$y, group2$y, group3$y, group4$y, group5$y)
g = c(rep(1,n1), rep(2,n2), rep(3,n3), rep(4,n4), rep(5,n5))
plot(y, type = "l")

A_start = rep(0,100)
# A_start[2:6] = c(0.5, 0.9, 1.5, 1.9, 2.5)
# AA = c(group1$A, group2$A, group3$A, group4$A)

cluster = rep(0, length(y))
# cluster[AA == A_start[2]] = 1
# cluster[AA == A_start[3]] = 2
# cluster[AA == A_start[4]] = 3
# cluster[AA == A_start[5]] = 4
# cluster[AA == A_start[6]] = 5


sum(length(spp))/length(y)

n = length(y)
J = length(unique(g))
nrep = 2500
set.seed(1234)

run = calcium_gibbs(Nrep = nrep, 
                    y = y, g = g,
                    cal = c(0,y),
                    clO = cluster, clD = c(1,2,3,4,5), 
                    A_start = A_start,
                    b_start = 0,
                    gamma_start = 0.6, 
                    sigma2_start = 0.004, 
                    tau2_start = 0.00001, 
                    p_start = 0.001, 
                    alpha = 3, beta_start = 1,
                    maxL_start = 20,
                    max_xiK = 100, 
                    kappa_D = 0.5, 
                    c0 = 0, varC0 = 0.1, 
                    hyp_A1 = 10, hyp_A2 = 10, 
                    hyp_b1 = 0, hyp_b2 = 1, 
                    hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                    hyp_tau21 = 1000, hyp_tau22 = 1, 
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    hyp_beta1 = 1, hyp_beta2 = 1,
                    hyp_maxL = 10,
                    eps_beta = 1,
                    eps_gamma = 0.007,
                    eps_A = 0.002,
                    eps_maxL = 10)

  


#plot(function(x) dgamma(x, 10, 7), 0, 3)
# burnin = 1:500
# run$calcium = run$calcium[,-burnin]
# run$cluster = run$cluster[,-burnin]
# run$b = run$b[-burnin]
# run$gamma = run$gamma[-burnin]
# run$sigma2 = run$sigma2[-burnin]
# run$tau2 = run$tau2[-burnin]
# run$A = run$A[,-burnin]
# run$p = run$p[-burnin]
# save(run, file = "res_sim_1709_par4_low.Rdata")

burnin = 1:2000
plot(1:length(run$p[-burnin]), run$p[-burnin], type = "l")
lines(1:length(run$p[-burnin]), cumsum(run$p[-burnin])/1:length(run$p[-burnin]), col =2)

plot(1:length(run$sigma2[-burnin]), run$sigma2[-burnin], type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(run$sigma2[-burnin]), cumsum(run$sigma2[-burnin])/1:length(run$sigma2[-burnin]), col =2)

plot(1:length(run$tau[-burnin]), run$tau[-burnin], type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(run$tau[-burnin]), cumsum(run$tau[-burnin])/1:length(run$tau[-burnin]), col =2)

plot(1:length(run$b[-burnin]), run$b[-burnin], type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run$b[-burnin]), cumsum(run$b[-burnin])/1:length(run$b[-burnin]), col =2)

plot(1:length(run$gamma[-burnin]), run$gamma[-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run$gamma[-burnin]), cumsum(run$gamma[-burnin])/1:length(run$gamma[-burnin]), col =2)

plot(1:length(run$beta[-burnin]), run$beta[-burnin], type = "l", xlab = "iterazioni", ylab = "beta")
lines(1:length(run$beta[-burnin]), cumsum(run$beta[-burnin])/1:length(run$beta[-burnin]), col =2)

plot(1:length(run$maxL[-burnin]), run$maxL[-burnin], type = "l", xlab = "iterazioni", ylab = "maxL")
lines(1:length(run$maxL[-burnin]), cumsum(run$maxL[-burnin])/1:length(run$maxL[-burnin]), col =2)



mean(run$sigma2[-burnin])
mean(run$tau2[-burnin])
mean(run$b[-burnin])
mean(run$gamma[-burnin])




#apply(run$A, 2, function(x) length(which(x>0)) ) # numero di clusters sulle osservazioni
barplot(table(apply(run$A[,-burnin], 2, function(x) length(which(x>0))+1 )))

burnin = 1:2000
AA = matrix(0,length(run$b[-burnin]),n)
for(i in 1:length(run$b[-burnin]))
{
  ii = i + max(burnin)
  AA[i, t(run$clusterO)[ii,] >0] = run$A[run$clusterO[run$clusterO[,ii] >0,ii]+1,ii]
}
est_spikes = colMeans(AA) 
est_spikes[which( apply(t(run$clusterO)[-burnin,], 2, function(x) mean(x != 0))<0.8)] = 0
times = which(est_spikes>0)

times
spp

#plot(1:length(AA[129,-burnin]), AA[129,-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")

sum(sapply(spp, function(x) !(x %in% times))) / length(spp)  ### spikes non identificati: falsi negativi
sum(sapply(times, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi

AA[,which(est_spikes == 0)] = 0
barplot(table( apply(AA, 1, function(x) length(unique(x))) ))

A_ind = AA[apply(AA, 1, function(x) length(unique( x )))==6,]
dataa = data.frame(A = A_ind[A_ind>0])
ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(bins = 35, aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() 
 # scale_x_continuous(limits = c(0.1,2.7))



barplot(table(apply(run$clusterD, 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni


# analizzo il caso = 3
mat_clusterD = matrix(NA, J, J)
ind3 = which( apply(run$clusterD, 2, function(x) length(unique(x)) ) ==3 ) # <-- qui metti quanti cluster vuoi vedere
mat_heatmap = expand.grid(J1 = unique(g),
                          J2 = unique(g))
for(i in 1:J)
{
  for(j in 1:i)
  {
    mat_clusterD[i,j] = sum( apply(run$clusterD[,ind3], 2, function(x) x[i] == x[j] ) )
    mat_clusterD[j,i] = mat_clusterD[i,j] 
  }
}
df_heat = data.frame(J1 = as.factor(mat_heatmap[,1]),
                     J2 = as.factor(mat_heatmap[,2]),
                     val = c(mat_clusterD)/mat_clusterD[1,1],
                     lab = round( c(mat_clusterD)/mat_clusterD[1,1] , 3) )
df_heat$lab[df_heat$J1==df_heat$J2] = 1:J
df_heat$val[df_heat$J1==df_heat$J2] = NA

ggplot(data = df_heat) +
  geom_tile( aes(x = J1, y = J2, fill = val)) +
  geom_text(aes(x = J1, y = J2, label = lab), size=2.5) +
  scale_fill_gradient(low = magma(3)[3], high = inferno(3)[2]) +
  ylim(rev(levels(df_heat$J2))) 




