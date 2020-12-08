library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)



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


#----------# parametri
sigma2 = 0.003
tau2 = 0.0001
n = 20000
gamma = 0.7
b = 0

# prob di spike "distanti"
pp = 0.012

# per quanti istanti successivi ho spikes
m = 7


# probabilità di spike per m istanti successivi a uno spike
pm = 0.17


set.seed(1234)
spp = times_spike(n, m, pp, pm)


spp = which(spp == 1)
length(spp)

prob = c(0.2, 0.2, 0.15, 0.15, 0.15, 0.15)
par = c(0.3, 0.50, 0.7, 0.9, 1.1, 1.5)

length(par)

data <- sim_data(n = n, sigma2 = sigma2, tau2 = tau2, time_spike = spp,
                   gamma = gamma, b = b,
                   prob = prob, par = par)

y = data$y
plot(y, type = "l")

clus = kmeans(diff(y,1)[(diff(y,1))>0.3], centers = 8)
A_start = rep(0,50)
A_start[2:9] = c(clus$centers)
cluster = numeric(length(y))
cluster[ (diff(y,1))>0.3] = clus$cluster

save(y, file="y_inner201120.Rdata")
save(data, file="data_inner201120.Rdata")

save(A_start, file="A_start_inner201120.Rdata")
save(cluster, file="cluster_inner201120.Rdata")

n = length(y)

nrep = 2000
burnin = 1:500
set.seed(1234)
sourceCpp('calcium_DP_innermixture.cpp')
run = calcium_gibbs(Nrep = nrep, y = y, 
                            cal = c(0,y),
                            cl = cluster, 
                            A_start = A_start,
                            b_start = 0,
                            gamma_start = 0.6, 
                            sigma2_start = 0.002, 
                            tau2_start = 0.001, 
                            p_start = 0.01, 
                            c0 = 0, varC0 = 0.1, 
                            alpha = 1, m = 1,
                            hyp_A1 = 4, hyp_A2 = 4, 
                            hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                            hyp_tau21 = 1000, hyp_tau22 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 1,
                            hyp_p1 = 1, hyp_p2 = 999,
                            eps_gamma = 0.01,
                            eps_A = 0.002)


# burnin = 1:500
# run$calcium = run$calcium[,-burnin]
# run$cluster = run$cluster[,-burnin]
# run$b = run$b[-burnin]
# run$gamma = run$gamma[-burnin]
# run$sigma2 = run$sigma2[-burnin]
# run$tau2 = run$tau2[-burnin]
# run$A = run$A[,-burnin]
# run$p = run$p[-burnin]
# save(run, file = "res_inner_201120.Rdata")

burnin = 1:500
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



mean(run$sigma2[-burnin])
mean(run$tau2[-burnin])
mean(run$b[-burnin])
mean(run$gamma[-burnin])



library(ggplot2)
int = 1:n
ggplot(data = data.frame(x = rep(int,2), 
                         y = c( y[int], rowMeans(run$calcium[,-burnin])[-1][int] ), 
                         col = as.factor( c(rep(1, length(int)), rep(2, length(int)) ))), 
       aes(x = x, y = y ) ) + 
  geom_line(aes(color = col)) +
  theme_bw() +
  scale_x_continuous(name = "t") +
  scale_y_continuous(name = expression(y[t])) +
  scale_color_manual(values = c("#000000", "#00BFC4", "#F8766D"),
                     labels = c("Observed trace", "Estimate"))+
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(t=-0.09, r=0, b=0, l=0, unit="cm"),
        legend.text=element_text(size=12)) 


burnin = 1:1700
AA = matrix(0, length(run$b[-burnin]), n)
for(i in 1:length(run$b[-burnin]))
{
  ii = i 
  AA[i, t(run$clus)[ii,] >0] = run$A[run$clus[run$clus[,ii] >0,ii]+1,ii]
}

est_spikes = colMeans(AA) ## run
est_spikes[which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))<0.8)] = 0
times = which(est_spikes>0)

times
spp

sum( est_spikes>0 )
length(spp)

sum(sapply(spp, function(x) !(x %in% times))) / length(spp)  ### spikes non identificati falsi negativi
sum(sapply(times, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi


