library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)
# install.packages("RcppProgress")
sourceCpp('calcium_DP_innermixture.cpp')
sourceCpp('calcium_DP_mixture_nonconjugate_2708.cpp')

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


# set.seed(123)
# data <- sim_data(n = 10000, sigma2 = 0.002, tau2 = 0.0005, time_spike = c(380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500),
#                  gamma = 0.9, b = 0,
#                  prob = c(0.23, 0.44, 0.33), par = c(2.6, 1.1, 0.6))
# 
# 
# y = data$y
# plot(y, type = "l")
# 
# 
# clus = data$s
# clus[clus>0] = data$k
# data$k
# c(2.6, 1, 0.5)[data$k]
# 
# 1-length(data$k)/length(data$y)
# 
# 
# A_start = rep(0,50)
# #A_start[2:4] = c(2.6, 1.1, 0.6)
# n = length(y)
# 
# plot(function(x) dgamma(x, 10,5),0,4)
# 
# nrep = 500
# set.seed(123)
# run = calcium_gibbs(Nrep = nrep, y = y, 
#                     cal = c(0,y),
#                     cl = rep(0,n), 
#                     A_start = A_start,
#                     b_start = 0,
#                     gamma_start = 0.8, 
#                     sigma2_start = 0.002, 
#                     tau2_start = 0.001, 
#                     p_start = 0.99, 
#                     c0 = 0, varC0 = 0.1, 
#                     alpha = 1, m = 5,
#                     hyp_A1 = 10, hyp_A2 = 5, 
#                     kappa = 0,
#                     hyp_b1 = 0, hyp_b2 = 1, 
#                     hyp_sigma21 = 1000, hyp_sigma22 = 1, 
#                     hyp_tau21 = 1000, hyp_tau22 = 1, 
#                     hyp_gamma1 = 1, hyp_gamma2 = 1,
#                     hyp_p1 = 99, hyp_p2 = 1,
#                     eps_gamma = 0.005,
#                     eps_A = 0.002)
# 
# 
# str(run) 
# burnin = 1:300
# which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# # (380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500)
# 
# plot(1:length(run$p), run$p, type = "l")
# lines(1:length(run$p), cumsum(run$p)/1:length(run$p), col =2)
# 
# plot(1:length(run$sigma2), run$sigma2, type = "l", xlab = "iterazioni", ylab = "sigma2")
# lines(1:length(run$sigma2), cumsum(run$sigma2)/1:length(run$sigma2), col =2)
# 
# plot(1:length(run$tau), run$tau, type = "l", xlab = "iterazioni", ylab = "tau2")
# lines(1:length(run$tau), cumsum(run$tau)/1:length(run$tau), col =2)
# 
# plot(1:length(run$b), run$b, type = "l", xlab = "iterazioni", ylab = "b")
# lines(1:length(run$b), cumsum(run$b)/1:length(run$b), col =2)
# 
# plot(1:length(run$gamma), run$gamma, type = "l", xlab = "iterazioni", ylab = "gamma")
# lines(1:length(run$gamma), cumsum(run$gamma)/1:length(run$gamma), col =2)
# 
# 
# 
# 
# iter = 1:length(run$gamma)
# obs = 1:500
# image(1:(length(iter)), 1:(length(obs)), t(run$cluster[obs,iter]), 
#       axes = F,
#       xlab = "iterazioni", ylab = "osservazione",
#       col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
# axis(1, at = seq(1, max(iter), by = 10))
# axis(2, at = seq(1,length(obs), by = 1), labels = obs)
# 
# obs = 501:1000
# image(1:(length(iter)), 1:(length(obs)), t(run$cluster[obs,iter]), 
#       axes = F,
#       xlab = "iterazioni", ylab = "osservazione",
#       col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = F)))
# axis(1, at = seq(1, max(iter), by = 10))
# axis(2, at = seq(1,length(obs), by = 1), labels = obs)
# 
# 
# burnin = 1:300
# 
# mean(run$sigma2[-burnin])
# mean(run$tau2[-burnin])
# mean(run$b[-burnin])
# mean(run$gamma[-burnin])
# 
# plot(0:n, c(0,data$c), type = "l")
# lines(0:n, rowMeans(run$calcium[,-burnin]), col = "turquoise")
# 
# 
# 
# plot(1:n, y, type = "l")
# AA = matrix(0,nrep,n)
# for(i in 1:nrep)
# {
#   AA[i, t(run$clus)[i,] >0] = run$A[run$clus[run$clus[,i] >0,i]+1,i]
# }
# 
# calcium = matrix(0, nrep, n+1)
# for(i in 1:nrep)
# {
#   for(j in 2:(n+1))
#   {
#     calcium[i,j] = calcium[i, j-1] * mean(run$gamma[-burnin]) + AA[i,j-1]
#   }
# }
# 
# plot(1:n, y, type = "l")
# lines(1:n, mean(run$b[-burnin]) + colMeans(calcium)[2:(n+1)], col = "turquoise3", lwd = 1.5)
# 
# abline(v = which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8), lty = 3, col = "salmon")
# which( apply(t(run$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# # (380, 1000, 1002, 1300, 2000, 2990, 4000, 4700, 4701, 5500)
# 
# hist(apply(run$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")
# 
# 
# ### number of clusters
# # head(t(run$cluster[375:383,-burnin]))
# # barplot( table(unlist(apply(run$cluster[,-burnin], 2, function(x) unique(x[x>0])))), main="number of cluster")
# plot(1:nrep, apply(run$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)
# 
# 
# ### cluster parameter
# minA = min( which(apply(run$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
# minA -1
# run$A = run$A[1:(minA +1),]
# str(run$A)
# 
# out_A = t(run$A)
# out_A[400:500,1:(minA)]







spp = c(380,385,
        985,
        1000, 1001, 1003, 1004, 1005, 1007, 
        1012,
        1021, 1022, 1024,
        1100, 
        2800, 
        2990, 2991, 2993, 2995, 2999, 3001, 3002,
        3101,3104,
        4000, 4004,
        4700, 4701, 4703,
        5500,
        6000,6003,
        6233, 
        6250,
        7100,7120,
        8700,
        8800, 8802, 8803, 8804)

set.seed(1234)
data2 <- sim_data(n = 10000, sigma2 = 0.004, tau2 = 0.00002, time_spike = spp,
                 gamma = 0.6, b = 0,
                 prob = c(0.43, 0.3, 0.13), par = c(0.3, 0.41, 0.65))

y2 = data2$y
plot(y2, type = "l")

A_start = rep(0,50)
A_start[2:4] = c(0.3, 0.45, 0.65)
n = length(y2)
cluster = rep(0,length(y2))
cluster[data2$s == 1] = data2$k

nrep = 900
set.seed(1234)
# run7 = calcium_gibbs(Nrep = nrep, y = y2,
#                      cal = c(0,y2),
#                      cl = cluster,
#                      A_start = A_start,
#                      b_start = 0,
#                      gamma_start = 0.9,
#                      sigma2_start = 0.002,
#                      tau2_start = 0.001,
#                      p_start = 0.99,
#                      c0 = 0, varC0 = 0.1,
#                      alpha = 1, m = 5,
#                      hyp_A1 = 10, hyp_A2 = 10,
#                      hyp_b1 = 0, hyp_b2 = 1,
#                      hyp_sigma21 = 1000, hyp_sigma22 = 1,
#                      hyp_tau21 = 1000, hyp_tau22 = 1,
#                      hyp_gamma1 = 1, hyp_gamma2 = 1,
#                      hyp_p1 = 99, hyp_p2 = 1,
#                      eps_gamma = 0.005,
#                      eps_A = 0.002)

run7 = calcium_gibbs(Nrep = nrep, y = y2, 
                            cal = c(0,y2),
                            cl = cluster, 
                            A_start = A_start,
                            b_start = 0,
                            gamma_start = 0.9, 
                            sigma2_start = 0.002, 
                            tau2_start = 0.001, 
                            p_start = 0.01, 
                            c0 = 0, varC0 = 0.1, 
                            alpha = 1, m = 1,
                            hyp_A1 = 10, hyp_A2 = 10, 
                            hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                            hyp_tau21 = 1000, hyp_tau22 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 1,
                            hyp_p1 = 1, hyp_p2 = 99,
                            eps_gamma = 0.01,
                            eps_A = 0.002)


# burnin = 1:500
# run7$calcium = run7$calcium[,-burnin]
# run7$cluster = run7$cluster[,-burnin]
# run7$b = run7$b[-burnin]
# run7$gamma = run7$gamma[-burnin]
# run7$sigma2 = run7$sigma2[-burnin]
# run7$tau2 = run7$tau2[-burnin]
# run7$A = run7$A[,-burnin]
# run7$p = run7$p[-burnin]
# save(run7, file = "res_sim_1709_par4_low.Rdata")

burnin = 1
plot(1:length(run7$p[-burnin]), run7$p[-burnin], type = "l")
lines(1:length(run7$p[-burnin]), cumsum(run7$p[-burnin])/1:length(run7$p[-burnin]), col =2)

plot(1:length(run7$sigma2[-burnin]), run7$sigma2[-burnin], type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(run7$sigma2[-burnin]), cumsum(run7$sigma2[-burnin])/1:length(run7$sigma2[-burnin]), col =2)

plot(1:length(run7$tau[-burnin]), run7$tau[-burnin], type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(run7$tau[-burnin]), cumsum(run7$tau[-burnin])/1:length(run7$tau[-burnin]), col =2)

plot(1:length(run7$b[-burnin]), run7$b[-burnin], type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run7$b[-burnin]), cumsum(run7$b[-burnin])/1:length(run7$b[-burnin]), col =2)

plot(1:length(run7$gamma[-burnin]), run7$gamma[-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run7$gamma[-burnin]), cumsum(run7$gamma[-burnin])/1:length(run7$gamma[-burnin]), col =2)



mean(run7$sigma2[-burnin])
mean(run7$tau2[-burnin])
mean(run7$b[-burnin])
mean(run7$gamma[-burnin])



AA = matrix(0,(length(run7$b)-max(burnin)),n)
for(i in 1:(length(run7$b)-max(burnin)) )
{
  ii = i + (max(burnin))
  AA[i, t(run7$clus)[ii,] >0] = run7$A[run7$clus[run7$clus[,ii] >0,ii]+1,ii]
}

est_spikes = colMeans(AA)
est_spikes[which( apply(t(run7$clus)[-burnin,], 2, function(x) mean(x != 0))<0.8)] = 0

int = 1:2000
plot(int, y2[int], type = "l")
lines(int, rowMeans(run7$calcium[,-burnin])[-1][int], col="turquoise", lwd = 1.2)



library(ggplot2)
int = 1:n
ggplot(data = data.frame(x = rep(int,2), 
                         y = c( y2[int], rowMeans(run7$calcium[,-burnin])[-1][int] ), 
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



AA = matrix(0,100,n)
for(i in 1:100)
{
  ii = i 
  AA[i, t(run7$clus)[ii,] >0] = run7$A[run7$clus[run7$clus[,ii] >0,ii]+1,ii]
}

est_spikes1 = colMeans(AA) ## run7
est_spikes1[which( apply(t(run7$clus)[-burnin,], 2, function(x) mean(x != 0))<0.8)] = 0
times1 = which(est_spikes1>0)
sum( est_spikes1>0)

sum(sapply(spp, function(x) !(x %in% times1))) / length(spp)  ### spikes non identificati falsi negativi
sum(sapply(times1, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi

AA = matrix(0,500,n)
for(i in 1:500)
{
  ii = i 
  AA[i, t(run7$clus)[ii,] >0] = run7$A[run7$clus[run7$clus[,ii] >0,ii]+1,ii]
}

est_spikes2 = colMeans(AA) ## run7
est_spikes2[which( apply(t(run7$clus)[-burnin,], 2, function(x) mean(x != 0))<0.9)] = 0
times2 = which(est_spikes2>0)
sum( est_spikes2>0)

sum(sapply(spp, function(x) !(x %in% times2)))/ length(spp)  ### spikes non identificati
sum(sapply(times2, function(x) !(x %in% spp))) / (n-length(spp))  ### falsi positivi


AA = matrix(0,500,n)
for(i in 1:500 )
{
  ii = i 
  AA[i, t(run8$clus)[ii,] >0] = run8$A[run8$clus[run8$clus[,ii] >0,ii]+1,ii]
}

est_spikes3 = colMeans(AA) ## run8
est_spikes3[which( apply(t(run8$clus)[-burnin,], 2, function(x) mean(x != 0))<0.9)] = 0
times3 = which(est_spikes3>0)
sum( est_spikes3>0)

sum(sapply(spp, function(x) !(x %in% times3))) / length(spp)### spikes non identificati
sum(sapply(times3, function(x) !(x %in% spp)))  / (n-length(spp))### falsi positivi

stimes1[times1 == 4700] = 4698
times2[times2 == 4700] = 4698
times3[times3 == 4700] = 4698
spp[spp == 4700] = 4698

times1[times1 == 4703] = 4704
times2[times2 == 4703] = 4704
times3[times3 == 4703] = 4704
spp[spp == 4703] = 4704

length(spp)


int = 4600:4800
ggplot(data = data.frame(x = rep(int,2), 
                         y = c( y2[int], rowMeans(run7$calcium[,-burnin])[-1][int] ), 
                         col = as.factor( c(rep(1, length(int)), rep(2, length(int)) ))), 
       aes(x = x, y = y ) ) + 
  geom_line(aes(color = col)) +
  theme_bw() +
  scale_x_continuous(name = "t", limits = c(min(int), max(int))) +
#  scale_y_continuous(name = expression(y[t]), limits = c(-0.35,0.67)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(t=-0.09, r=0, b=0, l=0, unit="cm"),
        legend.text.align = 0,
        legend.text=element_text(size=11,
                                 margin = margin(l = 0.003, unit = "cm")) ) +
  geom_point(data = data.frame(points = c(spp, times1, times2, times3),
                               h = c(rep(-0.2,length(spp)), 
                                     rep(-0.25, length(times1)), 
                                     rep(-0.28, length(times2)), 
                                     rep(-0.31, length(times3))),
                               set = as.factor( c(rep(3,length(spp)), 
                                                  rep(4, length(times1)), 
                                                  rep(5, length(times2)), 
                                                  rep(6, length(times3))) ) ),
             aes(x = points, y = h, color = set), shape = 3, size = 2 ) +
  scale_color_manual(values = c("#000000", "#00BFC4", c("black", viridis(4)[1:3]) ) ,
                     labels = c(expression(y[t]), expression(c[t]), "True spike", "4", "5", "6"))+
  guides(colour = guide_legend(byrow = TRUE, nrow = 1,
                               override.aes = list(linetype = c("solid", "solid", rep("blank", 4)),
                                                   shape = c(NA,NA,rep(3, 4)))),
         label.position = "right")
 


hist(apply(run7$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


A3 = AA[apply(AA, 1, function(x) length(unique( x[x!=0] )))==3,]
dataa = data.frame(A = A3[A3>0])
ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(bins = 35, aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(limits = c(0.15,0.8))


### number of clusters
plot(1:(nrep-max(burnin)), apply(run7$cluster[,-burnin], 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)
table(apply(run7$cluster[,-burnin], 2, function(x) length(unique(x[x>0]))))/1500


### cluster parameter
minA = min( which(apply(run7$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA -1
run7$A = run7$A[1:(minA +1),]
str(run7$A)

out_A = t(run7$A)
out_A[400:500,1:(minA)]

