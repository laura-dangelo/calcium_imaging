library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)
sourceCpp('./SourceCPP/calcium_DP_innermixture.cpp')

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
                 prob = c(0.5, 0.3, 0.2), par = c(0.3, 0.5, 0.8))

y2 = data2$y
plot(y2, type = "l")

A_start = rep(0,50)
#  A_start[2:4] = c(0.3, 0.45, 0.65)
n = length(y2)
cluster = rep(0,length(y2))
#  cluster[data2$s == 1] = data2$k

nrep = 1500
set.seed(1234)

run = calcium_gibbs(Nrep = nrep, y = y2, 
                            cal = c(0,y2),
                            cl = cluster, 
                            A_start = A_start,
                            b_start = 0,
                            gamma_start = 0.6, 
                            sigma2_start = 0.002, 
                            tau2_start = 0.001, 
                            p_start = 0.01, 
                            c0 = 0, varC0 = 0.1, 
                            alpha = 1, m = 1,
                            hyp_A1 = 7, hyp_A2 = 10, 
                            hyp_b1 = 0, hyp_b2 = 1, 
                            hyp_sigma21 = 1000, hyp_sigma22 = 1, 
                            hyp_tau21 = 1000, hyp_tau22 = 1, 
                            hyp_gamma1 = 1, hyp_gamma2 = 1,
                            hyp_p1 = 1, hyp_p2 = 99,
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
# save(run, file = "res_inner_181020.Rdata")

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
                         y = c( y2[int], rowMeans(run$calcium[,-burnin])[-1][int] ), 
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




# grafico se fai simulazioni con diversi valori del parametro
# int = 4600:4800
# ggplot(data = data.frame(x = rep(int,2), 
#                          y = c( y2[int], rowMeans(run$calcium[,-burnin])[-1][int] ), 
#                          col = as.factor( c(rep(1, length(int)), rep(2, length(int)) ))), 
#        aes(x = x, y = y ) ) + 
#   geom_line(aes(color = col)) +
#   theme_bw() +
#   scale_x_continuous(name = "t", limits = c(min(int), max(int))) +
# #  scale_y_continuous(name = expression(y[t]), limits = c(-0.35,0.67)) +
#   theme(legend.title = element_blank(),
#         legend.position = "bottom",
#         legend.margin = margin(t=-0.09, r=0, b=0, l=0, unit="cm"),
#         legend.text.align = 0,
#         legend.text=element_text(size=11,
#                                  margin = margin(l = 0.003, unit = "cm")) ) +
#   geom_point(data = data.frame(points = c(spp, times1, times2, times3),
#                                h = c(rep(-0.2,length(spp)), 
#                                      rep(-0.25, length(times1)), 
#                                      rep(-0.28, length(times2)), 
#                                      rep(-0.31, length(times3))),
#                                set = as.factor( c(rep(3,length(spp)), 
#                                                   rep(4, length(times1)), 
#                                                   rep(5, length(times2)), 
#                                                   rep(6, length(times3))) ) ),
#              aes(x = points, y = h, color = set), shape = 3, size = 2 ) +
#   scale_color_manual(values = c("#000000", "#00BFC4", c("black", viridis(4)[1:3]) ) ,
#                      labels = c(expression(y[t]), expression(c[t]), "True spike", "4", "5", "6"))+
#   guides(colour = guide_legend(byrow = TRUE, nrow = 1,
#                                override.aes = list(linetype = c("solid", "solid", rep("blank", 4)),
#                                                    shape = c(NA,NA,rep(3, 4)))),
#          label.position = "right")
 


hist(apply(run$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")

### number of clusters
plot(1:(nrep-max(burnin)), apply(run$cluster[,-burnin], 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)
barplot(table(apply(run$cluster[,-burnin], 2, function(x) length(unique(x[x>0])))))


A3 = AA[apply(AA, 1, function(x) length(unique( x[x!=0] )))==3,]
dataa = data.frame(A = c(A3[c(A3)>0]) )
ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(limits = c(0.1,0.99))



### cluster parameter
minA = min( which(apply(run$A[-1,-burnin], 1, function(x) sum(x == 0)) == nrep-(max(burnin))) )
minA -1
run$A = run$A[1:(minA +1),]

out_A = t(run$A)
out_A[400:500,1:(minA)]

