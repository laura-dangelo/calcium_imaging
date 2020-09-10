library(Rcpp)
library(RcppDist)
sourceCpp('calcium_DP_mixture_nonconjugate_2708.cpp')

data <- read.csv("data/sessionB.csv", header = FALSE)


y_real = c(data$V1)
length(y_real)

rm(list = ("data"))
plot(1:length(y_real), y_real, type = "l")
str(y_real)

n = length(y_real)
# nrep = 100
# 
# 
# out = list()
# out$calcium = matrix(c(0,y_real),length(y_real)+1, 1)
# out$cluster = matrix(0, length(y_real), 1)
# out$A = matrix(0, 100, 1)
# out$AA = matrix(NA, length(y_real), 1)
# out$b = 0
# out$gamma = 0.5
# out$sigma2 = 0.001
# out$tau2 = 0.0003
# out$p = 0.9
# 
# plot(function(x) dgamma(x, 12, 7), 0, 5)
# 
# start <- Sys.time()
# debug = calcium_gibbs(Nrep = nrep, y = y_real, 
#                       cal = c(out$calcium[,length(out$b)]),
#                       cl = c(out$cluster[,length(out$b)]), 
#                       A_start = c(out$A[,length(out$b)]),
#                       b_start = c(out$b[length(out$b)]),
#                       gamma_start = c(out$gamma[length(out$b)]), 
#                       sigma2_start = c(out$sigma2[length(out$b)]), 
#                       tau2_start = c(out$tau2[length(out$b)]), 
#                       p_start = c(out$p[length(out$b)]), 
#                       c0 = 0, varC0 = 0.1, 
#                       alpha = 1, m = 5,
#                       hyp_A1 = 12, hyp_A2 = 7, 
#                       hyp_b1 = 0, hyp_b2 = 1, 
#                       hyp_sigma21 = 1000, hyp_sigma22 = 1, 
#                       hyp_tau21 = 1000, hyp_tau22 = 1, 
#                       hyp_gamma1 = 1, hyp_gamma2 = 1,
#                       hyp_p1 = 99, hyp_p2 = 1,
#                       eps_gamma = 0.003,
#                       eps_A = 0.0009)
# end <- Sys.time()
# end - start

# 
# out$calcium = cbind(out$calcium, debug$calcium)
# out$cluster = cbind(out$cluster, debug$cluster)
# out$b = c(out$b, debug$b)
# out$gamma = c(out$gamma, debug$gamma)
# out$sigma2 = c(out$sigma2, debug$sigma2)
# out$tau2 = c(out$tau2, debug$tau2)
# out$A = cbind(out$A, debug$A)
# out$p = c(out$p, debug$p)

# burnin = 1:400
# out$calcium = out$calcium[,-burnin]
# out$cluster = out$cluster[,-burnin]
# out$b = out$b[-burnin]
# out$gamma = out$gamma[-burnin]
# out$sigma2 = out$sigma2[-burnin]
# out$tau2 = out$tau2[-burnin]
# out$A = out$A[,-burnin]
# out$p = out$p[-burnin]


# save(out, file = "out_0809.Rdata")

burnin = c(1:100)

plot(1:length(out$p[-burnin]), out$p[-burnin], type = "l", xlab = "iterazioni", ylab = "p")
lines(1:length(out$p[-burnin]), cumsum(out$p[-burnin])/1:length(out$p[-burnin]), col =2)

plot(1:length(out$sigma2[-burnin]), out$sigma2[-burnin], type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(out$sigma2[-burnin]), cumsum(out$sigma2[-burnin])/1:length(out$sigma2[-burnin]), col =2)

plot(1:length(out$tau2[-burnin]), out$tau2[-burnin], type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(out$tau2[-burnin]), cumsum(out$tau2[-burnin])/1:length(out$tau2[-burnin]), col =2)

plot(1:length(out$b[-burnin]), out$b[-burnin], type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(out$b[-burnin]), cumsum(out$b[-burnin])/1:length(out$b[-burnin]), col =2)

plot(1:length(out$gamma[-burnin]), out$gamma[-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$gamma[-burnin]), cumsum(out$gamma[-burnin])/1:length(out$gamma[-burnin]), col =2)


mean(out$sigma2[-burnin])
mean(out$tau2[-burnin])
mean(out$b[-burnin])
mean(out$gamma[-burnin])

# plot(1:n, y_real, col="turquoise", type = "l")
# lines(0:n, rowMeans(out$calcium[,-burnin]))


obs = 501:1000
iter = 1:length(out$gamma)
image(1:(length(iter)), 1:(length(obs)), t(out$cluster[obs,iter]), 
      axes = F,
      xlab = "iterazioni", ylab = "osservazione",
      col = c("#ffffff", hcl.colors(10, "YlOrRd", rev = FALSE)))
axis(1, at = seq(1, max(iter), by = 10))
axis(2, at = seq(1,length(obs), by = 1), labels = obs)



# plot(1:n, y_real, type = "l")
AA = matrix(0,(length(out$b)-max(burnin)),n)
for(i in 1:(length(out$b)-max(burnin)))
{
  ii = i + (max(burnin))
  AA[i, t(out$clus)[ii,] >0] = out$A[out$clus[out$clus[,ii] >0,ii]+1,ii]
}

est_spikes = colMeans(AA)
est_spikes[which( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))<0.6)] = 0

plot(1:n, y_real, col = "turquoise", type = "l")
lines(1:n, rowMeans(out$calcium[,-burnin])[-1], col="salmon")
lines(1:n, est_spikes + mean(out$b))


n_spikes = sum( apply(t(out$clus)[-burnin,], 2, function(x) mean(x != 0))>0.8)
# hist(apply(out$clus[,-burnin], 1, function(x) mean(x != 0)), main = "Distr. of spike probabilities")


### number of clusters
plot(1:length(out$b), apply(out$cluster, 2, function(x) length(unique(x[x>0]))), pch=19, cex=0.2)


### cluster parameter
minA = min( which(apply(out$A[-1,-burnin], 1, function(x) sum(x == 0)) == length(out$b)-(max(burnin))) )
minA -1
out$A = out$A[1:(minA +1),]
str(out$A)

out_A = t(out$A)
out_A[300:500,1:(minA)]

t(apply(out_A[450:500,1:(minA)], 1, sort, decreasing = T))



### analisi degli stimoli ###
## static gratings ##
stat_grat <- read.csv("data/static_grating.csv", header = T)
head(stat_grat)
apply(stat_grat[,1:3], 2, unique)

# stat_grat$orientation==90) & (stat_grat$spatial_frequency == 0.02) & (stat_grat$phase == 0.00)
stat_grat1 = stat_grat[(stat_grat$orientation==90) & (stat_grat$spatial_frequency == 0.02) &
                         stat_grat$phase == 0.00,]
stat_grat1 = stat_grat1[!is.na(stat_grat1$orientation), ]
str(stat_grat1)


tmp = apply(stat_grat1[,4:5], 1, function(x) est_spikes[x[1]:x[2]])
est1 = matrix(0, length(tmp), 9, byrow = T)
for(i in 1:nrow(est1))
{
  est1[i,1:length(tmp[[i]])] = tmp[[i]]
}
plot(1:length(est1), c(est1), type = "l") # no spikes


# stat_grat$orientation==30 & (stat_grat$spatial_frequency == 0.02) 
## & (stat_grat$phase == 0.50)
stat_grat1 = stat_grat[(stat_grat$orientation==30) & (stat_grat$spatial_frequency == 0.02) &
                       (stat_grat$phase == 0.50),]

## ori = 30, freq = 0.02
### pahse = 0.00 -> niente
### phase = 0.25 -> niente
### phase = 0.50 -> 2 spikes
stat_grat1 = data.frame(stat_grat1[!is.na(stat_grat1$orientation), ], row.names = NULL)
head(stat_grat1)

est1 = t(apply(stat_grat1[,4:5], 1, function(x) est_spikes[x[1]:x[2]]))
str(est1)
est1
plot(1:length(est1), c(est1), type = "l")

which(apply(est1, 1, function(x) sum(x>0))>0)
stat_grat1[which(apply(est1, 1, function(x) sum(x>0))>0),]

int = 1:n
plot(int, y_real[int], type = "l", col = "turquoise3")
for(i in 1:nrow(stat_grat1))
{
  rect(xleft = stat_grat1[i,4], ybottom = -1, xright = stat_grat1[i,5], ytop = 10, col = "#f5c8bc", border = NA)
  #lines(stat_grat1[i,4]:stat_grat1[i,5], est_spikes[stat_grat1[i,4]:stat_grat1[i,5]], lwd = 0.7, col = 2)
}
lines(int, y_real[int], type = "l", col = "turquoise3")
segments(x0 = int, x1 = int, y0 = mean(out$b[-burnin]), y1 = est_spikes[int], lwd = 1.3, col = 1)
for(i in 1:nrow(stat_grat1))
{
  points((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0], 
         rep(-0.1, length((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0]) ), 
         pch = 3, col = 2, cex = 1.2)
}


par(mfrow=c(2,1))
int = 11857+c(-100:100)
plot(int, y_real[int], type = "l", col = "turquoise3")
for(i in 1:nrow(stat_grat1))
{
  rect(xleft = stat_grat1[i,4], ybottom = -1, xright = stat_grat1[i,5], ytop = 3, col = "#f5c8bc", border = NA)
  #lines(stat_grat1[i,4]:stat_grat1[i,5], est_spikes[stat_grat1[i,4]:stat_grat1[i,5]], lwd = 0.7, col = 2)
}
lines(int, y_real[int], type = "l", col = "turquoise3")
segments(x0 = int, x1 = int, y0 = mean(out$b[-burnin]), y1 = est_spikes[int], lwd = 1.3, col = 1)
for(i in 1:nrow(stat_grat1))
{
  points((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0], 
         rep(-0.1, length((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0]) ), 
         pch = 3, col = 2, cex = 1.2)
}

int = 106023+c(-100:100)
plot(int, y_real[int], type = "l", col = "turquoise3")
for(i in 1:nrow(stat_grat1))
{
  rect(xleft = stat_grat1[i,4], ybottom = -1, xright = stat_grat1[i,5], ytop = 3, col = "#f5c8bc", border = NA)
  #lines(stat_grat1[i,4]:stat_grat1[i,5], est_spikes[stat_grat1[i,4]:stat_grat1[i,5]], lwd = 0.7, col = 2)
}
lines(int, y_real[int], type = "l", col = "turquoise3")
segments(x0 = int, x1 = int, y0 = mean(out$b[-burnin]), y1 = est_spikes[int], lwd = 1.3, col = 1)
for(i in 1:nrow(stat_grat1))
{
  points((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0], 
         rep(-0.1, length((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0]) ), 
         pch = 3, col = 2, cex = 1.2)
}
par(mfrow=c(1,1))



# stat_grat$orientation==30 & (stat_grat$spatial_frequency == 0.02)
##& (stat_grat$phase == 0.75)
stat_grat1 = stat_grat[(stat_grat$orientation==30) & (stat_grat$spatial_frequency == 0.02) &
                         (stat_grat$phase == 0.75),]

## ori = 30, freq = 0.02
### pahse = 0.00 -> niente
### phase = 0.25 -> niente
### phase = 0.50 -> 2 spikes
### phase = 0.75 -> 6 spikes
stat_grat1 = data.frame(stat_grat1[!is.na(stat_grat1$orientation), ], row.names = NULL)
stat_grat1

est1 = apply(stat_grat1[,4:5], 1, function(x) est_spikes[x[1]:x[2]] )
str(est1)
est1
plot(1:length(unlist(est1)), unlist(est1), type = "l")

which(lapply(est1, function(x) sum(x>0))>0)
stat_grat1[which(lapply(est1, function(x) sum(x>0))>0),]


par(mfrow=c(3,1))
int = 768+c(-100:100)
plot(int, y_real[int], type = "l", col = "turquoise3")
for(i in 1:nrow(stat_grat1))
{
  rect(xleft = stat_grat1[i,4], ybottom = -1, xright = stat_grat1[i,5], ytop = 3, col = "#f5c8bc", border = NA)
  #lines(stat_grat1[i,4]:stat_grat1[i,5], est_spikes[stat_grat1[i,4]:stat_grat1[i,5]], lwd = 0.7, col = 2)
}
lines(int, y_real[int], type = "l", col = "turquoise3")
segments(x0 = int, x1 = int, y0 = mean(out$b[-burnin]), y1 = est_spikes[int], lwd = 1.3, col = 1)
for(i in 1:nrow(stat_grat1))
{
  points((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0], 
         rep(-0.1, length((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0]) ), 
         pch = 3, col = 2, cex = 1.2)
}

int = 6156+c(-100:100)
plot(int, y_real[int], type = "l", col = "turquoise3")
for(i in 1:nrow(stat_grat1))
{
  rect(xleft = stat_grat1[i,4], ybottom = -1, xright = stat_grat1[i,5], ytop = 3, col = "#f5c8bc", border = NA)
  #lines(stat_grat1[i,4]:stat_grat1[i,5], est_spikes[stat_grat1[i,4]:stat_grat1[i,5]], lwd = 0.7, col = 2)
}
lines(int, y_real[int], type = "l", col = "turquoise3")
segments(x0 = int, x1 = int, y0 = mean(out$b[-burnin]), y1 = est_spikes[int], lwd = 1.3, col = 1)
for(i in 1:nrow(stat_grat1))
{
  points((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0], 
         rep(-0.1, length((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0]) ), 
         pch = 3, col = 2, cex = 1.2)
}

int = 59665+c(-100:100)
plot(int, y_real[int], type = "l", col = "turquoise3")
for(i in 1:nrow(stat_grat1))
{
  rect(xleft = stat_grat1[i,4], ybottom = -1, xright = stat_grat1[i,5], ytop = 3, col = "#f5c8bc", border = NA)
  #lines(stat_grat1[i,4]:stat_grat1[i,5], est_spikes[stat_grat1[i,4]:stat_grat1[i,5]], lwd = 0.7, col = 2)
}
lines(int, y_real[int], type = "l", col = "turquoise3")
segments(x0 = int, x1 = int, y0 = mean(out$b[-burnin]), y1 = est_spikes[int], lwd = 1.3, col = 1)
for(i in 1:nrow(stat_grat1))
{
  points((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0], 
         rep(-0.1, length((stat_grat1[i,4]:stat_grat1[i,5])[est_spikes[stat_grat1[i,4]:stat_grat1[i,5]]>0]) ), 
         pch = 3, col = 2, cex = 1.2)
}


par(mfrow=c(1,1))
