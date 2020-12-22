library(Rcpp)
library(RcppDist)
library(ggplot2)
library(viridis)
sourceCpp('./SourceCPP/calcium_nested_gMFM2.cpp')


data <- read.csv("../data/cellula2.csv", header = FALSE)

y_real = c(data$V1)
length(y_real)
rm(list = ("data"))
# plot(1:length(y_real), y_real, type = "l")
# str(y_real)
n = length(y_real)


g = rep(4, length(y_real))
J = 4
int = g
stat_grat <- read.csv("../data/static_grating.csv", header = TRUE)

g[min(stat_grat$start[stat_grat$start < 30000]):max(stat_grat$end[stat_grat$end < 30000])] = 1
int[min(stat_grat$start[stat_grat$start < 30000]):max(stat_grat$end[stat_grat$end < 30000])] = 11

g[min(stat_grat$start[(stat_grat$start > 30000) & 
                        (stat_grat$start < 90000)]):max(stat_grat$end[(stat_grat$end > 30000) & 
                                                                        (stat_grat$end < 90000)])] = 1
int[min(stat_grat$start[(stat_grat$start > 30000) & 
                          (stat_grat$start < 90000)]):max(stat_grat$end[(stat_grat$end > 30000) & 
                                                                          (stat_grat$end < 90000)])] = 12
g[min(stat_grat$start[stat_grat$start > 90000]):max(stat_grat$end[stat_grat$end > 90000])] = 1
int[min(stat_grat$start[stat_grat$start > 90000]):max(stat_grat$end[stat_grat$end > 90000])] = 13


nat_scene <- read.csv("../data/natural_scene.csv", header = TRUE)
g[min(nat_scene$start[nat_scene$start < 35000]):max(nat_scene$end[nat_scene$end < 35000])] = 2
int[min(nat_scene$start[nat_scene$start < 35000]):max(nat_scene$end[nat_scene$end < 35000])] = 21
g[min(nat_scene$start[(nat_scene$start > 35000) & 
                        (nat_scene$start < 60000)]):max(nat_scene$end[(nat_scene$end > 35000) & 
                                                                        (nat_scene$end < 60000)])] = 2
int[min(nat_scene$start[(nat_scene$start > 35000) & 
                          (nat_scene$start < 60000)]):max(nat_scene$end[(nat_scene$end > 35000) & 
                                                                          (nat_scene$end < 60000)])] = 22
g[min(nat_scene$start[nat_scene$start > 60000]):max(nat_scene$end[nat_scene$end > 60000])] = 2
int[min(nat_scene$start[nat_scene$start > 60000]):max(nat_scene$end[nat_scene$end > 60000])] = 23


nat_movie <- read.csv("../data/natural_movie_one.csv", header = TRUE)
g[min(nat_movie$start):max(nat_movie$end)] = 3
int[min(nat_movie$start):max(nat_movie$end)] = 3





df = data.frame(x = 1:n, y = y_real, g = g, interval = int)

df_rect = data.frame(start = sapply(unique(df$interval)[-1], function(x) min(df$x[df$interval==x])),
                     end = sapply(unique(df$interval)[-1], function(x) max(df$x[df$interval==x])),
                     Stimulus = as.factor(c("Static grating","Natural scene","Natural scene",
                                            "Static grating","Natural movie","Natural scene",
                                            "Static grating")) )


cols = c("#91ff00", "#00fffb","#ff3700")

ggplot(data = df) +
  geom_rect(data = df_rect, inherit.aes = FALSE,
            aes(xmin = start, 
                xmax = end, 
                ymin = -Inf, 
                ymax = Inf, fill = Stimulus), alpha = 0.12 ) +
  scale_fill_manual(values = cols) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level") +
  annotate("text", x = 8000, y = 2.8, label = expression( "~G"[1]), size = 8 ) +
  annotate("text", x = 23500, y = 2.8, label = expression( "~G"[2]), size = 8 ) +
  annotate("text", x = 35000, y = 2.8, label = expression( "~G"[3]), size = 8 ) +
  annotate("text", x = 46500, y = 2.8, label = expression( "~G"[2]), size = 8 ) +
  annotate("text", x = 61500, y = 2.8, label = expression( "~G"[1]), size = 8 )+
  annotate("text", x = 74500, y = 2.8, label = expression( "~G"[4]), size = 8 )+
  annotate("text", x = 86500, y = 2.8, label = expression( "~G"[2]), size = 8 )+
  annotate("text", x = 104500, y = 2.8, label = expression( "~G"[1]), size = 8 )


# plot solo di un intervallo
int = 61000:61750

int = 20600:26000
ggplot(data = df[int,]) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")



out = list()
out$calcium = matrix(c(0,y_real),length(y_real)+1, 1)
out$clusterO = matrix(0, length(y_real), 1)
out$clusterD = matrix(c(1,2,2,3), length(unique(g)), 1)
out$A = matrix(0, 100, 1)
out$b = 0
out$gamma = 0.6
out$sigma2 = 0.001
out$tau2 = 0.0003
out$p = 0.001
out$alpha = 0.5
out$beta = 0.5
out$maxK = 6
out$maxL = 15

clus = kmeans(y_real[y_real>0.4], centers = 5)
out$A[2:6,1] = c(clus$centers)
out$clusterO[y_real>0.4,1] = clus$cluster


rm(list=c("stat_grat","nat_movie","nat_scene","int", "clus"))

nrep = 500
start <- Sys.time()
run = calcium_gibbs(Nrep = nrep, 
                    y = y_real,
                    g = g,                      
                    cal = c(out$calcium[,length(out$b)]),
                    clO = c(out$clusterO[,length(out$b)]), 
                    clD = c(out$clusterD[,length(out$b)]),
                    A_start = c(out$A[,length(out$b)]),
                    b_start = c(out$b[length(out$b)]),
                    gamma_start = c(out$gamma[length(out$b)]),
                    sigma2_start = c(out$sigma2[length(out$b)]),
                    tau2_start = c(out$tau2[length(out$b)]),
                    p_start = c(out$p[length(out$b)]),
                    alpha_start = c(out$alpha[length(out$b)]), 
                    beta_start = c(out$beta[length(out$b)]),
                    maxK_start = c(out$maxK[length(out$b)]),
                    maxL_start = c(out$maxL[length(out$b)]),
                    c0 = 0, varC0 = 0.1,
                    hyp_A1 = 5, hyp_A2 = 7,
                    hyp_b1 = 0, hyp_b2 = 1,
                    hyp_sigma21 = 1000, hyp_sigma22 = 1,
                    hyp_tau21 = 1000, hyp_tau22 = 1,
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    hyp_alpha1 = 3, hyp_alpha2 = 3,
                    hyp_beta1 = 3, hyp_beta2 = 6,
                    hyp_maxK1 = 2, hyp_maxK2 = 4, hyp_maxK3 = 3,
                    hyp_maxL1 = 2, hyp_maxL2 = 4, hyp_maxL3 = 3,
                    eps_alpha = 0.5, eps_beta = 0.7,
                    eps_gamma = 0.005,
                    eps_A = 0.02,
                    eps_maxK = 4, eps_maxL = 5)
end <- Sys.time()
end - start

str(run)

out$calcium = cbind(out$calcium, run$calcium)
out$clusterO = cbind(out$clusterO, run$clusterO)
out$clusterD = cbind(out$clusterD, run$clusterD)
out$b = c(out$b, run$b)
out$gamma = c(out$gamma, run$gamma)
out$sigma2 = c(out$sigma2, run$sigma2)
out$tau2 = c(out$tau2, run$tau2)
out$A = cbind(out$A, run$A)
out$alpha = c(out$alpha, run$alpha)
out$beta = c(out$beta, run$beta)
out$maxK = c(out$maxK, run$maxK)
out$maxL = c(out$maxL, run$maxL)
out$p = c(out$p, run$p)

rm(list=("run"))

burnin = 900:1100
out$calcium = out$calcium[,-burnin]
out$clusterO = out$clusterO[,-burnin]
out$clusterD = out$clusterD[,-burnin]
out$b = out$b[-burnin]
out$gamma = out$gamma[-burnin]
out$sigma2 = out$sigma2[-burnin]
out$tau2 = out$tau2[-burnin]
out$A = out$A[,-burnin]
out$p = out$p[-burnin]
out$alpha = out$alpha[-burnin]
out$beta = out$beta[-burnin]
out$maxK = out$maxK[-burnin]
out$maxL = out$maxL[-burnin]
str(out)

# save(out, file = "res_realdata_181220c.Rdata")
out$calcium = NULL
# burnin 2000


burnin = 1:10
plot(1:length(out$p[-burnin]), out$p[-burnin], type = "l")
lines(1:length(out$p[-burnin]), cumsum(out$p[-burnin])/1:length(out$p[-burnin]), col =2)

plot(1:length(out$sigma2[-burnin]), out$sigma2[-burnin], type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(out$sigma2[-burnin]), cumsum(out$sigma2[-burnin])/1:length(out$sigma2[-burnin]), col =2)

plot(1:length(out$tau[-burnin]), out$tau[-burnin], type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(out$tau[-burnin]), cumsum(out$tau[-burnin])/1:length(out$tau[-burnin]), col =2)

plot(1:length(out$b[-burnin]), out$b[-burnin], type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(out$b[-burnin]), cumsum(out$b[-burnin])/1:length(out$b[-burnin]), col =2)

plot(1:length(out$gamma[-burnin]), out$gamma[-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(out$gamma[-burnin]), cumsum(out$gamma[-burnin])/1:length(out$gamma[-burnin]), col =2)

plot(1:length(out$alpha[-burnin]), out$alpha[-burnin], type = "l", xlab = "iterazioni", ylab = "alpha")
lines(1:length(out$alpha[-burnin]), cumsum(out$alpha[-burnin])/1:length(out$alpha[-burnin]), col =2)

plot(1:length(out$beta[-burnin]), out$beta[-burnin], type = "l", xlab = "iterazioni", ylab = "beta")
lines(1:length(out$beta[-burnin]), cumsum(out$beta[-burnin])/1:length(out$beta[-burnin]), col =2)

plot(1:length(out$maxL[-burnin]), out$maxL[-burnin], type = "l", xlab = "iterazioni", ylab = "maxL")
lines(1:length(out$maxL[-burnin]), cumsum(out$maxL[-burnin])/1:length(out$maxL[-burnin]), col =2)

plot(1:length(out$maxK[-burnin]), out$maxK[-burnin], type = "l", xlab = "iterazioni", ylab = "maxK")
lines(1:length(out$maxK[-burnin]), cumsum(out$maxK[-burnin])/1:length(out$maxL[-burnin]), col =2)

out$clusterD[,100]


burnin = 1:1200
barplot(table(apply(out$clusterO[,-burnin], 2, function(x) length(unique(x)) )))

AA_gMFM = matrix(0,length(out$b[-burnin]),n)
n = nrow(out$clusterO)
for(i in 1:length(out$b[-burnin]))
{
  ii = i + max(burnin)
  AA_gMFM[i, t(out$clusterO)[ii,] >0] = out$A[out$clusterO[out$clusterO[,ii] >0,ii]+1,ii]
}



# save(AA_gMFM, file = "AA_realdata_181220.Rdata")
est_spikes = colMeans(AA_gMFM) 
est_spikes[which( apply(t(out$clusterO)[-burnin,], 2, function(x) mean(x != 0))<0.5)] = 0
times = which(est_spikes>0)
length(times)

AA_gMFM[,which(est_spikes == 0)] = 0
barplot(table( apply(AA_gMFM, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(AA_gMFM, 1, function(x) length(unique(x))) )), "names"))

A_ind = AA_gMFM[apply(AA_gMFM, 1, function(x) length(unique( x )))==moda,]
datAA_gMFM = data.frame(A = A_ind[A_ind>0])
ggplot(data = datAA_gMFM, aes(x = A)) + 
  geom_histogram(bins = 30, aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  xlim(c(0.1,1.5))





#--------------------------------------------#
int = which(g==1)
int = which(g==2)
int = which(g==3)
int = which(g==4)
int = which((g==3)|(g==2))

length(times[times %in% int])
length(times[times %in% int])/length(int)


intt = which(g==1)
subsetAA = AA_gMFM[,intt]
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
dataa = data.frame(A = A_ind[A_ind>0])
plot1 <- ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(aes(x = A[A<1.11], y = ..density..), bins = 19, col = "#00AFBB", fill = "#00AFBB", alpha = 0.3, size = 0.2) +   
#  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(limits=c(0.1,1.75), breaks = seq(0.15,1.75, by=0.25), name = "") + 
  scale_y_continuous(name = "Static gratings") +
  theme(panel.grid.minor = element_blank())



intt = which((g==2))
subsetAA = AA_gMFM[,intt]
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
dataa = data.frame(A = A_ind[A_ind>0])
plot2 <- ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(aes(y = ..density..), bins = 19, col = "#00AFBB", fill = "#00AFBB", alpha = 0.3, size = 0.2) +   
#  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(limits=c(0.1,1.75), breaks = seq(0.15,1.75,by=0.25), name = "") + 
  scale_y_continuous(name = "Natural scene")+
  theme(panel.grid.minor = element_blank())



intt = which((g==3))
subsetAA = AA_gMFM[,intt]
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
dataa = data.frame(A = A_ind[A_ind>0])
plot3 <- ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(aes(y = ..density..), bins = 19, col = "#00AFBB", fill = "#00AFBB", alpha = 0.3, size = 0.2) +   
#  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(limits=c(0.1,1.75), breaks = seq(0.15,1.75,by=0.25), name = "A") + 
  scale_y_continuous(name = "Natural movie")+
  theme(panel.grid.minor = element_blank())


require(gridExtra)
grid.arrange(plot1, plot2, plot3, nrow=3)

# tot spikes = 1476
# group 1: 460 - 0.01018398
# group 2: 829 - 0.01851397 
# group 3: 186 - 0.02059801 -- media 0.1886
# group 4: 1
#--------------------------------------------#



burnin = 1:10
barplot(table(apply(out$clusterD[,-burnin], 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
moda = as.numeric(attr(which.max(table(apply(out$clusterD[,-burnin], 2, function(x) length(unique(x)) ))), "names"))

mat_clusterD = matrix(NA, J, J)
ind3 = which( apply(out$clusterD[,-burnin], 2, function(x) length(unique(x)) ) ==moda ) 
mat_heatmap = expand.grid(J1 = 1:4,
                          J2 = 1:4)
for(i in 1:J)
{
  for(j in 1:i)
  {
    mat_clusterD[i,j] = sum( apply(out$clusterD[,ind3], 2, function(x) x[i] == x[j] ) )
    mat_clusterD[j,i] = mat_clusterD[i,j] 
  }
}
df_heat = data.frame(J1 = as.factor(mat_heatmap[,1]),
                     Group = as.factor(mat_heatmap[,2]),
                     val = c(mat_clusterD)/mat_clusterD[1,1],
                     lab = round( c(mat_clusterD)/mat_clusterD[1,1] , 3) )
df_heat$lab[df_heat$J1==df_heat$Group] = 1:J
df_heat$val[df_heat$J1==df_heat$Group] = NA

ggplot(data = df_heat) +
  geom_tile( aes(x = J1, y = Group, fill = val)) +
  geom_text(aes(x = J1, y = Group, label = lab), size=2.5) +
  scale_fill_gradient(low = magma(3)[3], high = inferno(3)[2]) +
  ylim(rev(levels(df_heat$Group))) +
  scale_x_discrete(name = "Group") +
  theme(legend.title = element_blank())



interval = 1:n
interval = which(g==1)[which(g==1)<30000]
interval = which(g==2)
interval = which(g==3)
interval = which(g==4)

dff = data.frame(x = (1:n)[interval], y = y_real[interval], AA = colMeans(AA_gMFM)[interval])
ggplot(data = dff) +
  geom_line(aes(x = x, y = y), col = "turquoise") +
  geom_line(aes(x = x, y = AA)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")



ggplot(data = dff) +
  geom_rect(data = df_rect, inherit.aes = FALSE,
            aes(xmin = start, 
                xmax = end, 
                ymin = -Inf, 
                ymax = Inf, fill = Stimulus), alpha = 0.12 ) +
  scale_fill_manual(values = cols) +
  geom_line(aes(x = x, y = y), size = 0.8) +
  geom_line(aes(x = x, y = AA), col = "gold") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")





post_clusterO = out$clusterO
rm(list=c("out"))
str(post_clusterO)
post_clusterO_positive = 