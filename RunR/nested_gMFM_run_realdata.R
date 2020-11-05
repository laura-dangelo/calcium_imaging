library(Rcpp)
library(RcppDist)
library(ggplot2)
sourceCpp('./SourceCPP/calcium_nested_gMFM2.cpp')
library(viridis)

data <- read.csv("../data/cellula2.csv", header = FALSE)

y_real = c(data$V1)
length(y_real)
rm(list = ("data"))
plot(1:length(y_real), y_real, type = "l")
str(y_real)



g = rep(4, length(y_real))
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




plotAllLayers <- function(df)
{
  p = ggplot(data=df, aes(x=x, y=y)) +
    geom_line() +
    scale_y_continuous(limits = range(y_real))
  for(i in unique(int)[-1])
  { 
    p = p + 
      geom_rect(data = df, aes(xmin = min(x[interval == i]), 
                               xmax = max(x[interval == i]), 
                               ymin = -5, 
                               ymax = 5), fill = "yellow", alpha = 0.2) 
  }
  return(p)
}

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
  theme(legend.position = "bottom")


n = length(y_real)
nrep = 100

out = list()
out$calcium = matrix(c(0,y_real),length(y_real)+1, 1)
out$clusterO = matrix(0, length(y_real), 1)
out$clusterD = matrix(0, length(unique(g)), 1)
out$A = matrix(0, 100, 1)
out$AA = matrix(NA, length(y_real), 1)
out$b = 0
out$gamma = 0.5
out$sigma2 = 0.001
out$tau2 = 0.0003
out$p = 0.001
out$alpha = 0.5
out$beta = 0.5
out$maxK = 6
out$maxL = 15

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
                    alpha = 1, m = 1,
                    hyp_A1 = 5, hyp_A2 = 7,
                    hyp_b1 = 0, hyp_b2 = 1,
                    hyp_sigma21 = 1000, hyp_sigma22 = 1,
                    hyp_tau21 = 1000, hyp_tau22 = 1,
                    hyp_gamma1 = 1, hyp_gamma2 = 1,
                    hyp_p1 = 1, hyp_p2 = 999,
                    eps_gamma = 0.005,
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


# burnin = 1:400
# out$calcium = out$calcium[,-burnin]
# out$clusterO = out$clusterO[,-burnin]
# out$b = out$b[-burnin]
# out$gamma = out$gamma[-burnin]
# out$sigma2 = out$sigma2[-burnin]
# out$tau2 = out$tau2[-burnin]
# out$A = out$A[,-burnin]
# out$p = out$p[-burnin]
# out$alpha = out$alpha[-burnin]
# out$beta = out$beta[-burnin]
# out$maxK = out$maxK[-burnin]
# out$maxL = out$maxL[-burnin]


# save(out, file = "res_realdata_051120.Rdata")

