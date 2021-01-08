# library(ggplot2)
# library(viridis)
library(mclust)

gammapar = 8
nsim = 31
filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par", 
                  gammapar, "/scen2_run_DP_gammapar", gammapar, "_sim", nsim, ".Rdata")
load(file = filename)
filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/data/data_scen2_seed", 
                  nsim, ".Rdata")
load(file = filename)

y = out$y
g = out$g
A = out$A
s = out$s
spp = which(s>0)
k = out$k
plot(y[2000:2600], type = "l")

n = length(y)
J = length(unique(g))
n1 = n2 = n3 = n4 = n5 = 5000
prob1 = c(0.2, 0.2, 0.15, 0.15, 0.15, 0.15)
par1 = c(0.3, 0.50, 0.7, 0.9, 1.1, 1.5)

prob2 = rep(0.25, 4)
par2 = c(0.3, 0.90, 1.5, 1.8)

prob3 = c(0.4, 0.2, 0.4)
par3 = c(0.5, 0.90, 1.5)

### parametri: 
unip <- sort(unique(c(par1,par2,par3)))
unip
length(unip)


plot(1:length(run_DP$p), run_DP$p, type = "l")
lines(1:length(run_DP$p), cumsum(run_DP$p)/1:length(run_DP$p), col =2)

run_DP$b
run_DP$sigma2
run_DP$tau2
run_DP$gamma
plot(1:length(run_DP$A[2, ]), run_DP$A[2, ], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[3, ]), run_DP$A[3, ], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[4, ]), run_DP$A[4, ], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5, ]), run_DP$A[5, ], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5, ]), run_DP$A[6, ], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5, ]), run_DP$A[7, ], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5, ]), run_DP$A[8, ], type = "l", xlab = "iterazioni", ylab = "A")

burnin = 1:100
barplot(table(apply(run_DP$clusterO[,-burnin], 2, function(x) length(unique(x)) )))


burnin = 1:800
false_negatives = numeric(50)
false_positives = numeric(50)
# number_clustersO = numeric(50)
# number_clustersD = numeric(50)
# number_clustersO_j = matrix(NA, 50, J)

rand_indexO = matrix(0, 50, 201)
rand_indexD = matrix(0, 50, 201)

for(nsim in 2:50)
{
  filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par", 
                    gammapar, "/scen2_run_DP_gammapar", gammapar, "_sim", nsim, ".Rdata")
  load(file = filename)
  filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/data/data_scen2_seed", 
                    nsim, ".Rdata")
  load(file = filename) 
  y = out$y
  g = out$g
  A = out$A
  s = out$s
  spp = which(s>0)
  k = out$k
  n = length(y)
  
  burnin = 1:800
  AA_DP = matrix(0,length(run_DP$p[-burnin]),n)
  for(i in 1:length(run_DP$p[-burnin]))
  {
    ii = i + max(burnin)
    AA_DP[i, t(run_DP$clusterO)[ii,] >0] = run_DP$A[run_DP$clusterO[run_DP$clusterO[,ii] >0,ii]+1,ii]
  }
  filenameAA = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par",
                      gammapar, "/resAA_DP_gammapar", gammapar, "_", nsim, ".Rdata")
  # load(file = filenameAA)
  est_spikes = colMeans(AA_DP) 
  est_spikes[which( apply(t(run_DP$clusterO)[-burnin,], 2, function(x) mean(x != 0))<0.5)] = 0
  times = which(est_spikes>0)
  
  if(length(times)>0)
  {
    save(AA_DP, file = filenameAA)
    
    AA_cluster = apply(AA_DP, 1, rank )
    rand_indexO[nsim,] = apply(AA_cluster, 2, function(x) adjustedRandIndex(x,rank(A)) )
    rand_indexD[nsim,] = apply(run_DP$clusterD[,-burnin], 2, function(x) adjustedRandIndex(x, c(1,2,3,1) ) )
    
    false_negatives[nsim] = sum(sapply(spp, function(x) !(x %in% times))) / length(spp)  ### spikes non identificati: falsi negativi
    false_positives[nsim] = sum(sapply(times, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
    
  }
  if(nsim%%5 == 0) print(nsim)
}


TP = numeric(50)
TN = numeric(50)
FN = numeric(50)
FP = numeric(50)
for(nsim in 1:50)
{
  filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/data/data_scen2_seed", 
                    nsim, ".Rdata")
  load(file = filename) 
  
  s = out$s
  spp = which(s>0)
  n = length(out$y)
  
  FN[nsim] = round(false_negatives[nsim] * length(spp))
  FP[nsim] = round(false_positives[nsim] * (n-length(spp)))
  
  TN[nsim] = (n-length(spp)) - FP[nsim]
  TP[nsim] = length(spp) - FN[nsim]
}
save(FN, file=paste0("scen2DP_par", gammapar, "_FN.Rdata") )
save(FP, file=paste0("scen2DP_par", gammapar, "_FP.Rdata") )
save(TN, file=paste0("scen2DP_par", gammapar, "_TN.Rdata") )
save(TP, file=paste0("scen2DP_par", gammapar, "_TP.Rdata") )





save(rand_indexO, file=paste0("scen2DP_par", gammapar, "_randindexO.Rdata") )
save(rand_indexD, file=paste0("scen2DP_par", gammapar, "_randindexD.Rdata") )

rowMeans(rand_indexO)
mean(rowMeans(rand_indexO))
sd(rowMeans(rand_indexO))
summary(rowMeans(rand_indexO))

mean(rowMeans(rand_indexD))
sd(rowMeans(rand_indexD))



save(false_negatives, file=paste0("scen2DP_par", gammapar, "_false_neg.Rdata") )
save(false_positives, file=paste0("scen2DP_par", gammapar, "_false_pos.Rdata") )
# save(number_clustersD, file=paste0("scen2_par", gammapar, "_nclusD.Rdata") )
# save(number_clustersO, file=paste0("scen2_par", gammapar, "_nclusO.Rdata") )
# save(number_clustersO_j, file=paste0("scen2_par", gammapar, "_nclusOj.Rdata") )

false_negatives
summary(false_negatives)
sd(false_negatives)

false_positives
summary(false_positives)
sd(false_positives)

number_clustersO
table(number_clustersO)

number_clustersD

number_clustersO_j




AA_DP[,which(est_spikes == 0)] = 0
barplot(table( apply(AA_DP, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(AA_DP, 1, function(x) length(unique(x))) )), "names"))

A_ind = AA_DP[apply(AA_DP, 1, function(x) length(unique( x )))==moda,]
datAA_DP = data.frame(A = A_ind[A_ind>0])
ggplot(data = datAA_DP, aes(x = A)) + 
  geom_histogram(bins = 30, aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(breaks = unip)






#--------------------------------------------#
int = 1:n1
int = n1:(n1+n2)
int = (n1+n2):(n1+n2+n3)
int = (n1+n2+n3):(n1+n2+n3+n4)
int = (n1+n2+n3+n4):(n1+n2+n3+n4+n5)


subsetAA = AA_DP[,int]

plot(1:nrow(subsetAA), apply(subsetAA, 1, function(x) length(unique(x))) , type = "l", xlab = "iterazioni", ylab = "maxK")
lines(1:length(run_DP$maxK[-burnin]), cumsum(run_DP$maxK[-burnin])/1:length(run_DP$maxL[-burnin]), col =2)

barplot(table( apply(subsetAA, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))

A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
dataa = data.frame(A = A_ind[A_ind>0])
ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(bins = 35, aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(breaks = unip)
#--------------------------------------------#

val = matrix(NA, 50, J*J)

for(nsim in 1:50)
{
  filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen", 1,"/Res/par", 
                    gammapar, "/scen", 1, "_run_DP_gammapar", gammapar, "_sim", nsim, ".Rdata")
  load(file = filename)
  
  #barplot(table(apply(run_DP$clusterD, 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
  moda = as.numeric(attr(which.max(table(apply(run_DP$clusterD, 2, function(x) length(unique(x)) ))), "names"))
  
  mat_clusterD = matrix(NA, J, J)
  ind3 = which( apply(run_DP$clusterD, 2, function(x) length(unique(x)) ) ==moda ) 
  mat_heatmap = expand.grid(J1 = unique(g),
                            J2 = unique(g))
  for(i in 1:J)
  {
    for(j in 1:i)
    {
      mat_clusterD[i,j] = sum( apply(run_DP$clusterD[,ind3], 2, function(x) x[i] == x[j] ) )
      mat_clusterD[j,i] = mat_clusterD[i,j] 
    }
  }
  val[nsim,] = c(mat_clusterD)/mat_clusterD[1,1]
}

df_heat = data.frame(J1 = as.factor(mat_heatmap[,1]),
                     j = as.factor(mat_heatmap[,2]),
                     val = colMeans(val),
                     lab = round( colMeans(val) , 3) )
df_heat$lab[df_heat$J1==df_heat$j] = 1:J
df_heat$val[df_heat$J1==df_heat$j] = NA
  
ggplot(data = df_heat) +
    geom_tile( aes(x = J1, y = j, fill = val)) +
    geom_text(aes(x = J1, y = j, label = lab), size=2.5) +
    scale_fill_gradient(low = magma(3)[3], high = inferno(3)[2]) +
    ylim(rev(levels(df_heat$j))) +
  scale_x_discrete(name="j")
  


















