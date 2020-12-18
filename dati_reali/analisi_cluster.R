library(mcclust)
library(mcclust.ext)

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/dati_reali/res_realdata_161220d.Rdata")
str(out)
clusterO = out$clusterO
rm(list=("out"))
str(clusterO)

g = rep(4, nrow(clusterO))
J = 4
stat_grat <- read.csv("../data/static_grating.csv", header = TRUE)

g[min(stat_grat$start[stat_grat$start < 30000]):max(stat_grat$end[stat_grat$end < 30000])] = 1
g[min(stat_grat$start[(stat_grat$start > 30000) & 
                        (stat_grat$start < 90000)]):max(stat_grat$end[(stat_grat$end > 30000) & 
                                                                        (stat_grat$end < 90000)])] = 1
g[min(stat_grat$start[stat_grat$start > 90000]):max(stat_grat$end[stat_grat$end > 90000])] = 1


nat_scene <- read.csv("../data/natural_scene.csv", header = TRUE)
g[min(nat_scene$start[nat_scene$start < 35000]):max(nat_scene$end[nat_scene$end < 35000])] = 2
g[min(nat_scene$start[(nat_scene$start > 35000) & 
                        (nat_scene$start < 60000)]):max(nat_scene$end[(nat_scene$end > 35000) & 
                                                                        (nat_scene$end < 60000)])] = 2
g[min(nat_scene$start[nat_scene$start > 60000]):max(nat_scene$end[nat_scene$end > 60000])] = 2


nat_movie <- read.csv("../data/natural_movie_one.csv", header = TRUE)
g[min(nat_movie$start):max(nat_movie$end)] = 3
rm(list=c("nat_movie", "nat_scene", "stat_grat"))


group = 1
clusterO = clusterO[g==group,]
str(clusterO)

spike_yes = which(apply(clusterO, 1, function(x) sum(x>0) )>600)
clus_mat <- comp.psm(t(clusterO[spike_yes,])+1)


rm(clusterO)



# seq_BINDER_loss <- function(clust_mat){
#   res <- rep(0, nrow(clust_mat))
#   for(j in 1:ncol(clust_mat)){
#     for(l in j:ncol(clust_mat)){
#       temp <- mean(clust_mat[,l] == clust_mat[,j])
#       for(rep in 1:nrow(clust_mat)){
#         res[rep] <- res[rep] + abs(ifelse(clust_mat[rep,l] == clust_mat[rep,j], 1 , 0) - temp)
#       }
#     }
#   }
#   return(clust_mat[which.min(res),])
# }
# 
# seq_BINDER_loss(clust_mat = clus_mat)


minv <- minVI(clus_mat)
str(minv)
spike_cl <- cbind(spike_yes, minv$cl)

data <- read.csv("../data/cellula2.csv", header = FALSE)
y_real = c(data$V1)
length(y_real)
rm(list = ("data"))

int1 = 24300
int2 = 24800
plot(int1:int2,y_real[g==group][int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[g==group][int1:int2], type="l")


int1 = 25000
int2 = 26000
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")



which(diff(g)!=0)
#  744  15198  16100  30550  39580  54029  54931  69391  70293  79323  80226  96104  97382 113637
