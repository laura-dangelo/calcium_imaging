library(mcclust)
library(mcclust.ext)

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/dati_reali/res_realdata_181220c.Rdata")
str(out)
clusterO = out$clusterO
A_par = out$A
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

spike_yes = which(apply(clusterO, 1, function(x) sum(x>0) )>600)
clus_mat <- comp.psm(t(clusterO[spike_yes,])+1)

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

# per gruppo
# int1 = 24300
# int2 = 24800
# plot(int1:int2,y_real[g==group][int1:int2], type="l")
# rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
# abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
# lines(int1:int2,y_real[g==group][int1:int2], type="l")



which(diff(g)!=0)
#   744  15198  16100  30550  39580  54029  54931  69391  70293  79323  80226  96104  97382 113637
# 4 -   1  -   4   -  2   - 4  -   2   -  4  -   1  -   4   -  3  -   4   -  2   -  4   -  1  -   4

# senza gruppo
int1 = 1
int2 = length(y_real)
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")


### static grating
par(mfrow = c(1,1))
int1 = 744
int2 = 15198
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")

int1 = 54931
int2 = 69391
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")

int1 = 97382
int2 = 113637
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")





### natural scene
int1 = 16100
int2 = 30550
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")

int1 = 39580
int2 = 54029
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")

int1 = 80226
int2 = 96104
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")


### natural movie
int1 = 70293
int2 = 79323
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")






### zoom di stat grat
# int1 = 11700
# int2 = 12000
int1 = 64800
int2 = int1+500
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")

AA_gMFM[,spike_cl[rows,1]]

### zoom di nat scene
int1 = 87500
int2 = int1+500

int1 = 83340
int2 = int1+500
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")

### zoom di nat movie
int1 = 70400
int2 = 70900
plot(int1:int2,y_real[int1:int2], type="l")
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
abline(v=spike_cl[rows,1], col=spike_cl[rows,2]+1, lwd=0.4)
lines(int1:int2,y_real[int1:int2], type="l")





AA_gMFM = AA_gMFM[,spike_yes]
max(minv$cl)
z = density(AA_gMFM[,minv$cl==1])
z$x[z$y==max(z$y)]

z = density(AA_gMFM[,minv$cl==2])
z$x[z$y==max(z$y)]

z = density(AA_gMFM[,minv$cl==3])
z$x[z$y==max(z$y)]

z = density(AA_gMFM[,minv$cl==4])
z$x[z$y==max(z$y)]

z = density(AA_gMFM[,minv$cl==5])
z$x[z$y==max(z$y)]

z = density(AA_gMFM[,minv$cl==6])
z$x[z$y==max(z$y)]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(6)





int1 = 64750
int2 = int1+500
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
df = data.frame("x" = int1:int2, "y" = y_real[int1:int2],
                "cluster" = rep(0,length(int1:int2)))
df$cluster[df$x %in% spike_cl[rows,1]] = spike_cl[rows,2]+1
df$A = NA
df$A[df$cluster == 1] = "0.17"
df$A[df$cluster == 2] = "0.47"
df$A[df$cluster == 3] = "0.84"
df$A[df$cluster == 4] = "1.07"
df$A[df$cluster == 5] = "1.40"
df$A[df$cluster == 6] = "1.65"
df$A = as.factor(df$A)
levels(df$A) = c("0.17", "0.47", "0.84", "1.07", "1.40","1.65")
df$cluster[df$cluster==0] = NA


g1 = ggplot(data = df) +
  geom_line(aes(x = x, y = y)) +
  geom_vline(data = df[!is.na(df$A),], aes(xintercept = x, color = A), alpha = 0.5) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Static grating") +
  scale_color_manual(limits = c("0.17", "0.47", "0.84", "1.07", "1.40","1.65"), values = cols) +
  guides(colour = guide_legend(nrow=1, override.aes = list(size=3)))




int1 = 87500
int2 = int1+500
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
df = data.frame("x" = int1:int2, "y" = y_real[int1:int2],
                "cluster" = rep(0,length(int1:int2)))
df$cluster[df$x %in% spike_cl[rows,1]] = spike_cl[rows,2]+1
df$A = NA
df$A[df$cluster == 1] = "0.17"
df$A[df$cluster == 2] = "0.47"
df$A[df$cluster == 3] = "0.84"
df$A[df$cluster == 4] = "1.07"
df$A[df$cluster == 5] = "1.40"
df$A[df$cluster == 6] = "1.65"
df$A = as.factor(df$A)
df$cluster[df$cluster==0] = NA
levels(df$A) = c("0.17", "0.47", "0.84", "1.07", "1.40","1.65")

g2 = ggplot(data = df) +
  geom_line(aes(x = x, y = y)) +
  geom_vline(data = df[!is.na(df$A),], aes(xintercept = x, color = A), alpha = 0.5) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  scale_color_manual(values = cols[1:4]) +
  theme(legend.position = "none") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Natural scene") 




int1 = 70400
int2 = 70900
rows = which((spike_cl[,1]>int1)&(spike_cl[,1]<int2))
df = data.frame("x" = int1:int2, "y" = y_real[int1:int2],
                "cluster" = rep(0,length(int1:int2)))
df$cluster[df$x %in% spike_cl[rows,1]] = spike_cl[rows,2]+1
df$A = NA
df$A[df$cluster == 1] = "0.17"
df$A[df$cluster == 2] = "0.47"
df$A[df$cluster == 3] = "0.84"
df$A[df$cluster == 4] = "1.07"
df$A[df$cluster == 5] = "1.40"
df$A[df$cluster == 6] = "1.65"
df$cluster[df$cluster==0] = NA
levels(df$A) = c("0.17", "0.47", "0.84", "1.07", "1.40","1.65")

df

g3 = ggplot(data = df) +
  geom_line(aes(x = x, y = y)) +
  geom_vline(data = df[!is.na(df$A),], aes(xintercept = x, color = A), alpha = 0.5) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  scale_color_manual(values = cols) +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Natural movie") 



library(ggpubr)


ggarrange(g1, g2, g3, ncol=1, nrow=3, common.legend = TRUE, legend="bottom")

