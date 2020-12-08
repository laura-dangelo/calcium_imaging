library(ggplot2)
library(viridis)

gammapar = 4
nsim = 1

load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/data/y_scen1.Rdata")
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/data/g_scen1.Rdata")
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/data/spp_scen1.Rdata")

n = length(y)
J = length(unique(g))
n1 = n2 = n3 = n4 = n5 = n6 = 5000
prob1 = c(0.25, 0.25, 0.2, 0.15, 0.15)
par1 = c(0.5, 1.2, 1.7, 2, 2.3)

prob2 = rep(0.25, 4)
par2 = c(0.9, 1.2, 1.7, 2.3)

prob3 = c(0.4, 0.2, 0.4)
par3 = c(0.5, 1.2, 1.7)

prob4 = rep(1, 3)/3
par4 = c(0.5, 0.9, 2)

### parametri: 
unip <- sort(unique(c(par1,par2,par3,par4)))
unip
length(unip)

filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen", 1,"/Res/par", 
                  gammapar, "/scen", 1, "_run_gMFM_gammapar", gammapar, "_sim", nsim, ".Rdata")
load(file = filename)

str(run_gMFM)
filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen", 1,"/Res/par", 
                  gammapar, "/resAA", nsim, ".Rdata")
load(file = filename)
filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen", 1,"/Res/DP/resAA", nsim, ".Rdata")
load(file = filename)

#### gruppo 1
gen_mix <- function(N, prob, par)
{
  out <- numeric(N)
  clus = sample(1:length(prob), N, prob = prob, replace = T)
  out <- rnorm(N, par[clus], 0.00005)
}

int = 1:n1
prob1 = c(0.25, 0.25, 0.2, 0.15, 0.15)
par1 = c(0.5, 1.2, 1.7, 2, 2.3)
plot(density(gen_mix(10000, prob1, par1)))
subsetAA = AA_gMFM[,int]
subsetAAdp = AA_DP[,int]

barplot(table( apply(subsetAA, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
moda
barplot(table( apply(subsetAAdp, 1, function(x) length(unique(x))) ))
modadp = as.numeric(attr(which.max(table( apply(subsetAAdp, 1, function(x) length(unique(x))) )), "names"))
modadp

A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
A_inddp = subsetAA[apply(subsetAAdp, 1, function(x) length(unique( x )))== modadp,]
dataa = data.frame(A = A_ind[A_ind>0])
dataadp = data.frame(Adp = A_inddp[A_inddp>0])
ggplot(data = dataa) + 
  geom_histogram(data = dataadp, aes(x=Adp, y = ..density..), bins = 35, col = "salmon", fill = "salmon", alpha = 0.3) + 
  geom_histogram(bins = 35, aes(x=A, y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  geom_density(data = data.frame(true = gen_mix(10000, prob1, par1)), aes(x=true, y = ..density..) )+
  theme_bw() +
  scale_x_continuous(breaks = unip, name = "") +
  scale_y_continuous(name = "")



int = n1:(n1+n2)
prob2 = rep(0.25, 4)
par2 = c(0.9, 1.2, 1.7, 2.3)

subsetAA = AA_gMFM[,int]
subsetAAdp = AA_DP[,int]

barplot(table( apply(subsetAA, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
moda
barplot(table( apply(subsetAAdp, 1, function(x) length(unique(x))) ))
modadp = as.numeric(attr(which.max(table( apply(subsetAAdp, 1, function(x) length(unique(x))) )), "names"))
modadp

A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
A_inddp = subsetAA[apply(subsetAAdp, 1, function(x) length(unique( x )))== modadp,]
dataa = data.frame(A = A_ind[A_ind>0])
dataadp = data.frame(Adp = A_inddp[A_inddp>0])
ggplot(data = dataa) + 
  geom_histogram(data = dataadp, aes(x=Adp, y = ..density..), bins = 35, col = "salmon", fill = "salmon", alpha = 0.3) + 
  geom_histogram(bins = 35, aes(x=A, y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  geom_density(data = data.frame(true = gen_mix(10000, prob2, par2)), aes(x=true, y = ..density..) )+
  theme_bw() +
  scale_x_continuous(breaks = unip, name = "") +
  scale_y_continuous(name = "")







int = (n1+n2):(n1+n2+n3)
prob3 = c(0.4, 0.2, 0.4)
par3 = c(0.5, 1.2, 1.7)

subsetAA = AA_gMFM[,int]
subsetAAdp = AA_DP[,int]

barplot(table( apply(subsetAA, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
moda
barplot(table( apply(subsetAAdp, 1, function(x) length(unique(x))) ))
modadp = as.numeric(attr(which.max(table( apply(subsetAAdp, 1, function(x) length(unique(x))) )), "names"))
modadp

A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
A_inddp = subsetAA[apply(subsetAAdp, 1, function(x) length(unique( x )))== modadp,]
dataa = data.frame(A = A_ind[A_ind>0])
dataadp = data.frame(Adp = A_inddp[A_inddp>0])
ggplot(data = dataa) + 
  geom_histogram(data = dataadp, aes(x=Adp, y = ..density..), bins = 35, col = "salmon", fill = "salmon", alpha = 0.3) + 
  geom_histogram(bins = 35, aes(x=A, y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  geom_density(data = data.frame(true = gen_mix(10000, prob3, par3)), aes(x=true, y = ..density..) )+
  theme_bw() +
  scale_x_continuous(breaks = unip, name = "") +
  scale_y_continuous(name = "")





int = (n1+n2+n3):(n1+n2+n3+n4)
prob4 = rep(1, 3)/3
par4 = c(0.5, 0.9, 2)


subsetAA = AA_gMFM[,int]
subsetAAdp = AA_DP[,int]

barplot(table( apply(subsetAA, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(subsetAA, 1, function(x) length(unique(x))) )), "names"))
moda
barplot(table( apply(subsetAAdp, 1, function(x) length(unique(x))) ))
modadp = as.numeric(attr(which.max(table( apply(subsetAAdp, 1, function(x) length(unique(x))) )), "names"))
modadp

A_ind = subsetAA[apply(subsetAA, 1, function(x) length(unique( x )))== moda,]
A_inddp = subsetAA[apply(subsetAAdp, 1, function(x) length(unique( x )))== modadp,]
dataa = data.frame(A = A_ind[A_ind>0])
dataadp = data.frame(Adp = A_inddp[A_inddp>0])
ggplot(data = dataa) + 
  geom_histogram(data = dataadp, aes(x=Adp, y = ..density..), bins = 35, col = "salmon", fill = "salmon", alpha = 0.3) + 
  geom_histogram(bins = 35, aes(x=A, y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  geom_density(data = data.frame(true = gen_mix(10000, prob4, par4)), aes(x=true, y = ..density..) )+
  theme_bw() +
  scale_x_continuous(breaks = unip, name = "") +
  scale_y_continuous(name = "")




int = (n1+n2+n3+n4):(n1+n2+n3+n4+n5)
int = (n1+n2+n3+n4+n5):(n1+n2+n3+n4+n5+n6)






truen = c(6,5,4,4,6,5)

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/Res/par4/scen1_par4_nclusOj.Rdata")

diffn = t(apply(number_clustersO_j, 1, function(x) x - truen))
range(diffn)
matg = data.frame("True" = apply(diffn, 2, function(x) sum(x== 0)/nrow(diffn)),
                  "p1" = apply(diffn, 2, function(x) sum(x== 1)/nrow(diffn)),
                  "p2" = apply(diffn, 2, function(x) sum(x== 2)/nrow(diffn)),
                  "p3" = apply(diffn, 2, function(x) sum(x== 3)/nrow(diffn)),
                  "p4" = apply(diffn, 2, function(x) sum(x== 4)/nrow(diffn)))
datag = data.frame( "met" = rep("gMFM", 30),
                    "value" = c(rep(0,6), rep(1,6), rep(2,6), rep(3,6), rep(4,6)),
                    "group" = rep(1:6,5),
                    "dim" = c(as.matrix(matg)) )
datag

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/Res/DP/scen1_DP_nclusOj.Rdata")
diffndp = t(apply(number_clustersO_j, 1, function(x) x - truen))
range(diffndp)
matgdp = data.frame("True" = apply(diffndp, 2, function(x) sum(x== 0)/nrow(diffndp)),
                    "p1" = apply(diffndp, 2, function(x) sum(x== 1)/nrow(diffndp)),
                    "p2" = apply(diffndp, 2, function(x) sum(x== 2)/nrow(diffndp)),
                    "p3" = apply(diffndp, 2, function(x) sum(x== 3)/nrow(diffndp)),
                    "p4" = apply(diffndp, 2, function(x) sum(x== 4)/nrow(diffndp)))
datagdp = data.frame( "met" = rep("CAM", 30),
                      "value" = c(rep(0,6), rep(1,6), rep(2,6), rep(3,6), rep(4,6)),
                      "group" = rep(1:6,5),
                      "dim" = c(as.matrix(matgdp)) )
datagdp



g1 = ggplot(data = rbind(datag,datagdp), aes(x = value, y = group, size = dim, color=met)) +
  geom_point(alpha = 0.5) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(-1, 10)) +
  theme_bw() +
  scale_colour_manual(values = c("#00AFBB","salmon")) + 
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank() ,
        panel.grid.minor.y = element_blank() ,
        panel.grid.major.x = element_blank() ) +
  scale_x_continuous(labels = c("True", "+1", "+2", "+3", "+4"), name="") +
  scale_y_continuous(trans = "reverse", breaks = 1:6, labels = 1:6, name = "Group")
















