burnin = 1:100
plot(1:length(run_DP$p[-burnin]), run_DP$p[-burnin], type = "l")
lines(1:length(run_DP$p[-burnin]), cumsum(run_DP$p[-burnin])/1:length(run_DP$p[-burnin]), col =2)

plot(1:length(run_DP$sigma2[-burnin]), run_DP$sigma2[-burnin], type = "l", xlab = "iterazioni", ylab = "sigma2")
lines(1:length(run_DP$sigma2[-burnin]), cumsum(run_DP$sigma2[-burnin])/1:length(run_DP$sigma2[-burnin]), col =2)

plot(1:length(run_DP$tau[-burnin]), run_DP$tau[-burnin], type = "l", xlab = "iterazioni", ylab = "tau2")
lines(1:length(run_DP$tau[-burnin]), cumsum(run_DP$tau[-burnin])/1:length(run_DP$tau[-burnin]), col =2)

plot(1:length(run_DP$b[-burnin]), run_DP$b[-burnin], type = "l", xlab = "iterazioni", ylab = "b")
lines(1:length(run_DP$b[-burnin]), cumsum(run_DP$b[-burnin])/1:length(run_DP$b[-burnin]), col =2)

plot(1:length(run_DP$gamma[-burnin]), run_DP$gamma[-burnin], type = "l", xlab = "iterazioni", ylab = "gamma")
lines(1:length(run_DP$gamma[-burnin]), cumsum(run_DP$gamma[-burnin])/1:length(run_DP$gamma[-burnin]), col =2)


plot(1:length(run_DP$A[2,-burnin]), run_DP$A[2,-burnin], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[3,-burnin]), run_DP$A[3,-burnin], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[4,-burnin]), run_DP$A[4,-burnin], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5,-burnin]), run_DP$A[5,-burnin], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5,-burnin]), run_DP$A[6,-burnin], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5,-burnin]), run_DP$A[7,-burnin], type = "l", xlab = "iterazioni", ylab = "A")
plot(1:length(run_DP$A[5,-burnin]), run_DP$A[8,-burnin], type = "l", xlab = "iterazioni", ylab = "A")


barplot(table(apply(run_DP$clusterO[,-burnin], 2, function(x) length(unique(x)) )))


burnin = 1:1800
AA_DP = matrix(0,length(run_DP$b[-burnin]),n)
for(i in 1:length(run_DP$b[-burnin]))
{
  ii = i + max(burnin)
  AA_DP[i, t(run_DP$clusterO)[ii,] >0] = run_DP$A[run_DP$clusterO[run_DP$clusterO[,ii] >0,ii]+1,ii]
}
save(AA_DP, file = "scen3_AA_DP.Rdata")
est_spikes = colMeans(AA) 
est_spikes[which( apply(t(run_DP$clusterO)[-burnin,], 2, function(x) mean(x != 0))<0.5)] = 0
times = which(est_spikes>0)

sum(sapply(spp, function(x) !(x %in% times))) / length(spp)  ### spikes non identificati: falsi negativi
sum(sapply(times, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi

AA[,which(est_spikes == 0)] = 0
barplot(table( apply(AA, 1, function(x) length(unique(x))) ))
moda = as.numeric(attr(which.max(table( apply(AA, 1, function(x) length(unique(x))) )), "names"))

A_ind = AA[apply(AA, 1, function(x) length(unique( x )))==moda,]
dataa = data.frame(A = A_ind[A_ind>0])
ggplot(data = dataa, aes(x = A)) + 
  geom_histogram(bins = 30, aes(y = ..density..), col = "#00AFBB", fill = "#00AFBB", alpha = 0.3) +   
  stat_density(aes(y = ..density..), fill = 1, alpha = 0, col = 1) + 
  theme_bw() +
  scale_x_continuous(breaks = unip)






#--------------------------------------------#
int = 1:n1
int = n1:(n1+n2)
int = (n1+n2):(n1+n2+n3)
int = (n1+n2+n3):(n1+n2+n3+n4)

subsetAA = AA[,int]

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




barplot(table(apply(run_DP$clusterD, 2, function(x) length(unique(x)) ))) # quanti cluster di distribuzioni
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
df_heat = data.frame(J1 = as.factor(mat_heatmap[,1]),
                     J2 = as.factor(mat_heatmap[,2]),
                     val = c(mat_clusterD)/mat_clusterD[1,1],
                     lab = round( c(mat_clusterD)/mat_clusterD[1,1] , 3) )
df_heat$lab[df_heat$J1==df_heat$J2] = 1:J
df_heat$val[df_heat$J1==df_heat$J2] = NA

ggplot(data = df_heat) +
  geom_tile( aes(x = J1, y = J2, fill = val)) +
  geom_text(aes(x = J1, y = J2, label = lab), size=2.5) +
  scale_fill_gradient(low = magma(3)[3], high = inferno(3)[2]) +
  ylim(rev(levels(df_heat$J2))) 


