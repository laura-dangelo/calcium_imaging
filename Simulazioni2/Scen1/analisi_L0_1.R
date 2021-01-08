# devtools::install_github("jewellsean/FastLZeroSpikeInference")
library(FastLZeroSpikeInference)




### scen1
lambdas = seq(0.01, 0.08, length.out = 20)
false_negatives = matrix(NA, 50, length(lambdas))
false_positives = matrix(NA, 50, length(lambdas))
TP = matrix(NA, 50, length(lambdas))
TN = matrix(NA, 50, length(lambdas))
FN = matrix(NA, 50, length(lambdas))
FP = matrix(NA, 50, length(lambdas))

for(nsim in 1:50)
{
  filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/data/data_scen1_seed", 
                    nsim, ".Rdata")
  load(file = filename) 
  y = out$y
  g = out$g
  A = out$A
  s = out$s
  spp = which(s>0)
  k = out$k
  n = length(y)
  
  for(i in 1:length(lambdas))
  {
    fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.6, lambda = lambdas[i])
    FN[nsim,i] = sum(sapply(spp, function(x) !(x %in% fit$spikes)))
    FP[nsim,i] = sum(sapply(fit$spikes, function(x) !(x %in% spp)))
    TN[nsim,i] = (n-length(spp)) - FP[nsim,i]
    TP[nsim,i] = length(spp) - FN[nsim,i]
    false_negatives[nsim,i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
    false_positives[nsim,i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
  }
  
}



lambda_best = matrix(NA, 50, 1)
for(nsim in 1:50)
{
  means = apply( cbind(FP[nsim,],FN[nsim,]), 1, mean  )
  lambda_best[nsim] = lambdas[which.min(means)]
}

table(lambda_best)

lambda = 0.028

false_negatives = NULL
false_positives = NULL

for(nsim in 1:50)
{
  filename = paste0("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/data/data_scen1_seed", 
                    nsim, ".Rdata")
  load(file = filename) 
  y = out$y
  g = out$g
  A = out$A
  s = out$s
  spp = which(s>0)
  k = out$k
  n = length(y)
  
  fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.6, lambda = lambda)
  false_negatives[nsim] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
  false_positives[nsim] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
  
  FN[nsim] = sum(sapply(spp, function(x) !(x %in% fit$spikes)))
  FP[nsim] = sum(sapply(fit$spikes, function(x) !(x %in% spp)))
  TN[nsim] = (n-length(spp)) - FP[nsim]
  TP[nsim] = length(spp) - FN[nsim]
}


save(false_negatives, file=paste0("L0_1_false_neg2.Rdata") )
save(false_positives, file=paste0("L0_1_false_pos2.Rdata") )
 
save(FN, file=paste0("L0_1_FN.Rdata") )
save(FP, file=paste0("L0_1_FP.Rdata") )
save(TN, file=paste0("L0_1_TN.Rdata") )
save(TP, file=paste0("L0_1_TP.Rdata") )
