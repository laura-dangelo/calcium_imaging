# devtools::install_github("jewellsean/FastLZeroSpikeInference")
library(FastLZeroSpikeInference)




### scen1
lambdas = seq(0.01, 0.08, length.out = 10)
false_negatives = matrix(NA, 50, length(lambdas))
false_positives = matrix(NA, 50, length(lambdas))

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
    false_negatives[nsim,i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
    false_positives[nsim,i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
  }
  
}



lambda_best = matrix(NA, 50, 1)
for(nsim in 1:50)
{
  means = apply( cbind(false_positives[nsim,],false_negatives[nsim,]), 1, mean  )
  lambda_best[nsim] = lambdas[which.min(means)]
}

table(lambda_best)

lambda = 0.018

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
  
}


save(false_negatives, file=paste0("L0_1_false_neg.Rdata") )
save(false_positives, file=paste0("L0_1_false_pos.Rdata") )
 