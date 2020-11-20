# devtools::install_github("jewellsean/FastLZeroSpikeInference")
library(FastLZeroSpikeInference)




### scen1
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/data/y_scen1.Rdata")
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/data/g_scen1.Rdata")
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen1/data/spp_scen1.Rdata")
n = length(y)
false_negatives = NULL
false_positives = NULL
lambdas = seq(0.005, 0.1, length.out = 10)
for(i in 1:10)
{
  fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.6, lambda = lambdas[i])
  false_negatives[i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
  false_positives[i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
}
save(false_negatives, file = "scen1_false_neg_L0.Rdata")
save(false_positives, file = "scen1_false_pos_L0.Rdata")
cbind(lambdas, false_positives, false_negatives)

### scen2
load("home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen2/data/y_scen2.Rdata")
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen2/data/g_scen2.Rdata")
load("~/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen2/data/spp_scen2.Rdata")
n = length(y)
false_negatives2 = NULL
false_positives2 = NULL
lambdas = seq(0.005, 0.1, length.out = 10)
for(i in 1:10)
{
  fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.7, lambda = lambdas[i])
  false_negatives2[i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
  false_positives2[i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
}
cbind(false_positives2, false_negatives2)


false_negatives2 = NULL
false_positives2 = NULL
for(i in 1:50)
{
  fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.7, lambda = 0.026)
  false_negatives2[i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
  false_positives2[i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
}

save(false_negatives2, file = "scen2_false_neg_L0.Rdata")
save(false_positives2, file = "scen2_false_pos_L0.Rdata")
mean(false_positives2)
mean(false_negatives2)



### scen3
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen3/Data/y_scen3.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen3/Data/g_scen3.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen3/Data/spp_scen3.Rdata")
n = length(y)
false_negatives3 = NULL
false_positives3 = NULL
lambdas = seq(0.005, 0.05, length.out = 10)
for(i in 1:10)
{
  fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.6, lambda = lambdas[i])
  false_negatives3[i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
  false_positives3[i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
}
save(false_negatives3, file = "scen3_false_neg_L0.Rdata")
save(false_positives3, file = "scen3_false_pos_L0.Rdata")
cbind(false_positives3, false_negatives3)


### scen4
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen4/data/y_scen1.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen4/data/g_scen1.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni/Scen4/data/spp_scen1.Rdata")
n = length(y)
false_negatives3 = NULL
false_positives3 = NULL
lambdas = seq(0.005, 0.05, length.out = 10)
for(i in 1:10)
{
  fit <- FastLZeroSpikeInference::estimate_spikes(dat = y, gam = 0.6, lambda = lambdas[i])
  false_negatives3[i] = sum(sapply(spp, function(x) !(x %in% fit$spikes))) / length(spp)  ### spikes non identificati: falsi negativi
  false_positives3[i] = sum(sapply(fit$spikes, function(x) !(x %in% spp))) / (n-length(spp)) ### falsi positivi
}
save(false_negatives3, file = "scen3_false_neg_L0.Rdata")
save(false_positives3, file = "scen3_false_pos_L0.Rdata")
cbind(false_positives3, false_negatives3)
