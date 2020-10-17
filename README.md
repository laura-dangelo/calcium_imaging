# calcium_imaging
La cartella SourceCPP contiene i codici Rcpp dei diversi MCMC:
* calcium_DP_mixture_nonconjugate.cpp : outer spike and slab, Gamma base measure
* calcium_DP_innermixture.cpp: inner : inner spike and slab
* calcium_DP_inner_MFM.cpp : prova con mixture of finite mixtures (problema numerico sui coeff.)
* calcium_CAM.cpp : modello con CAM

La cartella RunR sono i codici R per ottenere le catene + simulazione dei dati e analisi.
