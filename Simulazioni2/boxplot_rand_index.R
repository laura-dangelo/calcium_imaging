library(ggplot2)

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/par8/scen1gMFM_par8_randindexD.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/par8/scen1gMFM_par8_randindexO.Rdata")

rand_indexO_gMFM1 = rand_indexO
rand_indexD_gMFM1 = rand_indexD



load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/DP_par8/scen1DP_par8_randindexO.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/DP_par8/scen1DP_par8_randindexD.Rdata")

rand_indexO_DP1 = rand_indexO
rand_indexD_DP1 = rand_indexD




load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/par8/scen2_par8_randindexD.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/par8/scen2_par8_randindexO.Rdata")

rand_indexO_gMFM2 = rand_indexO
rand_indexD_gMFM2 = rand_indexD


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par8/scen2DP_par8_randindexO.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par8/scen2DP_par8_randindexD.Rdata")

rand_indexO_DP2 = rand_indexO
rand_indexD_DP2 = rand_indexD




load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_randindexD.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_randindexO.Rdata")

rand_indexO_gMFM3 = rand_indexO
rand_indexD_gMFM3 = rand_indexD


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_randindexO.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_randindexD.Rdata")

rand_indexO_DP3 = rand_indexO
rand_indexD_DP3 = rand_indexD



data_rand = data.frame("value" = c( rowMeans(rand_indexO_gMFM1), rowMeans(rand_indexO_DP1),
                                    rowMeans(rand_indexO_gMFM2), rowMeans(rand_indexO_DP2),
                                    rowMeans(rand_indexO_gMFM3), rowMeans(rand_indexO_DP3) ),
                       "model" = rep( c( rep("fCAM", 50), rep("CAM", 50) ), 3),
                       "scenario" = c( rep("Scenario 1", 100), rep("Scenario 2", 100), rep("Scenario 3", 100) )
)

ggplot(data_rand, aes(x = model, y = value, fill = scenario)) + 
  geom_boxplot() +
  facet_wrap(~ scenario, scales = "free") +
  xlab("") +
  ylab("Observational clusters") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 



data_randD = data.frame("value" = c( rowMeans(rand_indexD_gMFM1), rowMeans(rand_indexD_DP1),
                                    rowMeans(rand_indexD_gMFM2), rowMeans(rand_indexD_DP2),
                                    rowMeans(rand_indexD_gMFM3), rowMeans(rand_indexD_DP3) ),
                       "model" = rep( c( rep("fCAM", 50), rep("CAM", 50) ), 3),
                       "scenario" = c( rep("Scenario 1", 100), rep("Scenario 2", 100), rep("Scenario 3", 100) )
)

ggplot(data_randD, aes(x = model, y = value, fill = scenario)) + 
  geom_boxplot() +
  facet_wrap(~ scenario, scales = "free") +
  xlab("") +
  ylab("Distributional clusters") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 







data_rand1 = data.frame("value" = c( rowMeans(rand_indexO_gMFM1), rowMeans(rand_indexO_DP1),
                                    rowMeans(rand_indexD_gMFM1), rowMeans(rand_indexD_DP1) ),
                       "model" = rep( c( rep("fCAM", 50), rep("CAM", 50) ), 2),
                       "cluster" = c( rep("Observational clusters", 100), rep("Distributional clusters", 100) )
                      )

g1 = ggplot(data_rand1, aes(x = model, y = value, fill = cluster)) + 
  geom_boxplot() +
  facet_wrap(~ cluster, scales = "free") +
  xlab("") +
  ylab("Scenario 1") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = NA)) 



data_rand2 = data.frame("value" = c( rowMeans(rand_indexO_gMFM2), rowMeans(rand_indexO_DP2),
                                     rowMeans(rand_indexD_gMFM2), rowMeans(rand_indexD_DP2) ),
                        "model" = rep( c( rep("fCAM", 50), rep("CAM", 50) ), 2),
                        "cluster" = c( rep("Observational clusters", 100), rep("Distributional clusters", 100) )
)

g2 = ggplot(data_rand2, aes(x = model, y = value, fill = cluster)) + 
  geom_boxplot() +
  facet_wrap(~ cluster, scales = "free") +
  xlab("") +
  ylab("Scenario 2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = NA)) 


data_rand3 = data.frame("value" = c( rowMeans(rand_indexO_gMFM3), rowMeans(rand_indexO_DP3),
                                     rowMeans(rand_indexD_gMFM3), rowMeans(rand_indexD_DP3) ),
                        "model" = rep( c( rep("fCAM", 50), rep("CAM", 50) ), 2),
                        "cluster" = c( rep("Observational clusters", 100), rep("Distributional clusters", 100) )
)

g3 = ggplot(data_rand3, aes(x = model, y = value, fill = cluster)) + 
  geom_boxplot() +
  facet_wrap(~ cluster, scales = "free") +
  xlab("") +
  ylab("Scenario 3") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = NA)) 


require(gridExtra)
grid.arrange(g1, g2, g3, nrow = 3)
