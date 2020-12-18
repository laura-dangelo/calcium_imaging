library(ggplot2)
require(gridExtra)

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/par8/scen1_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/par8/scen1_par8_false_pos.Rdata")

false_positives_gMFM = false_positives
false_negatives_gMFM = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/DP_par8/scen1DP_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/DP_par8/scen1DP_par8_false_pos.Rdata")

false_positives_DP = false_positives
false_negatives_DP = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/L0_1_false_pos.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/L0_1_false_neg.Rdata")

false_positives_L0 = false_positives
false_negatives_L0 = false_negatives


data_scen1 <- data.frame("value" = c( false_positives_gMFM, false_negatives_gMFM, 
                                     false_positives_DP, false_negatives_DP,
                                     false_positives_L0, false_negatives_L0),
                         "rate" = rep(c(rep("False positive rate",50), rep("False negative rate", 50)),3),
                         "model" = c(rep("gMFM", 100), rep("DP", 100), rep("L0", 100))
                         )

p1 = ggplot(data_scen1, aes(x = model, y = value, fill = rate)) + 
  geom_boxplot() +
  facet_wrap(~ rate, scales = "free") +
  xlab("") +
  ylab("Scenario 1") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 





load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/par8/scen2_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/par8/scen2_par8_false_pos.Rdata")

false_positives_gMFM = false_positives
false_negatives_gMFM = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par8/scen2DP_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/DP_par8/scen2DP_par8_false_pos.Rdata")

false_positives_DP = false_positives
false_negatives_DP = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/L0_2_false_pos.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen2/L0_2_false_neg.Rdata")

false_positives_L0 = false_positives
false_negatives_L0 = false_negatives


data_scen2 <- data.frame("value" = c( false_positives_gMFM, false_negatives_gMFM, 
                                      false_positives_DP, false_negatives_DP,
                                      false_positives_L0, false_negatives_L0),
                         "rate" = rep(c(rep("False positive rate",50), rep("False negative rate", 50)),3),
                         "model" = c(rep("gMFM", 100), rep("DP", 100), rep("L0", 100))
)

p2 = ggplot(data_scen2, aes(x = model, y = value, fill = rate)) + 
  geom_boxplot() +
  facet_wrap(~ rate, scales = "free") +
  xlab("") +
  ylab("Scenario 2") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 





load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_false_pos.Rdata")

false_positives_gMFM = false_positives
false_negatives_gMFM = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_false_pos.Rdata")

false_positives_DP = false_positives
false_negatives_DP = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/L0_3_false_pos.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/L0_3_false_neg.Rdata")

false_positives_L0 = false_positives
false_negatives_L0 = false_negatives


data_scen3 <- data.frame("value" = c( false_positives_gMFM, false_negatives_gMFM, 
                                      false_positives_DP, false_negatives_DP,
                                      false_positives_L0, false_negatives_L0),
                         "rate" = rep(c(rep("False positive rate",50), rep("False negative rate", 50)),3),
                         "model" = c(rep("gMFM", 100), rep("DP", 100), rep("L0", 100))
)

p3 = ggplot(data_scen3, aes(x = model, y = value, fill = rate)) + 
  geom_boxplot() +
  facet_wrap(~ rate, scales = "free") +
  xlab("") +
  ylab("Scenario 3") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13)) 




grid.arrange(p1, p2, p3, nrow = 3)
