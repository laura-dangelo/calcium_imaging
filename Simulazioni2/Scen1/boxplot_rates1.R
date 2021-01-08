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


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/L0_1_false_pos2.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/L0_1_false_neg2.Rdata")

false_positives_L0 = false_positives
false_negatives_L0 = false_negatives


data_scen1 <- data.frame("value" = c( false_positives_gMFM, false_negatives_gMFM, 
                                      false_positives_DP, false_negatives_DP,
                                      false_positives_L0, false_negatives_L0),
                         "rate" = rep(c(rep("False positive rate",50), rep("False negative rate", 50)),3),
                         "model" = c(rep("fCAM", 100), rep("CAM", 100), rep("L0", 100))
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
        axis.title.y = element_text(size = 13),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = NA))




load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/par8/scen1_par8_FN.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/par8/scen1_par8_FP.Rdata")

rate_gMFM = (FP+FN)/30000


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/DP_par8/scen1DP_par8_FN.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/DP_par8/scen1DP_par8_FP.Rdata")

rate_DP = (FP+FN)/30000

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/L0_1_FP.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/L0_1_FN.Rdata")

rate_L0 = (FP+FN)/30000

data_scen1 <- data.frame("value" = c( false_positives_gMFM, false_negatives_gMFM, rate_gMFM,
                                      false_positives_DP, false_negatives_DP, rate_DP,
                                      false_positives_L0, false_negatives_L0, rate_L0),
                         "rate" = rep(c(rep("False positive rate",50), rep("False negative rate", 50), rep("Misclassification rate",50)),3),
                         "model" = c(rep("gMFM", 150), rep("DP", 150), rep("L0", 150))
)

p11 = ggplot(data_scen1, aes(x = model, y = value, fill = rate)) + 
  geom_boxplot() +
  facet_wrap(~ rate, scales = "free") +
  xlab("") +
  ylab("Scenario 1") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        panel.background = element_rect(fill = "transparent", colour = NA))


