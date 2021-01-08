library(ggplot2)
require(gridExtra)

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_false_pos.Rdata")

false_positives_gMFM = false_positives
false_negatives_gMFM = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_false_neg.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_false_pos.Rdata")

false_positives_DP = false_positives
false_negatives_DP = false_negatives


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/L0_1_false_pos2.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/L0_1_false_neg2.Rdata")

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
  ylab("Scenario 1") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(size = 10), axis.title = element_text(size = 10), text = element_text(size = 13),
        axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        rect = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent", colour = NA))




load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_FN.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/par8/scen3_par8_FP.Rdata")

rate_gMFM = (FP+FN)/20000


load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_FN.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/DP_par8/scen3DP_par8_FP.Rdata")

rate_DP = (FP+FN)/20000

load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/L0_3_FP.Rdata")
load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/L0_3_FN.Rdata")

rate_L0 = (FP+FN)/20000

data_scen3 <- data.frame("value" = c( false_positives_gMFM, false_negatives_gMFM, rate_gMFM,
                                      false_positives_DP, false_negatives_DP, rate_DP,
                                      false_positives_L0, false_negatives_L0, rate_L0),
                         "rate" = rep(c(rep("False positive rate",50), rep("False negative rate", 50), rep("Misclassification rate",50)),3),
                         "model" = c(rep("gMFM", 150), rep("DP", 150), rep("L0", 150))
)

p31 = ggplot(data_scen3, aes(x = model, y = value, fill = rate)) + 
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


grid.arrange(p1, p2, p3, nrow = 3)
grid.arrange(p11, p21, p31, nrow = 3)
