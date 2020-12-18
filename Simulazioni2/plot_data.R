load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen1/data/data_scen1_seed5.Rdata")

y = out$y
g = out$g
n = length(y)
s = out$s

df = data.frame(x = 1:n, y = y, g = g, interval = 1:n)

df_rect = data.frame(start = seq(1, n, by = 5000),
                     end = seq(5000, n, by = 5000),
                     Stimulus = as.factor(c(1,2,3,4,1,2)) )


ggplot(data = df) +
  geom_rect(data = df_rect, inherit.aes = FALSE,
            aes(xmin = start, 
                xmax = end, 
                ymin = -Inf, 
                ymax = Inf, fill = Stimulus), alpha = 0.12 ) +
 # scale_fill_manual(values = cols) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")


15400:17400
ggplot(data = df[22500:24000,]) +
  # geom_rect(data = df_rect, inherit.aes = FALSE,
  #           aes(xmin = start, 
  #               xmax = end, 
  #               ymin = -Inf, 
  #               ymax = Inf, fill = Stimulus), alpha = 0.12 ) +
  # scale_fill_manual(values = cols) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")














load("/home/laura/Documents/Dottorato/2.06 Calcium imaging/calcium_imaging/calcium_imaging/Simulazioni2/Scen3/data/data_scen3_seed5.Rdata")
y = out$y
g = out$g
n = length(y)

df = data.frame(x = 1:n, y = y, g = g, interval = 1:n)

df_rect = data.frame(start = seq(1, n, by = 5000),
                     end = seq(5000, n, by = 5000),
                     Stimulus = as.factor(c(1,2,3,1,2)) )

cols = c("#91ff00", "#00fffb","#ff3700")

ggplot(data = df) +
  geom_rect(data = df_rect, inherit.aes = FALSE,
            aes(xmin = start, 
                xmax = end, 
                ymin = -Inf, 
                ymax = Inf, fill = Stimulus), alpha = 0.12 ) +
  scale_fill_manual(values = cols) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")


