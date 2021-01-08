plot(out$y, type="l")
plot(out$y[1:3000], type="l")

plot(out$y[20500:23000], type="l")


data <- read.csv("../data/cellula2.csv", header = FALSE)

y_real = c(data$V1)
length(y_real)
rm(list = ("data"))
plot(y_real[20000:25000], type="l")


library(ggplot2)

df= data.frame(x=1:length(y_real), y = y_real)
int = 20600:26000
#int=85500:90000
ggplot(data = df[int,]) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")+
  theme(rect = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))



df= data.frame(x=1:length(out$y), y = out$y)
int = 1:25000
int=5500:7000
ggplot(data = df[int,]) +
  geom_line(aes(x = x, y = y)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_continuous(name = "Time") +
  scale_y_continuous(name = "Calcium level")+
  theme(rect = element_rect(fill="transparent", colour=NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA))
