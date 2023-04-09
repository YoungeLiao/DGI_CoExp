# === 1. load data ===
library(readxl)
superoxide <- data.frame(read_excel('./04_ROS/data/superoxide.xlsx'))
head(superoxide)
# Samples Time     Rate Light
# 1      B1   12 5.603586  Blue
# 2      B2   12 4.789975  Blue
# 3      B3   12 3.705544  Blue
# 4      D1   12 3.372503  Dark
# 5      D2   12 3.206403  Dark
# 6      D3   12 4.195470  Dark

library(ggplot2)
library(ggpubr)
results <- compare_means(superoxide~Light, data=superoxide, group.by = 'Stage',# ref.group = 'Dark', 
                         method = 't.test')

# === 2. Visualization ===
Cust_palette <- c('#40A9FE', '#252525', '#FEA940') 
# Cust_palette <- c('#65D4CF', '#587DF7', '#FE7875', '#FED565', '#585858') # green, blue, red, yellow, solar
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme<-theme_classic()+
  theme(# plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 14),
    # legend.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.position = c(0.15, 0.9), # 调整legend位置
    legend.background = element_blank(),
    plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                     units="inches"))

p <- ggboxplot(superoxide, x="Stage", y="superoxide", color = "Light",
               palette = Cust_palette, width = 0.5,
               add = "jitter", add.params = list(size = 2, alpha = 0.7),
               bxp.errorbar = T,
               bxp.errorbar.width = 0.4)

p + labs(x = 'Time (h)', y = 'Total ROS (abs410 s-1)') + mytheme 


