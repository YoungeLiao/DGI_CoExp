# ===== default setting =====
library(readxl)
library(ggplot2)
library(ggpubr)
setwd('/Users/yangliao/Documents/GitHub/OSProtein/04_ROS_Enzyme')
# Cust_palette <- c('#40A9FE', '#252525', '#FEA940') 
# Cust_palette <- c('#65D4CF', '#587DF7', '#FE7875', '#FED565', '#585858') # green, blue, red, yellow, solar
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme<-theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        legend.title = element_blank(),
        # legend.position = c(0.9, 0.9),
        legend.position = 'none',
        legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
mytheme1 <- 
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        legend.title = element_blank(),
        # legend.position = c(0.9, 0.9),
        legend.position = 'right',
        legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
# ===== NO3 =====
# === 1. load data ===
NO3 <- data.frame(read_excel('./data/EPS_SelfCatalysis.xlsx', sheet = 'Removal'))
head(NO3)
NO3


NO3[,6:7] <- apply(NO3[,6:7], 2, as.numeric)*100
class(NO3$Removal)

# === 2. Visualization: nitrate removal ===
Cust_palette <- c('#424242','#82B0FE',  '#FED180') 
p <- ggboxplot(NO3, x="Group", y="Removal", color = "Group",
               palette = Cust_palette, width = 0.5,
               add = "jitter", add.params = list(size = 3, alpha = 0.7),
               # facet.by = "Catalyst",
               bxp.errorbar = T,
               bxp.errorbar.width = 0.4) + 
  ylim(0, 100) + 
  labs(y = 'Removal rate (%)') + mytheme

p <- p + ylim(50, 100) + labs(y = 'Removal rate (%)', x = '') + mytheme
p
