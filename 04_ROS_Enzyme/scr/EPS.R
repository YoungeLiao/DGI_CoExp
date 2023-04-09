# ===== default setting =====
library(readxl)
library(ggplot2)
library(ggpubr)

# Cust_palette <- c('#40A9FE', '#252525', '#FEA940') 
# Cust_palette <- c('#65D4CF', '#587DF7', '#FE7875', '#FED565', '#585858') # green, blue, red, yellow, solar
top.mar=0.1
right.mar=0.1
bottom.mar=0.1
left.mar=0.1
mytheme<-theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
       # legend.title = element_blank(),
        # legend.position = c(0.9, 0.9), 
        legend.position = 'none',
        legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
mytheme1 <- theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.text = element_text(size = 14),
        # legend.title = element_text(size = 16),
        legend.title = element_blank(),
        # legend.position = c(0.9, 0.9), # 调整legend位置
        legend.position = 'right',
        legend.background = element_blank(),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
# ===== EPS =====
# === 1. load data ===
EPS_raw <- data.frame(read_excel('./data/EPS_SelfCatalysis.xlsx', sheet = 'EPS'))
library(reshape)
data_melt <- melt(EPS_raw)

custo_color <- c('#CE92D7', '#80CAC4', "#FFED6FFF")

p <- ggplot(data_melt, aes(x=light, y=value/3, fill=variable)) +
  geom_bar(stat="identity",
           width = 0.6) + 
  # labs(title = 'DEGs of yellow light') + 
  scale_fill_manual(values = custo_color) + 
  theme_bw() + 
  mytheme1
p
p +  scale_y_continuous(expand = c(0,0)) + 
  labs(x = '', y = 'Mean concentration (mg L-1)')

## triplicate
# ggbarplot(data_melt, x="light", y="value", add = "mean_se", color = "variable",
#           palette = "jco", position = position_dodge(0.8))
