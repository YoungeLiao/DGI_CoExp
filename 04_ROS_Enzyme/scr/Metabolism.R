# === 1. load data ===
library(readxl)
library(ggplot2)
setwd('/Users/yangliao/Documents/GitHub/OSProtein')
NO3 <- data.frame(read_excel('./04_ROS_Enzyme/data/Metabolism.xlsx', sheet = 'NO3'))
NO2 <- data.frame(read_excel('./04_ROS_Enzyme/data/Metabolism.xlsx', sheet = 'NO2'))
C <- data.frame(read_excel('./04_ROS_Enzyme/data/Metabolism.xlsx', sheet = 'C'))
Biomass <- data.frame(read_excel('./04_ROS_Enzyme/data/Metabolism.xlsx', sheet = 'Biomass'))

# define function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# define theme 
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
    legend.position = c(0.92, 0.9), # 调整legend位置
    legend.background = element_blank(),
    plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                     units="inches"))
mytheme1 <- theme_classic()+
  theme(# plot.title = element_text(hjust = 0.5,size = 16, face = 'bold'),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = 'bold'),
    legend.text = element_text(size = 14),
    # legend.title = element_text(size = 16),
    legend.title = element_blank(),
    legend.position = c(0.1, 0.9), # 调整legend位置
    legend.background = element_blank(),
    plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                     units="inches"))

# === NO3- === 
df2 <- data_summary(NO3, varname = 'value', groupnames = c('Time', 'Light'))
Cust_palette <- c('#587DF7', '#585858', '#FED565') 

p <- ggplot() +
  geom_bar(data=df2,
           mapping = aes(x = Time, y = value, 
                         fill = Light), 
           # color = Light,
           size = 1.2,
           alpha = 0.8,
           position="dodge", stat="identity", width = 0.6
           ) +
  scale_fill_manual(values = Cust_palette) +
  # scale_color_manual(values = Cust_palette) # 设置颜色
  geom_jitter(data=NO3,
                mapping=aes(x=Time, y=value, 
                            color = Light), 
                # shape = 21,
                size = 2,
                # height = 0.02, 
                width = 0.2, 
                alpha = 1) +
  scale_color_manual(values = Cust_palette) + 
  geom_errorbar(data=df2,
                  mapping=aes(x = Time, 
                              ymin = value, 
                              ymax = value + sd, color = Light), 
                  width = 0.6, size = 0.6, alpha = 0.8,
                  position="dodge", stat="identity")
# beautify 
p + scale_y_continuous(expand = c(0, 0.3, 0.1, 0)) + 
  labs(y="Nitrate (mg NO3-N L-1)", x="Time") + 
  mytheme

# === NO2- ===
df2 <- data_summary(NO2, varname = 'value', groupnames = c('Time', 'Light'))
Cust_palette <- c('#587DF7', '#585858', '#FED565') 

p <- ggplot() +
  geom_bar(data=df2,
           mapping = aes(x = Time, y = value, 
                         fill = Light), 
           # color = Light,
           size = 1.2,
           alpha = 0.8,
           position="dodge", stat="identity", width = 0.6
  ) +
  scale_fill_manual(values = Cust_palette) +
  # scale_color_manual(values = Cust_palette) 
  geom_jitter(data=NO2,
              mapping=aes(x=Time, y=value, 
                          color = Light), 
              # shape = 21,
              size = 2,
              # height = 0.02, 
              width = 0.2, 
              alpha = 1) +
  scale_color_manual(values = Cust_palette) + 
  geom_errorbar(data=df2,
                mapping=aes(x = Time, 
                            ymin = value, 
                            ymax = value + sd, color = Light), 
                width = 0.6, size = 0.6, alpha = 0.8,
                position="dodge", stat="identity")
# beautify 
p + scale_y_continuous(expand = c(0, 0.3, 0.1, 0)) + 
  labs(y="Nitrite (mg NO2-N L-1)", x="Time") + 
  mytheme

# === Acetate ===
df2 <- data_summary(C, varname = 'value', groupnames = c('Time', 'Light'))
Cust_palette <- c('#587DF7', '#585858', '#FED565') 

p <- ggplot() +
  geom_bar(data=df2,
           mapping = aes(x = Time, y = value, 
                         fill = Light), 
           # color = Light,
           size = 1.2,
           alpha = 0.8,
           position="dodge", stat="identity", width = 0.6
  ) +
  scale_fill_manual(values = Cust_palette) +
  # scale_color_manual(values = Cust_palette) # 设置颜色
  geom_jitter(data=C,
              mapping=aes(x=Time, y=value, 
                          color = Light), 
              # shape = 21,
              size = 2,
              # height = 0.02, 
              width = 0.2, 
              alpha = 1) +
  scale_color_manual(values = Cust_palette) + 
  geom_errorbar(data=df2,
                mapping=aes(x = Time, 
                            ymin = value, 
                            ymax = value + sd, 
                            color = Light), 
                width = 0.6, size = 0.6,
                alpha = 0.8,
                position="dodge", stat="identity")
# beautify 
p + scale_y_continuous(expand = c(0, 0.3, 0.1, 0)) + 
  labs(y="Acetate (mg C L-1)", x="Time") + 
  mytheme


# Biomass
df2 <- data_summary(Biomass, varname = 'value', groupnames = c('Time', 'Light'))
Cust_palette <- c('#587DF7', '#585858', '#FED565') 

p <- ggplot() +
  geom_bar(data=df2,
           mapping = aes(x = Time, y = value, 
                         fill = Light), 
           # color = Light,
           size = 1.2,
           alpha = 0.8,
           position="dodge", stat="identity", width = 0.6
  ) +
  scale_fill_manual(values = Cust_palette) +
  # scale_color_manual(values = Cust_palette) 
  geom_jitter(data=Biomass,
              mapping=aes(x=Time, y=value, 
                          color = Light), 
              # shape = 21,
              size = 2,
              # height = 0.02, 
              width = 0.2, 
              alpha = 1) +
  scale_color_manual(values = Cust_palette) + 
  geom_errorbar(data=df2,
                mapping=aes(x = Time, 
                            ymin = value, 
                            ymax = value + sd, 
                            color = Light), 
                width = 0.6, size = 0.6,
                alpha = 0.8,
                position="dodge", stat="identity")
# beautify 
p + scale_y_continuous(expand = c(0, 0, 0.1, 0)) + 
  labs(y="Nitrite (mg NO2-N L-1)", x="Time") + 
  mytheme1
