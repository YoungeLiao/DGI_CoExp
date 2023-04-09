# Line plot for system performance 
library(readxl)
library(ggplot2)
library(ggpubr)
rawdata <- read_excel('./04_ROS_Enzyme/data/AddO2.xlsx')
df2 <- as.data.frame(rawdata)
head(df2)
# Time         Group  value
# 1    2         Solar 648.8937
# 2    2         Solar 525.3622
# 3    2         Solar 502.5021
# 4    2 Solar_Co_C3N4 519.5859
# 5    2 Solar_Co_C3N4 522.8411
# 6    2 Solar_Co_C3N4 522.3102

# stat analysis: mean, sd —— error bar
attach(df2)
aov.mean<-aggregate(value,by=list(Group, Time),FUN=mean) # mean
aov.sd<-aggregate(value,by=list(Group, Time),FUN=sd) # sd
detach(df2)

plot_data <- data.frame(aov.mean, sd=aov.sd$x)
col_name <- c('Group', 'Time', 'value', 'sd')
colnames(plot_data) <- col_name
class(plot_data)
plot_data[,2:4] <- as.data.frame(apply(plot_data[,2:4], 2, as.numeric)) 

pale_4 <- c('#999999','#131313', '#ACC5FE', '#587DF7')
pale_2 <- c('#131313', '#587DF7')
p <- ggplot(plot_data, aes(x=Time, y=value, Group = Group)) + 
  geom_line(aes(linetype=Group, color=Group)) +
  geom_point(aes(color=Group),size=4, alpha = 0.7) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=2,
                position=position_dodge(0.05)) +
  # 手动设置line的颜色
  scale_color_manual(values=pale_2) +
  ylab('Nitrate (mg NO3-N·L-1)') + 
  xlab('Time (h)') 
p

p <- p +   
  # theme_bw() + 
  theme_classic() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face='bold'),
        legend.text = element_text(size = 14),
        legend.title = element_blank(),
        legend.position=c(0.2,0.15)) # the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
p

# === significance analysis ===
# data
control <- rep(c(41.07, 25.49, 36.14),1)
test <- rep(c(65.51, 68.57, 71.07), 1)
rawdata <- c(control, test)
group <- c(rep('control', 3), rep('·O2-', 3))
data <- data.frame(rawdata, group)
results <- compare_means(rawdata~group, data=data, # group.by = 'group',# ref.group = 'Dark', 
                         method = 't.test')
results
# .y.     group1  group2      p p.adj p.format p.signif method
# <chr>   <chr>   <chr>   <dbl> <dbl> <chr>    <chr>    <chr> 
#   1 rawdata control ·O2-   0.0108 0.011 0.011    *        T-test
