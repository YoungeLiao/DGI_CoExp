# === 1. load data ===
library(readxl)
superoxide <- data.frame(read_excel('./04_ROS_Enzyme/data/gradient_superoxide.xlsx'))
head(superoxide)

# === 2. Visualization ===
library(ggplot2)
library(ggpubr)
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
    # legend.position = c(0.9, 0.9), # 调整legend位置
    legend.position = 'none',
    legend.background = element_blank(),
    plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                     units="inches"))

p <- ggboxplot(superoxide, x="c", y="Superoxide", color = "c",
               palette = "jco", width = 0.5,
               add = "jitter", add.params = list(size = 3, alpha = 0.7),
               bxp.errorbar = T,
               bxp.errorbar.width = 0.4)


# p + stat_compare_means(method = "anova", label.y = 20) + # Add global p-value
#   stat_compare_means(label = "p.signif", method = "t.test", ref.group = "0 mg L-1") 
# my_comparisons <- list(c("0 mg L-1", "100 mg L-1"), c("0 mg L-1", "500 mg L-1"), c("0 mg L-1", "1000 mg L-1"))
my_comparisons <- list(c("0", "100"), c("0", "500"), c("0", "1000"))
p <- p + stat_compare_means(comparisons = my_comparisons, 
                       label = "p.signif", method = "t.test",
                       label.y = c(17000, 19000, 21000))


p + labs(x = 'Nitrate (mg NO3-N L-1)', y = 'Superoxide') + mytheme 
